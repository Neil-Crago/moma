//! Real-data check on 3D+t human sperm flagellum centrelines (3D-SpermFlagella,
//! Zenodo 10.5281/zenodo.15299846, CC BY 4.0; see data/README.md).
//!
//! Layout of each cell's `X.csv`/`Y.csv`/`Z.csv`: one row per point along the tail
//! (head to tip, about 0.1 um apart), one column per time point. The track of any
//! one point is therefore a 3D curve in time that oscillates at the flagellar beat
//! around a drifting axis, which is the moving-axis problem.
//!
//! There is an internal ground truth: every point of one tail beats at the *same*
//! frequency, and the amplitude grows from the head toward the tip. Neither is
//! given to the tools. This reports, per cell, the period each position finds on
//! its own, the period found by pooling all positions (`shared_period_3d`), and
//! how well they agree.
//!
//! Usage: `cargo run --release --example sperm_beat -- <dir containing cells> [-v]`

use moma::helix::{
    dominant_period_3d, fill_gaps_3d, moving_axis_period_3d, shared_period_3d, to_cylindrical, V3,
};
use std::path::{Path, PathBuf};

fn matrix(path: &Path) -> Vec<Vec<f64>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap_or(f64::NAN)).collect())
        .collect()
}

fn cells(dir: &Path, out: &mut Vec<PathBuf>) {
    if dir.join("X.csv").is_file() && dir.join("Y.csv").is_file() && dir.join("Z.csv").is_file() {
        out.push(dir.to_path_buf());
    }
    if let Ok(rd) = std::fs::read_dir(dir) {
        let mut subs: Vec<_> = rd.flatten().map(|e| e.path()).filter(|p| p.is_dir()).collect();
        subs.sort();
        for s in subs {
            cells(&s, out);
        }
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    v[v.len() / 2]
}

struct Cell {
    name: String,
    frames: usize,
    shared: f64,
    own: Vec<f64>,
    amplitudes: Vec<f64>,
}

fn analyse(dir: &Path) -> Option<Cell> {
    let m = [matrix(&dir.join("X.csv")), matrix(&dir.join("Y.csv")), matrix(&dir.join("Z.csv"))];
    let (rows, frames) = (m[0].len(), m[0].first()?.len());
    if rows < 20 || frames < 24 || m.iter().any(|a| a.len() != rows || a.iter().any(|r| r.len() != frames)) {
        return None;
    }
    // Probe positions along the tail, away from the poorly resolved tip.
    let probes = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75];
    let max_period = (frames as f64 / 3.0).min(40.0);
    let mut tracks: Vec<Vec<V3>> = Vec::new();
    for f in probes {
        let s = ((rows - 1) as f64 * f) as usize;
        let track: Vec<V3> = (0..frames).map(|t| [m[0][s][t], m[1][s][t], m[2][s][t]]).collect();
        let missing = track.iter().filter(|p| p.iter().any(|v| !v.is_finite())).count();
        if missing * 5 > frames {
            continue;
        }
        let pts = fill_gaps_3d(&track);
        if pts.iter().all(|p| p.iter().all(|v| v.is_finite())) {
            tracks.push(pts);
        }
    }
    if tracks.len() < 4 {
        return None;
    }
    let shared = shared_period_3d(&tracks, 3.0, max_period)?;
    let own: Vec<f64> = tracks.iter().filter_map(|t| dominant_period_3d(t, 3.0, max_period)).collect();
    let trim = (shared / 2.0).ceil() as usize + 1;
    if frames <= 2 * trim + 2 {
        return None;
    }
    let amplitudes = tracks
        .iter()
        .map(|t| {
            let axis = moving_axis_period_3d(t, shared);
            let (polar, _) = to_cylindrical(t, &axis);
            let mut r: Vec<f64> = polar[trim..frames - trim].iter().map(|p| p.radius).collect();
            median(&mut r)
        })
        .collect();
    let name = dir.file_name()?.to_string_lossy().into_owned();
    Some(Cell { name, frames, shared, own, amplitudes })
}

/// Relation of a position's own period to the shared one: 1 if within 10%,
/// 2/3 if it is the half/third, 4/5 if double/triple, 0 if unrelated.
fn family(own: f64, shared: f64) -> usize {
    let r = own / shared;
    [(1.0, 1), (0.5, 2), (1.0 / 3.0, 3), (2.0, 4), (3.0, 5)]
        .iter()
        .find(|(t, _)| (r / t - 1.0).abs() < 0.10)
        .map_or(0, |&(_, k)| k)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let root = PathBuf::from(args.next().unwrap_or_else(|| "data".into()));
    let verbose = args.any(|a| a == "-v");
    let mut dirs = Vec::new();
    cells(&root, &mut dirs);
    let results: Vec<Cell> = dirs.iter().filter_map(|d| analyse(d)).collect();
    println!("{} cell folders found, {} analysable\n", dirs.len(), results.len());
    if results.is_empty() {
        return;
    }
    let agree = |c: &Cell| c.own.iter().filter(|&&p| family(p, c.shared) == 1).count() as f64 / c.own.len() as f64;
    if verbose {
        println!("cell                              frames  shared period  own agree  amp head->tip (um)");
        for c in &results {
            println!(
                "{:<32}  {:>6}  {:>13.2}  {:>8.0}%  {:.2} -> {:.2}",
                c.name.chars().take(32).collect::<String>(),
                c.frames,
                c.shared,
                100.0 * agree(c),
                c.amplitudes[0],
                c.amplitudes[c.amplitudes.len() - 1]
            );
        }
        println!();
    }
    let n = results.len() as f64;
    let positions: usize = results.iter().map(|c| c.own.len()).sum();
    let mut counts = [0usize; 6];
    for c in &results {
        for &p in &c.own {
            counts[family(p, c.shared)] += 1;
        }
    }
    let pct = |k: usize| 100.0 * counts[k] as f64 / positions as f64;
    let share = |limit: f64| 100.0 * results.iter().filter(|c| agree(c) >= limit).count() as f64 / n;
    let growing = 100.0
        * results.iter().filter(|c| c.amplitudes[c.amplitudes.len() - 1] > 1.5 * c.amplitudes[0]).count() as f64
        / n;
    let mut shared: Vec<f64> = results.iter().map(|c| c.shared).collect();
    println!("positions whose own strongest period equals the shared period (+-10%): {:.0}%", pct(1));
    println!(
        "positions on a harmonic of it: half {:.0}%, third {:.0}%, double {:.0}%, triple {:.0}%; unrelated {:.0}%",
        pct(2), pct(3), pct(4), pct(5), pct(0)
    );
    println!("cells where >= 75% of positions agree with the shared period: {:.0}%", share(0.75));
    println!("cells whose tip-side amplitude exceeds 1.5x the head-side amplitude: {growing:.0}%");
    println!(
        "shared beat period over cells (samples): median {:.1}, range {:.1} to {:.1}",
        median(&mut shared),
        shared[0],
        shared[shared.len() - 1]
    );
}
