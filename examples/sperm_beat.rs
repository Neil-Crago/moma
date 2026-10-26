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

struct Rng(u64);
impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        (-2.0 * self.uniform().ln()).sqrt() * (std::f64::consts::TAU * self.uniform()).cos()
    }
}

fn detrend(v: &[f64]) -> Vec<f64> {
    let n = v.len() as f64;
    let tm = (n - 1.0) / 2.0;
    let mean = v.iter().sum::<f64>() / n;
    let var_t: f64 = (0..v.len()).map(|i| (i as f64 - tm).powi(2)).sum();
    let slope = v.iter().enumerate().map(|(i, x)| (i as f64 - tm) * (x - mean)).sum::<f64>() / var_t;
    v.iter().enumerate().map(|(i, x)| x - mean - slope * (i as f64 - tm)).collect()
}

/// Largest periodogram line at integer frequencies `kmin..=kmax` (cycles per
/// record), summed over components, as a share of the total variance.
fn band_share(cols: &[Vec<f64>], kmin: usize, kmax: usize) -> f64 {
    let n = cols[0].len();
    let total: f64 = cols.iter().map(|c| c.iter().map(|v| v * v).sum::<f64>()).sum();
    let mut best = 0.0f64;
    for k in kmin..=kmax {
        let w = std::f64::consts::TAU * k as f64 / n as f64;
        let p: f64 = cols
            .iter()
            .map(|c| {
                let (mut re, mut im) = (0.0, 0.0);
                for (i, v) in c.iter().enumerate() {
                    re += v * (w * i as f64).cos();
                    im += v * (w * i as f64).sin();
                }
                re * re + im * im
            })
            .sum();
        best = best.max(p / (n as f64 * total));
    }
    best
}

/// One-sided p-value: how often AR(1) red noise with each component's variance
/// and lag-1 autocorrelation puts as much of its variance into the band as the
/// track does. Linear drift is removed first; curved drift stays in the track and
/// is absorbed by the red-noise null, which makes the test conservative.
fn red_noise_p(track: &[V3], kmin: usize, kmax: usize, surrogates: usize, seed: u64) -> Option<f64> {
    let n = track.len();
    if kmin > kmax || kmax >= n / 2 {
        return None;
    }
    let cols: Vec<Vec<f64>> = (0..3).map(|k| detrend(&track.iter().map(|p| p[k]).collect::<Vec<_>>())).collect();
    let observed = band_share(&cols, kmin, kmax);
    let params: Vec<(f64, f64)> = cols
        .iter()
        .map(|c| {
            let var = c.iter().map(|v| v * v).sum::<f64>() / n as f64;
            let cov = c.windows(2).map(|w| w[0] * w[1]).sum::<f64>() / n as f64;
            ((cov / var.max(1e-300)).clamp(-0.99, 0.99), var.sqrt())
        })
        .collect();
    let mut rng = Rng(seed);
    let mut at_least = 0usize;
    for _ in 0..surrogates {
        let sur: Vec<Vec<f64>> = params
            .iter()
            .map(|&(phi, sd)| {
                let innov = sd * (1.0 - phi * phi).sqrt();
                let mut x = sd * rng.normal();
                let v: Vec<f64> = (0..n)
                    .map(|_| {
                        x = phi * x + innov * rng.normal();
                        x
                    })
                    .collect();
                detrend(&v)
            })
            .collect();
        if band_share(&sur, kmin, kmax) >= observed {
            at_least += 1;
        }
    }
    Some((at_least as f64 + 1.0) / (surrogates as f64 + 1.0))
}

struct Cell {
    name: String,
    frames: usize,
    shared: f64,
    own: Vec<f64>,
    amplitudes: Vec<f64>,
    capacitating: bool,
    head_period: Option<f64>,
    head_p: Option<f64>,
    beat_p: Option<f64>,
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
    let capacitating = !dir.to_string_lossy().contains("Non-Capacitating");
    // Head path: mean of the first rows, searched for a period slower than the beat.
    let head: Vec<V3> = (0..frames)
        .map(|t| {
            let mut s = [0.0; 3];
            for r in 0..5 {
                for k in 0..3 {
                    s[k] += m[k][r][t] / 5.0;
                }
            }
            s
        })
        .collect();
    let head_period = if head.iter().all(|p| p.iter().all(|v| v.is_finite())) {
        dominant_period_3d(&head, 1.5 * shared, frames as f64 / 2.0)
    } else {
        None
    };
    let seed = frames as u64 * 7919 + rows as u64;
    let nf = frames as f64;
    // Slow band: periods from 1.5x the beat up to half the record.
    let head_p = if head.iter().all(|p| p.iter().all(|v| v.is_finite())) {
        red_noise_p(&head, 2, (nf / (1.5 * shared)).floor() as usize, 300, seed)
    } else {
        None
    };
    // Positive control: a mid-tail track in a band around the beat itself.
    let beat_p = red_noise_p(
        &tracks[tracks.len() / 2],
        (nf / (1.2 * shared)).ceil() as usize,
        (nf / (0.8 * shared)).floor() as usize,
        300,
        seed + 1,
    );
    Some(Cell { name, frames, shared, own, amplitudes, capacitating, head_period, head_p, beat_p })
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
    // The piezo sweeps at 90 Hz and each up or down sweep is one volume (29-45
    // planes at 5000-8000 fps), giving 180 volumes per second.
    const VOLUMES_PER_SECOND: f64 = 180.0;
    println!("\nby condition (beat in Hz assumes {VOLUMES_PER_SECOND} volumes/s):");
    for (label, cap) in [("capacitating", true), ("non-capacitating", false)] {
        let g: Vec<&Cell> = results.iter().filter(|c| c.capacitating == cap).collect();
        if g.is_empty() {
            continue;
        }
        let mut p: Vec<f64> = g.iter().map(|c| c.shared).collect();
        let mut a: Vec<f64> = g.iter().map(|c| c.amplitudes[c.amplitudes.len() - 1]).collect();
        let mut ratio: Vec<f64> =
            g.iter().map(|c| c.amplitudes[c.amplitudes.len() - 1] / c.amplitudes[0]).collect();
        let (p, a, ratio) = (median(&mut p), median(&mut a), median(&mut ratio));
        println!(
            "  {label:<17} n={:<3} beat {:.1} samples = {:.1} Hz, tip amplitude {:.2} um, tip/head {:.1}x",
            g.len(), p, VOLUMES_PER_SECOND / p, a, ratio
        );
    }
    let report = |label: &str, ps: Vec<f64>| {
        if ps.is_empty() {
            return;
        }
        let sig = ps.iter().filter(|&&p| p < 0.05).count() as f64 / ps.len() as f64;
        println!("{label}: {} cells testable, {:.0}% significant at p<0.05 (5% expected from noise alone)", ps.len(), 100.0 * sig);
    };
    println!("\nred-noise (AR(1)) surrogate tests, 300 surrogates per cell:");
    report("  tail beat band, positive control", results.iter().filter_map(|c| c.beat_p).collect());
    report("  head slow band                   ", results.iter().filter_map(|c| c.head_p).collect());
    let mut heads: Vec<f64> =
        results.iter().filter_map(|c| c.head_period.map(|h| h / c.shared)).collect();
    if !heads.is_empty() {
        println!(
            "\nhead-path period (slower than the beat) found in {} cells; ratio to beat period: median {:.1}, range {:.1} to {:.1}",
            heads.len(),
            median(&mut heads),
            heads[0],
            heads[heads.len() - 1]
        );
    }
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
