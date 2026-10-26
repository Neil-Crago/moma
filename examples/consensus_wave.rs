//! Real-data check of circular consensus on sperm flagellum centrelines
//! (3D-SpermFlagella, Zenodo 10.5281/zenodo.15299846; see data/README.md).
//!
//! A beating flagellum carries a travelling wave: the beat phase at consecutive
//! points along the tail advances by a roughly constant step. So the phase steps
//! between neighbouring positions are votes for one number, which is exactly what
//! `circular_mean` combines, and its concentration says how much they agree.
//! Three checks, none of which tell the tools the answer:
//!
//! 1. Real tails against the same positions in shuffled order (negative control):
//!    shuffling destroys the wave, so concentration must collapse.
//! 2. Leave-one-out: predict a held-out position's phase from its neighbour plus
//!    the consensus step, against (a) no step and (b) the single previous step.
//! 3. Calibration: does a low concentration warn of bad predictions?
//!
//! Usage: `cargo run --release --example consensus_wave -- <dir containing cells>`

use moma::consensus::circular_mean;
use moma::helix::{fill_gaps_3d, shared_period_3d, V3};
use std::f64::consts::{PI, TAU};
use std::path::{Path, PathBuf};

const PROBES: usize = 48;

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

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn wrap(a: f64) -> f64 {
    (a + PI).rem_euclid(TAU) - PI
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.total_cmp(b));
    v[v.len() / 2]
}

/// Complex Fourier coefficient at angular frequency `w` of each coordinate,
/// after a linear detrend and a Hann taper: `(re, im)` per coordinate.
fn coefficient(track: &[V3], w: f64) -> [(f64, f64); 3] {
    let n = track.len();
    let tm = (n as f64 - 1.0) / 2.0;
    let var_t: f64 = (0..n).map(|i| (i as f64 - tm).powi(2)).sum();
    let mut out = [(0.0, 0.0); 3];
    for k in 0..3 {
        let mean = track.iter().map(|p| p[k]).sum::<f64>() / n as f64;
        let slope = (0..n).map(|i| (i as f64 - tm) * (track[i][k] - mean)).sum::<f64>() / var_t;
        for (i, p) in track.iter().enumerate() {
            let x = p[k] - mean - slope * (i as f64 - tm);
            let hann = 0.5 - 0.5 * (TAU * (i as f64 + 0.5) / n as f64).cos();
            out[k].0 += hann * x * (w * i as f64).cos();
            out[k].1 -= hann * x * (w * i as f64).sin();
        }
    }
    out
}

/// Phase and amplitude of the beat along a fixed direction `u`.
fn phase_amp(c: &[(f64, f64); 3], u: [f64; 3]) -> (f64, f64) {
    let re: f64 = (0..3).map(|k| u[k] * c[k].0).sum();
    let im: f64 = (0..3).map(|k| u[k] * c[k].1).sum();
    (im.atan2(re), re.hypot(im))
}

/// The direction in which the combined coefficients are strongest.
fn main_direction(all: &[[(f64, f64); 3]]) -> [f64; 3] {
    let mut m = [[0.0; 3]; 3];
    for c in all {
        for a in 0..3 {
            for b in 0..3 {
                m[a][b] += c[a].0 * c[b].0 + c[a].1 * c[b].1;
            }
        }
    }
    let mut u = [1.0, 0.7, 0.3];
    for _ in 0..100 {
        let v: Vec<f64> = (0..3).map(|a| (0..3).map(|b| m[a][b] * u[b]).sum()).collect();
        let n = v.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-300);
        u = [v[0] / n, v[1] / n, v[2] / n];
    }
    u
}

struct Cell {
    concentration: f64,
    shuffled: f64,
    err_none: f64,
    err_single: f64,
    err_consensus: f64,
    err_local: f64,
}

/// Steps between consecutive entries of `order`, with weights from amplitudes.
fn steps(phase: &[f64], amp: &[f64], order: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let mut s = Vec::new();
    let mut w = Vec::new();
    for pair in order.windows(2) {
        s.push(wrap(phase[pair[1]] - phase[pair[0]]));
        w.push(amp[pair[0]] * amp[pair[1]]);
    }
    (s, w)
}

fn analyse(dir: &Path, rng: &mut Rng) -> Option<Cell> {
    let m = [matrix(&dir.join("X.csv")), matrix(&dir.join("Y.csv")), matrix(&dir.join("Z.csv"))];
    let (rows, frames) = (m[0].len(), m[0].first()?.len());
    if rows < 40 || frames < 24 || m.iter().any(|a| a.len() != rows || a.iter().any(|r| r.len() != frames)) {
        return None;
    }
    let mut tracks: Vec<Vec<V3>> = Vec::new();
    for j in 0..PROBES {
        let s = ((rows - 1) as f64 * (0.05 + 0.55 * j as f64 / (PROBES - 1) as f64)) as usize;
        let t: Vec<V3> = (0..frames).map(|t| [m[0][s][t], m[1][s][t], m[2][s][t]]).collect();
        if t.iter().filter(|p| p.iter().any(|v| !v.is_finite())).count() * 5 > frames {
            return None;
        }
        let t = fill_gaps_3d(&t);
        if t.iter().any(|p| p.iter().any(|v| !v.is_finite())) {
            return None;
        }
        tracks.push(t);
    }
    let period = shared_period_3d(&tracks, 3.0, (frames as f64 / 3.0).min(40.0))?;
    let w = TAU / period;
    let coeffs: Vec<_> = tracks.iter().map(|t| coefficient(t, w)).collect();
    let u = main_direction(&coeffs);
    let (phase, amp): (Vec<f64>, Vec<f64>) = coeffs.iter().map(|c| phase_amp(c, u)).unzip();

    let order: Vec<usize> = (0..PROBES).collect();
    let (st, wt) = steps(&phase, &amp, &order);
    let concentration = circular_mean(&st, &wt)?.concentration;

    // Negative control: the same positions in a random order, averaged over draws.
    let mut shuffled = 0.0;
    let draws = 20;
    for _ in 0..draws {
        let mut p = order.clone();
        for i in (1..p.len()).rev() {
            p.swap(i, (rng.next() % (i as u64 + 1)) as usize);
        }
        let (s, w) = steps(&phase, &amp, &p);
        shuffled += circular_mean(&s, &w).map_or(0.0, |c| c.concentration) / draws as f64;
    }

    // Leave one interior position out and predict it from its predecessor.
    let (mut e0, mut e1, mut e2) = (Vec::new(), Vec::new(), Vec::new());
    let mut e3 = Vec::new();
    for j in 3..PROBES - 3 {
        let rest: Vec<usize> = (0..PROBES).filter(|&i| i != j).collect();
        // Steps that skip over `j` span two spacings and are excluded.
        let (mut s, mut wgt) = (Vec::new(), Vec::new());
        for pair in rest.windows(2) {
            if pair[1] - pair[0] == 1 {
                s.push(wrap(phase[pair[1]] - phase[pair[0]]));
                wgt.push(amp[pair[0]] * amp[pair[1]]);
            }
        }
        let consensus = circular_mean(&s, &wgt)?.angle;
        let consensus = wrap(consensus);
        let single = wrap(phase[j - 1] - phase[j - 2]);
        // Local consensus: the two steps on each side of the held-out position.
        let near = [(j - 3, j - 2), (j - 2, j - 1), (j + 1, j + 2), (j + 2, j + 3)];
        let ls: Vec<f64> = near.iter().map(|&(a, b)| wrap(phase[b] - phase[a])).collect();
        let lw: Vec<f64> = near.iter().map(|&(a, b)| amp[a] * amp[b]).collect();
        let local = wrap(circular_mean(&ls, &lw)?.angle);
        e3.push(wrap(phase[j] - (phase[j - 1] + local)).abs());
        e0.push(wrap(phase[j] - phase[j - 1]).abs());
        e1.push(wrap(phase[j] - (phase[j - 1] + single)).abs());
        e2.push(wrap(phase[j] - (phase[j - 1] + consensus)).abs());
    }
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    Some(Cell {
        concentration,
        shuffled,
        err_none: mean(&e0),
        err_single: mean(&e1),
        err_consensus: mean(&e2),
        err_local: mean(&e3),
    })
}

fn main() {
    let root = PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| "data".into()));
    let mut dirs = Vec::new();
    cells(&root, &mut dirs);
    let mut rng = Rng(12345);
    let mut results: Vec<Cell> = dirs.iter().filter_map(|d| analyse(d, &mut rng)).collect();
    println!("{} cell folders, {} analysable ({PROBES} positions each, 5% to 60% of the tail)\n", dirs.len(), results.len());
    if results.is_empty() {
        return;
    }
    let n = results.len() as f64;
    let mut real: Vec<f64> = results.iter().map(|c| c.concentration).collect();
    let mut shuf: Vec<f64> = results.iter().map(|c| c.shuffled).collect();
    println!("1. concentration of the phase steps (1 = every step the same)");
    println!("   real order:     median {:.2}", median(&mut real));
    println!("   shuffled order: median {:.2}", median(&mut shuf));
    let above = results.iter().filter(|c| c.concentration > c.shuffled).count() as f64 / n;
    println!("   cells where real order beats its own shuffle: {:.0}%\n", 100.0 * above);

    let avg = |f: fn(&Cell) -> f64| results.iter().map(f).sum::<f64>() / n;
    println!("2. leave-one-out prediction of a held-out position's phase (mean abs error, radians)");
    println!("   no step (copy neighbour):      {:.4}", avg(|c| c.err_none));
    println!("   single previous step:          {:.4}", avg(|c| c.err_single));
    println!("   consensus of all other steps:  {:.4}", avg(|c| c.err_consensus));
    println!("   local consensus (4 nearest):   {:.4}\n", avg(|c| c.err_local));

    results.sort_by(|a, b| a.concentration.total_cmp(&b.concentration));
    println!("3. calibration: consensus prediction error by concentration third");
    let third = results.len() / 3;
    for (label, range) in [("low", 0..third), ("middle", third..2 * third), ("high", 2 * third..results.len())] {
        let g = &results[range];
        let conc = g.iter().map(|c| c.concentration).sum::<f64>() / g.len() as f64;
        let err = g.iter().map(|c| c.err_consensus).sum::<f64>() / g.len() as f64;
        let single = g.iter().map(|c| c.err_single).sum::<f64>() / g.len() as f64;
        println!("   {label:<7} concentration {conc:.4}: consensus error {err:.4}, single-step error {single:.4} ({} cells)", g.len());
    }
}
