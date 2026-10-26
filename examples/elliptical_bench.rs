//! Synthetic benchmark for elliptical, irregular rotation (the pea-shoot regime),
//! scored against known ground truth: true axis, true period and true total turns.
//!
//! The orbit is an ellipse (minor/major = 0.5) around a drifting, snaking axis; the
//! rotation rate wobbles and, in the "reversing" case, the sense flips back and forth.
//!
//! Run with: `cargo run --release --example elliptical_bench`

use moma::helix::{
    auto_axis_3d, coverage_axis_3d, moving_axis_3d, moving_midrange_axis_3d, to_cylindrical, V3,
};
use moma::spiral::{total_turns_smoothed, unwrap_phase, window_for_period};
use std::f64::consts::TAU;

const N: usize = 4000;
const PERIOD: f64 = 800.0;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

fn axis(t: f64) -> V3 {
    [1.5 * (TAU * t / 5000.0).sin(), 0.02 * t, 1.0 * (TAU * t / 7000.0).cos()]
}

/// Orbit angle: mean rate 2π/PERIOD plus a sinusoidal wobble of `wobble` radians.
fn angle(t: f64, wobble: f64) -> f64 {
    TAU * t / PERIOD + wobble * (TAU * t / 2300.0).sin()
}

fn track(wobble: f64, noise: f64, rng: &mut Lcg) -> (Vec<V3>, f64) {
    let mut true_turns = 0.0;
    let mut prev = angle(0.0, wobble);
    let pts = (0..N)
        .map(|i| {
            let t = i as f64;
            let th = angle(t, wobble);
            true_turns += (th - prev).abs() / TAU;
            prev = th;
            let a = axis(t);
            [
                a[0] + 3.0 * th.cos() + noise * rng.next(),
                a[1] + noise * rng.next(),
                a[2] + 1.5 * th.sin() + noise * rng.next(),
            ]
        })
        .collect();
    (pts, true_turns)
}

fn rms(v: impl Iterator<Item = f64>) -> f64 {
    let (s, n) = v.fold((0.0, 0), |(s, n), x| (s + x * x, n + 1));
    (s / n.max(1) as f64).sqrt()
}

fn axis_err(est: &[V3], lo: usize, hi: usize) -> f64 {
    rms((lo..hi).map(|i| {
        let t = axis(i as f64);
        (0..3).map(|k| (est[i][k] - t[k]).powi(2)).sum::<f64>().sqrt()
    }))
}

fn main() {
    println!("true period {PERIOD}; ellipse semi-axes 3.0 x 1.5; N = {N}\n");
    println!("case        noise  method            window  axis err  turns (true)   turns err");
    // wobble 0 is steady; 1.0 is irregular but monotonic; 4.0 reverses direction.
    for (case, wobble) in [("steady", 0.0), ("irregular", 1.0), ("reversing", 4.0)] {
        for noise in [0.05, 0.3] {
            let mut rng = Lcg(7);
            let (pts, true_turns) = track(wobble, noise, &mut rng);
            let (lo, hi) = (600, N - 600);
            let auto = auto_axis_3d(&pts, 100.0, 2000.0);
            let mut runs: Vec<(String, usize, Vec<V3>)> = [101, 301]
                .iter()
                .map(|&w| ("fixed".to_string(), w, moving_axis_3d(&pts, w)))
                .collect();
            if let Some((ax, p)) = auto {
                let w = window_for_period(p);
                runs.push((format!("midrange p={p:.0}"), w, moving_midrange_axis_3d(&pts, w, 0.02)));
                runs.push((
                    format!("coverage p={p:.0}"),
                    w,
                    coverage_axis_3d(&pts, w, 3 * w, 0.02),
                ));
                runs.push((format!("auto p={p:.0}"), w, ax));
            }
            for (name, w, ax) in runs {
                let (polar, _) = to_cylindrical(&pts, &ax);
                let u = unwrap_phase(&polar);
                let turns = total_turns_smoothed(&u[lo..hi], 101);
                let truth = true_turns * (hi - lo) as f64 / N as f64;
                println!(
                    "{case:<10}  {noise:<5}  {name:<16}  {w:<6}  {:<8.3}  {turns:>6.1} ({truth:>5.1})  {:>+6.0}%",
                    axis_err(&ax, lo, hi),
                    100.0 * (turns - truth) / truth
                );
            }
        }
    }
}
