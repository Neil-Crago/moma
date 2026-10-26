//! Does consensus help? Several reasonable axis estimators are run on the synthetic
//! elliptical benchmark (known true axis), then merged with `consensus_axis`.
//!
//! The honest comparison is against the single best estimator (which you would
//! not know in advance on real data), the average estimator, and the plain
//! average of all estimates. Consensus is expected to match or beat the typical
//! estimator and the plain mean, and to track the best one only when the others
//! are not all biased the same way.
//!
//! Run with: `cargo run --release --example consensus_bench`

use moma::consensus::consensus_axis;
use moma::helix::{
    auto_axis_3d, coverage_axis_3d, moving_axis_3d, moving_axis_period_3d, moving_midrange_axis_3d, V3,
};
use moma::spiral::window_for_period;
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

fn track(wobble: f64, noise: f64, rng: &mut Lcg) -> Vec<V3> {
    (0..N)
        .map(|i| {
            let t = i as f64;
            let th = TAU * t / PERIOD + wobble * (TAU * t / 2300.0).sin();
            let a = axis(t);
            [
                a[0] + 3.0 * th.cos() + noise * rng.next(),
                a[1] + noise * rng.next(),
                a[2] + 1.5 * th.sin() + noise * rng.next(),
            ]
        })
        .collect()
}

fn err(est: &[V3], lo: usize, hi: usize) -> f64 {
    let s: f64 = (lo..hi)
        .map(|i| {
            let t = axis(i as f64);
            (0..3).map(|k| (est[i][k] - t[k]).powi(2)).sum::<f64>()
        })
        .sum();
    (s / (hi - lo) as f64).sqrt()
}

fn main() {
    println!("axis rms error against the true axis (lower is better); N = {N}, true period {PERIOD}\n");
    println!("case       noise  best    median  worst   plain mean  consensus  weights (per estimator)");
    let (lo, hi) = (600, N - 600);
    for (case, wobble) in [("steady", 0.0), ("irregular", 1.0), ("reversing", 4.0)] {
        for noise in [0.05, 0.3] {
            let pts = track(wobble, noise, &mut Lcg(7));
            let Some((auto, p)) = auto_axis_3d(&pts, 100.0, 2000.0) else {
                println!("{case:<9}  {noise:<5}  period not found");
                continue;
            };
            let w = window_for_period(p);
            let estimates: Vec<Vec<V3>> = vec![
                moving_axis_3d(&pts, 101),
                moving_axis_3d(&pts, 301),
                moving_axis_period_3d(&pts, p),
                moving_midrange_axis_3d(&pts, w, 0.02),
                coverage_axis_3d(&pts, w, 3 * w, 0.02),
                auto,
            ];
            let mut errs: Vec<f64> = estimates.iter().map(|e| err(e, lo, hi)).collect();
            let plain: Vec<V3> = (0..N)
                .map(|t| {
                    let mut s = [0.0; 3];
                    for e in &estimates {
                        for k in 0..3 {
                            s[k] += e[t][k] / estimates.len() as f64;
                        }
                    }
                    s
                })
                .collect();
            let c = consensus_axis(&estimates, None).expect("consensus");
            let weights: Vec<String> = c.weights.iter().map(|w| format!("{w:.2}")).collect();
            let (plain_err, cons_err) = (err(&plain, lo, hi), err(&c.axis, lo, hi));
            errs.sort_by(|a, b| a.total_cmp(b));
            println!(
                "{case:<9}  {noise:<5}  {:<6.3}  {:<6.3}  {:<6.3}  {plain_err:<10.3}  {cons_err:<9.3}  {}",
                errs[0],
                errs[errs.len() / 2],
                errs[errs.len() - 1],
                weights.join(" ")
            );
        }
    }
    println!("\nestimators, in weight order: fixed 101, fixed 301, period-matched mean, midrange, coverage, auto");
}
