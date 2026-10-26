//! Simulated swimming flagellum: a helical wave on a bending, snaking axis whose
//! wave speed changes along its length, sampled with noise.
//!
//! Everything is built from a known ground truth (axis, local wave rate, radius),
//! so each estimator can be scored on three things a biologist would care about:
//! where the axis runs, how fast the wave turns, and how many turns it makes.
//!
//! Run with: `cargo run --example flagellum`

use moma::helix::{adaptive_axis_3d, moving_axis_3d, to_cylindrical, V3};
use moma::spiral::{backward_steps, quantise_phase, unwrap_phase};
use std::f64::consts::TAU;

const N: usize = 1600;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

/// Ground truth, as a function of position `s` along the flagellum.
fn axis(s: f64) -> V3 {
    // A bending body: amplitude grows towards the tail, like a swimming wave.
    let amp = 1.0 + 4.0 * s / N as f64;
    [0.12 * s, amp * (0.012 * s).sin(), amp * 0.6 * (0.009 * s).cos()]
}
/// Local wave rate (rad/sample): the wave speeds up by 2.5x along the length.
fn rate(s: f64) -> f64 {
    0.14 + 0.21 * s / N as f64
}
fn radius(s: f64) -> f64 {
    0.8 + 0.2 * (0.004 * s).sin()
}

fn unit(a: V3) -> V3 {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn simulate(noise: f64, rng: &mut Lcg) -> (Vec<V3>, Vec<f64>) {
    let mut phase = 0.0;
    let mut truth_phase = Vec::with_capacity(N);
    let pts = (0..N)
        .map(|i| {
            let s = i as f64;
            phase += rate(s);
            truth_phase.push(phase);
            let (a, b) = (axis((s - 1.0).max(0.0)), axis(s + 1.0));
            let t = unit([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
            let e1 = unit(cross(t, [0.0, 0.0, 1.0]));
            let e2 = cross(t, e1);
            let c = axis(s);
            let r = radius(s);
            let mut p = [0.0; 3];
            for k in 0..3 {
                p[k] = c[k] + r * (phase.cos() * e1[k] + phase.sin() * e2[k]) + noise * rng.next();
            }
            p
        })
        .collect();
    (pts, truth_phase)
}

fn rms(v: impl Iterator<Item = f64>) -> f64 {
    let (s, n) = v.fold((0.0, 0), |(s, n), x| (s + x * x, n + 1));
    (s / n.max(1) as f64).sqrt()
}

fn main() {
    let (lo, hi) = (130, N - 130);
    let true_turns = (lo..hi).map(|i| rate(i as f64)).sum::<f64>() / TAU;
    println!("true turns over the analysed span: {true_turns:.2}\n");
    println!("noise  method        axis err  rate err  turns   turns err  backward (mod 8)");

    for noise in [0.0, 0.05, 0.15] {
        let (pts, _) = simulate(noise, &mut Lcg(11));
        let methods: [(&str, Vec<V3>); 4] = [
            ("fixed w=15", moving_axis_3d(&pts, 15)),
            ("fixed w=31", moving_axis_3d(&pts, 31)),
            ("fixed w=63", moving_axis_3d(&pts, 63)),
            ("adaptive 31", adaptive_axis_3d(&pts, 31, 3)),
        ];
        for (name, ax) in methods {
            let (polar, _) = to_cylindrical(&pts, &ax);
            let u = unwrap_phase(&polar);
            let axis_err = rms((lo..hi).map(|i| {
                let t = axis(i as f64);
                ((ax[i][0] - t[0]).powi(2) + (ax[i][1] - t[1]).powi(2) + (ax[i][2] - t[2]).powi(2)).sqrt()
            }));
            // Compare local rate over a 20-sample span; the sign depends on the frame.
            let rate_err = rms((lo..hi).step_by(5).map(|i| {
                ((u[i + 10] - u[i - 10]).abs() / 20.0) - rate(i as f64)
            }));
            let turns = ((u[hi - 1] - u[lo]) / TAU).abs();
            let q = quantise_phase(&polar[lo..hi], 8);
            println!(
                "{noise:<5}  {name:<12}  {axis_err:<8.3}  {rate_err:<8.3}  {turns:<6.2}  {:<9.2}  {}",
                turns - true_turns,
                backward_steps(&q, 8)
            );
        }
        println!();
    }
}
