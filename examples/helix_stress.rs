//! Stress test against ground truth: how close is each method's estimated axis to
//! the true (known) axis of a synthetic spiral, under noise, dropouts and spikes?
//!
//! Run with: `cargo run --example helix_stress`

use moma::helix::{moving_axis_3d, repair_3d, to_cylindrical, V3};
use moma::spiral::{
    adaptive_axis, backward_steps, moving_axis, quantise_phase, refine_axis, to_polar,
    unwrap_phase, winding_number, Polar,
};
use std::f64::consts::TAU;

const N: usize = 1500;
const RATE: f64 = 0.2;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

fn axis2d(t: f64) -> (f64, f64) {
    (0.08 * t, 3.0 * (0.01 * t).sin())
}
fn axis3d(t: f64) -> V3 {
    [0.1 * t, 2.0 * (0.005 * t).sin(), 2.0 * (0.004 * t).cos()]
}
fn radius(t: f64) -> f64 {
    1.0 + 0.0004 * t
}

fn snake2d(noise: f64, rng: &mut Lcg) -> Vec<(f64, f64)> {
    (0..N)
        .map(|i| {
            let t = i as f64;
            let a = axis2d(t);
            (
                a.0 + radius(t) * (RATE * t).cos() + noise * rng.next(),
                a.1 + radius(t) * (RATE * t).sin() + noise * rng.next(),
            )
        })
        .collect()
}

fn helix(noise: f64, rng: &mut Lcg) -> Vec<V3> {
    (0..N)
        .map(|i| {
            let t = i as f64;
            let a = axis3d(t);
            [
                a[0] + noise * rng.next(),
                a[1] + radius(t) * (RATE * t).cos() + noise * rng.next(),
                a[2] + radius(t) * (RATE * t).sin() + noise * rng.next(),
            ]
        })
        .collect()
}

fn rms(v: impl Iterator<Item = f64>) -> f64 {
    let (s, n) = v.fold((0.0, 0), |(s, n), x| (s + x * x, n + 1));
    (s / n.max(1) as f64).sqrt()
}

fn axis_err_2d(axis: &[(f64, f64)], lo: usize, hi: usize) -> f64 {
    rms((lo..hi).map(|i| {
        let t = axis2d(i as f64);
        (axis[i].0 - t.0).hypot(axis[i].1 - t.1)
    }))
}

fn radius_err(polar: &[Polar], lo: usize, hi: usize) -> f64 {
    rms((lo..hi).map(|i| polar[i].radius - radius(i as f64)))
}

fn main() {
    let (lo, hi) = (80, N - 80);

    println!("== 2D, error against the true axis (rms, lower is better) ==");
    println!("noise  method       axis err  radius err  fit rms (self-reported)");
    for noise in [0.0, 0.05, 0.15, 0.3] {
        let pts = snake2d(noise, &mut Lcg(1));
        let row = |name: &str, axis: &[(f64, f64)], fit_rms: Option<f64>| {
            println!(
                "{noise:<5}  {name:<11}  {:<8.3}  {:<10.3}  {}",
                axis_err_2d(axis, lo, hi),
                radius_err(&to_polar(&pts, axis), lo, hi),
                fit_rms.map_or("-".into(), |r| format!("{r:.3}"))
            );
        };
        row("fixed w=31", &moving_axis(&pts, 31), None);
        row("adaptive", &adaptive_axis(&pts, 31, 3), None);
        let fit = refine_axis(&pts, 31, 3);
        row("joint fit", &fit.axis, Some(fit.rms_residual));
        println!("       (true noise rms: {:.3})", noise * (2.0f64 / 3.0).sqrt());
    }

    let expected = RATE * (hi - lo - 1) as f64 / TAU;
    println!("\n== 3D helix, expected turns over the window: {expected:.2} ==");
    println!("case                          axis err  turns   backward steps (mod 8)");
    for (label, noise, damage, repair) in [
        ("clean", 0.0, false, false),
        ("noise 0.1", 0.1, false, false),
        ("noise 0.1 + damage, no repair", 0.1, true, false),
        ("noise 0.1 + damage, repaired", 0.1, true, true),
    ] {
        let mut pts = helix(noise, &mut Lcg(3));
        if damage {
            for i in (300..1300).step_by(97) {
                pts[i] = [f64::NAN; 3];
                pts[i + 40] = [50.0, -50.0, 50.0];
            }
        }
        let used = if repair { repair_3d(&pts, 31, 4.0) } else { pts };
        // Without repair, NaNs would poison everything; substitute zeros so the
        // naive pipeline at least runs, as a careless user might.
        let used: Vec<V3> = used.iter().map(|p| p.map(|v| if v.is_finite() { v } else { 0.0 })).collect();
        let axis = moving_axis_3d(&used, 31);
        let (polar, _) = to_cylindrical(&used, &axis);
        let axis_err = rms((lo..hi).map(|i| {
            let t = axis3d(i as f64);
            let d = [axis[i][0] - t[0], axis[i][1] - t[1], axis[i][2] - t[2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        }));
        let u = unwrap_phase(&polar);
        let q = quantise_phase(&polar[lo..hi], 8);
        println!(
            "{label:<29} {axis_err:<8.3}  {:<6.2}  {}",
            winding_number(&u[lo..hi]).abs(),
            backward_steps(&q, 8)
        );
    }
}
