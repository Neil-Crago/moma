//! Synthetic snaking spiral with measurement noise.
//!
//! Generates a spiral (slowly growing radius) around an axis that snakes along a
//! sine path, adds noise, then compares a single fixed centre with moving axes of
//! different window lengths. Reports radius regularity and how well the true
//! rotation rate (0.2 rad/sample) is recovered.
//!
//! Run with: `cargo run --example snaking_spiral`

use moma::spiral::{adaptive_axis, moving_axis, quantise_phase, radius_spread, to_polar, unwrap_phase};

const TRUE_RATE: f64 = 0.2;
const N: usize = 1200;

/// Small deterministic generator, so results are reproducible without a seed API.
struct Lcg(u64);
impl Lcg {
    /// Uniform in [-1, 1).
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

fn generate(noise: f64, rng: &mut Lcg) -> Vec<(f64, f64)> {
    (0..N)
        .map(|i| {
            let t = i as f64;
            let axis = (0.08 * t, 3.0 * (0.01 * t).sin());
            let r = 1.0 + 0.0003 * t; // slowly widening spiral
            (
                axis.0 + r * (TRUE_RATE * t).cos() + noise * rng.next(),
                axis.1 + r * (TRUE_RATE * t).sin() + noise * rng.next(),
            )
        })
        .collect()
}

/// Median absolute deviation of the per-sample phase rate from the truth.
fn rate_error(unwrapped: &[f64]) -> f64 {
    let mut e: Vec<f64> = unwrapped.windows(2).map(|w| ((w[1] - w[0]) - TRUE_RATE).abs()).collect();
    e.sort_by(|a, b| a.partial_cmp(b).unwrap());
    e[e.len() / 2]
}

fn main() {
    let (lo, hi) = (60, N - 60); // trim edge bias
    println!("noise  axis          radius spread   median rate error");
    for noise in [0.0, 0.05, 0.15] {
        let pts = generate(noise, &mut Lcg(42));

        let n = pts.len() as f64;
        let c = pts.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0 / n, a.1 + p.1 / n));
        let fixed = to_polar(&pts, &vec![c; pts.len()]);
        println!(
            "{noise:<5}  fixed         {:<14.3}  {:.3}",
            radius_spread(&fixed[lo..hi]),
            rate_error(&unwrap_phase(&fixed)[lo..hi])
        );

        for window in [11, 31, 63, 127] {
            let polar = to_polar(&pts, &moving_axis(&pts, window));
            println!(
                "{noise:<5}  moving w={window:<4} {:<14.3}  {:.3}",
                radius_spread(&polar[lo..hi]),
                rate_error(&unwrap_phase(&polar)[lo..hi])
            );
        }
        println!();
    }

    let pts = generate(0.05, &mut Lcg(42));
    let polar = to_polar(&pts, &moving_axis(&pts, 31));
    let q = quantise_phase(&polar[lo..hi], 8);
    println!("phase residues mod 8 (first 24): {:?}", &q[..24]);

    // Adaptive window on a chirp whose rotation speed changes 4x along the path.
    println!("\nchirp (rate 0.1 -> 0.4), noise 0.05:");
    let mut rng = Lcg(7);
    let mut phase = 0.0;
    let chirp: Vec<(f64, f64)> = (0..N)
        .map(|i| {
            let t = i as f64;
            phase += 0.1 + 0.3 * t / N as f64;
            (
                0.08 * t + phase.cos() + 0.05 * rng.next(),
                3.0 * (0.01 * t).sin() + phase.sin() + 0.05 * rng.next(),
            )
        })
        .collect();
    for (name, axis) in [
        ("fixed w=31   ", moving_axis(&chirp, 31)),
        ("fixed w=63   ", moving_axis(&chirp, 63)),
        ("adaptive (31)", adaptive_axis(&chirp, 31, 3)),
    ] {
        println!("{name}  radius spread {:.3}", radius_spread(&to_polar(&chirp, &axis)[lo..hi]));
    }
}
