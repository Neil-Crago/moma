//! Run the moving-axis pipeline on a real 3D track (x,y,z per row).
//!
//! There is no ground truth for real data, so this reports *consistency*: how much
//! the radius, turn count and rotation period change when the method or window
//! changes. A trustworthy result should not depend strongly on those choices.
//!
//! Usage: `cargo run --release --example real_track -- data/pea_exp3box8AT1.csv`
//! (any file readable by `moma::helix::parse_xyz`).

use moma::helix::{auto_axis_3d, dominant_period_3d, moving_axis_3d, parse_xyz, repair_3d, to_cylindrical};
use moma::spiral::{total_turns, unwrap_phase, window_for_period};

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| "data/pea_exp3box8AT1.csv".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let mut pts = parse_xyz(&text).unwrap_or_else(|e| panic!("cannot parse {path}: {e}"));
    let n = pts.len();
    if pts.iter().any(|p| p.iter().any(|v| !v.is_finite())) {
        pts = repair_3d(&pts, 101, 20.0);
    }
    let Some((_, period)) = auto_axis_3d(&pts, 100.0, 3000.0) else {
        println!("{path}: {n} samples: no period found");
        return;
    };
    // Windows of 0.8, 1.0 and 1.25 periods: a method that depends on this choice is not trustworthy.
    let mut rows = Vec::new();
    for scale in [0.8, 1.0, 1.25] {
        let w = window_for_period(period * scale);
        let axis = moving_axis_3d(&pts, w);
        let (polar, _) = to_cylindrical(&pts, &axis);
        let trim = (w / 2).min(n / 4);
        let mut radii: Vec<f64> = polar[trim..n - trim].iter().map(|p| p.radius).collect();
        radii.sort_by(|a, b| a.total_cmp(b));
        let turns = total_turns(&unwrap_phase(&polar)[trim..n - trim]);
        let p2 = dominant_period_3d(&pts[trim..n - trim], 100.0, 3000.0).unwrap_or(f64::NAN);
        rows.push((radii[radii.len() / 2], turns, p2));
    }
    let spread = |f: fn(&(f64, f64, f64)) -> f64| {
        let v: Vec<f64> = rows.iter().map(f).collect();
        let (lo, hi) = v.iter().fold((f64::MAX, f64::MIN), |(l, h), x| (l.min(*x), h.max(*x)));
        100.0 * (hi - lo) / (0.5 * (hi + lo)).max(1e-12)
    };
    println!(
        "{path}: n={n} period={period:.0} radius={:.3} total_turns={:.2} (spread across windows: radius {:.0}%, turns {:.0}%)",
        rows[1].0,
        rows[1].1,
        spread(|r| r.0),
        spread(|r| r.1)
    );
}
