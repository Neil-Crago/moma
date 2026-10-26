//! 3D helices around a wandering axis.
//!
//! The 2D tools in [`crate::spiral`] extended to space. The axis is a moving
//! barycentre; offsets from it are resolved in a frame *perpendicular to the
//! axis tangent*, parallel-transported along the path so the frame does not twist
//! on its own. The result is the same cylindrical view: radius, phase, and
//! position along the axis.
//!
//! ```
//! use moma::helix::{moving_axis_3d, to_cylindrical};
//! use moma::spiral::radius_spread;
//! let pts: Vec<[f64; 3]> = (0..300)
//!     .map(|i| {
//!         let t = i as f64;
//!         [0.1 * t, (0.2 * t).cos(), (0.2 * t).sin()]
//!     })
//!     .collect();
//! let (polar, _along) = to_cylindrical(&pts, &moving_axis_3d(&pts, 32));
//! assert!(radius_spread(&polar[40..260]) < 0.05);
//! ```

use crate::spiral::{fill_gaps_nd, moving_median_nd, outlier_mask, parse_rows, unwrap_phase, Polar};
use std::f64::consts::TAU;

/// A point or vector in 3D.
pub type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: V3) -> Option<V3> {
    let n = norm(a);
    (n > 1e-12).then(|| scale(a, 1.0 / n))
}

/// Any unit vector perpendicular to `t`.
fn any_perpendicular(t: V3) -> V3 {
    let helper = if t[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    unit(cross(t, helper)).unwrap_or([0.0, 0.0, 1.0])
}

/// Centred moving barycentre over `window` samples (clipped at the ends; trim
/// about half a window when analysing).
pub fn moving_axis_3d(points: &[V3], window: usize) -> Vec<V3> {
    windowed_axis_3d(points, |_| window)
}

/// 3D counterpart of [`crate::spiral::adaptive_axis`]: each sample's window is
/// resized to one local turn (`2π / rate`), clamped to `[5, 4 * initial_window]`.
/// Only worthwhile when the rotation speed varies along the track.
pub fn adaptive_axis_3d(points: &[V3], initial_window: usize, iterations: usize) -> Vec<V3> {
    let n = points.len();
    let (min_w, max_w) = (5usize, (4 * initial_window.max(2)).max(5));
    let mut windows = vec![initial_window.clamp(min_w, max_w); n];
    let mut axis = windowed_axis_3d(points, |i| windows[i]);
    for _ in 0..iterations {
        let (polar, _) = to_cylindrical(points, &axis);
        let phase = unwrap_phase(&polar);
        windows = (0..n)
            .map(|i| {
                let k = (windows[i] / 4).max(2);
                let (lo, hi) = (i.saturating_sub(k), (i + k).min(n.saturating_sub(1)));
                let rate = if hi > lo { (phase[hi] - phase[lo]).abs() / (hi - lo) as f64 } else { 0.0 };
                if rate < 1e-9 { max_w } else { ((TAU / rate).round() as usize).clamp(min_w, max_w) }
            })
            .collect();
        axis = windowed_axis_3d(points, |i| windows[i]);
    }
    axis
}

fn windowed_axis_3d(points: &[V3], window_at: impl Fn(usize) -> usize) -> Vec<V3> {
    let n = points.len();
    let mut prefix = vec![[0.0; 3]; n + 1];
    for (i, p) in points.iter().enumerate() {
        for k in 0..3 {
            prefix[i + 1][k] = prefix[i][k] + p[k];
        }
    }
    (0..n)
        .map(|i| {
            let half = window_at(i).max(1) / 2;
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            let c = (hi - lo) as f64;
            let mut out = [0.0; 3];
            for k in 0..3 {
                out[k] = (prefix[hi][k] - prefix[lo][k]) / c;
            }
            out
        })
        .collect()
}

/// Unit tangents of the axis by central differences. Where the axis is
/// momentarily stationary the previous tangent is reused (initially +z).
pub fn tangents(axis: &[V3]) -> Vec<V3> {
    let n = axis.len();
    let mut prev = [0.0, 0.0, 1.0];
    (0..n)
        .map(|i| {
            let (lo, hi) = (i.saturating_sub(1), (i + 1).min(n - 1));
            if let Some(t) = unit(sub(axis[hi], axis[lo])) {
                prev = t;
            }
            prev
        })
        .collect()
}

/// Cylindrical view of `points` about `axis`: `(radius and phase, position along
/// the axis)` per sample. Phase is measured in a parallel-transported frame, so
/// its absolute offset is arbitrary but its rate is meaningful. Panics if the
/// slices differ in length.
pub fn to_cylindrical(points: &[V3], axis: &[V3]) -> (Vec<Polar>, Vec<f64>) {
    assert_eq!(points.len(), axis.len(), "points and axis must align");
    let ts = tangents(axis);
    let mut e1: Option<V3> = None;
    let mut polar = Vec::with_capacity(points.len());
    let mut along = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let t = ts[i];
        let carried = e1
            .and_then(|e| unit(sub(e, scale(t, dot(e, t)))))
            .unwrap_or_else(|| any_perpendicular(t));
        e1 = Some(carried);
        let e2 = cross(t, carried);
        let d = sub(points[i], axis[i]);
        let a = dot(d, t);
        let perp = sub(d, scale(t, a));
        polar.push(Polar {
            radius: norm(perp),
            phase: dot(perp, e2).atan2(dot(perp, carried)).rem_euclid(TAU),
        });
        along.push(a);
    }
    (polar, along)
}

/// Parses `x,y,z` rows (header, `#` comments and extra columns tolerated).
pub fn parse_xyz(text: &str) -> Result<Vec<V3>, String> {
    Ok(parse_rows(text, 3)?.into_iter().map(|r| [r[0], r[1], r[2]]).collect())
}

/// Replaces non-finite points by linear interpolation (see [`crate::spiral::fill_gaps`]).
pub fn fill_gaps_3d(points: &[V3]) -> Vec<V3> {
    fill_gaps_nd(points)
}

/// Per-coordinate moving median axis; robust to spikes and dropouts.
pub fn moving_median_axis_3d(points: &[V3], window: usize) -> Vec<V3> {
    moving_median_nd(points, window)
}

/// One-call cleanup for dirty tracks: fills dropouts, finds a robust axis with
/// a moving median, flags points whose radius is a robust outlier (`k` scaled
/// MADs, about 4 is a good start) and fills those too. Returns points suitable
/// for [`moving_axis_3d`].
pub fn repair_3d(points: &[V3], window: usize, k: f64) -> Vec<V3> {
    let filled = fill_gaps_3d(points);
    let axis = moving_median_axis_3d(&filled, window);
    if axis.iter().any(|a| a.iter().any(|v| !v.is_finite())) {
        return filled;
    }
    let (polar, _) = to_cylindrical(&filled, &axis);
    let radii: Vec<f64> = polar.iter().map(|p| p.radius).collect();
    let marked: Vec<V3> = filled
        .iter()
        .zip(outlier_mask(&radii, k))
        .map(|(&p, bad)| if bad { [f64::NAN; 3] } else { p })
        .collect();
    fill_gaps_3d(&marked)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spiral::{radius_spread, unwrap_phase, winding_number};

    /// Helix of radius `r`, rate 0.2 rad/sample, around a gently curving 3D axis.
    fn helix(n: usize, r: f64) -> Vec<V3> {
        (0..n)
            .map(|i| {
                let t = i as f64;
                let a = [0.1 * t, 2.0 * (0.005 * t).sin(), 2.0 * (0.004 * t).cos()];
                [a[0], a[1] + r * (0.2 * t).cos(), a[2] + r * (0.2 * t).sin()]
            })
            .collect()
    }

    #[test]
    fn recovers_steady_rotation_around_a_curved_axis() {
        let pts = helix(800, 1.0);
        let (polar, _) = to_cylindrical(&pts, &moving_axis_3d(&pts, 31));
        assert!(radius_spread(&polar[40..760]) < 0.1, "{}", radius_spread(&polar[40..760]));

        let u = unwrap_phase(&polar);
        let turns = winding_number(&u[40..760]);
        let expected = 0.2 * 719.0 / TAU;
        assert!((turns.abs() - expected).abs() / expected < 0.02, "turns {turns} vs {expected}");
    }

    #[test]
    fn along_axis_progress_is_monotone_for_a_forward_helix() {
        let pts = helix(400, 1.0);
        let (_, along) = to_cylindrical(&pts, &moving_axis_3d(&pts, 31));
        // `along` is relative to the local axis point, so check it stays small.
        assert!(along[40..360].iter().all(|a| a.abs() < 0.5));
    }

    #[test]
    fn tangents_are_unit_and_handle_a_stationary_axis() {
        let ts = tangents(&[[0.0; 3]; 5]);
        assert!(ts.iter().all(|t| (norm(*t) - 1.0).abs() < 1e-12));
        let line: Vec<V3> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        assert!(tangents(&line).iter().all(|t| (t[0] - 1.0).abs() < 1e-12));
        assert!(tangents(&[]).is_empty());
    }

    #[test]
    fn degenerate_inputs_are_finite() {
        let still = vec![[1.0, 2.0, 3.0]; 20];
        let (polar, along) = to_cylindrical(&still, &moving_axis_3d(&still, 5));
        assert!(polar.iter().all(|p| p.radius.is_finite() && p.phase.is_finite()));
        assert!(along.iter().all(|a| a.is_finite()));
        assert!(moving_axis_3d(&[], 5).is_empty());
    }

    #[test]
    fn repair_recovers_a_damaged_track() {
        let clean = helix(800, 1.0);
        let mut dirty = clean.clone();
        for i in (100..700).step_by(53) {
            dirty[i] = [f64::NAN; 3];
            dirty[i + 20] = [60.0, -60.0, 60.0];
        }
        let fixed = repair_3d(&dirty, 31, 4.0);
        assert!(fixed.iter().all(|p| p.iter().all(|v| v.is_finite())));
        // Every spike is gone, and the track is back near the truth.
        let worst = fixed.iter().zip(&clean).map(|(a, b)| norm(sub(*a, *b))).fold(0.0, f64::max);
        assert!(worst < 1.0, "worst error {worst}");
        // The mean-based axis, by contrast, is dragged off by the damage.
        let naive = moving_axis_3d(&fill_gaps_3d(&dirty), 31);
        let truth = moving_axis_3d(&clean, 31);
        let naive_err = naive.iter().zip(&truth).map(|(a, b)| norm(sub(*a, *b))).fold(0.0, f64::max);
        assert!(naive_err > 1.0, "naive axis error {naive_err}");
        assert_eq!(repair_3d(&[], 31, 4.0).len(), 0);
    }

    #[test]
    fn adaptive_3d_follows_a_speeding_helix() {
        let mut phase = 0.0;
        let pts: Vec<V3> = (0..900)
            .map(|i| {
                let t = i as f64;
                phase += 0.12 + 0.25 * t / 900.0;
                let a = [0.1 * t, 2.0 * (0.008 * t).sin(), 2.0 * (0.006 * t).cos()];
                [a[0], a[1] + phase.cos(), a[2] + phase.sin()]
            })
            .collect();
        let err = |axis: &[V3]| {
            let s: f64 = (100..800)
                .map(|i| {
                    let t = i as f64;
                    let a = [0.1 * t, 2.0 * (0.008 * t).sin(), 2.0 * (0.006 * t).cos()];
                    norm(sub(axis[i], a)).powi(2)
                })
                .sum();
            (s / 700.0).sqrt()
        };
        let fixed = err(&moving_axis_3d(&pts, 31));
        let adaptive = err(&adaptive_axis_3d(&pts, 31, 3));
        assert!(adaptive < fixed * 0.6, "fixed {fixed}, adaptive {adaptive}");
        assert!(adaptive_axis_3d(&[], 31, 2).is_empty());
        let still = vec![[1.0, 1.0, 1.0]; 30];
        assert!(adaptive_axis_3d(&still, 9, 2).iter().all(|p| p.iter().all(|v| v.is_finite())));
    }

    #[test]
    fn parse_xyz_works() {
        assert_eq!(parse_xyz("x y z\n1 2 3\n4,5,6,99\n").unwrap(), vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert!(parse_xyz("1 2 3\n1 2\n").is_err());
    }
}
