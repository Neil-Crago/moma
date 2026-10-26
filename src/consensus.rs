//! Barycentric consensus for moving origins.
//!
//! Several participants each propose a position. The consensus is their weighted
//! barycentre. Two flavours are provided:
//!
//! * [`ring_barycentre`]: positions live on a ring `Z_m`. The arithmetic mean is
//!   wrong there (the mean of 1 and m-1 is not m/2), so positions are mapped to
//!   angles and the weighted *circular* mean is used.
//! * [`plane_barycentre`]: positions are 2D points; the result is an
//!   [`OriginShift`] relative to a reference point.
//!
//! ```
//! use moma::consensus::ring_barycentre;
//! // 1 and 9 on a ring of 10 are neighbours around 0, not far apart.
//! assert_eq!(ring_barycentre(&[1, 9], &[1.0, 1.0], 10), Some(0));
//! ```

use crate::barycentric::OriginShift;
use std::f64::consts::TAU;

/// Weighted circular mean of ring positions, rounded to the nearest residue.
///
/// Returns `None` when there is nothing to average: empty input, `modulus < 2`,
/// mismatched lengths, non-positive total weight, or when votes cancel exactly
/// (e.g. two equal-weight opposite points), where no barycentre exists.
pub fn ring_barycentre(values: &[u64], weights: &[f64], modulus: u64) -> Option<u64> {
    if values.is_empty() || values.len() != weights.len() || modulus < 2 {
        return None;
    }
    if weights.iter().any(|w| !w.is_finite() || *w < 0.0) {
        return None;
    }
    let total: f64 = weights.iter().sum();
    if total <= 0.0 {
        return None;
    }
    let m = modulus as f64;
    let (mut s, mut c) = (0.0, 0.0);
    for (v, w) in values.iter().zip(weights) {
        let a = ((v % modulus) as f64) / m * TAU;
        s += w * a.sin();
        c += w * a.cos();
    }
    if (s * s + c * c).sqrt() < 1e-9 * total {
        return None;
    }
    let angle = s.atan2(c).rem_euclid(TAU);
    Some(((angle / TAU * m).round() as u64) % modulus)
}

/// Weighted barycentre of 2D points, returned as the shift from `reference`.
/// Returns `None` for empty/mismatched input or non-positive total weight.
pub fn plane_barycentre(
    points: &[(f64, f64)],
    weights: &[f64],
    reference: (f64, f64),
) -> Option<OriginShift> {
    if points.is_empty() || points.len() != weights.len() {
        return None;
    }
    let total: f64 = weights.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return None;
    }
    let (x, y) = points
        .iter()
        .zip(weights)
        .fold((0.0, 0.0), |(x, y), ((px, py), w)| (x + px * w, y + py * w));
    Some(OriginShift { dx: x / total - reference.0, dy: y / total - reference.1 })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wraps_around_zero() {
        assert_eq!(ring_barycentre(&[1, 9], &[1.0, 1.0], 10), Some(0));
        assert_eq!(ring_barycentre(&[2, 4], &[1.0, 1.0], 10), Some(3));
    }

    #[test]
    fn weights_pull_the_centre() {
        let c = ring_barycentre(&[0, 4], &[1.0, 3.0], 16).unwrap();
        assert!(c > 2 && c <= 4, "got {c}");
    }

    #[test]
    fn degenerate_inputs() {
        assert_eq!(ring_barycentre(&[], &[], 10), None);
        assert_eq!(ring_barycentre(&[1], &[1.0, 2.0], 10), None);
        assert_eq!(ring_barycentre(&[1], &[1.0], 1), None);
        assert_eq!(ring_barycentre(&[1], &[0.0], 10), None);
        assert_eq!(ring_barycentre(&[1], &[f64::NAN], 10), None);
        assert_eq!(ring_barycentre(&[0, 5], &[1.0, 1.0], 10), None); // exact cancellation
    }

    #[test]
    fn large_values_reduce_mod() {
        assert_eq!(ring_barycentre(&[u64::MAX, u64::MAX], &[1.0, 1.0], 10), Some(5));
    }

    #[test]
    fn plane_shift() {
        let s = plane_barycentre(&[(0.0, 0.0), (4.0, 2.0)], &[1.0, 1.0], (0.0, 0.0)).unwrap();
        assert_eq!((s.dx, s.dy), (2.0, 1.0));
        assert!(plane_barycentre(&[], &[], (0.0, 0.0)).is_none());
    }
}
