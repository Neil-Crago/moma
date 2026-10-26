//! Barycentric consensus: one answer from several competing positions.
//!
//! Several participants (people, sensors, estimators) each propose a position.
//! The consensus is a weighted barycentre, and how it is computed depends on where
//! the positions live:
//!
//! * **On a ring `Z_m`** ([`ring_barycentre`], [`ring_consensus`], [`circular_mean`]):
//!   the arithmetic mean is wrong when positions wrap (the mean of 1 and m-1 is
//!   not m/2), so the weighted *circular* mean is used. Its concentration says
//!   how strongly the participants agree.
//! * **In space** ([`weighted_barycentre`], [`geometric_median`],
//!   [`plane_barycentre`]): the mean is the least-squares answer and is dragged
//!   away by a single wild participant; the geometric median is the robust one.
//! * **For a whole path** ([`consensus_axis`]): several estimates of the same
//!   moving axis (different windows, different methods, different tracks) are
//!   merged sample by sample, and each estimator is weighted by how well it
//!   agrees with the rest, so a method that fails is outvoted instead of
//!   averaged in.
//!
//! ```
//! use moma::consensus::ring_barycentre;
//! // 1 and 9 on a ring of 10 are neighbours around 0, not far apart.
//! assert_eq!(ring_barycentre(&[1, 9], &[1.0, 1.0], 10), Some(0));
//! ```

use crate::barycentric::OriginShift;
use std::f64::consts::TAU;

fn valid_weights(n: usize, weights: &[f64]) -> Option<f64> {
    if n == 0 || weights.len() != n || weights.iter().any(|w| !w.is_finite() || *w < 0.0) {
        return None;
    }
    let total: f64 = weights.iter().sum();
    (total > 0.0).then_some(total)
}

/// A circular mean and how concentrated the votes were around it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circular {
    /// Mean direction in radians, in `[0, TAU)`.
    pub angle: f64,
    /// Length of the weighted mean vector, from 0 (votes cancel) to 1 (all agree).
    pub concentration: f64,
}

/// Weighted circular mean of angles (radians).
///
/// Returns `None` for empty or mismatched input, negative or non-finite weights,
/// non-positive total weight, or votes that cancel exactly (no direction exists).
pub fn circular_mean(angles: &[f64], weights: &[f64]) -> Option<Circular> {
    let total = valid_weights(angles.len(), weights)?;
    if angles.iter().any(|a| !a.is_finite()) {
        return None;
    }
    let (mut s, mut c) = (0.0, 0.0);
    for (a, w) in angles.iter().zip(weights) {
        s += w * a.sin();
        c += w * a.cos();
    }
    let concentration = (s * s + c * c).sqrt() / total;
    if concentration < 1e-9 {
        return None;
    }
    Some(Circular { angle: s.atan2(c).rem_euclid(TAU), concentration })
}

/// A ring consensus and the strength of agreement behind it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RingConsensus {
    /// The consensus residue.
    pub value: u64,
    /// 0 to 1; see [`Circular::concentration`]. Low values mean the participants
    /// are spread around the ring and the consensus is fragile.
    pub concentration: f64,
}

/// Weighted circular mean of ring positions, with its concentration.
///
/// Same `None` cases as [`ring_barycentre`].
pub fn ring_consensus(values: &[u64], weights: &[f64], modulus: u64) -> Option<RingConsensus> {
    if modulus < 2 {
        return None;
    }
    let m = modulus as f64;
    let angles: Vec<f64> = values.iter().map(|v| ((v % modulus) as f64) / m * TAU).collect();
    let c = circular_mean(&angles, weights)?;
    Some(RingConsensus {
        value: ((c.angle / TAU * m).round() as u64) % modulus,
        concentration: c.concentration,
    })
}

/// Weighted circular mean of ring positions, rounded to the nearest residue.
///
/// Returns `None` when there is nothing to average: empty input, `modulus < 2`,
/// mismatched lengths, non-positive total weight, or when votes cancel exactly
/// (e.g. two equal-weight opposite points), where no barycentre exists.
pub fn ring_barycentre(values: &[u64], weights: &[f64], modulus: u64) -> Option<u64> {
    ring_consensus(values, weights, modulus).map(|c| c.value)
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

fn dist<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>().sqrt()
}

/// Weighted mean of points in any dimension. `None` for empty or mismatched
/// input, bad weights, or non-finite coordinates.
pub fn weighted_barycentre<const N: usize>(points: &[[f64; N]], weights: &[f64]) -> Option<[f64; N]> {
    let total = valid_weights(points.len(), weights)?;
    if points.iter().any(|p| p.iter().any(|v| !v.is_finite())) {
        return None;
    }
    let mut out = [0.0; N];
    for (p, w) in points.iter().zip(weights) {
        for k in 0..N {
            out[k] += w * p[k];
        }
    }
    Some(out.map(|v| v / total))
}

/// Weighted geometric median: the point minimising the weighted sum of distances.
///
/// Unlike the mean it is barely moved by a minority of wild points (a participant
/// that is arbitrarily far away has a bounded effect). Uses Weiszfeld's
/// iteration with the Vardi-Zhang correction for iterates that land on a data
/// point. Same `None` cases as [`weighted_barycentre`].
pub fn geometric_median<const N: usize>(points: &[[f64; N]], weights: &[f64]) -> Option<[f64; N]> {
    let mut y = weighted_barycentre(points, weights)?;
    for _ in 0..500 {
        let mut num = [0.0; N];
        let mut den = 0.0;
        let mut pull = [0.0; N];
        let mut on_point = 0.0;
        for (p, w) in points.iter().zip(weights) {
            let d = dist(p, &y);
            if d < 1e-12 {
                on_point += w;
                continue;
            }
            den += w / d;
            for k in 0..N {
                num[k] += w / d * p[k];
                pull[k] += w * (p[k] - y[k]) / d;
            }
        }
        if den == 0.0 {
            return Some(y);
        }
        let t = num.map(|v| v / den);
        let next = if on_point > 0.0 {
            let r = pull.iter().map(|v| v * v).sum::<f64>().sqrt();
            if r <= on_point {
                return Some(y);
            }
            let g = on_point / r;
            let mut n = [0.0; N];
            for k in 0..N {
                n[k] = (1.0 - g) * t[k] + g * y[k];
            }
            n
        } else {
            t
        };
        let step = dist(&next, &y);
        y = next;
        if step < 1e-12 * (1.0 + y.iter().map(|v| v.abs()).sum::<f64>()) {
            break;
        }
    }
    Some(y)
}

/// The merged axis, and which estimators it trusted.
#[derive(Debug, Clone)]
pub struct AxisConsensus<const N: usize> {
    /// The consensus axis, one point per sample.
    pub axis: Vec<[f64; N]>,
    /// Final weight of each estimator, summing to 1. An estimator that disagrees
    /// with the others everywhere ends up with a small weight.
    pub weights: Vec<f64>,
    /// Per sample, the weighted rms distance of the estimators from the
    /// consensus. Large where the estimators disagree, so it shows where the
    /// axis is poorly determined.
    pub spread: Vec<f64>,
}

/// Merges several estimates of the same moving axis into one.
///
/// Each element of `estimates` is a full axis (one point per sample, all the
/// same length), for instance from different windows or methods, or from
/// different tracks of the same structure. At each sample the consensus is the
/// weighted geometric median of the estimates. Estimator weights are then
/// re-fitted by iterative reweighting: an estimator's weight is its `prior`
/// divided by its rms distance to the consensus, so one that is wrong throughout
/// is outvoted rather than averaged in. `prior` defaults to equal weights.
///
/// Returns `None` for no estimates, estimates of different or zero length,
/// non-finite points, or bad priors.
///
/// This cannot help when most estimators share the same bias (for example all
/// use windows that are too short); it protects against a minority failing.
pub fn consensus_axis<const N: usize>(
    estimates: &[Vec<[f64; N]>],
    prior: Option<&[f64]>,
) -> Option<AxisConsensus<N>> {
    let k = estimates.len();
    let len = estimates.first()?.len();
    if len == 0
        || estimates.iter().any(|e| {
            e.len() != len || e.iter().any(|p| p.iter().any(|v| !v.is_finite()))
        })
    {
        return None;
    }
    let ones = vec![1.0; k];
    let prior = prior.unwrap_or(&ones);
    let prior_total = valid_weights(k, prior)?;
    let prior: Vec<f64> = prior.iter().map(|w| w / prior_total).collect();
    let mut weights = prior.clone();
    let mut axis: Vec<[f64; N]> = Vec::new();
    let mut column = vec![[0.0; N]; k];
    for _ in 0..30 {
        axis = (0..len)
            .map(|t| {
                for (c, e) in column.iter_mut().zip(estimates) {
                    *c = e[t];
                }
                geometric_median(&column, &weights).unwrap_or(column[0])
            })
            .collect();
        let d: Vec<f64> = estimates
            .iter()
            .map(|e| (e.iter().zip(&axis).map(|(p, a)| dist(p, a).powi(2)).sum::<f64>() / len as f64).sqrt())
            .collect();
        let floor = (0.05 * d.iter().sum::<f64>() / k as f64).max(1e-12);
        let raw: Vec<f64> = prior.iter().zip(&d).map(|(p, d)| p / d.max(floor)).collect();
        let total: f64 = raw.iter().sum();
        let next: Vec<f64> = raw.iter().map(|w| w / total).collect();
        let change = next.iter().zip(&weights).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        weights = next;
        if change < 1e-6 {
            break;
        }
    }
    let spread = (0..len)
        .map(|t| {
            estimates
                .iter()
                .zip(&weights)
                .map(|(e, w)| w * dist(&e[t], &axis[t]).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    Some(AxisConsensus { axis, weights, spread })
}

/// [`consensus_axis`] for 2D axes given as `(x, y)` tuples, as used by
/// [`crate::spiral`]. The result uses `[x, y]` arrays.
pub fn consensus_axis_2d(
    estimates: &[Vec<(f64, f64)>],
    prior: Option<&[f64]>,
) -> Option<AxisConsensus<2>> {
    let arrays: Vec<Vec<[f64; 2]>> =
        estimates.iter().map(|e| e.iter().map(|&(x, y)| [x, y]).collect()).collect();
    consensus_axis(&arrays, prior)
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

    #[test]
    fn circular_mean_reports_agreement() {
        let c = circular_mean(&[0.1, TAU - 0.1], &[1.0, 1.0]).unwrap();
        assert!(c.angle < 1e-9 || (TAU - c.angle) < 1e-9);
        assert!((c.concentration - 0.1f64.cos()).abs() < 1e-9);
        let tight = ring_consensus(&[3, 3, 3], &[1.0; 3], 12).unwrap();
        assert_eq!(tight.value, 3);
        assert!((tight.concentration - 1.0).abs() < 1e-12);
        let loose = ring_consensus(&[0, 4, 8], &[1.0, 1.0, 1.0], 12);
        assert!(loose.is_none_or(|c| c.concentration < 1e-6));
        let spread = ring_consensus(&[0, 3, 6], &[1.0; 3], 12).unwrap();
        assert!(spread.concentration < 0.5);
        assert!(circular_mean(&[f64::NAN], &[1.0]).is_none());
    }

    #[test]
    fn median_ignores_a_wild_participant() {
        let pts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1000.0, 1000.0]];
        let w = [1.0; 5];
        let mean = weighted_barycentre(&pts, &w).unwrap();
        let med = geometric_median(&pts, &w).unwrap();
        assert!(mean[0] > 100.0);
        assert!(med[0] < 1.5 && med[1] < 1.5, "median {med:?}");
        // Landing on a data point (majority at the origin) must not stall or diverge.
        let heavy = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let m = geometric_median(&heavy, &[1.0; 5]).unwrap();
        assert!(dist(&m, &[0.0, 0.0]) < 1e-9, "median {m:?}");
        assert_eq!(geometric_median::<2>(&[], &[]), None);
        assert_eq!(geometric_median(&[[1.0, 2.0]], &[0.0]), None);
        assert_eq!(geometric_median(&[[3.0, 4.0]], &[2.0]), Some([3.0, 4.0]));
    }

    fn line(offset: f64, wobble: f64, n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let s = i as f64;
                [s * 0.1 + offset + wobble * (s * 1.7).sin(), (s * 0.05).sin() + offset, offset]
            })
            .collect()
    }

    #[test]
    fn consensus_outvotes_a_failed_estimator() {
        let n = 100;
        let truth = line(0.0, 0.0, n);
        let good = [line(0.0, 0.02, n), line(0.0, -0.03, n), line(0.0, 0.04, n), line(0.0, 0.01, n)];
        let bad = line(5.0, 0.0, n);
        let mut all: Vec<Vec<[f64; 3]>> = good.to_vec();
        all.push(bad);
        let c = consensus_axis(&all, None).unwrap();
        let rms = |a: &[[f64; 3]]| {
            (a.iter().zip(&truth).map(|(p, t)| dist(p, t).powi(2)).sum::<f64>() / n as f64).sqrt()
        };
        assert!(rms(&c.axis) < 0.05, "consensus error {}", rms(&c.axis));
        let sum: f64 = c.weights.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9);
        assert!(c.weights[4] < 0.05, "failed estimator weight {}", c.weights[4]);
        assert!(c.weights[..4].iter().all(|w| *w > 0.15));
        assert_eq!(c.spread.len(), n);
        // The plain mean is dragged by the failure.
        let mean: Vec<[f64; 3]> = (0..n)
            .map(|t| weighted_barycentre(&all.iter().map(|e| e[t]).collect::<Vec<_>>(), &[1.0; 5]).unwrap())
            .collect();
        assert!(rms(&mean) > 5.0 * rms(&c.axis));
    }

    #[test]
    fn consensus_axis_validates_input() {
        let a = vec![[0.0, 0.0]; 4];
        assert!(consensus_axis::<2>(&[], None).is_none());
        assert!(consensus_axis(&[a.clone(), vec![[0.0, 0.0]; 3]], None).is_none());
        assert!(consensus_axis(&[vec![[f64::NAN, 0.0]; 4]], None).is_none());
        assert!(consensus_axis::<2>(&[Vec::new()], None).is_none());
        assert!(consensus_axis(&[a.clone(), a.clone()], Some(&[1.0])).is_none());
        assert!(consensus_axis(&[a.clone(), a.clone()], Some(&[0.0, 0.0])).is_none());
        // Identical estimators: consensus is that estimate, equal weights, zero spread.
        let c = consensus_axis(&[a.clone(), a.clone()], None).unwrap();
        assert!(c.spread.iter().all(|s| *s < 1e-12));
        assert!((c.weights[0] - 0.5).abs() < 1e-9);
        let two = consensus_axis_2d(&[vec![(1.0, 2.0); 3], vec![(1.0, 2.0); 3]], None).unwrap();
        assert_eq!(two.axis[0], [1.0, 2.0]);
    }

    #[test]
    fn prior_weights_are_respected() {
        let n = 20;
        let a = line(0.0, 0.0, n);
        let b = line(1.0, 0.0, n);
        let heavy_a = consensus_axis(&[a.clone(), b.clone()], Some(&[10.0, 1.0])).unwrap();
        assert!(heavy_a.weights[0] > heavy_a.weights[1]);
    }
}
