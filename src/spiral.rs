//! Moving-axis analysis of spiralling signals.
//!
//! A point tracing a spiral around an axis that itself wanders ("snaking") looks
//! irregular in fixed coordinates. Re-expressed relative to a *moving origin*
//! (the local barycentre of the path), the same data collapses to a steady
//! rotation: near-constant radius and a smoothly advancing phase.
//!
//! ```
//! use moma::spiral::{moving_axis, to_polar, radius_spread};
//! // A clean circle around a drifting centre.
//! let pts: Vec<(f64, f64)> = (0..400)
//!     .map(|i| {
//!         let t = i as f64;
//!         (0.05 * t + (0.2 * t).cos(), (0.2 * t).sin())
//!     })
//!     .collect();
//! let axis = moving_axis(&pts, 32);
//! let polar = to_polar(&pts, &axis);
//! assert!(radius_spread(&polar[40..360]) < 0.05);
//! ```

use std::f64::consts::TAU;

/// A point expressed relative to a local axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Polar {
    pub radius: f64,
    /// Angle in `[0, TAU)`.
    pub phase: f64,
}

/// Centred moving barycentre of `points` over `window` samples (clipped at the
/// ends). For best results `window` should span about one full turn; near the
/// ends the clipped window is biased, so trim a half-window when analysing.
pub fn moving_axis(points: &[(f64, f64)], window: usize) -> Vec<(f64, f64)> {
    windowed_axis(points, |_| window)
}

/// Axis whose window adapts to the local rotation period.
///
/// Starts from a fixed `initial_window`, measures the local angular rate of the
/// path about that axis, then resizes each sample's window to one turn
/// (`2π / |rate|`), clamped to `[5, 4 * initial_window]`. Repeats `iterations`
/// times (2-3 is plenty). Where the rotation stalls the window hits the cap.
pub fn adaptive_axis(
    points: &[(f64, f64)],
    initial_window: usize,
    iterations: usize,
) -> Vec<(f64, f64)> {
    let n = points.len();
    let (min_w, max_w) = (5usize, (4 * initial_window.max(2)).max(5));
    let mut windows = vec![initial_window.clamp(min_w, max_w); n];
    let mut axis = windowed_axis(points, |i| windows[i]);
    for _ in 0..iterations {
        let phase = unwrap_phase(&to_polar(points, &axis));
        let new: Vec<usize> = (0..n)
            .map(|i| {
                // Rate over a span of about half the current window, so noise
                // averages out without smearing real changes in speed.
                let k = (windows[i] / 4).max(2);
                let (lo, hi) = (i.saturating_sub(k), (i + k).min(n.saturating_sub(1)));
                let rate = if hi > lo { (phase[hi] - phase[lo]).abs() / (hi - lo) as f64 } else { 0.0 };
                if rate < 1e-9 {
                    max_w
                } else {
                    ((TAU / rate).round() as usize).clamp(min_w, max_w)
                }
            })
            .collect();
        windows = new;
        axis = windowed_axis(points, |i| windows[i]);
    }
    axis
}

fn windowed_axis(points: &[(f64, f64)], window_at: impl Fn(usize) -> usize) -> Vec<(f64, f64)> {
    let n = points.len();
    let mut px = Vec::with_capacity(n + 1);
    let mut py = Vec::with_capacity(n + 1);
    px.push(0.0);
    py.push(0.0);
    for &(x, y) in points {
        px.push(px.last().unwrap() + x);
        py.push(py.last().unwrap() + y);
    }
    (0..n)
        .map(|i| {
            let half = window_at(i).max(1) / 2;
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            let k = (hi - lo) as f64;
            ((px[hi] - px[lo]) / k, (py[hi] - py[lo]) / k)
        })
        .collect()
}

/// Expresses each point in polar form around the matching axis point.
/// Panics if the slices differ in length.
pub fn to_polar(points: &[(f64, f64)], axis: &[(f64, f64)]) -> Vec<Polar> {
    assert_eq!(points.len(), axis.len(), "points and axis must align");
    points
        .iter()
        .zip(axis)
        .map(|(&(x, y), &(ax, ay))| {
            let (dx, dy) = (x - ax, y - ay);
            Polar { radius: dx.hypot(dy), phase: dy.atan2(dx).rem_euclid(TAU) }
        })
        .collect()
}

/// Unwraps phases so consecutive samples differ by less than half a turn,
/// giving a continuous, monotone-if-rotating angle.
pub fn unwrap_phase(polar: &[Polar]) -> Vec<f64> {
    let mut out = Vec::with_capacity(polar.len());
    let mut offset = 0.0;
    let mut prev: Option<f64> = None;
    for p in polar {
        if let Some(q) = prev {
            let d = p.phase - q;
            if d > TAU / 2.0 {
                offset -= TAU;
            } else if d < -TAU / 2.0 {
                offset += TAU;
            }
        }
        out.push(p.phase + offset);
        prev = Some(p.phase);
    }
    out
}

/// Quantises each phase to a residue in `Z_modulus`. With the moving axis as the
/// origin, this is the MOMA view of the spiral. Returns an empty vec if
/// `modulus == 0`.
pub fn quantise_phase(polar: &[Polar], modulus: u64) -> Vec<u64> {
    if modulus == 0 {
        return Vec::new();
    }
    polar
        .iter()
        .map(|p| ((p.phase / TAU * modulus as f64).floor() as u64) % modulus)
        .collect()
}

/// Coefficient of variation of the radius (std dev / mean); 0 for a perfect
/// circle. Returns `f64::INFINITY` for empty input or zero mean radius.
pub fn radius_spread(polar: &[Polar]) -> f64 {
    if polar.is_empty() {
        return f64::INFINITY;
    }
    let n = polar.len() as f64;
    let mean = polar.iter().map(|p| p.radius).sum::<f64>() / n;
    if mean <= 0.0 {
        return f64::INFINITY;
    }
    let var = polar.iter().map(|p| (p.radius - mean).powi(2)).sum::<f64>() / n;
    var.sqrt() / mean
}

/// Nearest odd window length to `period` (at least 3). A centred moving mean over
/// exactly one period removes that period and all its harmonics, however
/// non-circular the orbit, which makes this the right window for the axis.
pub fn window_for_period(period: f64) -> usize {
    if !period.is_finite() {
        return 3;
    }
    ((2.0 * ((period - 1.0) / 2.0).round() + 1.0).max(3.0)) as usize
}

/// Dominant oscillation period (in samples) of a 2D track, searched between
/// `min_period` and `max_period`. A linear drift is removed and a Hann taper
/// applied before a fine periodogram; the strongest interior spectral peak wins.
/// Returns `None` for short or non-finite input, or an empty search range
/// (`max_period` is capped at half the track length so two cycles are seen).
pub fn dominant_period(points: &[(f64, f64)], min_period: f64, max_period: f64) -> Option<f64> {
    let xs: Vec<f64> = points.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = points.iter().map(|p| p.1).collect();
    dominant_period_cols(&[xs, ys], min_period, max_period)
}

/// Axis from a moving mean whose window is set from the measured dominant
/// period. Returns the axis and the period used.
pub fn auto_axis(
    points: &[(f64, f64)],
    min_period: f64,
    max_period: f64,
) -> Option<(Vec<(f64, f64)>, f64)> {
    let p = dominant_period(points, min_period, max_period)?;
    Some((moving_axis(points, window_for_period(p)), p))
}

pub(crate) fn dominant_period_cols(
    cols: &[Vec<f64>],
    min_period: f64,
    max_period: f64,
) -> Option<f64> {
    let n = cols.first()?.len();
    if n < 8 || cols.iter().any(|c| c.len() != n || c.iter().any(|v| !v.is_finite())) {
        return None;
    }
    let max_p = max_period.min(n as f64 / 2.0);
    if min_period.is_nan() || max_p.is_nan() || min_period < 2.0 || max_p <= min_period {
        return None;
    }
    let prepared: Vec<Vec<f64>> = cols.iter().map(|c| detrend_and_taper(c)).collect();
    let (f_lo, f_hi) = (1.0 / max_p, 1.0 / min_period);
    let step = 1.0 / (4.0 * n as f64);
    let count = ((f_hi - f_lo) / step).ceil() as usize + 1;
    let power: Vec<f64> = (0..count)
        .map(|s| {
            let w = TAU * (f_lo + s as f64 * step);
            prepared
                .iter()
                .map(|c| {
                    let (mut re, mut im) = (0.0, 0.0);
                    for (i, v) in c.iter().enumerate() {
                        let (sn, cs) = (w * i as f64).sin_cos();
                        re += v * cs;
                        im += v * sn;
                    }
                    re * re + im * im
                })
                .sum()
        })
        .collect();
    let interior_peak = (1..count.saturating_sub(1))
        .filter(|&s| power[s] > power[s - 1] && power[s] >= power[s + 1])
        .max_by(|&a, &b| power[a].total_cmp(&power[b]));
    let best = interior_peak.or_else(|| (0..count).max_by(|&a, &b| power[a].total_cmp(&power[b])))?;
    (power[best] > 0.0).then(|| 1.0 / (f_lo + best as f64 * step))
}

fn detrend_and_taper(c: &[f64]) -> Vec<f64> {
    let n = c.len();
    let tm = (n as f64 - 1.0) / 2.0;
    let mean = c.iter().sum::<f64>() / n as f64;
    let var_t: f64 = (0..n).map(|i| (i as f64 - tm).powi(2)).sum();
    let slope = (0..n).map(|i| (i as f64 - tm) * (c[i] - mean)).sum::<f64>() / var_t.max(1e-12);
    (0..n)
        .map(|i| {
            let hann = 0.5 * (1.0 - (TAU * i as f64 / (n as f64 - 1.0)).cos());
            (c[i] - mean - slope * (i as f64 - tm)) * hann
        })
        .collect()
}

/// Total rotation travelled, in turns, ignoring direction (sum of |Δphase| / 2π
/// over an unwrapped phase). Unlike [`winding_number`] it is not cancelled by
/// reversals, so it suits motion whose rotation sense changes; it also counts
/// jitter, so smooth noisy data first.
pub fn total_turns(unwrapped: &[f64]) -> f64 {
    unwrapped.windows(2).map(|w| (w[1] - w[0]).abs()).sum::<f64>() / TAU
}

/// [`total_turns`] after smoothing the unwrapped phase with a centred moving mean
/// of `window` samples, so noise-driven back-and-forth steps do not inflate the
/// count. Choose a window well below the shortest real reversal timescale.
pub fn total_turns_smoothed(unwrapped: &[f64], window: usize) -> f64 {
    total_turns(&smooth_symmetric(unwrapped, window))
}

/// Net number of full turns (signed) in an unwrapped phase series.
pub fn winding_number(unwrapped: &[f64]) -> f64 {
    match (unwrapped.first(), unwrapped.last()) {
        (Some(a), Some(b)) => (b - a) / TAU,
        _ => 0.0,
    }
}

/// Signed shortest steps between consecutive residues on `Z_modulus`
/// (range `(-m/2, m/2]`). Empty if `modulus < 2`. A steady rotation gives a
/// steady run of small positive steps.
pub fn quantised_steps(residues: &[u64], modulus: u64) -> Vec<i64> {
    if modulus < 2 {
        return Vec::new();
    }
    let m = modulus as i64;
    residues
        .windows(2)
        .map(|w| {
            let d = ((w[1] % modulus) as i64 - (w[0] % modulus) as i64).rem_euclid(m);
            if d > m / 2 { d - m } else { d }
        })
        .collect()
}

/// Number of ring steps that move backwards: reversals of rotation, or jitter
/// from noise. A cheap slip signature.
pub fn backward_steps(residues: &[u64], modulus: u64) -> usize {
    quantised_steps(residues, modulus).iter().filter(|&&d| d < 0).count()
}

/// Replaces non-finite points (gaps) by linear interpolation between the
/// nearest valid neighbours, holding the nearest value at the ends. Intended
/// for axis estimation; it does not restore the lost rotation.
pub fn fill_gaps(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let rows: Vec<[f64; 2]> = points.iter().map(|p| [p.0, p.1]).collect();
    fill_gaps_nd(&rows).into_iter().map(|r| (r[0], r[1])).collect()
}

/// Dimension-generic [`fill_gaps`]; a row is a gap if any component is non-finite.
pub(crate) fn fill_gaps_nd<const D: usize>(points: &[[f64; D]]) -> Vec<[f64; D]> {
    let valid = |p: &[f64; D]| p.iter().all(|v| v.is_finite());
    let idx: Vec<usize> = (0..points.len()).filter(|&i| valid(&points[i])).collect();
    if idx.is_empty() {
        return points.to_vec();
    }
    (0..points.len())
        .map(|i| {
            if valid(&points[i]) {
                return points[i];
            }
            let pos = idx.partition_point(|&j| j < i);
            match (pos.checked_sub(1).map(|p| idx[p]), idx.get(pos).copied()) {
                (Some(a), Some(b)) => {
                    let f = (i - a) as f64 / (b - a) as f64;
                    let mut out = [0.0; D];
                    for (k, o) in out.iter_mut().enumerate() {
                        *o = points[a][k] + f * (points[b][k] - points[a][k]);
                    }
                    out
                }
                (Some(a), None) => points[a],
                (None, Some(b)) => points[b],
                (None, None) => unreachable!("idx is non-empty"),
            }
        })
        .collect()
}

/// Moving *median* axis: per-coordinate median over a centred window. A full
/// turn of a circle has median equal to its centre, and unlike
/// [`moving_axis`] a handful of spikes cannot drag it away, so use it as the
/// first pass before [`reject_outliers`]. Non-finite values are ignored.
pub fn moving_median_axis(points: &[(f64, f64)], window: usize) -> Vec<(f64, f64)> {
    let rows: Vec<[f64; 2]> = points.iter().map(|p| [p.0, p.1]).collect();
    moving_median_nd(&rows, window).into_iter().map(|r| (r[0], r[1])).collect()
}

pub(crate) fn moving_median_nd<const D: usize>(points: &[[f64; D]], window: usize) -> Vec<[f64; D]> {
    let n = points.len();
    let half = window.max(1) / 2;
    (0..n)
        .map(|i| {
            let (lo, hi) = (i.saturating_sub(half), (i + half + 1).min(n));
            let mut out = [f64::NAN; D];
            for (k, o) in out.iter_mut().enumerate() {
                let vals: Vec<f64> =
                    points[lo..hi].iter().map(|p| p[k]).filter(|v| v.is_finite()).collect();
                if !vals.is_empty() {
                    *o = median(&vals);
                }
            }
            out
        })
        .collect()
}

/// Marks points whose radius about `axis` is a robust outlier (more than `k`
/// scaled MADs from the median radius) as gaps and fills them. Use an axis from
/// a first pass; `k` of about 4 is a sensible start.
pub fn reject_outliers(points: &[(f64, f64)], axis: &[(f64, f64)], k: f64) -> Vec<(f64, f64)> {
    let r: Vec<f64> = to_polar(points, axis).iter().map(|p| p.radius).collect();
    let bad = outlier_mask(&r, k);
    let marked: Vec<(f64, f64)> = points
        .iter()
        .zip(&bad)
        .map(|(&p, &b)| if b { (f64::NAN, f64::NAN) } else { p })
        .collect();
    fill_gaps(&marked)
}

/// Flags radii more than `k` scaled MADs from the median; non-finite radii are
/// flagged too.
pub(crate) fn outlier_mask(radii: &[f64], k: f64) -> Vec<bool> {
    let finite: Vec<f64> = radii.iter().copied().filter(|x| x.is_finite()).collect();
    if finite.is_empty() {
        return vec![true; radii.len()];
    }
    let med = median(&finite);
    let dev: Vec<f64> = finite.iter().map(|x| (x - med).abs()).collect();
    let mad = (1.4826 * median(&dev)).max(1e-12);
    radii.iter().map(|x| !x.is_finite() || (x - med).abs() > k * mad).collect()
}

pub(crate) fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.total_cmp(b));
    let n = s.len();
    if n == 0 { 0.0 } else if n % 2 == 1 { s[n / 2] } else { 0.5 * (s[n / 2 - 1] + s[n / 2]) }
}

/// Result of a joint axis/phase fit.
#[derive(Debug, Clone)]
pub struct Fit {
    pub axis: Vec<(f64, f64)>,
    pub polar: Vec<Polar>,
    /// Root-mean-square distance between the data and the fitted
    /// axis-plus-smooth-rotation model: an honest per-sample uncertainty.
    pub rms_residual: f64,
}

/// Joint fit of axis and rotation.
///
/// Starting from [`adaptive_axis`], alternately (1) smooths radius and unwrapped
/// phase into a model rotation and (2) re-estimates the axis from the data with
/// that rotation subtracted. Because the rotation is removed, the axis can be
/// estimated over a shorter window than one turn, so it follows a snaking axis
/// more closely and is less biased by a growing radius.
pub fn refine_axis(points: &[(f64, f64)], window: usize, iterations: usize) -> Fit {
    let n = points.len();
    let mut axis = adaptive_axis(points, window, 2);
    let model = |axis: &[(f64, f64)]| {
        let polar = to_polar(points, axis);
        let phase = smooth_symmetric(&unwrap_phase(&polar), window);
        let radius = smooth_symmetric(&polar.iter().map(|p| p.radius).collect::<Vec<_>>(), window);
        (0..n).map(|i| (radius[i] * phase[i].cos(), radius[i] * phase[i].sin())).collect::<Vec<_>>()
    };
    for _ in 0..iterations {
        let offset = model(&axis);
        let resid: Vec<(f64, f64)> =
            (0..n).map(|i| (points[i].0 - offset[i].0, points[i].1 - offset[i].1)).collect();
        axis = windowed_axis(&resid, |_| (window / 2).max(3));
    }
    let offset = model(&axis);
    let rms_residual = if n == 0 {
        0.0
    } else {
        ((0..n)
            .map(|i| {
                (points[i].0 - axis[i].0 - offset[i].0).powi(2)
                    + (points[i].1 - axis[i].1 - offset[i].1).powi(2)
            })
            .sum::<f64>()
            / n as f64)
            .sqrt()
    };
    Fit { polar: to_polar(points, &axis), axis, rms_residual }
}

/// Centred moving average whose window shrinks symmetrically at the ends, so a
/// linear trend is preserved right up to the boundary.
pub(crate) fn smooth_symmetric(values: &[f64], window: usize) -> Vec<f64> {
    let n = values.len();
    let half = window.max(1) / 2;
    let mut prefix = Vec::with_capacity(n + 1);
    prefix.push(0.0);
    for v in values {
        prefix.push(prefix.last().unwrap() + v);
    }
    (0..n)
        .map(|i| {
            let h = half.min(i).min(n - 1 - i);
            (prefix[i + h + 1] - prefix[i - h]) / (2 * h + 1) as f64
        })
        .collect()
}

/// Parses whitespace- or comma-separated numeric rows with at least `cols`
/// columns. `#` comments and blank lines are skipped; a non-numeric first data
/// line is treated as a header. Extra columns (e.g. time) are ignored.
pub(crate) fn parse_rows(text: &str, cols: usize) -> Result<Vec<Vec<f64>>, String> {
    let mut rows = Vec::new();
    let mut seen_data = false;
    for (ln, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parsed: Result<Vec<f64>, _> = line
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .take(cols)
            .map(str::parse::<f64>)
            .collect();
        match parsed {
            Ok(v) if v.len() == cols => rows.push(v),
            _ if !seen_data => {}
            _ => return Err(format!("line {}: expected {cols} numeric columns", ln + 1)),
        }
        seen_data = true;
    }
    Ok(rows)
}

/// Parses `x,y` rows (see [`parse_rows`] rules) into 2D points.
pub fn parse_xy(text: &str) -> Result<Vec<(f64, f64)>, String> {
    Ok(parse_rows(text, 2)?.into_iter().map(|r| (r[0], r[1])).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Elliptical, drifting orbit with a known period.
    fn ellipse(n: usize, period: f64) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let t = TAU * i as f64 / period;
                (0.01 * i as f64 + 3.0 * t.cos(), 1.5 * t.sin())
            })
            .collect()
    }

    #[test]
    fn finds_period_of_elliptical_orbit() {
        let p = dominant_period(&ellipse(4000, 900.0), 100.0, 2000.0).unwrap();
        assert!((p - 900.0).abs() < 30.0, "period {p}");
    }

    #[test]
    fn auto_axis_removes_the_orbit() {
        let pts = ellipse(4000, 900.0);
        let (axis, _) = auto_axis(&pts, 100.0, 2000.0).unwrap();
        let mid = 2000;
        let r = ((pts[mid].0 - axis[mid].0).powi(2) + (pts[mid].1 - axis[mid].1).powi(2)).sqrt();
        assert!(r > 1.0 && (axis[mid].0 - 0.01 * mid as f64).abs() < 0.5);
    }

    #[test]
    fn smoothed_turns_ignore_jitter() {
        let phase: Vec<f64> = (0..400)
            .map(|i| TAU * i as f64 / 100.0 + if i % 2 == 0 { 0.2 } else { -0.2 })
            .collect();
        let raw = total_turns(&phase);
        let smooth = total_turns_smoothed(&phase, 11);
        assert!(raw > 5.0 && (smooth - 3.99).abs() < 0.2, "{raw} {smooth}");
    }

    #[test]
    fn total_turns_ignores_direction() {
        let phase = [0.0, TAU, 0.0];
        assert!((total_turns(&phase) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn period_rejects_bad_input() {
        assert!(dominant_period(&[(0.0, 0.0); 4], 2.0, 10.0).is_none());
        assert!(dominant_period(&ellipse(100, 20.0), f64::NAN, 50.0).is_none());
        assert_eq!(window_for_period(f64::NAN), 3);
        assert_eq!(window_for_period(900.0) % 2, 1);
    }

    /// Spiral of radius `r` around a centre snaking along a sine path.
    fn snake(n: usize, r: f64) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let t = i as f64;
                let c = (0.08 * t, 3.0 * (0.01 * t).sin());
                (c.0 + r * (0.2 * t).cos(), c.1 + r * (0.2 * t).sin())
            })
            .collect()
    }

    #[test]
    fn moving_axis_straightens_a_snaking_spiral() {
        let pts = snake(600, 1.0);
        let fixed_centre = {
            let n = pts.len() as f64;
            let (sx, sy) = pts.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0, a.1 + p.1));
            vec![(sx / n, sy / n); pts.len()]
        };
        let fixed = radius_spread(&to_polar(&pts, &fixed_centre)[40..560]);
        let moving = radius_spread(&to_polar(&pts, &moving_axis(&pts, 31))[40..560]);
        assert!(fixed > 0.5, "fixed origin should look irregular, got {fixed}");
        assert!(moving < 0.05, "moving axis should look steady, got {moving}");
    }

    #[test]
    fn unwrapped_phase_advances_steadily() {
        let pts = snake(600, 1.0);
        let polar = to_polar(&pts, &moving_axis(&pts, 31));
        let u = unwrap_phase(&polar);
        let rates: Vec<f64> = u[40..560].windows(2).map(|w| w[1] - w[0]).collect();
        assert!(rates.iter().all(|r| (r - 0.2).abs() < 0.03), "rate drifted");
    }

    #[test]
    fn quantised_phase_cycles_through_ring() {
        let pts = snake(600, 1.0);
        let polar = to_polar(&pts, &moving_axis(&pts, 31));
        let q = quantise_phase(&polar[40..560], 8);
        assert!(q.iter().all(|&r| r < 8));
        assert_eq!((0..8).filter(|r| q.contains(r)).count(), 8);
        assert!(quantise_phase(&polar, 0).is_empty());
    }

    /// Spiral whose angular rate ramps from 0.1 to 0.4 rad/sample.
    fn chirp(n: usize) -> Vec<(f64, f64)> {
        let mut phase = 0.0;
        (0..n)
            .map(|i| {
                let t = i as f64;
                phase += 0.1 + 0.3 * t / n as f64;
                (0.08 * t + phase.cos(), 3.0 * (0.01 * t).sin() + phase.sin())
            })
            .collect()
    }

    #[test]
    fn adaptive_window_beats_fixed_on_a_chirp() {
        let pts = chirp(900);
        let spread = |axis: &[(f64, f64)]| radius_spread(&to_polar(&pts, axis)[100..800]);
        let fixed = spread(&moving_axis(&pts, 31));
        let adaptive = spread(&adaptive_axis(&pts, 31, 3));
        assert!(adaptive < fixed * 0.5, "fixed {fixed}, adaptive {adaptive}");
    }

    #[test]
    fn adaptive_axis_edge_cases() {
        assert!(adaptive_axis(&[], 31, 2).is_empty());
        assert_eq!(adaptive_axis(&[(1.0, 2.0)], 31, 2), vec![(1.0, 2.0)]);
        // A stationary point never rotates; this must not panic or produce NaN.
        let still = vec![(1.0, 1.0); 50];
        assert!(adaptive_axis(&still, 11, 2).iter().all(|p| p.0.is_finite() && p.1.is_finite()));
    }

    #[test]
    fn winding_and_ring_steps() {
        let pts = snake(600, 1.0);
        let polar = to_polar(&pts, &moving_axis(&pts, 31));
        let turns = winding_number(&unwrap_phase(&polar)[40..560]);
        assert!((turns - 0.2 * 519.0 / TAU).abs() < 0.2, "turns {turns}");

        let q = quantise_phase(&polar[40..560], 8);
        assert_eq!(backward_steps(&q, 8), 0);
        assert_eq!(quantised_steps(&[7, 0, 1], 8), vec![1, 1]);
        assert_eq!(quantised_steps(&[1, 0], 8), vec![-1]);
        assert!(quantised_steps(&[1, 0], 1).is_empty());
        assert_eq!(winding_number(&[]), 0.0);
    }

    #[test]
    fn gaps_and_outliers_are_repaired() {
        let mut pts = snake(600, 1.0);
        let clean = pts.clone();
        for p in &mut pts[200..215] {
            *p = (f64::NAN, f64::NAN);
        }
        pts[300] = (500.0, -500.0);
        let filled = fill_gaps(&pts);
        assert!(filled.iter().all(|p| p.0.is_finite() && p.1.is_finite()));

        let axis = moving_axis(&clean, 31);
        let cleaned = reject_outliers(&filled, &axis, 4.0);
        assert!(cleaned[300].0.abs() < 100.0, "outlier survived");
        assert_eq!(fill_gaps(&[(f64::NAN, 0.0)]).len(), 1);
        assert_eq!(fill_gaps(&[(f64::NAN, 0.0), (2.0, 4.0)])[0], (2.0, 4.0));
        assert!(reject_outliers(&[], &[], 4.0).is_empty());
    }

    #[test]
    fn refine_axis_reports_honest_residual() {
        let pts = snake(600, 1.0);
        let fit = refine_axis(&pts, 31, 3);
        assert!(fit.rms_residual < 0.05, "rms {}", fit.rms_residual);
        assert!(radius_spread(&fit.polar[40..560]) < 0.05);
        assert_eq!(refine_axis(&[], 31, 2).rms_residual, 0.0);
    }

    #[test]
    fn parse_xy_handles_headers_comments_and_errors() {
        let text = "# track\nx,y,t\n1.0, 2.0, 0\n3 4 1\n\n";
        assert_eq!(parse_xy(text).unwrap(), vec![(1.0, 2.0), (3.0, 4.0)]);
        assert!(parse_xy("1,2\nfoo,3\n").is_err());
        assert!(parse_xy("").unwrap().is_empty());
    }

    #[test]
    fn edge_cases() {
        assert!(moving_axis(&[], 5).is_empty());
        assert_eq!(radius_spread(&[]), f64::INFINITY);
        let one = moving_axis(&[(1.0, 2.0)], 9);
        assert_eq!(one, vec![(1.0, 2.0)]);
    }
}
