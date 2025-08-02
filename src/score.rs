/// Scores resonance strength from a spectrum or autocorrelation series
pub fn score_signal_to_noise(data: &[f64]) -> f64 {
    if data.is_empty() { return 0.0; }
    let max = data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    max / mean
}

/// Scores peak sharpness using normalized kurtosis
pub fn score_kurtosis(data: &[f64]) -> f64 {
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    let variance = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / data.len() as f64;
    let fourth_moment = data.iter().map(|x| (x - mean).powi(4)).sum::<f64>() / data.len() as f64;
    if variance == 0.0 { return 0.0; }
    fourth_moment / variance.powi(2)
}