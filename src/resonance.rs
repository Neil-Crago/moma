use rustfft::{FftPlanner, num_complex::Complex};
use arima::acf::acf;

/// Trait for resonance detection
pub trait ResonanceDetector {
    fn detect(&self, data: &[f64]) -> Vec<f64>;
}

/// FFT-based detector
pub struct FFTResonance;

impl ResonanceDetector for FFTResonance {
    fn detect(&self, data: &[f64]) -> Vec<f64> {
        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(data.len());

        let mut buffer: Vec<Complex<f64>> = data.iter()
            .map(|&x| Complex { re: x, im: 0.0 })
            .collect();

        fft.process(&mut buffer);
        buffer.iter().map(|c| c.norm()).collect()
    }

}

/// Autocorrelation-based detector
pub struct AutoCorrelationResonance {
    pub max_lag: usize,
}

impl ResonanceDetector for AutoCorrelationResonance {
    fn detect(&self, data: &[f64]) -> Vec<f64> {
        // Pass Some(self.max_lag) as the second argument, and use unwrap() to get Vec<f64>
        acf(data, Some(self.max_lag), false).unwrap()
    }
}

pub fn score_resonance_strength(spectrum: &[f64]) -> f64 {
    let max = spectrum.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean = spectrum.iter().sum::<f64>() / spectrum.len() as f64;
    max / mean // signal-to-noise ratio
}


pub fn detect_modular_resonance(data: &[f64], max_lag: usize) -> Vec<f64> {
    acf(data, Some(max_lag), false).unwrap()
}

pub fn detect_fft_resonance(data: &[f64]) -> Vec<f64> {
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(data.len());

    let mut buffer: Vec<Complex<f64>> = data.iter()
        .map(|&x| Complex { re: x, im: 0.0 })
        .collect();

    fft.process(&mut buffer);

    buffer.iter().map(|c| c.norm()).collect()
}
