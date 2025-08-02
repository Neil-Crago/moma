use moma::cosmo::{Planet, BarycenterSimulator};
use moma::resonance::{FFTResonance, AutoCorrelationResonance};
use moma::score::{score_signal_to_noise, score_kurtosis};
use moma::utils;

pub fn run_simulation() {
    let planets = vec![
        Planet { name: "Jupiter", mass: 317.8, orbital_radius: 5.2, phase: 0.0, entropy: None },
        Planet { name: "Saturn",  mass: 95.2,  orbital_radius: 9.5, phase: 1.0, entropy: None },
        Planet { name: "Earth",   mass: 1.0,   orbital_radius: 1.0, phase: 2.0, entropy: None },
    ];

    let mut sim = BarycenterSimulator::new(planets);

    for _ in 0..100 {
        sim.step(0.1);
    }

    sim.export_history("origin_shift.csv");

    

let fft = FFTResonance;
let spectrum = sim.detect_resonance(fft);
let snr = score_signal_to_noise(&spectrum);
let sharpness = score_kurtosis(&spectrum);
println!("FFT Resonance SNR: {:.3}, Sharpness: {:.3}", snr, sharpness);

    let _ = crate::utils::write_csv("fft_spectrum.csv", &spectrum);

    let ac = AutoCorrelationResonance { max_lag: 30 };
    let autocorr = sim.detect_resonance(ac);
    let _ = moma::utils::write_csv("autocorr.csv", &autocorr);
}

fn main() {
    run_simulation();
}   