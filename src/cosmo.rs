use crate::barycentric::OriginShift;
use crate::resonance::ResonanceDetector;
use crate::utils::write_csv;


#[derive(Debug, Clone)]
pub struct CosmoEntropyPulse {
    pub amplitude: f64,
    pub frequency: f64,
    pub phase: f64,
}

impl CosmoEntropyPulse {
    pub fn value_at(&self, time: f64) -> f64 {
        self.amplitude * (self.frequency * time + self.phase).sin()
    }
}


#[derive(Debug, Clone)]
pub struct Planet {
    pub name: &'static str,
    pub mass: f64,
    pub orbital_radius: f64,
    pub phase: f64,
    pub entropy: Option<CosmoEntropyPulse>, // Optional entropy modulation
}


pub struct BarycenterSimulator {
    pub planets: Vec<Planet>,
    pub time: f64,
    pub history: Vec<f64>, // magnitude of origin shifts
}

impl BarycenterSimulator {
    pub fn new(planets: Vec<Planet>) -> Self {
        Self { planets, time: 0.0, history: Vec::new() }
    }

    pub fn step(&mut self, dt: f64) -> OriginShift {
        self.time += dt;
        let mut x = 0.0;
        let mut y = 0.0;
        let total_mass: f64 = self.planets.iter().map(|p| p.mass).sum();

        for p in &self.planets {
    let angle = self.time / p.orbital_radius + p.phase;
    let entropy_mod = p.entropy.as_ref().map_or(1.0, |e| 1.0 + e.value_at(self.time));
    let effective_mass = p.mass * entropy_mod;

    x += effective_mass * p.orbital_radius * angle.cos();
    y += effective_mass * p.orbital_radius * angle.sin();
}

        let origin = OriginShift {
            dx: x / total_mass,
            dy: y / total_mass,
        };

        self.history.push(origin.magnitude());
        origin
    }

    pub fn export_history(&self, path: &str) {
        let _ = write_csv(path, &self.history);
    }

    pub fn detect_resonance<D: ResonanceDetector>(&self, detector: D) -> Vec<f64> {
        detector.detect(&self.history)
    }
}

