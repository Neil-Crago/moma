use crate::codon::{CodonRing, EntropyPulse};
use crate::biosig::BioSignature;

pub struct MutationSnapshot {
    pub step: usize,
    pub spread: f64,
    pub entropy: f64,
    pub resonance: f64,
    pub drift: isize,
}

pub struct MutationTimeline {
    pub ring: CodonRing,
    pub pulses: Vec<EntropyPulse>, // Multiple pulses
    pub steps: usize,
    pub snapshots: Vec<MutationSnapshot>,
}

impl MutationTimeline {
    pub fn new(ring: CodonRing, pulses: Vec<EntropyPulse>, steps: usize) -> Self {
        Self {
            ring,
            pulses,
            steps,
            snapshots: Vec::new(),
        }
    }

    pub fn run(&mut self) {
        for step in 0..self.steps {
            let pulse = &self.pulses[step % self.pulses.len()];
            self.ring.apply_entropy_pulse(pulse);

            let snapshot = MutationSnapshot {
                step,
                spread: self.ring.modular_spread(),
                entropy: self.ring.entropy_footprint(),
                resonance: self.ring.resonance_score(),
                drift: self.ring.drift_vector(),
            };
            self.snapshots.push(snapshot);
        }
    }
}
