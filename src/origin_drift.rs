// src/origin_drifter.rs
// This module provides the OriginDrifter struct and related functionality.
// It allows for simulating a modular ring with a drifting origin and computing residue distributions.
// The OriginDrifter can step through the ring, recording residue distributions and their entropies
// at each step.
// The compute_entropy function calculates the entropy of a given residue distribution.
// The ResidueSnapshot struct captures the state of the origin and residue distribution at each step.
// The module is designed to be used in a modular arithmetic context, where the origin can drift
// and residues are computed relative to the current origin.
// The module is intended for use in applications that require tracking residue distributions
// in a modular arithmetic setting, such as cryptography, number theory, or simulations involving
// modular systems.
// This file is part of the MOMA project.

pub use crate::origin_drift::drifter::{
    OriginDrifter, ResidueSnapshot
};

pub mod drifter {

use std::collections::HashMap;

/// Represents a modular ring with a drifting origin.
pub struct OriginDrifter {
    modulus: usize,
    origin: usize,
    velocity: isize,
    history: Vec<ResidueSnapshot>,
}

/// Captures the residue distribution at a given origin.
#[derive(Debug, Clone)]
pub struct ResidueSnapshot {
    pub origin: usize,
    pub residues: Vec<usize>,
    pub entropy: f64,
}

impl OriginDrifter {
    /// Create a new OriginDrifter with given modulus, initial origin, and velocity.
    pub fn new(modulus: usize, origin: usize, velocity: isize) -> Self {
        OriginDrifter {
            modulus,
            origin: origin % modulus,
            velocity,
            history: Vec::new(),
        }
    }

    /// Advance the origin by one step and record the residue snapshot.
    pub fn step(&mut self) {
        self.origin = ((self.origin as isize + self.velocity).rem_euclid(self.modulus as isize)) as usize;
        let residues = self.compute_residues();
        let entropy = compute_entropy(&residues, self.modulus);
        self.history.push(ResidueSnapshot {
            origin: self.origin,
            residues,
            entropy,
        });
    }

    /// Run for `n` steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// Compute residues relative to current origin.
    fn compute_residues(&self) -> Vec<usize> {
        (0..self.modulus)
            .map(|x| (x + self.modulus - self.origin) % self.modulus)
            .collect()
    }

    /// Access the history of snapshots.
    pub fn history(&self) -> &[ResidueSnapshot] {
        &self.history
    }
}

/// Compute entropy of a residue distribution.
fn compute_entropy(residues: &[usize], _modulus: usize) -> f64 {
    let mut counts = HashMap::new();
    for &r in residues {
        *counts.entry(r).or_insert(0) += 1;
    }

    let total = residues.len() as f64;
    counts
        .values()
        .map(|&count| {
            let p = count as f64 / total;
            -p * p.log2()
        })
        .sum()
}
}

