//! Provides tools for modeling the "influence" of composite numbers.

use crate::primes::primes;
use std::collections::HashMap;

/// A tool to model the gravitational-like "influence" of composite numbers.
///
/// This struct calculates a "mass" for each composite number in a given range
/// based on its number of prime factors. It can then be used to calculate the
/// total influence exerted by these masses at a specific point in the number line.
#[derive(Debug)]
pub struct CompositeInfluence {
    /// A map from a composite number to its calculated mass.
    pub composite_masses: HashMap<u64, f64>,
}

impl CompositeInfluence {
    /// Creates a new `CompositeInfluence` field for a given number range.
    /// The "mass" of each composite is calculated using `primes::prime_factor_mass`.
    pub fn new(range_start: u64, range_end: u64) -> Self {
        let composite_masses = (range_start..=range_end)
            .filter(|&n| !primes::is_prime(n))
            .map(|n| (n, primes::prime_factor_mass(n) as f64))
            .collect();
        Self { composite_masses }
    }

    /// Calculates the total influence exerted by all composite masses at a given point.
    ///
    /// The influence of each composite number is weighted by the inverse square of its
    /// distance to the target point, simulating a gravitational field.
    ///
    /// # Parameters
    /// - `point`: The number line coordinate to measure the influence at.
    pub fn influence_at_point(&self, point: f64) -> f64 {
        self.composite_masses
            .iter()
            .map(|(&composite, &mass)| {
                // Use inverse square law for influence falloff
                let dist_sq = (point - composite as f64).powi(2);
                mass / dist_sq.max(1.0) // Avoid division by zero
            })
            .sum()
    }
}