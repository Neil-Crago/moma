//! # Resonance Module
//!
//! Provides tools to find "resonance" events, where a MOMA signature aligns
//! with another mathematical property of its prime context.

use crate::core::{MomaRing, OriginStrategy};
use crate::primes;
use std::marker::PhantomData;

/// A function pointer type that defines a property of a prime number.
/// This is used as the target for resonance checks.
pub type PrimePropertyFn = fn(u64) -> u64;

/// An analyzer that finds primes where the MOMA signature "resonates" with
/// another property of the prime.
///
/// Resonance occurs when `ring.signature(p) % property_fn(p) == 0`.
pub struct ResonanceFinder<S: OriginStrategy> {
    ring: MomaRing<S>,
    property_fn: PrimePropertyFn,
    _strategy: PhantomData<S>,
}

impl<S: OriginStrategy> ResonanceFinder<S> {
    /// Creates a new `ResonanceFinder`.
    ///
    /// # Arguments
    /// * `modulus` - The modulus for the internal `MomaRing`.
    /// * `strategy` - The `OriginStrategy` to use for generating signatures.
    /// * `property_fn` - A function that defines the property to check for resonance against.
    ///   For example, `primes::prime_factor_mass` could be used.
    pub fn new(modulus: u64, strategy: S, property_fn: PrimePropertyFn) -> Self {
        Self {
            ring: MomaRing::new(modulus, strategy),
            property_fn,
            _strategy: PhantomData,
        }
    }

    /// Checks a single prime for a resonance event.
    ///
    /// # Returns
    /// `Some(signature)` if resonance occurs, otherwise `None`.
    pub fn check_prime(&self, p: u64) -> Option<u64> {
        let signature = self.ring.signature(p);
        let property_value = (self.property_fn)(p);

        // Avoid division by zero and check for resonance.
        if property_value > 0 && signature % property_value == 0 {
            Some(signature)
        } else {
            None
        }
    }

    /// Finds all primes within a given range that exhibit resonance.
    ///
    /// # Arguments
    /// * `start_range` - The beginning of the range to search for primes.
    /// * `end_range` - The end of the range to search.
    ///
    /// # Returns
    /// A `Vec` of tuples `(prime, signature)` for each resonance event found.
    pub fn find_in_range(&self, start_range: u64, end_range: u64) -> Vec<(u64, u64)> {
        let mut p = primes::next_prime(start_range.saturating_sub(1));
        let mut resonances = Vec::new();

        while p < end_range {
            if let Some(signature) = self.check_prime(p) {
                resonances.push((p, signature));
            }
            p = primes::next_prime(p);
        }
        resonances
    }
}