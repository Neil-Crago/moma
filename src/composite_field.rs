//! # Composite Field Module
//!
//! Provides a simple utility for identifying composite numbers within a given range.

use crate::primes;

/// A utility struct that finds all composite numbers within a specified range.
#[derive(Debug)]
pub struct CompositeField {
    /// The start of the number range to analyze.
    pub range_start: u64,
    /// The end of the number range to analyze.
    pub range_end: u64,
}

impl CompositeField {
    /// Creates a new `CompositeField` for a given range.
    pub fn new(range_start: u64, range_end: u64) -> Self {
        Self {
            range_start,
            range_end,
        }
    }

    /// Generates a vector containing all composite numbers in the specified range.
    ///
    /// # Returns
    /// A `Vec<u64>` of the composite numbers.
    pub fn composites(&self) -> Vec<u64> {
        (self.range_start..=self.range_end)
            .filter(|&n| n > 1 && !primes::is_prime(n))
            .collect()
    }
}