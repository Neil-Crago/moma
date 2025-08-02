//! Provides tools for calculating entropy.

use std::collections::HashMap;
use std::hash::Hash;

/// A generic struct to calculate the Shannon entropy of a sequence of items.
///
/// Entropy is a measure of the uncertainty or randomness in a set of data.
/// A higher entropy score implies a more uniform and less predictable distribution.
#[derive(Debug, Default)]
pub struct Entropy<T> {
    frequencies: HashMap<T, u64>,
    count: u64,
}

impl<T: Eq + Hash> Entropy<T> {
    /// Creates a new, empty `Entropy` calculator.
    pub fn new() -> Self {
        Self {
            frequencies: HashMap::new(),
            count: 0,
        }
    }

    /// Adds an item to the sequence being analyzed.
    pub fn add(&mut self, item: T) {
        *self.frequencies.entry(item).or_insert(0) += 1;
        self.count += 1;
    }

    /// Adds multiple items from an iterator to the sequence.
    pub fn add_all<I>(&mut self, iter: I)
    where
        I: IntoIterator<Item = T>,
    {
        for item in iter {
            self.add(item);
        }
    }

    /// Calculates the total Shannon entropy of the distribution of items seen so far.
    ///
    /// The formula used is H(X) = -Σ [P(x) * log₂(P(x))] for all x in X.
    ///
    /// # Returns
    /// The total entropy as an `f64`. Returns `0.0` if no items have been added.
    pub fn total_entropy(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }

        self.frequencies
            .values()
            .map(|&count| {
                let probability = count as f64 / self.count as f64;
                if probability > 0.0 {
                    -probability * probability.log2()
                } else {
                    0.0
                }
            })
            .sum()
    }
}