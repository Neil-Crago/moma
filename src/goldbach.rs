//! Provides tools for exploring Goldbach's conjecture.

use crate::primes;
use std::collections::HashSet;

/// A tool to efficiently find Goldbach pairs for even numbers.
///
/// Goldbach's conjecture states that every even integer greater than 2 is the
/// sum of two prime numbers. This struct uses a pre-computed set of primes
/// for fast lookups to find these pairs.
#[derive(Debug)]
pub struct GoldbachProjector {
    prime_set: HashSet<u64>,
}

impl GoldbachProjector {
    /// Creates a new `GoldbachProjector` with a prime number database
    /// generated up to a specified limit.
    pub fn new(limit: u64) -> Self {
        let prime_set = (2..=limit).filter(|&n| primes::is_prime(n)).collect();
        Self { prime_set }
    }

    /// Finds all unique pairs of primes `(p1, p2)` that sum to a given even number `n`.
    ///
    /// The method ensures `p1 <= p2` to avoid duplicate pairs like `(3, 7)` and `(7, 3)`.
    ///
    /// # Parameters
    /// - `n`: The even number to find Goldbach pairs for.
    ///
    /// # Returns
    /// A `Vec` of tuples `(p1, p2)`. Returns an empty vector if `n` is odd or too small.
    pub fn project(&self, n: u64) -> Vec<(u64, u64)> {
        if n <= 2 || n % 2 != 0 {
            return Vec::new();
        }

        self.prime_set
            .iter()
            .filter(|&&p1| p1 <= n / 2) // Iterate up to n/2 to ensure unique pairs
            .filter_map(|&p1| {
                let p2 = n - p1;
                if self.prime_set.contains(&p2) {
                    Some((p1, p2))
                } else {
                    None
                }
            })
            .collect()
    }
}