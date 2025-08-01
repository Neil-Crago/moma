//! # MOMA: Moving Origin Modular Arithmetic
//!
//! A framework for exploring number theory concepts using a "moving origin"
//! in modular arithmetic.
//!
//! The core idea of MOMA is that for a given modulus `m`, the "zero point" or "origin"
//! for the calculation `x mod m` is not fixed. Instead, it shifts based on a
//! contextual value, typically a prime number `p`. This allows for the study of
//! properties of numbers (like primes) in a dynamic relational context.
//!
//! ## Core Concepts
//!
//! - **`MomaRing`**: The central object for calculations. It is defined by a `modulus` and an `OriginStrategy`.
//! - **`OriginStrategy`**: A trait that defines *how* the origin moves. The crate provides several
//!   common strategies, and users can easily implement their own.
//! - **Residue**: The result of a MOMA calculation, computed as `(value + origin) % modulus`.
//!
//! ## Example Usage
//!
//! ```
//! use moma::{MomaRing, strategy}; 
//!
//! // Create a MOMA ring with modulus 37.
//! // The origin will be the gap between a prime and the previous prime.
//! let ring = MomaRing::new(37, strategy::PrimeGap);
//!
//! // Let's analyze the prime p = 29.
//! let p = 29;
//!
//! // The value we want to find the residue of. In my original research,
//! // this was the sum of a prime and its predecessor.
//! let value_to_test = p + primes::prev_prime(p); // 29 + 23 = 52
//!
//! // Calculate the MOMA residue.
//! // The origin for p=29 is (29 - 23) = 6.
//! // The residue is (52 + 6) % 37 = 58 % 37 = 21.
//! let residue = ring.residue(value_to_test, p);
//!
//! println!("For p={}, the MOMA residue of {} is {}", p, value_to_test, residue);
//! assert_eq!(residue, 21);
//! ```

pub mod origin_drift;
pub mod massfield;

// Re-export key components for easier access.
pub use crate::analysis::{CompositeDampener, MassField};
pub use crate::core::{MomaRing, OriginStrategy};

/// Core MOMA structures and traits.
pub mod core {
    use crate::primes;

    /// Defines a strategy for calculating the moving origin for a given prime context.
    ///
    /// This trait is the cornerstone of MOMA's flexibility. By implementing this trait,
    /// users can create custom logic for how the modular origin should be determined.
    pub trait OriginStrategy {
        /// Calculates the origin based on a contextual prime `p`.
        ///
        /// # Parameters
        /// - `p`: The prime number providing the context for the origin calculation.
        fn calculate_origin(&self, p: u64) -> u64;
    }

    /// The central struct for performing Moving Origin Modular Arithmetic.
    ///
    /// A `MomaRing` is configured with a modulus and a chosen `OriginStrategy`.
    /// It then calculates residues by shifting the input value by the dynamically
    /// computed origin before applying the modulus.
    pub struct MomaRing<S: OriginStrategy> {
        pub modulus: u64,
        strategy: S,
    }

    impl<S: OriginStrategy> MomaRing<S> {
        /// Creates a new `MomaRing` with a given modulus and origin strategy.
        ///
        /// # Parameters
        /// - `modulus`: The modulus for the arithmetic operations.
        /// - `strategy`: An instance of a struct that implements `OriginStrategy`.
        pub fn new(modulus: u64, strategy: S) -> Self {
            Self { modulus, strategy }
        }

        /// Calculates the MOMA residue for a value within a prime context.
        ///
        /// This is the primary operation of the `MomaRing`. It first calculates the
        /// origin using the ring's strategy and the provided `prime_context`,
        /// then computes `(value + origin) % modulus`.
        ///
        /// # Parameters
        /// - `value`: The input value to map to the ring.
        /// - `prime_context`: The prime number used to determine the origin shift.
        pub fn residue(&self, value: u64, prime_context: u64) -> u64 {
            // Ensure modulus is not zero to prevent division by zero panic.
            if self.modulus == 0 {
                return value;
            }
            let origin = self.strategy.calculate_origin(prime_context);
            (value.wrapping_add(origin)) % self.modulus
        }

        /// A convenience method for calculating the "signature" of a prime.
        ///
        /// The signature is defined as the residue of the sum of a prime and its
        /// immediate predecessor. This is a common use case in MOMA-based analysis.
        ///
        /// # Parameters
        /// - `p`: The prime for which to calculate the signature.
        pub fn signature(&self, p: u64) -> u64 {
            if p < 3 { return 0; } // prev_prime(2) is problematic, handle edge case.
            let input = p.wrapping_add(primes::prev_prime(p));
            self.residue(input, p)
        }
    }
}

/// Implementations of various origin strategies.
pub mod strategy {
    use crate::core::OriginStrategy;
    use crate::primes;

    /// An origin strategy where the origin is fixed to a constant value.
    #[derive(Debug, Clone, Copy)]
    pub struct Fixed(pub u64);
    impl OriginStrategy for Fixed {
        fn calculate_origin(&self, _p: u64) -> u64 {
            self.0
        }
    }

    /// An origin strategy where the origin is the gap between a prime and its predecessor.
    /// `origin(p) = p - p_prev`
    #[derive(Debug, Clone, Copy)]
    pub struct PrimeGap;
    impl OriginStrategy for PrimeGap {
        fn calculate_origin(&self, p: u64) -> u64 {
            if p < 3 { return 0; }
            p - primes::prev_prime(p)
        }
    }

    /// An origin strategy where the origin is the sum of prime factors of all
    /// composite numbers in the gap between a prime and its successor.
    /// `origin(p) = Σ mass(c)` for `c` in `(p, p_next)`.
    #[derive(Debug, Clone, Copy)]
    pub struct CompositeMass;
    impl OriginStrategy for CompositeMass {
        fn calculate_origin(&self, p: u64) -> u64 {
            let p_next = primes::next_prime(p);
            (p + 1..p_next)
                .filter(|&n| !primes::is_prime(n))
                .map(primes::prime_factor_mass)
                .sum()
        }
    }
}

/// Tools for number-theoretic analysis related to MOMA.
pub mod analysis {
    use crate::primes;

    /// A tool to analyze the "dampening" of composite numbers within a range.
    ///
    /// The dampening score measures how many composites in a range [lower, upper]
    /// are divisible by a given set of small primes. A higher score means the
    /// composites in the range are "less random" and more likely to be multiples
    /// of small primes.
    pub struct CompositeDampener {
        pub lower: u64,
        pub upper: u64,
        pub small_primes: Vec<u64>,
    }

    impl CompositeDampener {
        /// Calculates the dampening score for the given range.
        ///
        /// The score is the ratio of composites hit by `small_primes` to the total
        /// number of composites in the range.
        pub fn score(&self) -> f64 {
            let composites: Vec<u64> = (self.lower + 1..self.upper)
                .filter(|&n| !primes::is_prime(n))
                .collect();

            if composites.is_empty() {
                return 0.0;
            }

            let hits = composites
                .iter()
                .filter(|&c| self.small_primes.iter().any(|sp| c % sp == 0))
                .count();

            hits as f64 / composites.len() as f64
        }
    }

    /// A tool to analyze the "mass" of composite numbers between consecutive primes.
    pub struct MassField {
        pub range_start: u64,
        pub range_end: u64,
    }

    impl MassField {
        /// Generates a map of `(prime, composite_mass_in_next_gap)`.
        ///
        /// The "mass" is the sum of the count of prime factors for each composite
        /// number between a prime `p` and its successor `p_next`.
        pub fn generate_mass_map(&self) -> Vec<(u64, u64)> {
            let mut map = Vec::new();
            let mut p = primes::next_prime(self.range_start.saturating_sub(1));

            while p < self.range_end {
                let p_next = primes::next_prime(p);
                let mass = (p + 1..p_next)
                    .filter(|&n| !primes::is_prime(n))
                    .map(primes::prime_factor_mass)
                    .sum();
                map.push((p, mass));
                p = p_next;
            }
            map
        }
    }
}

/// Utility functions for prime number operations.
///
/// NOTE: For a high-performance production crate, consider replacing these
/// with a dependency on a specialized library like `primal`.
pub mod primes {

    /// A basic primality test.
    pub fn is_prime(n: u64) -> bool {
        if n < 2 { return false; }
        if n == 2 || n == 3 { return true; }
        if n % 2 == 0 || n % 3 == 0 { return false; }
        let mut i = 5;
        while i * i <= n {
            if n % i == 0 || n % (i + 2) == 0 {
                return false;
            }
            i += 6;
        }
        true
    }

    /// Finds the next prime number strictly greater than `n`.
    pub fn next_prime(n: u64) -> u64 {
        if n < 2 { return 2; }
        // Start with the next odd number.
        let mut x = if n % 2 == 0 { n + 1 } else { n + 2 };
        loop {
            if is_prime(x) {
                return x;
            }
            x += 2; // Only check odd numbers.
        }
    }

    /// Finds the greatest prime number strictly less than `n`.
    /// Returns 0 if no such prime exists (e.g., for n <= 2).
    pub fn prev_prime(n: u64) -> u64 {
        if n <= 2 { return 0; }
        let mut x = n - 1;
        while x >= 2 {
            if is_prime(x) {
                return x;
            }
            x -= 1;
        }
        0
    }

    /// Calculates the "mass" of a number, defined as the count of its prime factors
    /// with multiplicity. For example, `prime_factor_mass(12) = mass(2*2*3) = 3`.
    pub fn prime_factor_mass(n: u64) -> u64 {
        if n < 2 { return 0; }
        let mut count = 0;
        let mut temp_n = n;
        let mut factor = 2;
        while factor * factor <= temp_n {
            while temp_n % factor == 0 {
                count += 1;
                temp_n /= factor;
            }
            factor += 1;
        }
        if temp_n > 1 {
            count += 1;
        }
        count
    }
}
