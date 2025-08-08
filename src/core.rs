//! Core MOMA structures and traits.

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

