/// Implementations of various origin strategies.
pub mod strategy {
    use crate::core::core::OriginStrategy;
    use crate::primes::primes;

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
