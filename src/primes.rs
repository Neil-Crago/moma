//! Utility functions for prime number operations.
///
/// NOTE: For a high-performance production crate, consider replacing these
/// with a dependency on a specialized library like `primal`
/// 
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