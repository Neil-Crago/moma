use crate::primes;

/// A tool to analyze the "mass" of composite numbers between consecutive primes.
///
/// This struct defines a range and can generate a map of each prime in that range
/// to the "composite mass" found in the gap immediately following it. The "mass"
/// is the sum of the count of prime factors for each composite number.
pub struct MassField {
    /// The start of the number range to analyze.
    pub range_start: u64,
    /// The end of the number range to analyze.
    pub range_end: u64,
}

impl MassField {
    /// Creates a new `MassField` for a given range.
    pub fn new(range_start: u64, range_end: u64) -> Self {
        Self {
            range_start,
            range_end,
        }
    }

    /// Generates a map of `(prime, composite_mass_in_next_gap)`.
    ///
    /// - It iterates through each prime `p` within the specified range.
    /// - For each `p`, it calculates the total composite mass in the interval `(p, p_next)`.
    /// - Returns a `Vec` of tuples, where each tuple contains the starting prime
    ///   and the calculated mass of the subsequent gap.
    ///
    /// # Example
    /// ```
    /// use crate::massfield::MassField;
    /// let field = MassField::new(1, 20);
    /// let mass_map = field.generate_mass_map();
    /// // For p=13, p_next=17. Composites are 14, 15, 16.
    /// // Mass = mass(14) + mass(15) + mass(16) = 2 + 2 + 4 = 8.
    /// assert!(mass_map.contains(&(13, 8)));
    /// ```
    pub fn generate_mass_map(&self) -> Vec<(u64, u64)> {
        let mut map = Vec::new();
        // Start with the first prime at or after the range_start.
        let mut p = primes::next_prime(self.range_start.saturating_sub(1));

        while p < self.range_end {
            let p_next = primes::next_prime(p);
            // Stop if the next prime goes beyond the desired range.
            if p_next > self.range_end {
                break;
            }

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
