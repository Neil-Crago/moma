
/// Tools for number-theoretic analysis related to MOMA.
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
    /// Creates a new `CompositeDampener`.
    pub fn new(lower: u64, upper: u64, small_primes: Vec<u64>) -> Self {
        Self {
            lower,
            upper,
            small_primes,
        }
    }

    /// Calculates the dampening score for the given range.
    ///
    /// The score is the ratio of composites hit by `small_primes` to the total
    /// number of composites in the range. It ranges from 0.0 to 1.0.
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

