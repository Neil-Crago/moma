//! # Prime Gap Field Analyzer
//!
//! This program demonstrates advanced analysis of prime numbers using the `moma_crate`.
//! It introduces the concept of a `PrimeGapField`, a structure for studying the
//! statistical properties of the gaps between consecutive primes.
//!
//! The main analyses performed are:
//! 1.  **Barycentric Offset**: Measuring how much each prime gap deviates from its local average.
//! 2.  **Modular Class Entropy**: Calculating the Shannon entropy of prime gaps modulo `n`.
//! 3.  **Composite Influence**: Modeling how nearby composite numbers "influence" prime gaps.
//! 4.  **Goldbach Projection**: Using the collected prime data to find Goldbach pairs.

use moma::primes::primes;
use std::collections::{HashMap, HashSet};

// --- Main Analysis Structures ---

/// Represents a single gap between two consecutive prime numbers.
#[derive(Debug, Clone)]
pub struct PrimeGap {
    /// The prime number at the start of the gap.
    pub start_prime: u64,
    /// The prime number at the end of the gap.
    pub end_prime: u64,
    /// The size of the gap (`end_prime - start_prime`).
    pub size: u64,
    /// The modular class of the gap size, i.e., `size % modulus`.
    pub mod_class: u64,
    /// The "barycentric offset," representing how much the gap's size
    /// deviates from the local average gap size. Can be modified by other analyses.
    pub bary_offset: f64,
}

/// A data structure for analyzing a sequence of prime gaps.
///
/// It holds a collection of `PrimeGap` instances and provides methods for
/// statistical analysis like filtering, entropy scoring, and more.
pub struct PrimeGapField {
    /// The collection of prime gaps in the field.
    pub gaps: Vec<PrimeGap>,
    /// The modulus used for calculating `mod_class` for each gap.
    pub modulus: u64,
    /// A map holding the calculated Shannon entropy for each modular class.
    pub entropy_scores: HashMap<u64, f64>,
}

impl PrimeGapField {
    /// Creates a new `PrimeGapField` from a slice of primes and a modulus.
    ///
    /// # Panics
    /// Panics if the provided `primes` slice has fewer than two elements.
    pub fn new(primes: &[u64], modulus: u64) -> Self {
        assert!(primes.len() >= 2, "Need at least two primes to form a gap.");

        let gaps = primes
            .windows(2)
            .enumerate()
            .map(|(i, window)| {
                let p1 = window[0];
                let p2 = window[1];
                let gap_size = p2 - p1;

                // Calculate the average of a small window of gaps around the current one.
                // The window includes the two preceding, the current, and the next gap.
                let local_avg = Self::calculate_local_avg(primes, i + 1);
                let bary_offset = gap_size as f64 - local_avg;

                PrimeGap {
                    start_prime: p1,
                    end_prime: p2,
                    size: gap_size,
                    mod_class: gap_size % modulus,
                    bary_offset,
                }
            })
            .collect();

        Self {
            gaps,
            modulus,
            entropy_scores: HashMap::new(),
        }
    }

    /// Filters gaps where the absolute barycentric offset exceeds a threshold.
    /// This is useful for finding unusually large or small gaps.
    pub fn filter_by_bary_offset(&self, threshold: f64) -> Vec<&PrimeGap> {
        self.gaps
            .iter()
            .filter(|gap| gap.bary_offset.abs() > threshold)
            .collect()
    }

    /// Filters gaps belonging to a specific modular class.
    pub fn filter_by_mod_class(&self, target_class: u64) -> Vec<&PrimeGap> {
        self.gaps
            .iter()
            .filter(|gap| gap.mod_class == target_class)
            .collect()
    }

    /// Calculates the Shannon entropy for the distribution of gap modular classes.
    /// The results are stored in the `entropy_scores` field.
    pub fn calculate_entropy(&mut self) {
        if self.gaps.is_empty() {
            return;
        }
        let mut frequencies = HashMap::new();
        for gap in &self.gaps {
            *frequencies.entry(gap.mod_class).or_insert(0) += 1;
        }

        let total_gaps = self.gaps.len() as f64;
        self.entropy_scores = frequencies
            .into_iter()
            .map(|(class, count)| {
                let p = count as f64 / total_gaps;
                let entropy = if p > 0.0 { -p * p.log2() } else { 0.0 };
                (class, entropy)
            })
            .collect();
    }

    /// Modifies the `bary_offset` of each gap based on the "influence" of nearby composites.
    /// This simulates a "gravitational" pull from numbers with high prime factor mass.
    pub fn apply_composite_influence(&mut self, influence_field: &CompositeInfluence) {
        for gap in &mut self.gaps {
            // Calculate the total influence on the midpoint of the gap.
            let gap_midpoint = gap.start_prime as f64 + (gap.size as f64 / 2.0);
            let total_influence: f64 = influence_field
                .composite_masses
                .iter()
                .map(|(&composite, &mass)| {
                    // Use inverse square law for influence falloff
                    let dist_sq = (gap_midpoint - composite as f64).powi(2);
                    mass / dist_sq.max(1.0) // Avoid division by zero
                })
                .sum();

            // Modulate the existing offset by this calculated influence.
            gap.bary_offset += total_influence;
        }
    }

    /// Suggests Goldbach pairs for an even number using the primes available in the field.
    /// A Goldbach pair `(p1, p2)` consists of two primes such that `p1 + p2 = even_n`.
    pub fn project_goldbach(&self, even_n: u64) -> Vec<(u64, u64)> {
        if even_n % 2 != 0 {
            return Vec::new(); // Goldbach conjecture is for even numbers
        }
        // For efficient lookups, put all primes from the field into a HashSet.
        let prime_set: HashSet<u64> = self
            .gaps
            .iter()
            .flat_map(|gap| [gap.start_prime, gap.end_prime])
            .collect();

        prime_set
            .iter()
            .filter_map(|&p1| {
                if p1 > even_n / 2 { return None; } // Avoid duplicate pairs like (5, 7) and (7, 5)
                let p2 = even_n - p1;
                if prime_set.contains(&p2) {
                    Some((p1, p2))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Private helper to calculate the local average gap size around a given index.
    fn calculate_local_avg(primes: &[u64], index: usize) -> f64 {
        // Defines a window of 2 gaps before and 1 after the current one.
        let start = index.saturating_sub(2);
        let end = (index + 1).min(primes.len() - 2);

        if start >= end { return 0.0; }

        let window = &primes[start..=end + 1];
        let total_gap_size: u64 = window.windows(2).map(|w| w[1] - w[0]).sum();
        let count = window.len() - 1;

        total_gap_size as f64 / count.max(1) as f64
    }
}

/// Models the "influence" of composite numbers within a given range.
/// Each composite is assigned a "mass" based on its number of prime factors.
#[derive(Debug)]
pub struct CompositeInfluence {
    /// A map from a composite number to its calculated mass.
    pub composite_masses: HashMap<u64, f64>,
}

impl CompositeInfluence {
    /// Creates a new `CompositeInfluence` field for a given number range.
    pub fn new(range_start: u64, range_end: u64) -> Self {
        let composite_masses = (range_start..=range_end)
            .filter(|&n| !primes::is_prime(n))
            .map(|n| {
                // The "mass" is the count of prime factors (from moma_crate).
                let mass = primes::prime_factor_mass(n) as f64;
                (n, mass)
            })
            .collect();
        Self { composite_masses }
    }
}


// --- Main Application Logic ---

fn main() {
    println!("\n--- Prime Gap Field Analysis --- 🌌");

    // 1. Generate a list of primes to analyze.
    let primes: Vec<u64> = (1..=100).filter(|&n| primes::is_prime(n)).collect();

    // 2. Create a PrimeGapField with a modulus of 6.
    //    The choice of 6 is interesting because all primes > 3 are of the form 6k ± 1.
    let mut field = PrimeGapField::new(&primes, 6);

    // 3. Calculate and display the entropy of the gap classes.
    field.calculate_entropy();
    println!("\n📊 Entropy Scores for Prime Gaps mod 6:");
    let mut sorted_entropy: Vec<_> = field.entropy_scores.iter().collect();
    sorted_entropy.sort_by_key(|(class, _)| *class);
    for (class, score) in sorted_entropy {
        println!("   Class {:>2}: {:.4}", class, score);
    }

    // 4. Find gaps with a significant barycentric offset.
    let outliers = field.filter_by_bary_offset(3.0);
    println!("\n🌠 Outlier Gaps (Barycentric Offset > 3.0):");
    for gap in outliers {
        println!(
            "   Gap [{:>3}, {:>3}] | Size: {:>2} | Offset: {:>+.2}",
            gap.start_prime, gap.end_prime, gap.size, gap.bary_offset
        );
    }
    
    // 5. Project Goldbach pairs for an even number.
    let even_n = 96;
    let goldbach_pairs = field.project_goldbach(even_n);
    println!("\n✨ Goldbach Projections for {}:", even_n);
    // Format the output nicely
    let pair_strings: Vec<String> = goldbach_pairs.iter().map(|(p1, p2)| format!("{}+{}", p1, p2)).collect();
    println!("   {} = {}", even_n, pair_strings.join(" = "));

}


// --- Unit Tests ---

#[cfg(test)]
mod tests {
    use super::*;

    fn get_test_primes() -> Vec<u64> {
        vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    }

    #[test]
    fn test_field_creation() {
        let primes = get_test_primes();
        let field = PrimeGapField::new(&primes, 6);
        // We have 15 primes, so we expect 14 gaps.
        assert_eq!(field.gaps.len(), 14);
        assert_eq!(field.modulus, 6);
        // The first gap is 3-2=1. Its mod 6 class should be 1.
        assert_eq!(field.gaps[0].size, 1);
        assert_eq!(field.gaps[0].mod_class, 1);
        // The second gap is 5-3=2. Its mod 6 class should be 2.
        assert_eq!(field.gaps[1].size, 2);
        assert_eq!(field.gaps[1].mod_class, 2);
    }

    #[test]
    fn test_mod_class_filter() {
        let primes = get_test_primes();
        let field = PrimeGapField::new(&primes, 6);
        // Gaps of size 2, 4, 6. mod 6 classes are 2, 4, 0.
        // Gaps of size 2: (3,5), (11,13), (17,19), (29,31), (41,43) -> 5 gaps
        // Gaps of size 4: (7,11), (13,17), (19,23), (43,47) -> 4 gaps
        let class_2_gaps = field.filter_by_mod_class(2);
        assert_eq!(class_2_gaps.len(), 6);
        let class_4_gaps = field.filter_by_mod_class(4);
        assert_eq!(class_4_gaps.len(), 5);
    }
    
    #[test]
    fn test_goldbach_projection() {
        let primes = get_test_primes();
        let field = PrimeGapField::new(&primes, 48); // Even number must be <= sum of largest two primes
        let pairs = field.project_goldbach(48);
        // Expected pairs for 48: (5, 43), (7, 41), (11, 37), (17, 31), (19, 29)
        let mut expected = vec![(5, 43), (7, 41), (11, 37), (17, 31), (19, 29)];
        // The result might be in a different order, so we sort both to compare.
        let mut sorted_pairs = pairs;
        sorted_pairs.sort();
        expected.sort();
        assert_eq!(sorted_pairs, expected);
    }
    
    #[test]
    fn test_entropy_calculation() {
        let primes = get_test_primes();
        let mut field = PrimeGapField::new(&primes, 6);
        field.calculate_entropy();
        
        assert!(field.entropy_scores.contains_key(&0)); // Gaps of size 6 (e.g., 23->29)
        assert!(field.entropy_scores.contains_key(&1)); // Gap of size 1 (2->3)
        assert!(field.entropy_scores.contains_key(&2));
        assert!(field.entropy_scores.contains_key(&4));
        
        // Total entropy should be the sum of individual scores
        let total_entropy: f64 = field.entropy_scores.values().sum();
        assert!(total_entropy > 0.0);
    }
}
