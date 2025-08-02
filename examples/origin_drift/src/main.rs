//! # OriginDrift Analysis Example
//!
//! This program demonstrates how to use the generic `OriginDrift` struct
//! from the `moma` crate to compare the volatility of different `OriginStrategy`
//! implementations.
//!
//! It creates two separate drift analyzers—one for `PrimeGap` and one for
//! `CompositeMass`—feeds them the same sequence of primes, and then compares
//! their final "drift magnitude" scores.

use moma::origin_drift::OriginDrift;
use moma::primes::primes;
use moma::strategy::strategy;

fn main() {
    println!("\n--- OriginDrift Analysis: Comparing Strategy Volatility ---\n");

    // --- Setup ---
    let modulus = 100; // A common modulus for both analyzers.
    let num_primes_to_test = 20; // The number of primes to feed into each analyzer.

    // 1. Create two drift analyzers, one for each strategy we want to compare.
    let mut gap_drift = OriginDrift::new(modulus, strategy::PrimeGap);
    let mut mass_drift = OriginDrift::new(modulus, strategy::CompositeMass);

    println!(
        "Feeding the first {} primes into two drift analyzers (Modulus: {})...",
        num_primes_to_test, modulus
    );

    // --- Analysis Loop ---
    // 2. Iterate through a sequence of primes and feed each one into both analyzers.
    let mut p = 3; // Start with the first odd prime.
    for i in 0..num_primes_to_test {
        let sig1 = gap_drift.next(p);
        let sig2 = mass_drift.next(p);
        println!(
            "  {:>2}. p={:<3} | PrimeGap Sig: {:<3} | CompositeMass Sig: {:<3}",
            i + 1,
            p,
            sig1,
            sig2
        );
        p = primes::next_prime(p);
    }

    // --- Results ---
    // 3. Calculate and print the final drift magnitude for each strategy.
    let gap_drift_magnitude = gap_drift.drift_magnitude();
    let mass_drift_magnitude = mass_drift.drift_magnitude();

    println!("\n-------------------------------------------------");
    println!("Final Drift Magnitudes:");
    println!(
        "  - PrimeGap Strategy Volatility:      {:.2}",
        gap_drift_magnitude
    );
    println!(
        "  - CompositeMass Strategy Volatility: {:.2}",
        mass_drift_magnitude
    );
    println!("-------------------------------------------------");

    // 4. Compare the results to determine which strategy was more volatile.
    if mass_drift_magnitude > gap_drift_magnitude {
        println!("\n🏆 Conclusion: The CompositeMass strategy was more volatile.");
    } else {
        println!("\n🏆 Conclusion: The PrimeGap strategy was more volatile.");
    }
}
