//! # MOMA Cosmology Simulation Example
//!
//! This program demonstrates a speculative and creative application of the MOMA
//! framework to cosmology. It models a simplified "universe" where prime
//! numbers represent fundamental events or epochs in time.
//!
//! This example is a conceptual counterpart to the `bioinformatics` example.
//! Both simulations explore the interplay between order and chaos using the
//! same core MOMA tools, but apply them to different domains.
//!
//! The simulation models:
//! 1.  **Structure Formation (Resonance)**: It uses a `ResonanceFinder` to identify
//!     primes where the MOMA signature aligns with another property. These are
//!     treated as significant, structure-forming events in the universe's timeline.
//! 2.  **Chaotic Eras (Entropy Pulses)**: It calculates the entropy of recent MOMA
//!     signatures. A sudden spike in entropy is modeled as a "pulse" that
//!     represents a chaotic era, potentially disrupting the formation of structure.

use moma::entropy::Entropy;
use moma::primes;
use moma::resonance::{PrimePropertyFn, ResonanceFinder};
use moma::strategy;
use std::collections::VecDeque;

/// A function that returns the number of prime factors of a number.
/// This property is used as the target for our resonance search.
fn factor_mass_property(p: u64) -> u64 {
    primes::prime_factor_mass(p)
}

fn main() {
    println!("--- MOMA Cosmology Simulation ---");

    // --- Simulation Parameters ---
    let modulus = 120;
    let end_of_time = 5000; // Simulate events for primes up to this value.
    let history_window = 10;   // Calculate entropy over the last 10 resonant events.
    let chaos_threshold = 2.5; // An entropy value that triggers a "chaotic era".

    // --- Setup ---
    let resonance_finder = ResonanceFinder::new(
        modulus,
        strategy::CompositeMass,
        factor_mass_property as PrimePropertyFn,
    );

    let mut signature_history: VecDeque<u64> = VecDeque::with_capacity(history_window);
    let mut total_structures_formed = 0;
    let mut chaotic_eras_triggered = 0;

    println!("Simulating universe evolution up to p = {}...", end_of_time);
    println!("A structure forms when a resonance event occurs.");
    println!("A chaotic era begins when entropy H(signatures) > {}\n", chaos_threshold);

    // --- Simulation Loop ---
    let mut p = 2;
    while p < end_of_time {
        // Check for a resonance event at the current prime 'p'.
        if let Some(signature) = resonance_finder.check_prime(p) {
            total_structures_formed += 1;
            println!(
                "✨ Structure formed at epoch p={:<4} (Resonance Signature: {})",
                p, signature
            );

            // Add the new signature to our history for entropy calculation.
            if signature_history.len() >= history_window {
                signature_history.pop_front();
            }
            signature_history.push_back(signature);

            // Recalculate entropy with the new history.
            let mut entropy_calculator = Entropy::new();
            entropy_calculator.add_all(signature_history.iter().copied());
            let current_entropy = entropy_calculator.total_entropy();

            // Check if the new state triggers an entropy pulse.
            if current_entropy > chaos_threshold {
                chaotic_eras_triggered += 1;
                println!(
                    "💥 A Chaotic Era begins! Entropy spiked to H={:.2}.",
                    current_entropy
                );
            }
        }
        p = primes::next_prime(p);
    }

    // --- Final Report ---
    println!("\n--------------------");
    println!("Simulation Complete.");
    println!("Total Structures Formed: {}", total_structures_formed);
    println!("Total Chaotic Eras:      {}", chaotic_eras_triggered);
    println!("--------------------");
}