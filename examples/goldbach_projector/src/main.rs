//! # GoldbachProjector Example
//!
//! This program demonstrates how to use the `GoldbachProjector` struct from the
//! `moma` crate to efficiently find all prime number pairs that sum up to a
//! given set of even numbers.

use moma::goldbach::GoldbachProjector;

fn main() {
    println!("--- Goldbach Projector Demonstration ---");

    // 1. Create a projector with a database of primes up to 1000.
    //    This is efficient because the prime set is generated only once.
    let limit = 1000;
    let projector = GoldbachProjector::new(limit);
    println!("Prime database generated up to {}.\n", limit);

    // 2. Define a list of even numbers to test.
    let even_numbers = [48, 96, 100, 300, 774];

    // 3. For each even number, find and print its Goldbach pairs.
    for &n in &even_numbers {
        let pairs = projector.project(n);

        if pairs.is_empty() {
            println!("No Goldbach pairs found for {} (is it even and > 2?).", n);
        } else {
            // Format the output nicely for display.
            let pairs_str: Vec<String> =
                pairs.iter().map(|(p1, p2)| format!("{}+{}", p1, p2)).collect();
            println!(
                "✨ Found {} Goldbach pairs for {}:",
                pairs.len(),
                n
            );
            println!("   {} = {}\n", n, pairs_str.join(" = "));
        }
    }
}