//! # MassField Analysis Example
//!
//! This program demonstrates how to use the `MassField` struct from the `moma`
//! crate to analyze the distribution of "composite mass" between prime numbers.
//!
//! It calculates the mass for a given range and prints the results, highlighting
//! the prime gap with the highest concentration of composite matter.

use moma::massfield::MassField;

fn main() {
    println!("\n--- MassField Analysis ---");
    println!("Analyzing composite mass in prime gaps from 1 to 200.\n");

    // 1. Define the analysis range and create a new MassField.
    let start_range = 1;
    let end_range = 200;
    let field = MassField::new(start_range, end_range);

    // 2. Generate the map of (prime, mass_in_next_gap).
    let mass_map = field.generate_mass_map();

    // 3. Find the prime followed by the gap with the most composite mass.
    //    We use `max_by_key` to find the tuple with the largest second element (`mass`).
    let heaviest_gap = mass_map.iter().max_by_key(|&(_, mass)| mass);

    // 4. Print the results.
    println!("Calculated Composite Mass for each Prime Gap:");
    for (prime, mass) in &mass_map {
        println!("  - Gap after p={:<3}: Mass = {}", prime, mass);
    }

    if let Some((prime, mass)) = heaviest_gap {
        println!("\n-------------------------------------------------");
        println!(
            "🏆 The gap with the most composite mass follows p = {}.",
            prime
        );
        println!("   This gap has a total composite mass of {}.", mass);
        println!("-------------------------------------------------");
    } else {
        println!("\nNo prime gaps were found in the specified range.");
    }
}
