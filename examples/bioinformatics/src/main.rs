//! # MOMA Bioinformatics Example
//!
//! This program demonstrates a creative application of the MOMA framework to the
//! field of bioinformatics. It simulates how abstract mathematical patterns,
//! derived from prime numbers, could be mapped to concrete biological events
//! like genetic mutations within a chaotic environment.
//!
//! The core idea is to:
//! 1. Use a `MomaRing` to generate a "signature" for each prime in a sequence.
//! 2. Interpret this numeric signature as a location for a simulated point mutation.
//! 3. Track the **Shannon entropy** of the recent signatures.
//! 4. If entropy spikes, trigger a **"chaotic event"** or **"entropy pulse"**,
//!    simulating an environmental stressor that could increase mutation rates.
//! 5. Collect the data and visualize the signatures as a heatmap.

use moma::biosig::BioSigAnalyzer;
use moma::entropy::Entropy;
use moma::primes;
use moma::strategy;
use plotters::prelude::*;
use std::collections::VecDeque;

fn main() {
    println!("--- MOMA Bioinformatics Signature Analysis ---");

    // --- Configuration ---
    let modulus = 60;
    let dna_sequence = "AGCTGCGATCGTACGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";
    let end_time = 3000;
    let entropy_window_size = 15;
    let entropy_pulse_threshold = 3.5;

    let analyzer = BioSigAnalyzer::new(modulus, strategy::CompositeMass);

    println!("DNA Sequence Length: {}", dna_sequence.len());
    println!("MOMA Modulus: {}", modulus);
    println!("Analyzing primes up to {}...", end_time);
    println!("An entropy pulse occurs when H(signatures) > {}\n", entropy_pulse_threshold);


    // --- Simulation State ---
    let mut results = Vec::new();
    let mut signature_history: VecDeque<u64> = VecDeque::with_capacity(entropy_window_size);
    let mut mutation_events = 0;
    let mut chaotic_events = 0;


    // --- Analysis Loop ---
    let mut p = 2;
    while p < end_time {
        if let Some((signature, mutation)) = analyzer.analyze(p, dna_sequence) {
            mutation_events += 1;
            println!(
                "p={:<4} -> sig={:<2} -> Mutation at pos {:<2} | {:<10?} | {} -> {}",
                p,
                signature,
                signature % dna_sequence.len() as u64,
                mutation.mutation_type,
                mutation.original_codon,
                mutation.mutated_codon
            );
            results.push((p, signature));

            // --- Entropy Pulse Check ---
            // Add the new signature to our history for entropy calculation.
            if signature_history.len() >= entropy_window_size {
                signature_history.pop_front(); // Keep the window size fixed.
            }
            signature_history.push_back(signature);

            // Recalculate entropy with the new history.
            let mut entropy_calculator = Entropy::new();
            entropy_calculator.add_all(signature_history.iter().copied());
            let current_entropy = entropy_calculator.total_entropy();

            // Check if the new state triggers an entropy pulse.
            if current_entropy > entropy_pulse_threshold {
                chaotic_events += 1;
                println!(
                    "💥 CHAOTIC EVENT! High entropy (H={:.2}) suggests environmental stress.",
                    current_entropy
                );
            }
        }
        p = primes::next_prime(p);
    }

    // --- Final Report ---
    println!("\n--------------------");
    println!("Simulation Complete.");
    println!("Total Mutation Events: {}", mutation_events);
    println!("Total Chaotic Events:  {}", chaotic_events);
    println!("--------------------");


    // --- Visualization ---
    if !results.is_empty() {
        println!("\nGenerating heatmap of signatures...");
        if let Err(e) = plot_heatmap(&results, modulus) {
            eprintln!("Error plotting heatmap: {}", e);
        }
    } else {
        println!("\nNo results to plot.");
    }
}

/// Plots a heatmap of MOMA signatures and saves it as a PNG file.
fn plot_heatmap(data: &[(u64, u64)], modulus: u64) -> Result<(), Box<dyn std::error::Error>> {
    let out_file_name = "bio_heatmap.png";
    let root = BitMapBackend::new(out_file_name, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let min_prime = data.iter().map(|(p, _)| *p).min().unwrap_or(0);
    let max_prime = data.iter().map(|(p, _)| *p).max().unwrap_or(1);

    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .margin(20)
        .caption("MOMA Bio-Signature Heatmap", ("sans-serif", 40))
        .build_cartesian_2d(min_prime..max_prime, 0u64..modulus)?;

    chart
        .configure_mesh()
        .x_desc("Prime Number Context")
        .y_desc("MOMA Signature (Mutation Site)")
        .draw()?;

    chart.draw_series(
        data.iter().map(|&(prime, signature)| {
            let color = HSLColor(
                240.0 / 360.0 * (1.0 - signature as f64 / modulus as f64),
                0.7,
                0.5,
            );
            // Draw a single point for each signature
            Circle::new((prime, signature), 3, color.filled())
        }),
    )?;

    root.present()?;
    println!("Successfully saved heatmap to {}", out_file_name);
    Ok(())
}