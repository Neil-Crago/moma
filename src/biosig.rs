//! Provides a bridge between MOMA's numerical analysis and bioinformatics.
// This module defines a "biological signature" by mapping MOMA's numeric
// output to the effects of genetic mutations.

use crate::codon::CodonTable;
use crate::core::{MomaRing, OriginStrategy};
use crate::mutation::Mutation;
use std::marker::PhantomData;

/// An analyzer that generates MOMA signatures and interprets them as genetic mutations.
///
/// It uses a `MomaRing` to generate a numeric signature from a prime number, then
/// uses this signature to simulate a point mutation in a DNA sequence and analyzes
/// the resulting change in the amino acid sequence.
pub struct BioSigAnalyzer<S: OriginStrategy> {
    ring: MomaRing<S>,
    codon_table: CodonTable,
    _strategy: PhantomData<S>,
}

impl<S: OriginStrategy> BioSigAnalyzer<S> {
    /// Creates a new `BioSigAnalyzer`.
    ///
    /// # Arguments
    /// * `modulus` - The modulus for the internal `MomaRing`.
    /// * `strategy` - The `OriginStrategy` to use for generating signatures.
    pub fn new(modulus: u64, strategy: S) -> Self {
        Self {
            ring: MomaRing::new(modulus, strategy),
            codon_table: CodonTable::new(),
            _strategy: PhantomData,
        }
    }

    /// Generates a MOMA signature for a prime and analyzes its mutational effect.
    ///
    /// # Arguments
    /// * `p` - The prime number to use as the context for the MOMA signature.
    /// * `dna_sequence` - The DNA sequence to apply the simulated mutation to.
    ///
    /// # Returns
    /// An `Option<(u64, Mutation)>` containing the numeric signature and the
    /// resulting `Mutation` analysis. Returns `None` if the sequence is too short
    /// or the signature points to an invalid position.
    pub fn analyze(&self, p: u64, dna_sequence: &str) -> Option<(u64, Mutation)> {
        // 1. Generate the core MOMA signature.
        let signature = self.ring.signature(p);

        // 2. Use the signature to determine the mutation site.
        let mutation_pos = (signature % dna_sequence.len() as u64) as usize;

        // 3. Determine the codon affected by the mutation.
        let codon_start = (mutation_pos / 3) * 3;
        if codon_start + 3 > dna_sequence.len() {
            return None; // Not enough sequence left for a full codon.
        }
        let original_codon_str = &dna_sequence[codon_start..codon_start + 3];

        // 4. Translate the original codon.
        let original_aa = self.codon_table.translate(original_codon_str)?;

        // 5. Simulate the mutation by changing the base at the mutation position.
        let mut mutated_sequence = dna_sequence.to_string();
        let original_char = mutated_sequence.chars().nth(mutation_pos).unwrap();
        // Simple mutation: cycle through A -> C -> G -> T -> A
        let new_char = match original_char {
            'A' => 'C',
            'C' => 'G',
            'G' => 'T',
            'T' => 'A',
            _ => return None, // Invalid character in sequence
        };
        mutated_sequence.replace_range(mutation_pos..mutation_pos + 1, &new_char.to_string());

        // 6. Analyze the new, mutated codon.
        let mutated_codon_str = &mutated_sequence[codon_start..codon_start + 3];
        let mutated_aa = self.codon_table.translate(mutated_codon_str)?;

        // 7. Create and return the final analysis.
        let mutation = Mutation::new(
            original_codon_str.to_string(),
            mutated_codon_str.to_string(),
            original_aa,
            mutated_aa,
        );

        Some((signature, mutation))
    }
}