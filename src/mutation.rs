//! Provides structures for representing and analyzing genetic mutations.

use crate::codon::AminoAcid;

/// Represents the type of a point mutation's effect on the resulting amino acid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MutationType {
    /// The mutation does not change the amino acid.
    Silent,
    /// The mutation changes the amino acid to another.
    Missense,
    /// The mutation changes an amino acid codon to a stop codon.
    Nonsense,
}

/// Represents a single point mutation event.
#[derive(Debug, Clone)]
pub struct Mutation {
    /// The original codon before mutation.
    pub original_codon: String,
    /// The new codon after mutation.
    pub mutated_codon: String,
    /// The original amino acid.
    pub original_amino_acid: AminoAcid,
    /// The new amino acid (or Stop signal).
    pub mutated_amino_acid: AminoAcid,
    /// The classified type of the mutation.
    pub mutation_type: MutationType,
}

impl Mutation {
    /// Analyzes a mutation event and creates a new `Mutation` struct.
    ///
    /// # Arguments
    /// * `original_codon` - The 3-letter codon before the mutation.
    /// * `mutated_codon` - The 3-letter codon after the mutation.
    /// * `original_amino_acid` - The amino acid corresponding to the original codon.
    /// * `mutated_amino_acid` - The amino acid corresponding to the mutated codon.
    ///
    /// # Returns
    /// A new `Mutation` instance with the `mutation_type` correctly classified.
    pub fn new(
        original_codon: String,
        mutated_codon: String,
        original_amino_acid: AminoAcid,
        mutated_amino_acid: AminoAcid,
    ) -> Self {
        let mutation_type = if mutated_amino_acid == AminoAcid::Stop {
            MutationType::Nonsense
        } else if original_amino_acid == mutated_amino_acid {
            MutationType::Silent
        } else {
            MutationType::Missense
        };

        Self {
            original_codon,
            mutated_codon,
            original_amino_acid,
            mutated_amino_acid,
            mutation_type,
        }
    }
}