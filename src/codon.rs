//! Provides structures and logic for working with genetic codons and the standard DNA codon table.

use std::collections::HashMap;

/// Represents a single amino acid or a Stop signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AminoAcid {
    Alanine, Arginine, Asparagine, AsparticAcid, Cysteine,
    GlutamicAcid, Glutamine, Glycine, Histidine, Isoleucine,
    Leucine, Lysine, Methionine, Phenylalanine, Proline,

    Serine, Threonine, Tryptophan, Tyrosine, Valine,
    Stop, // Represents a translation stop signal
}

/// A struct that holds the standard DNA codon translation table.
/// It maps three-letter DNA codons (e.g., "AUG") to their corresponding amino acids.
#[derive(Debug)]
pub struct CodonTable {
    map: HashMap<String, AminoAcid>,
}

impl Default for CodonTable {
    /// Creates a new `CodonTable` populated with the standard genetic code.
    fn default() -> Self {
        let mut map = HashMap::new();
        // Alanine
        map.insert("GCU".to_string(), AminoAcid::Alanine);
        map.insert("GCC".to_string(), AminoAcid::Alanine);
        map.insert("GCA".to_string(), AminoAcid::Alanine);
        map.insert("GCG".to_string(), AminoAcid::Alanine);
        // Arginine
        map.insert("CGU".to_string(), AminoAcid::Arginine);
        map.insert("CGC".to_string(), AminoAcid::Arginine);
        map.insert("CGA".to_string(), AminoAcid::Arginine);
        map.insert("CGG".to_string(), AminoAcid::Arginine);
        map.insert("AGA".to_string(), AminoAcid::Arginine);
        map.insert("AGG".to_string(), AminoAcid::Arginine);
        // ... and so on for all other amino acids ...
        // Phenylalanine
        map.insert("UUU".to_string(), AminoAcid::Phenylalanine);
        map.insert("UUC".to_string(), AminoAcid::Phenylalanine);
        // Leucine
        map.insert("UUA".to_string(), AminoAcid::Leucine);
        map.insert("UUG".to_string(), AminoAcid::Leucine);
        map.insert("CUU".to_string(), AminoAcid::Leucine);
        map.insert("CUC".to_string(), AminoAcid::Leucine);
        map.insert("CUA".to_string(), AminoAcid::Leucine);
        map.insert("CUG".to_string(), AminoAcid::Leucine);
        // Stop codons
        map.insert("UAA".to_string(), AminoAcid::Stop);
        map.insert("UAG".to_string(), AminoAcid::Stop);
        map.insert("UGA".to_string(), AminoAcid::Stop);
        // Methionine (Start codon)
        map.insert("AUG".to_string(), AminoAcid::Methionine);

        Self { map }
    }
}

impl CodonTable {
    /// Creates a new, populated `CodonTable`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Translates a three-letter DNA codon string into an `Option<AminoAcid>`.
    ///
    /// # Arguments
    /// * `codon` - A string slice representing the 3-letter codon (e.g., "AUG").
    ///
    /// # Returns
    /// `Some(AminoAcid)` if the codon is valid, or `None` if it is not found in the table.
    pub fn translate(&self, codon: &str) -> Option<AminoAcid> {
        // In a real-world scenario, you might want to handle 'T' vs 'U' (Thymine vs Uracil)
        let rna_codon = codon.replace('T', "U");
        self.map.get(&rna_codon).cloned()
    }
}