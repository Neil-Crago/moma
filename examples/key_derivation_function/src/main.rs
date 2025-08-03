//! # MOMA Key Derivation Function (KDF) Demonstration
//!
//! This program implements a KDF using the moma crate. A KDF takes a
//! low-entropy secret (like a password) and derives a strong, fixed-length
//! cryptographic key from it.

use moma::core::MomaRing;

use moma::primes;
use moma::strategy;

use sha2::{Digest, Sha256};

/// A Key Derivation Function based on Moving Origin Modular Arithmetic.
pub struct MomaKdf<'a> {
    password: &'a [u8],
    salt: &'a [u8],
    iterations: u32,
    output_len: usize,
}

impl<'a> MomaKdf<'a> {
    /// Configures a new MOMA KDF derivation.
    ///
    /// # Parameters
    /// - `password`: The secret input, typically a user's password.
    /// - `salt`: A public, random value unique to each password.
    /// - `iterations`: The work factor. Higher values are more secure but slower.
    /// - `output_len`: The desired length of the final key in bytes.
    pub fn new(
        password: &'a [u8],
        salt: &'a [u8],
        iterations: u32,
        output_len: usize,
    ) -> Self {
        Self {
            password,
            salt,
            iterations,
            output_len,
        }
    }

    /// Derives the cryptographic key.
    ///
    /// This is the core function that performs the key stretching. It uses MOMA's
    /// signature function in a chained loop to create a computationally-intensive
    /// process that is difficult to brute-force.
    pub fn derive_key(&self) -> Vec<u8> {
        // --- 1. Seeding Phase ---
        // Use SHA-256 to deterministically generate MOMA parameters from the password and salt.
        // This ensures that the same inputs always produce the same key.
        let modulus_seed = Sha256::digest([self.password, self.salt].concat());
        let prime_seed = Sha256::digest([self.salt, self.password].concat());

        let modulus_u32 = u32::from_le_bytes(modulus_seed[..4].try_into().unwrap());
        let prime_u32 = u32::from_le_bytes(prime_seed[..4].try_into().unwrap());

        // The cast to u64 is necessary to match the type used in the moma crate.
        let modulus = primes::next_prime(modulus_u32 as u64);
        let mut current_prime = primes::next_prime(prime_u32 as u64);
        
        // --- 2. Iteration Phase ---
        // We use CompositeMass as it's a more complex strategy, making for a better KDF demo.
        let ring = MomaRing::new(modulus, strategy::PrimeGap);
        // let ring = MomaRing::new(modulus, strategy::CompositeMass);
        let mut derived_bytes = Vec::with_capacity(self.iterations as usize * 8);

        for _ in 0..self.iterations {
            // Calculate the MOMA signature.
            let residue = ring.signature(current_prime);
            derived_bytes.extend_from_slice(&residue.to_le_bytes());

            // Update the state in a dependent way. The next prime depends on the
            // previous residue, creating a chain that cannot be parallelized.
            current_prime = primes::next_prime(current_prime.wrapping_add(residue));
        }

        // --- 3. Finalization Phase ---
        // Hash the entire sequence of derived bytes to produce the final key.
        // This mixes all the results together into a uniformly distributed output.
        let final_hash = Sha256::digest(&derived_bytes);
        final_hash[..self.output_len.min(32)].to_vec()
    }
}


fn main() {
    println!("\n--- MOMA KDF Demonstration --- 🔑");

    let password = "a_very_secure_password_123";
    let salt = "a_unique_random_salt_for_each_user"; // In a real app, this MUST be random.
    
    // A low iteration count for a quick demo. Real applications would use >100,000.
    let iterations = 1_000; 
    let key_length = 32; // For a 256-bit key (e.g., for AES-256)

    println!("Password: \"{}\"", password);
    println!("Salt:     \"{}\"", salt);
    println!("Iterations: {}", iterations);
    println!("------------------------------------");

    // --- Derivation 1 ---
    let kdf = MomaKdf::new(password.as_bytes(), salt.as_bytes(), iterations, key_length);
    let derived_key1 = kdf.derive_key();
    println!("✅ Derived Key 1: {}", hex::encode(&derived_key1));

    // --- Derivation 2 (showing a different password yields a different key) ---
    let other_password = "a_different_password";
    let kdf_other = MomaKdf::new(other_password.as_bytes(), salt.as_bytes(), iterations, key_length);
    let derived_key_other = kdf_other.derive_key();
    println!("🔑 Derived Key (other pass): {}", hex::encode(&derived_key_other));

    // --- Derivation 3 (showing it is deterministic) ---
    let kdf_same = MomaKdf::new(password.as_bytes(), salt.as_bytes(), iterations, key_length);
    let derived_key_same = kdf_same.derive_key();
    println!("✅ Derived Key (same input): {}", hex::encode(&derived_key_same));

    assert_eq!(derived_key1, derived_key_same);
    assert_ne!(derived_key1, derived_key_other);
}
