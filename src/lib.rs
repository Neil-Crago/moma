//! # MOMA: Moving Origin Modular Arithmetic
//!
//! A framework for exploring number theory concepts using a "moving origin"
//! in modular arithmetic.
//!
//! The core idea of MOMA is that for a given modulus `m`, the "zero point" or "origin"
//! for the calculation `x mod m` is not fixed. Instead, it shifts based on a
//! contextual value, typically a prime number `p`. This allows for the study of
//! properties of numbers (like primes) in a dynamic relational context.
//!
//! ## Core Concepts
//!
//! - **`MomaRing`**: The central object for calculations. It is defined by a `modulus` and an `OriginStrategy`.
//! - **`OriginStrategy`**: A trait that defines *how* the origin moves. The crate provides several
//!   common strategies, and users can easily implement their own.
//! - **Residue**: The result of a MOMA calculation, computed as `(value + origin) % modulus`.
//!
//! ## Example Usage
//!
//! ```
//! use moma::core::core::{MomaRing};
//! use mom::strategy::strategy; 
//! 
//!
//! // Create a MOMA ring with modulus 37.
//! // The origin will be the gap between a prime and the previous prime.
//! let ring = MomaRing::new(37, strategy::PrimeGap);
//!
//! // Let's analyze the prime p = 29.
//! let p = 29;
//!
//! // The value we want to find the residue of. In my original research,
//! // this was the sum of a prime and its predecessor.
//! let value_to_test = p + primes::prev_prime(p); // 29 + 23 = 52
//!
//! // Calculate the MOMA residue.
//! // The origin for p=29 is (29 - 23) = 6.
//! // The residue is (52 + 6) % 37 = 58 % 37 = 21.
//! let residue = ring.residue(value_to_test, p);
//!
//! println!("For p={}, the MOMA residue of {} is {}", p, value_to_test, residue);
//! assert_eq!(residue, 21);
//! ```


// --- Module Declarations ---
pub mod analysis;
pub mod core;
pub mod massfield;
pub mod entropy;      
pub mod goldbach;    
pub mod influence;   
pub mod origin_drift;
pub mod primes;
pub mod strategy;
pub mod utils;
pub mod cosmo;
pub mod barycentric;
pub mod resonance;
pub mod score;
pub mod composite;

// --- Public API Re-exports ---
// This makes the most important structs directly accessible to users.
pub use crate::core::core::{MomaRing, OriginStrategy};
pub use crate::massfield::MassField;
pub use crate::origin_drift::OriginDrift;
pub use crate::primes::primes::{is_prime, next_prime, prev_prime};
pub use crate::strategy::strategy::PrimeGap;
pub use crate::cosmo::{BarycenterSimulator, EntropyPulse, Planet};
pub use crate::barycentric::OriginShift;
pub use crate::resonance::{AutoCorrelationResonance, FFTResonance, ResonanceDetector, score_resonance_strength, detect_fft_resonance, detect_modular_resonance};
pub use crate::score::{score_kurtosis, score_signal_to_noise};
pub use crate::composite::{CompositeDampener, MassDampener};
