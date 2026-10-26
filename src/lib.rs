//! # MOMA: Moving Origin Modular Arithmetic
//!
//! A toolkit built on one idea: measure things relative to an origin that *moves*.
//!
//! * **Moving-axis analysis** ([`spiral`], [`helix`]): a spiral around a wandering
//!   axis looks irregular in fixed coordinates but is a steady rotation relative
//!   to the axis. Estimate the axis (fixed, adaptive or jointly fitted), repair
//!   dirty tracks, and read off radius, phase, winding and ring signatures.
//! * **Ring barycentres** ([`consensus`]): weighted circular means on `Z_m` and
//!   planar barycentres, the right "average" when positions wrap.
//! * **The ring itself** ([`core`], [`strategy`]): `MomaRing` with a pluggable
//!   [`OriginStrategy`]. The prime-based strategies and the number-theory and
//!   bioinformatics modules are the crate's original, more speculative uses.
//!
//! The spiral and helix tools are validated on synthetic data only so far.

// --- Module Declarations ---
pub mod analysis;
pub mod biosig;
pub mod codon;
pub mod composite_field; // Renamed from composite.rs
pub mod core;
pub mod entropy;
pub mod goldbach;
pub mod influence;
pub mod massfield;
pub mod mutation;
pub mod origin_drift;
pub mod primes;
pub mod resonance; // New
pub mod strategy;
pub mod score;
pub mod barycentric; // New
pub mod consensus;
pub mod spiral;
pub mod helix;
pub mod utils;

// --- Public API Re-exports ---
// This makes the most important structs directly accessible to users.
pub use crate::core::{MomaRing, OriginStrategy};
pub use crate::analysis::CompositeDampener;
pub use crate::biosig::BioSigAnalyzer;
pub use crate::composite_field::CompositeField;
pub use crate::entropy::{Entropy, calculate_path_entropy, format_float_to_string};
pub use crate::goldbach::GoldbachProjector;
pub use crate::influence::CompositeInfluence;
pub use crate::massfield::MassField;
pub use crate::origin_drift::OriginDrift;
pub use crate::resonance::ResonanceFinder;
pub use crate::score::{score_signal_to_noise, score_kurtosis};
pub use crate::strategy::{Fixed, PrimeGap, CompositeMass};
pub use crate::primes::{is_prime, next_prime, prev_prime, prime_factor_mass};
pub use crate::mutation::{Mutation, MutationType};
pub use crate::codon::{CodonTable};
pub use crate::barycentric::{OriginShift};
pub use crate::utils::write_csv;
