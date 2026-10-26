//! # MOMA: Moving Origin Modular Arithmetic
//!
//! A toolkit built on one idea: measure things relative to an origin that *moves*.
//!
//! * **Moving-axis analysis** ([`spiral`], [`helix`]): a spiral around a wandering
//!   axis looks irregular in fixed coordinates but is a steady rotation relative
//!   to the axis. Estimate the axis (fixed, adaptive, period-matched or jointly
//!   fitted), repair dirty tracks, find the period, and read off radius, phase,
//!   winding and ring signatures.
//! * **Barycentric consensus** ([`consensus`]): weighted circular means on `Z_m`,
//!   robust (geometric-median) barycentres, and consensus of several axis
//!   estimates, where an estimator that disagrees with the rest is outvoted.
//!
//! See the README for what has been validated on real data and what has not.

pub mod barycentric;
pub mod consensus;
pub mod helix;
pub mod spiral;
pub mod utils;

pub use crate::barycentric::OriginShift;
pub use crate::utils::write_csv;
