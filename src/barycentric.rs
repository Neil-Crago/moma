//! Barycentric coordinate system utilities for origin shifts and transformations.
//! This module provides functionality to handle barycentric coordinates
//! and transformations in a 2D space, particularly useful for modular arithmetic
//! and number theory applications.
//! This module is designed to work with the MOMA framework, allowing for
//! efficient calculations of origin shifts and transformations in a modular context.
//! This module is part of the MOMA (Moving Origin Modular Arithmetic) framework,
//! which explores number theory, cryptography, and bioinformatics through the lens of a "moving
//! origin" in modular arithmetic.
//! 
#[derive(Debug, Clone, Copy)]
pub struct OriginShift {
    pub dx: f64,
    pub dy: f64,
}

impl OriginShift {
    pub fn zero() -> Self {
        OriginShift { dx: 0.0, dy: 0.0 }
    }

    pub fn magnitude(&self) -> f64 {
        (self.dx.powi(2) + self.dy.powi(2)).sqrt()
    }
}