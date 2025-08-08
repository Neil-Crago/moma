//! This module provides functionality to handle barycentric coordinates

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