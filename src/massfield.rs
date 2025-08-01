pub use crate::massfield::massfield_overlay::{
    MassFieldOverlay, ModularOverlay, OverlayRegistry
};

pub mod massfield_overlay {
use std::collections::HashMap;    
/// Represents a mass field over modular space
pub struct MassFieldOverlay {
    pub modulus: usize,
    pub masses: HashMap<usize, f64>, // composite → mass
    pub radius: usize,               // influence radius
}

/// Result of applying mass field to a residue ring
#[derive(Debug)]
pub struct DistortedResidueField {
    pub origin: usize,
    pub residues: Vec<usize>,
    pub distortions: Vec<f64>, // distortion per residue
}

impl MassFieldOverlay {
    /// Apply mass field to a residue ring at given origin
    pub fn apply(&self, origin: usize) -> DistortedResidueField {
        let residues: Vec<usize> = (0..self.modulus)
            .map(|x| (x + self.modulus - origin) % self.modulus)
            .collect();

        let distortions: Vec<f64> = residues
            .iter()
            .map(|&r| {
                self.masses.iter().map(|(&comp, &mass)| {
                    let dist = ((r as isize - comp as isize).abs().min(self.radius as isize)) as f64;
                    if dist == 0.0 { mass } else { mass / dist }
                }).sum()
            })
            .collect();

        DistortedResidueField {
            origin,
            residues,
            distortions,
        }
    }
}


pub trait ModularOverlay {
    fn apply(&self, origin: usize, modulus: usize) -> Vec<f64>;
    fn name(&self) -> &str;
}

impl ModularOverlay for MassFieldOverlay {
    fn apply(&self, origin: usize, modulus: usize) -> Vec<f64> {
        let residues: Vec<usize> = (0..modulus)
            .map(|x| (x + modulus - origin) % modulus)
            .collect();

        residues.iter().map(|&r| {
            self.masses.iter().map(|(&comp, &mass)| {
                let dist = ((r as isize - comp as isize).abs().min(self.radius as isize)) as f64;
                if dist == 0.0 { mass } else { mass / dist }
            }).sum()
        }).collect()
    }

    fn name(&self) -> &str {
        "MassFieldOverlay"
    }
}

pub struct OverlayRegistry {
    pub overlays: Vec<Box<dyn ModularOverlay>>,
}

impl Default for OverlayRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl OverlayRegistry {
    pub fn new() -> Self {
        Self { overlays: vec![] }
    }

    pub fn register(&mut self, overlay: Box<dyn ModularOverlay>) {
        self.overlays.push(overlay);
    }

    pub fn apply_all(&self, origin: usize, modulus: usize) -> Vec<(String, Vec<f64>)> {
        self.overlays.iter().map(|o| {
            (o.name().to_string(), o.apply(origin, modulus))
        }).collect()
    }
}
}

