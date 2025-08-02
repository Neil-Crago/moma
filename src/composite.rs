pub trait CompositeDampener {
    fn dampen(&self, value: f64) -> f64;
}

/// Example: mass-based dampening
pub struct MassDampener {
    pub mass: f64,
}

impl CompositeDampener for MassDampener {
    fn dampen(&self, value: f64) -> f64 {
        value / (1.0 + self.mass)
    }
}