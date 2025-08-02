use crate::core::core::{MomaRing, OriginStrategy};
use std::marker::PhantomData;

/// A tool to measure the "drift" or volatility of MOMA signatures over a sequence of primes.
///
/// This struct is generic over an `OriginStrategy`. It maintains a history of calculated
/// signatures and can compute the average change (drift) between consecutive signatures.
/// This can be used to analyze the stability or chaotic nature of a given strategy.
pub struct OriginDrift<S: OriginStrategy> {
    ring: MomaRing<S>,
    history: Vec<u64>,
    // PhantomData is used because S is part of the struct's logic but not a field.
    _strategy: PhantomData<S>,
}

impl<S: OriginStrategy> OriginDrift<S> {
    /// Creates a new `OriginDrift` analyzer for a given modulus and strategy.
    ///
    /// # Parameters
    /// - `modulus`: The modulus for the internal `MomaRing`.
    /// - `strategy`: An instance of a struct that implements `OriginStrategy`.
    pub fn new(modulus: u64, strategy: S) -> Self {
        Self {
            ring: MomaRing::new(modulus, strategy),
            history: Vec::new(),
            _strategy: PhantomData,
        }
    }

    /// Calculates the MOMA signature for the next prime in a sequence and records it.
    ///
    /// # Parameters
    /// - `p`: The prime number to analyze.
    ///
    /// # Returns
    /// The calculated signature for the prime `p`.
    pub fn next(&mut self, p: u64) -> u64 {
        let signature = self.ring.signature(p);
        self.history.push(signature);
        signature
    }

    /// Calculates the average absolute difference between consecutive signatures in the history.
    ///
    /// A higher value indicates greater "drift" or volatility for the chosen strategy.
    /// A value of 0.0 means the signatures have been stable or there's not enough history.
    ///
    /// # Returns
    /// The average drift magnitude as an `f64`.
    pub fn drift_magnitude(&self) -> f64 {
        if self.history.len() < 2 {
            return 0.0;
        }
        let deltas: Vec<f64> = self
            .history
            .windows(2)
            .map(|w| (w[1] as i64 - w[0] as i64).abs() as f64)
            .collect();

        deltas.iter().sum::<f64>() / deltas.len() as f64
    }

    /// Returns a slice of the recorded signature history.
    pub fn history(&self) -> &[u64] {
        &self.history
    }
}
