use crate::codon::CodonRing;
pub trait BioSignature {
    fn modular_spread(&self) -> f64;
    fn entropy_footprint(&self) -> f64;
    fn resonance_score(&self) -> f64;
    fn drift_vector(&self) -> isize;
}

impl BioSignature for CodonRing {
    fn modular_spread(&self) -> f64 {
        let mut counts = vec![0usize; self.modulus];
        for &r in &self.residues {
            counts[r % self.modulus] += 1;
        }
        let mean = self.residues.len() as f64 / self.modulus as f64;
        let variance = counts.iter()
            .map(|&c| (c as f64 - mean).powi(2))
            .sum::<f64>() / self.modulus as f64;
        variance.sqrt() // standard deviation
    }

    fn entropy_footprint(&self) -> f64 {
        let mut counts = vec![0usize; self.modulus];
        for &r in &self.residues {
            counts[r % self.modulus] += 1;
        }
        let total = self.residues.len() as f64;
        counts.iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / total;
                -p * p.log2()
            })
            .sum()
    }

    fn resonance_score(&self) -> f64 {
        // Placeholder: could use FFT or autocorrelation later
        let diffs: Vec<isize> = self.residues.windows(2)
            .map(|w| {
                let a = ((w[1] + self.origin) % self.modulus) as isize;
                let b = ((w[0] + self.origin) % self.modulus) as isize;
                a - b
            })
            .collect();
        let mean = diffs.iter().sum::<isize>() as f64 / diffs.len() as f64;
        let variance = diffs.iter()
            .map(|&d| (d as f64 - mean).powi(2))
            .sum::<f64>() / diffs.len() as f64;
        1.0 / (1.0 + variance) // lower variance = higher resonance
    }

    fn drift_vector(&self) -> isize {
        self.origin as isize
    }
}
