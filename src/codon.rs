use rand::Rng;

pub fn encode_codon(codon: &str) -> Option<usize> {
    let base_map = |c| match c {
        'A' => Some(0),
        'C' => Some(1),
        'G' => Some(2),
        'T' => Some(3),
        _ => None,
    };
    let chars: Vec<_> = codon.chars().collect();
    if chars.len() != 3 {
        return None;
    }
    Some(
        base_map(chars[0])? * 16 +
        base_map(chars[1])? * 4 +
        base_map(chars[2])?
    )
}

pub struct CodonRing {
    pub residues: Vec<usize>, // codon values
    pub origin: usize,
    pub modulus: usize,
}

impl CodonRing {
    pub fn from_codons(codons: &[&str], origin: usize) -> Self {
        let residues = codons.iter()
            .filter_map(|c| encode_codon(c))
            .collect();
        CodonRing {
            residues,
            origin,
            modulus: 64,
        }
    }
}

impl CodonRing {
    pub fn drift_origin(&mut self, delta: isize) {
        let new_origin = ((self.origin as isize + delta) % self.modulus as isize + self.modulus as isize) % self.modulus as isize;
        self.origin = new_origin as usize;
    }
}

impl CodonRing {
    pub fn shifted_residues(&self) -> Vec<usize> {
        self.residues.iter()
            .map(|r| (r + self.origin) % self.modulus)
            .collect()
    }
}

pub struct EntropyPulse {
    pub strength: f64,       // 0.0 to 1.0
    pub bias: Option<usize>, // Optional target residue
    pub mode: PulseMode,
}

pub enum PulseMode {
    OriginDrift,
    ResiduePerturb,
    TargetedDisruption,
}

impl CodonRing {
    pub fn apply_entropy_pulse(&mut self, pulse: &EntropyPulse) {
        match pulse.mode {
            PulseMode::OriginDrift => {
                let delta = ((pulse.strength * self.modulus as f64).round() as isize)
                    * if rand::random::<bool>() { 1 } else { -1 };
                self.drift_origin(delta);
            }
            PulseMode::ResiduePerturb => {
                for r in &mut self.residues {
                    *r = (*r + rand::rng().random_range(0u32..self.modulus as u32) as usize) % self.modulus;
                    if rand::random::<f64>() < pulse.strength {
                        *r = (*r + (rand::rng().random_range(0u32..self.modulus as u32) as usize)) % self.modulus;
                    }
                }
            }
            PulseMode::TargetedDisruption => {
                if let Some(target) = pulse.bias {
                    for r in &mut self.residues {
                        *r = rand::rng().random_range(0..self.modulus);
                        if *r == target && rand::random::<f64>() < pulse.strength {
                            *r = rand::rng().random_range(0..self.modulus);
                        }
                    }
                }
            }
        }
    }
}
