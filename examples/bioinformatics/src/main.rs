use moma::codon::{CodonRing, EntropyPulse, PulseMode};
use moma::biosig::BioSignature;
use moma::mutation::MutationTimeline;
use moma::codon::encode_codon;
use moma::utils::{
    plot_entropy_heatmap,
    plot_mutation_timeline,
};

fn plot_heatmap() -> Result<(), Box<dyn std::error::Error>> {
    let mut timeline = build_multi_pulse_timeline();
    timeline.run();

    plot_entropy_heatmap(&timeline, "entropy_heatmap.png")?;
    Ok(())
}

fn plot_timeline() -> Result<(), Box<dyn std::error::Error>> {
    let mut timeline = build_multi_pulse_timeline(); // your setup
    timeline.run();

    plot_mutation_timeline(&timeline, "timeline_plot.png")?;
    Ok(())
}

fn build_multi_pulse_timeline() -> MutationTimeline {
    let codons = vec!["ATG", "CGA", "TAA", "GGC", "TGC", "AAA"];
    let ring = CodonRing::from_codons(&codons, 0);

    let pulses = vec![
        EntropyPulse {
            strength: 0.1,
            bias: None,
            mode: PulseMode::OriginDrift,
        },
        EntropyPulse {
            strength: 0.3,
            bias: None,
            mode: PulseMode::ResiduePerturb,
        },
        EntropyPulse {
            strength: 0.2,
            bias: Some(encode_codon("CGA").unwrap()),
            mode: PulseMode::TargetedDisruption,
        },
    ];

    let mut timeline = MutationTimeline::new(ring, pulses, 60);
    timeline.run();
    timeline
   
    //let _ = timeline.export_csv("multi_pulse_timeline.csv");
}

fn mutation_timeline() -> MutationTimeline {
    let codons = vec!["ATG", "CGA", "TAA", "GGC", "TGC", "AAA"];
    let ring = CodonRing::from_codons(&codons, 0);

    let pulse = EntropyPulse {
        strength: 0.2,
        bias: None,
        mode: PulseMode::ResiduePerturb,
    };

    let mut timeline = MutationTimeline::new(ring, vec![pulse], 50);
    timeline.run();
    timeline
}



fn biosignature() {
    let codons = vec!["ATG", "CGA", "TAA", "GGC"];
    let ring = CodonRing::from_codons(&codons, 5);
    let signature = (
        ring.modular_spread(),
        ring.entropy_footprint(),
        ring.resonance_score(),
        ring.drift_vector(),
    );
    println!("BioSignature: {:?}", signature);
}

fn perturbed_residues() {
    let codons = vec!["ATG", "CGA", "TAA", "GGC"];
    let mut ring = CodonRing::from_codons(&codons, 0);

    let pulse = EntropyPulse {
        strength: 0.3,
        bias: Some(encode_codon("CGA").unwrap()),
        mode: PulseMode::TargetedDisruption,
    };

    ring.apply_entropy_pulse(&pulse);
    println!("Perturbed residues: {:?}", ring.residues);
}

fn original_residues() {
    let codons = vec!["ATG", "CGA", "TAA", "GGC"];
    let mut ring = CodonRing::from_codons(&codons, 0);
    println!("Original residues: {:?}", ring.residues);
    ring.drift_origin(5);
    println!("Shifted residues: {:?}", ring.shifted_residues());
}

fn main() {
    original_residues();
    perturbed_residues();
    biosignature();
    mutation_timeline();
    plot_timeline().unwrap();
    plot_heatmap().unwrap();
}


