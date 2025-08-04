use std::fs::File;
use std::io::{Write, BufWriter};
//use plotters::prelude::*;
use crate::mutation::MutationTimeline;

pub fn write_csv(path: &str, data: &[f64]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (i, value) in data.iter().enumerate() {
        writeln!(writer, "{},{}", i, value)?;
    }

    Ok(())
}

/* 
pub fn plot_mutation_timeline(timeline: &MutationTimeline, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (1280, 720)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_step = timeline.steps as f64;
    let max_entropy = timeline.snapshots.iter().map(|s| s.entropy).fold(0.0, f64::max);
    let max_resonance = timeline.snapshots.iter().map(|s| s.resonance).fold(0.0, f64::max);

    let mut chart = ChartBuilder::on(&root)
        .caption("Mutation Timeline", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0f64..max_step, 0f64..max_entropy.max(max_resonance))?;

    chart.configure_mesh()
        .x_desc("Step")
        .y_desc("Score")
        .draw()?;

    chart.draw_series(LineSeries::new(
        timeline.snapshots.iter().map(|s| (s.step as f64, s.entropy)),
        &RED,
    ))?.label("Entropy").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart.draw_series(LineSeries::new(
        timeline.snapshots.iter().map(|s| (s.step as f64, s.resonance)),
        &BLUE,
    ))?.label("Resonance").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}



pub fn plot_entropy_heatmap(timeline: &MutationTimeline, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let codon_count = 1;
    let steps = timeline.snapshots.len();

    let root = BitMapBackend::new(filename, (codon_count as u32 * 10, steps as u32 * 10)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_entropy = timeline.snapshots
        .iter()
        .map(|s| s.entropy)
        .fold(0.0, f64::max);

    for (step_idx, snapshot) in timeline.snapshots.iter().enumerate() {
        let entropy = snapshot.entropy;
        let intensity = (entropy / max_entropy).min(1.0);
        let color = RGBColor(
            (255.0 * intensity) as u8,
            (255.0 * (1.0 - intensity)) as u8,
            100,
        );
        root.draw(&Rectangle::new(
            [
                (0, step_idx as i32 * 10),
                (10, (step_idx + 1) as i32 * 10),
            ],
            color.filled(),
        ))?;
    }

    Ok(())
}

*/