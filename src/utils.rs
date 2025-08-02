use std::fs::File;
use std::io::{Write, BufWriter};

pub fn write_csv(path: &str, data: &[f64]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (i, value) in data.iter().enumerate() {
        writeln!(writer, "{},{}", i, value)?;
    }

    Ok(())
}