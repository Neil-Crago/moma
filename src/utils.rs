//! Utilities for writing data to CSV files and plotting mutation timelines.
 //! This module provides functions to write numerical data to CSV format and
 //! to plot mutation timelines using the Plotters library.
 //! It is designed to work with the MOMA framework, particularly for analyzing
 //! mutation data in bioinformatics applications.
 //! 
use std::fs::File;
use std::io::{Write, BufWriter};
//use plotters::prelude::*;

pub fn write_csv(path: &str, data: &[f64]) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (i, value) in data.iter().enumerate() {
        writeln!(writer, "{},{}", i, value)?;
    }

    Ok(())
}

