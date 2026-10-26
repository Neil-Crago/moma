//! Real-structure accuracy check: the alpha helices of a PDB entry.
//!
//! The C-alpha trace of an ideal alpha helix has known geometry: radius 2.3 Å,
//! 100° of rotation and 1.5 Å of rise per residue (3.6 residues per turn). Real
//! helices bend and distort, so values near those, not equal to them, are expected.
//! The recovered period, radius, rotation per residue and rise per residue are
//! printed for every long alpha-helix segment (class 1) listed in the file's HELIX
//! records. Uses the 4 Å flagellar filament entry 1UCU (a Salmonella flagellin).
//!
//! Run with: `cargo run --release --example pdb_helices -- data/1UCU.pdb`

use moma::helix::{auto_axis_3d, moving_axis_period_3d, to_cylindrical, V3};
use moma::spiral::unwrap_phase;
use std::f64::consts::TAU;

fn field(line: &str, a: usize, b: usize) -> &str {
    line.get(a..b.min(line.len())).unwrap_or("").trim()
}

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| "data/1UCU.pdb".into());
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));

    let mut ca: Vec<(char, i32, V3)> = Vec::new();
    let mut helices: Vec<(char, i32, i32)> = Vec::new();
    for line in text.lines() {
        if line.starts_with("ATOM") && field(line, 12, 16) == "CA" && !line[16..17].contains('B') {
            let num = |a, b| field(line, a, b).parse::<f64>();
            if let (Ok(x), Ok(y), Ok(z), Ok(r)) = (num(30, 38), num(38, 46), num(46, 54), num(22, 26)) {
                ca.push((line.as_bytes()[21] as char, r as i32, [x, y, z]));
            }
        } else if line.starts_with("HELIX")
            && field(line, 38, 40) == "1"
            && let (Ok(a), Ok(b)) = (field(line, 21, 25).parse(), field(line, 33, 37).parse())
        {
            helices.push((line.as_bytes()[19] as char, a, b));
        }
    }
    println!("{path}: {} C-alpha atoms, {} alpha-helix segments\n", ca.len(), helices.len());
    println!("expected for an ideal alpha helix: period 3.6, radius 2.3 A, 100 deg/residue, rise 1.5 A\n");
    println!("residues      n   period  radius  deg/res  rise/res  axis bend (deg)");

    let mut all = Vec::new();
    for (chain, a, b) in helices {
        let pts: Vec<V3> = ca.iter().filter(|c| c.0 == chain && c.1 >= a && c.1 <= b).map(|c| c.2).collect();
        let n = pts.len();
        if n < 16 {
            continue;
        }
        let Some((_, period)) = auto_axis_3d(&pts, 2.5, 6.0) else { continue };
        let axis = moving_axis_period_3d(&pts, period);
        let trim = (period / 2.0).ceil() as usize + 1;
        let (polar, _) = to_cylindrical(&pts, &axis);
        let u = unwrap_phase(&polar);
        let span = trim..n - trim;
        let steps = (span.len() - 1) as f64;
        let mut radii: Vec<f64> = polar[span.clone()].iter().map(|p| p.radius).collect();
        radii.sort_by(|x, y| x.total_cmp(y));
        let deg = ((u[n - trim - 1] - u[trim]) / steps).abs() * 360.0 / TAU;
        let rise = axis[span.clone()]
            .windows(2)
            .map(|p| (0..3).map(|k| (p[1][k] - p[0][k]).powi(2)).sum::<f64>().sqrt())
            .sum::<f64>()
            / steps;
        let end_dir = |i: usize, j: usize| {
            let d = [0, 1, 2].map(|k| axis[j][k] - axis[i][k]);
            let m = d.iter().map(|v| v * v).sum::<f64>().sqrt();
            d.map(|v| v / m)
        };
        let (h, q) = (trim + span.len() / 3, n - trim - span.len() / 3);
        let (d1, d2) = (end_dir(trim, h), end_dir(q, n - trim - 1));
        let bend = (0..3).map(|k| d1[k] * d2[k]).sum::<f64>().clamp(-1.0, 1.0).acos().to_degrees();
        println!(
            "{chain}{a:>4}-{b:<4}  {n:>4}   {period:>5.2}  {:>6.2}  {deg:>7.1}  {rise:>8.2}  {bend:>6.1}",
            radii[radii.len() / 2]
        );
        all.push((radii[radii.len() / 2], deg, rise, n));
    }
    if !all.is_empty() {
        let mean = |f: fn(&(f64, f64, f64, usize)) -> f64| all.iter().map(f).sum::<f64>() / all.len() as f64;
        println!(
            "\nmean over {} helices: radius {:.2} A, {:.1} deg/residue, rise {:.2} A/residue",
            all.len(),
            mean(|r| r.0),
            mean(|r| r.1),
            mean(|r| r.2)
        );
    }
}
