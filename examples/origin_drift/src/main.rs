use moma::origin_drift::{OriginDrifter};

fn main() {
    println!();
    println!("Origin Drifter Simulation");
    println!("===========================");
    let mut drifter = OriginDrifter::new(12, 0, 1);
    drifter.run(24);

    for snapshot in drifter.history() {
        println!(
            "Origin: {:2}, Entropy: {:.4}, Residues: {:?}",
            snapshot.origin, snapshot.entropy, snapshot.residues
        );
    }
    println!();
    println!("Simulation complete.\n");
    println!("Total steps: {}", drifter.history().len());
    println!("Final origin: {}", drifter.history().last().unwrap().origin);
    println!("Final entropy: {:.4}", drifter.history().last().unwrap().entropy);
    println!("Final residues: {:?}", drifter.history().last().unwrap().residues);
    println!("\nEnd of simulation.");
    println!();
}