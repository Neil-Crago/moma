use moma::massfield::{MassFieldOverlay, OverlayRegistry}; 
use std::collections::HashMap;

fn main() {
    println!();
    println!("Start of Mass Field Overlay Example");
    println!();
    let mut registry = OverlayRegistry::new();

    let masses = HashMap::from([
        (4, 2.0),
        (8, 1.5),
        (10, 3.0),
    ]);

    let mass_overlay = MassFieldOverlay {
        modulus: 12,
        masses,
        radius: 6,
    };

    registry.register(Box::new(mass_overlay));

    for origin in 0..12 {
        let results = registry.apply_all(origin, 12);
        for (name, distortions) in results {
            println!("Overlay: {} at Origin {}", name, origin);
            for (r, d) in (0..12).zip(distortions.iter()) {
                println!("Residue {} → Distortion {:.2}", r, d);
            }
            println!("---");
        }
    }
    println!();
    println!("End of Mass Field Overlay Example");
    println!();
}
