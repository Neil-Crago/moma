<p align="left">
  <img src="docs/moma_trim.png" alt="MOMA Logo"/>
</p>  


# MOMA: Moving Origin Modular Arithmetic

[![Crates.io](https://img.shields.io/crates/v/moma.svg?style=flat-square)](https://crates.io/crates/moma)
[![Docs.rs](https://img.shields.io/docsrs/moma?style=flat-square)](https://docs.rs/moma)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/neil-crago/moma/actions/workflows/rust.yml/badge.svg)](https://github.com/neil-crago/moma/actions/workflows/rust.yml)


**MOMA** is a Rust framework for exploring number theory, cryptography, and bioinformatics through the lens of **Moving Origin Modular Arithmetic**.

The crate is designed for researchers and developers who are interested in a novel, relational framework for analyzing  complex sequences.

---

## The Core Idea: A Barycenter for Numbers

The inspiration for this crate comes from the concept of a **barycenter** in astrophysics. Just as the Earth and Moon orbit a common center of mass that is not the exact center of the Earth, MOMA treats modular arithmetic as a system where the "zero point" or "origin" is not fixed.

This origin shifts dynamically based on a contextual value—typically a prime number `p`—and the chosen `OriginStrategy`. This provides a novel relational framework for analyzing complex systems.

> The original inspiration came from this NASA article: [What Is a Barycenter?](https://spaceplace.nasa.gov/barycenter/en/)


## Core Concepts

The MOMA framework is built on a few simple but powerful concepts:

* **`MomaRing`**: The primary object for all calculations. A ring is defined by a `modulus` and a chosen `OriginStrategy`.
* **`OriginStrategy`**: A trait that defines *how* the origin moves, making the framework highly extensible.
* **Analysis Tools**: A suite of modular tools for deeper analysis:
    * **Number Theory**: `MassField`, `OriginDrift`, `CompositeInfluence`, `CompositeDampener`, `GoldbachProjector`, `Entropy`, and `ResonanceFinder`.
    * **Bioinformatics**: `BioSigAnalyzer`, `CodonTable`, and `Mutation` for mapping numeric signatures to biological events.


## Features

* **Flexible Core**: A powerful and extensible system based on the `MomaRing` and `OriginStrategy` trait.
* **Advanced Analysis Tools**: A suite of high-level structs for statistical, number-theoretic, and bioinformatics analysis.
* **Cryptographic Primitives**: Demonstrates how MOMA can be used to build components like a Key Derivation Function.
* **Prime Number Utilities**: A helper `primes` module for primality testing and prime generation.
* **Pure Rust**: Built with safe, idiomatic Rust.

---

## Installation

Add MOMA to your `Cargo.toml`:

```toml
[dependencies]
moma = "0.3.6" # Replace with the latest version
````

or rather just run `cargo add moma` from the terminal


## Exploring MOMA

MOMA is a toolkit for exploration across different domains.

### Example 1: Finding "Resonance"

Use the `ResonanceFinder` to find primes where the MOMA signature is a multiple of another property of the prime, such as its number of prime factors (`prime_factor_mass`).

```rust
use moma::resonance::ResonanceFinder;
use moma::strategy;
use moma::primes;

// Find primes where the signature (using CompositeMass strategy)
// is divisible by the prime's factor mass.
let finder = ResonanceFinder::new(
    100,
    strategy::CompositeMass,
    primes::prime_factor_mass,
);

// Search for resonance events in the range 1 to 500.
let resonances = finder.find_in_range(1, 500);

println!("Found {} resonance events.", resonances.len());
for (prime, signature) in resonances {
    println!("  - Resonance at p={} with signature {}", prime, signature);
}
```

### Example 2: Bioinformatics Signature Analysis

Use the `BioSigAnalyzer` to map MOMA signatures to simulated genetic mutations.

```rust
use moma::biosig::BioSigAnalyzer;
use moma::strategy;

let analyzer = BioSigAnalyzer::new(60, strategy::CompositeMass);
let dna_sequence = "AGCTGCGATCGTACGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";

// Analyze the mutational effect of the signature derived from prime 13.
if let Some((signature, mutation)) = analyzer.analyze(13, dna_sequence) {
    println!("Signature for p=13 is {}", signature);
    println!("Resulting mutation type: {:?}", mutation.mutation_type);
}
```

## More Examples

There are more examples in the GitHub repository including:-

* **Cosmology** - For the Astrophysicists
* **BioInformatics** - For Biologists
* **Prime Gaps** - For the number theorists
* **Goldbach Projector** - For Goldbach Conjecture theorists
* **Origin Drift** - to understand MOMA's drifting Origin 
* **Key Derivation Function (KDF)** - For Cryptographers
* **Mass Field** - How to create a massfield for a specific range

## Author

Neil Crago — experimental mathematician


## Contributing

Contributions are welcome\! If you have an idea for a new `OriginStrategy`, an analysis tool, or find a bug, please feel free to open an issue or submit a pull request.

## License

This project is licensed under either of:

  * Apache License, Version 2.0, ([LICENSE-APACHE](https://www.google.com/search?q=LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
  * MIT license ([LICENSE-MIT](https://www.google.com/search?q=LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option

## Related Crates

This crate is part of a collection of crates by the same author:
These include:-
  * MOMA_simulation_engine
  * Fractal_Algebra
  * tma_engine
  * factorial_engine
  * fa_slow_ai
  * coheron
  * curvature