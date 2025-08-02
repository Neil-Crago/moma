# MOMA: Moving Origin Modular Arithmetic

[](https://www.google.com/search?q=https://crates.io/crates/moma)
[](https://www.google.com/search?q=https://docs.rs/moma)
[](https://opensource.org/licenses/MIT)
[](https://www.google.com/search?q=https://github.com/neil-crago/moma/actions)

**MOMA** is a Rust framework for exploring complex systems including cosmology, number theory, algorithmic data analysis and cryptography through the lens of **Moving Origin Modular Arithmetic**.

It should interest anyone who is interested in a novel, relational framework for analyzing the properties of complex systems.

-----

## The Core Idea: A Barycenter for Numbers

The inspiration for this crate comes from the concept of a **barycenter** in astrophysics. Just as the Earth and Moon orbit a common center of mass that is not the exact center of the Earth, MOMA treats modular arithmetic as a system where the "zero point" or "origin" is not fixed.

This origin shifts dynamically based on a contextual value—typically a prime number `p`—and the chosen `OriginStrategy`. This provides a novel relational framework for analyzing the properties of integers, much like understanding the true orbital mechanics of our solar system.

> The original inspiration came from this NASA article: [What Is a Barycenter?](https://spaceplace.nasa.gov/barycenter/en/)

-----

## Core Concepts

The MOMA framework is built on a few simple but powerful concepts:

  * **`MomaRing`**: The primary object for all calculations. A ring is defined by a `modulus` and a chosen `OriginStrategy`.

  * **`OriginStrategy`**: A trait that defines *how* the origin moves. This makes the framework highly extensible.

  * **Analysis Tools**: A suite of modular tools for deeper analysis:

      * **`MassField`**: Maps primes to the "composite mass" in the gap immediately following them.
      * **`OriginDrift`**: Measures the volatility or "drift" of MOMA signatures for any given strategy.
      * **`CompositeInfluence`**: Models the "gravitational" influence of nearby composite numbers.
      * **`GoldbachProjector`**: An efficient tool for finding prime pairs for even numbers.
      * **`Entropy`**: A generic calculator for the Shannon entropy of a sequence.

-----

## Features

  * **Flexible Core**: A powerful and extensible system based on the `MomaRing` and `OriginStrategy` trait.
  * **Advanced Analysis Tools**: A suite of high-level structs for statistical and number-theoretic analysis.
  * **Cryptographic Primitives**: Demonstrates how MOMA can be used to build components like a Key Derivation Function.
  * **Prime Number Utilities**: A helper `primes` module for primality testing and prime generation.
  * **Pure Rust**: Built with safe, idiomatic Rust.

-----

## Installation

Add MOMA to your `Cargo.toml`:

```toml
[dependencies]
moma = "0.1.9" # Replace with the latest version
```

-----

## Quick Start

The easiest way to get started is to create a `MomaRing` and calculate the "signature" of a prime.

```rust
use moma::core::core::MomaRing;
use moma::primes::primes;
use moma::strategy::strategy;

// 1. Create a MOMA ring with modulus 37 and the PrimeGap strategy.
let ring = MomaRing::new(37, strategy::PrimeGap);

// 2. Let's analyze the prime p = 29.
let p = 29;

// 3. Calculate the MOMA signature.
//    The origin for p=29 is (29 - 23) = 6.
//    The value is p + p_prev = 29 + 23 = 52.
//    The residue is (52 + 6) % 37 = 21.
let signature = ring.signature(p);

println!("For p={}, the MOMA signature is: {}", p, signature);
assert_eq!(signature, 21);
```

-----

## Exploring Further

MOMA is more than just a simple calculator; it's a toolkit for exploration.

### Example 1: Finding Goldbach Pairs

Use the `GoldbachProjector` to efficiently find prime pairs for an even number.

```rust
use moma::goldbach::GoldbachProjector;

// Create a projector with a database of primes up to 1000.
let projector = GoldbachProjector::new(1000);

// Find all pairs for the number 96.
let pairs = projector.project(96);

println!("Found {} Goldbach pairs for 96.", pairs.len());
// Example pair: (7, 89)
assert!(pairs.contains(&(7, 89)));
```

### Example 2: Measuring Strategy Volatility

Use `OriginDrift` to compare the stability of different `OriginStrategy` implementations.

```rust
use moma::origin_drift::OriginDrift;
use moma::strategy::strategy;
use moma::primes::primes;

// Create a drift analyzer for the PrimeGap strategy.
let mut drift_analyzer = OriginDrift::new(100, strategy::PrimeGap);

// Feed it a sequence of primes.
let mut p = 3;
for _ in 0..10 {
    drift_analyzer.next(p);
    p = primes::next_prime(p);
}

// A higher drift magnitude means the strategy is more volatile.
println!(
    "Drift magnitude for PrimeGap strategy: {:.2}",
    drift_analyzer.drift_magnitude()
);
```

-----

## Contributing

Contributions are welcome\! If you have an idea for a new `OriginStrategy`, an analysis tool, or find a bug, please feel free to open an issue or submit a pull request.

## License

This project is licensed under either of:

  * Apache License, Version 2.0, ([LICENSE-APACHE](https://www.google.com/search?q=LICENSE-APACHE) or [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0))
  * MIT license ([LICENSE-MIT](https://www.google.com/search?q=LICENSE-MIT) or [http://opensource.org/licenses/MIT](http://opensource.org/licenses/MIT))

at your option.