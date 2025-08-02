# MOMA: Moving Origin Modular Arithmetic

[](https://www.google.com/search?q=https://crates.io/crates/moma)
[](https://www.google.com/search?q=https://docs.rs/moma)
[](https://opensource.org/licenses/MIT)
[](https://www.google.com/search?q=https://github.com/neil-crago/moma/actions)

**MOMA** is a Rust framework for exploring number theory and cryptography through the lens of **Moving Origin Modular Arithmetic**.

The central idea is that the "zero point" or "origin" for modular arithmetic is not fixed. Instead, it shifts dynamically based on a contextual value—typically a prime number `p`. This provides a novel relational framework for analyzing the properties of integers.

The crate is designed for researchers and developers working in number theory, cryptography, and algorithmic data analysis.

-----

## Core Concepts

The MOMA framework is built on two simple but powerful concepts:

  * **`MomaRing`**: The primary object for all calculations. A ring is defined by a `modulus` and a chosen `OriginStrategy`.

  * **`OriginStrategy`**: A trait that defines *how* the origin moves. This makes the framework highly extensible. You can define custom logic for the origin's movement based on any property of the contextual prime.

The core calculation is the MOMA residue, defined as:

$residue = (value + origin) \\pmod{modulus}$

### Built-in Strategies

The crate includes several pre-built strategies to get you started:

  * `strategy::PrimeGap`: The origin is the gap between a prime and its predecessor ($p - p\_{prev}$).
  * `strategy::CompositeMass`: The origin is the sum of prime factors of all composite numbers in the gap between a prime and its successor.
  * `strategy::Fixed`: A simple strategy where the origin is a fixed constant.

-----

## Features

  * **Flexible Core**: A powerful and extensible system based on the `MomaRing` and `OriginStrategy` trait.
  * **Analysis Tools**: Includes high-level structs like `PrimeGapField` for statistical analysis of prime gaps.
  * **Cryptographic Primitives**: Demonstrates how MOMA can be used to build cryptographic components like the `MomaKdf` (Key Derivation Function).
  * **Prime Number Utilities**: A helper `primes` module for efficient primality testing and prime generation.
  * **Pure Rust**: Built with safe, idiomatic Rust.

-----

## Installation

Add MOMA to your `Cargo.toml`:

```toml
[dependencies]
moma = "0.1.7"
```

-----

## Quick Start

The easiest way to get started is to create a `MomaRing` and calculate the "signature" of a prime. The signature is the MOMA residue of the sum of a prime and its predecessor.

```rust
use moma::core::MomaRing;
use moma::primes;
use moma::strategy;

fn main() {
    // 1. Create a MOMA ring with modulus 37.
    //    The origin will be the gap between a prime and its predecessor.
    let ring = MomaRing::new(37, strategy::PrimeGap);

    // 2. Let's analyze the prime p = 29.
    let p = 29;
    
    // The previous prime is 23.
    let p_prev = primes::prev_prime(p);
    
    // The origin for p=29 is the prime gap: 29 - 23 = 6.
    let origin = strategy::PrimeGap.calculate_origin(p);
    assert_eq!(origin, 6);

    // 3. Calculate the MOMA signature.
    //    The value to test is p + p_prev = 29 + 23 = 52.
    //    The residue is (52 + origin) % 37 = (52 + 6) % 37 = 21.
    let signature = ring.signature(p);

    println!("For p={}, the MOMA signature is: {}", p, signature);
    assert_eq!(signature, 21);
}
```

-----

## Exploring Further

MOMA is more than just a simple calculator; it's a toolkit for exploration.

### Example 1: Prime Gap Analysis

You can use the `PrimeGapField` analyser to find statistical outliers in the sequence of primes.

```rust
// main.rs from the `moma-analysis` example
use moma::analysis::PrimeGapField;
use moma::primes;

// Generate primes up to 100
let primes: Vec<u64> = (1..=100).filter(|&n| primes::is_prime(n)).collect();

// Create a field to analyze gaps modulo 6
let field = PrimeGapField::new(&primes, 6);

// Find gaps that are significantly larger or smaller than their local average
let outliers = field.filter_by_bary_offset(3.0);

println!("Found {} outlier gaps:", outliers.len());
for gap in outliers {
    println!(
        "  Gap [{: >2}, {: >2}] of size {} deviates by {:.2}",
        gap.start_prime, gap.end_prime, gap.size, gap.bary_offset
    );
}
```

### Example 2: Cryptographic Key Derivation

MOMA can be used as a building block for cryptographic functions. The included `MomaKdf` example shows how to use it as a password-stretching function.

```rust
// main.rs from the `moma-kdf` example
// [Note: This requires the sha2 and hex crates]

// Create a KDF instance
let kdf = MomaKdf::new(
    "my-secret-password".as_bytes(),
    "random-public-salt".as_bytes(),
    1_000, // Iterations (work factor)
    32,    // Output key length in bytes
);

// Derive a strong 256-bit key
let key = kdf.derive_key();

println!("Derived AES-256 Key: {}", hex::encode(&key));
```

## Contributing

Contributions are welcome\! If you have an idea for a new `OriginStrategy`, an analysis tool, or find a bug, please feel free to open an issue or submit a pull request.

## License

This project is licensed under either of:

  * Apache License, Version 2.0, ([LICENSE-APACHE](https://www.google.com/search?q=LICENSE-APACHE) or [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0))
  * MIT license ([LICENSE-MIT](https://www.google.com/search?q=LICENSE-MIT) or [http://opensource.org/licenses/MIT](http://opensource.org/licenses/MIT))

at your option. 