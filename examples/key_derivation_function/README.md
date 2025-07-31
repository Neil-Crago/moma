# `Key_Derivation_Function`

#### Author: Neil Crago

## Description:

This is a practical implementation of a Key Derivation Function (KDF) using the moma_crate.

This code defines a MomaKdf struct that takes a password, a salt, and a work factor (iterations) to produce a secure key. It's designed to be a clear demonstration of how MOMA's unique properties can be used as a cryptographic primitive.

## MOMA KDF Implementation

First we add a hashing library to the Cargo.toml to handle the initial seeding. sha2 is a standard choice.

## How It Works

 * Seeding: The KDF takes the password and salt and hashes them to deterministically generate the two starting parameters for MOMA: the modulus and a starting_prime. Using a salt prevents attackers from using pre-computed "rainbow tables" to crack passwords.
 * Chained Iteration: The core of the KDF is a loop that runs for a specified number of iterations. In each step, it calculates a MOMA signature. Crucially, the input for the next step (current_prime) is derived from the output of the current step. This creates a sequential dependency chain, forcing any attacker to perform the work serially, which is the main principle behind modern KDFs.
 * Finalization: After all iterations, the long sequence of intermediate results is hashed one last time to produce the clean, fixed-size output key.

## ⚠️ Security Disclaimer

This code is a demonstration of the concept and is perfect for showcasing your MOMA framework. However, for use in a real-world, production security system, it would need to undergo extensive public peer review and cryptanalysis. Real-world KDFs contain many subtle protections against side-channel attacks and other advanced threats.
For production systems, you should always use established, audited libraries like argon2 or ring::pbkdf2.


## Major Revisions

  - First Commit => July 2025

## Performance and Accuracy

Not focused on Perf yet.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the
work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any
additional terms or conditions.
