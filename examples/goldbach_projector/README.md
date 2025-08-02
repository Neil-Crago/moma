# Example: Goldbach Projector

This example demonstrates how to use the `GoldbachProjector` struct from the `moma` crate to efficiently find Goldbach pairs for even numbers.

It first creates a `GoldbachProjector`, which pre-computes a set of primes up to a certain limit for fast lookups. It then iterates through a list of even numbers and, for each one, finds all unique pairs of prime numbers that sum up to it.

## Concepts Demonstrated

* Creating an efficient `GoldbachProjector` with a prime database.
* Finding all Goldbach pairs for multiple even numbers.
* The benefits of pre-computation for performance in number theory problems.

## How to Run

From the root of the `moma` crate repository, run the following command in your terminal:

```sh
cargo run --example goldbach_projector