# MOMA: Moving Origin Modular Arithmetic

![MOMA Logo](docs/moma_trim.png)

[![Crates.io](https://img.shields.io/crates/v/moma.svg?style=flat-square)](https://crates.io/crates/moma)
[![Docs.rs](https://img.shields.io/docsrs/moma?style=flat-square)](https://docs.rs/moma)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/neil-crago/moma/actions/workflows/rust.yml/badge.svg)](https://github.com/neil-crago/moma/actions/workflows/rust.yml)

**MOMA** is a Rust toolkit for analysing things relative to an origin that *moves*. Its main use is **spirals and helices around a wandering axis**: the kind of path you get from a flagellum, a tendril, a vortex track or a coiled filament. Seen in fixed coordinates these look irregular; relative to a moving axis they are a steady rotation you can measure.

The crate also keeps its original ring-based tools (modular arithmetic with a movable origin, prime-based signatures, a little bioinformatics), described under [Legacy and experimental](#legacy-and-experimental).

> **Status:** the spiral and helix tools are tested on synthetic data with known ground truth. On real data (20 pea-shoot circumnutation tracks, `examples/real_track.rs`) a fixed or adaptive window of 100–400 samples **gave inconsistent results**, because the motion is an elongated, irregular ellipse with a period of 350–1900 samples, far longer than those windows. Setting the window to the *measured* dominant period (`dominant_period_3d`, `auto_axis_3d`) fixes most of this: across windows of 0.8–1.25 periods the median radius varies by 7–31% and the total turns by 10–40% (versus about 85% before). That is consistency, not accuracy. On the synthetic elliptical benchmark (`examples/elliptical_bench.rs`, known truth) the period-matched window recovers axis (error about 0.05) and total turns (within 2%) for steady rotation, and total turns within 3% for a wobbling rate (axis error rises to about 0.6), at both noise levels tested, when turns are counted with `total_turns_smoothed`. When the rotation sense **reverses**, period detection locks onto a harmonic (400 instead of 800) and the result degrades (+29% turns at low noise, +315% at high noise). Raw `total_turns` counts jitter and is badly inflated by noise; always smooth. For the real tracks there is no ground truth, signed turns are meaningless when the rotation sense reverses (use `total_turns`), and total turns includes jitter. The tools remain best suited to roughly steady helices (flagella, tendrils).

---

## The Core Idea: A Barycenter for Paths

Just as the Earth and Moon orbit a common center of mass that is not the center of the Earth, a spiral rotates around an axis that is itself moving. Subtract the moving axis and what is left is radius, phase and progress along the axis, which are the quantities you actually care about. Phase lives on a circle, so it also maps naturally onto a ring `Z_m`, which is where the modular arithmetic comes in.

> The original inspiration came from this NASA article: [What Is a Barycenter?](https://spaceplace.nasa.gov/barycenter/en/)

## Core Concepts

* **`spiral` / `helix`**: moving-axis analysis in 2D and 3D.
* **`consensus`**: weighted barycentres, using a circular mean on `Z_m` (so 1 and 9 on `Z_10` average to 0, not 5).
* **`MomaRing` and `OriginStrategy`**: the underlying ring with a pluggable origin.

## Features

* **Axis estimation**: fixed window, adaptive window, and a joint axis-and-rotation fit that reports its own residual.
* **Dirty-data handling**: gap filling, outlier rejection and a median axis, with one-call repair for 3D tracks.
* **Ring signatures**: winding number, quantised phase, and a backward-step (slip) count.
* **Your own data**: `parse_xy` and `parse_xyz` read CSV-like text.
* **Pure Rust**, no heavy dependencies.

---

## Installation

Add MOMA to your `Cargo.toml`:

```toml
[dependencies]
moma = "0.5"
```

or simply run:

```bash
cargo add moma
```

## Quick start

```rust
use moma::helix::{adaptive_axis_3d, to_cylindrical, repair_3d};
use moma::spiral::{unwrap_phase, winding_number};

// `points` is your 3D track, e.g. from `moma::helix::parse_xyz(&text)`.
let points = repair_3d(&points, 31, 4.0);          // fill dropouts, drop spikes
let axis = adaptive_axis_3d(&points, 31, 3);        // moving axis
let (polar, along) = to_cylindrical(&points, &axis); // radius, phase, progress
let turns = winding_number(&unwrap_phase(&polar)).abs();
```

## Results on synthetic data

Run with `cargo run --example <name>`; each compares methods against the true axis.

* **`flagellum`**: a bending 3D helix whose wave speeds up 2.5x. The adaptive axis cuts axis error about 4x against the best fixed window (0.03 against 0.10) and recovers 52.05 of 52.24 true turns.
* **`helix_stress`**: noise, dropouts and spikes. Repairing a damaged track restores axis error from 1.9 to 0.03 and the winding count from 30 to 42.6 turns (true 42.6). The joint fit helps at low noise and gives no gain at high noise.
* **`snaking_spiral`**: 2D noise sweep and window-length comparison.

Caveats: the adaptive window only helps when rotation speed varies (it was worse than a fixed window on constant-speed data); the window must be roughly one turn; gap filling restores the axis, not the lost rotation; heavy noise degrades every method.

## Spirals around a moving axis

`moma::spiral` (2D) and `moma::helix` (3D) re-express a path relative to a moving origin (the local barycentre), with the 3D phase measured in a parallel-transported frame so it does not twist on its own.

- **Axis estimation:** `moving_axis` / `moving_axis_3d` (fixed window), `adaptive_axis` / `adaptive_axis_3d` (window sized to the local rotation period), `refine_axis` (2D joint fit).
- **Dirty data:** `fill_gaps`, `moving_median_axis`, `reject_outliers`, and `helix::repair_3d`.
- **Ring link:** `quantise_phase` maps phase to `Z_m`; `quantised_steps`, `backward_steps` and `winding_number` give slip and winding signatures.
- **Your data:** `parse_xy` / `helix::parse_xyz`.

### Barycentric consensus

`moma::consensus` computes weighted barycentres on the ring and in the plane (as an `OriginShift`).

## Legacy and experimental

The crate's original modules remain, but are more speculative and are not the focus: prime-based origin strategies (`PrimeGap`, `CompositeMass`), `ResonanceFinder`, `MassField`, `OriginDrift`, `GoldbachProjector`, and the bioinformatics helpers `BioSigAnalyzer`, `CodonTable` and `Mutation`.

```rust
use moma::resonance::ResonanceFinder;
use moma::strategy;
use moma::primes;

let finder = ResonanceFinder::new(100, strategy::CompositeMass, primes::prime_factor_mass);
let resonances = finder.find_in_range(1, 500);
println!("Found {} resonance events.", resonances.len());
```

Example crates for these (`prime_gaps`, `goldbach_projector`, `origin_drift`, `mass_field`) are in `examples/`.

## Author

Neil Crago

## Contributing

Contributions are welcome. If you have an idea for a new `OriginStrategy`, an analysis tool, or find a bug, please feel free to open an issue or submit a pull request.

## License

This project is licensed under either of:

* Apache License, Version 2.0 ([LICENSE-APACHE](https://www.apache.org/licenses/LICENSE-2.0))
* MIT license ([LICENSE-MIT](https://opensource.org/licenses/MIT))

at your option.

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
