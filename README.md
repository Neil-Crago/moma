# MOMA: Moving Origin Modular Arithmetic

![MOMA Logo](docs/moma_trim.png)

[![Crates.io](https://img.shields.io/crates/v/moma.svg?style=flat-square)](https://crates.io/crates/moma)
[![Docs.rs](https://img.shields.io/docsrs/moma?style=flat-square)](https://docs.rs/moma)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue?style=flat-square)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/neil-crago/moma/actions/workflows/rust.yml/badge.svg)](https://github.com/neil-crago/moma/actions/workflows/rust.yml)

**MOMA** is a Rust toolkit for analysing things relative to an origin that *moves*. Its main use is **spirals and helices around a wandering axis**: the kind of path you get from a flagellum, a tendril, a vortex track or a coiled filament. Seen in fixed coordinates these look irregular; relative to a moving axis they are a steady rotation you can measure.

The crate's name is historical: it began as modular arithmetic with a movable origin. That part, with the prime-based and bioinformatics modules, was removed in 1.0 (see the [CHANGELOG](CHANGELOG.md)); the ring survives as the place phase is quantised to, and as the setting for circular consensus.

> **Real-data validation:** see [Results on real data](#results-on-real-data) for the human sperm flagellum, flagellin PDB and pea-shoot checks, including what failed.

> **Real-structure check:** the alpha helices in a real flagellar filament structure (PDB 1UCU, `examples/pdb_helices.rs`) have known geometry (radius 2.3 A, 100 degrees and 1.5 A per residue). With a fractional window equal to the measured period (`moving_axis_period_3d`) the tools recover 2.20 A, 100.9 degrees and 1.49 A across five helices, and the measured period is 3.4 to 3.6 residues. An integer window overestimated radius and rise by 20%, because 3.6 is not an integer. This is a static structure, so it checks the geometry, not tracking over time.

> **Status:** the spiral and helix tools are tested on synthetic data with known ground truth. On real data (20 pea-shoot circumnutation tracks, `examples/real_track.rs`) a fixed or adaptive window of 100–400 samples **gave inconsistent results**, because the motion is an elongated, irregular ellipse with a period of 350–1900 samples, far longer than those windows. Setting the window to the *measured* dominant period (`dominant_period_3d`, `auto_axis_3d`) fixes most of this: across windows of 0.8–1.25 periods the median radius varies by 7–31% and the total turns by 10–40% (versus about 85% before). That is consistency, not accuracy. On the synthetic elliptical benchmark (`examples/elliptical_bench.rs`, known truth) the period-matched window recovers axis (error about 0.05) and total turns (within 2%) for steady rotation, and total turns within 3% for a wobbling rate (axis error rises to about 0.6), at both noise levels tested, when turns are counted with `total_turns_smoothed`. When the rotation sense **reverses**, the strongest spectral peak can be a harmonic (400 instead of 800), and even the right window fails because the orbit may not complete a loop inside it. `auto_axis` therefore uses the *longest significant* period (`longest_period`), and `coverage_axis` / `coverage_axis_3d` grow the window until one full turn is seen and take the midpoint of the coordinate range, which is not biased by dwelling. On the reversing benchmark this cut axis error from 1.6 to about 0.5 and kept total turns within 6% (noise 0.05) and 3% (noise 0.3); the plain mean axis was off by 22–99%. Raw `total_turns` counts jitter and is badly inflated by noise; always smooth. For the real tracks there is no ground truth, signed turns are meaningless when the rotation sense reverses (use `total_turns`), and total turns includes jitter. The tools remain best suited to roughly steady helices (flagella, tendrils).

---

## The Core Idea: A Barycenter for Paths

Just as the Earth and Moon orbit a common center of mass that is not the center of the Earth, a spiral rotates around an axis that is itself moving. Subtract the moving axis and what is left is radius, phase and progress along the axis, which are the quantities you actually care about. Phase lives on a circle, so it also maps naturally onto a ring `Z_m`, which is where the modular arithmetic comes in.

> The original inspiration came from this NASA article: [What Is a Barycenter?](https://spaceplace.nasa.gov/barycenter/en/)

## Core Concepts

* **`spiral` / `helix`**: moving-axis analysis in 2D and 3D.
* **`consensus`**: one answer from several competing positions or axis estimates: circular means on `Z_m` (so 1 and 9 on `Z_10` average to 0, not 5), robust geometric-median barycentres, and a path-level consensus that outvotes a failing estimator.

## Features

* **Axis estimation**: fixed window, adaptive window, period-matched window, midrange and coverage axes, and a joint axis-and-rotation fit that reports its own residual.
* **Period finding**: dominant and longest-significant periods, and one period pooled across simultaneous tracks.
* **Dirty-data handling**: gap filling, outlier rejection and a median axis, with one-call repair for 3D tracks.
* **Ring signatures**: winding number, quantised phase, and a backward-step (slip) count.
* **Consensus**: circular mean with concentration, geometric median, and `consensus_axis` with per-estimator weights and per-sample spread.
* **Your own data**: `parse_xy` and `parse_xyz` read CSV-like text.
* **Pure Rust**, no dependencies.

---

## Installation

Add MOMA to your `Cargo.toml`:

```toml
[dependencies]
moma = "1.0"
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

## Results on real data

Datasets and licences are listed in [`data/README.md`](data/README.md). Large files are not committed. Real data has no ground truth for the axis, so each check uses either known geometry or a property the data must satisfy that the tools are not told.

### Human sperm flagellum centrelines (3D+t)

`examples/sperm_beat.rs` on 135 cells from 3D-SpermFlagella (Zenodo 10.5281/zenodo.15299846, CC BY 4.0; 10 donors). Each point on a tail is a 3D track oscillating around a drifting axis. Volume rate 180 per second (90 Hz piezo, one volume per sweep, from the dataset paper), but the shortest records (65 volumes) do not match the stated 1 to 3.5 s, so Hz values are provisional.

Properties the data must satisfy, not given to the tools:

* **One beat frequency along a tail.** Each position choosing its own period agreed in only about 75% of cells; most of the rest locked on a harmonic that outweighed the fundamental. Pooling the positions with `shared_period_3d` (one vote per track) gives 94% of positions, and 94% of cells with at least 75% of positions, agreeing with the shared period. The remaining positions are unrelated to it, not harmonics, so they are probably noisy or short records.
* **Amplitude grows head to tip.** Tip amplitude exceeds 1.5x the head-side amplitude in 96% of cells.

Descriptive results (not independent: cells share donors, and no significance test was run):

| Condition | Cells | Beat | Tip amplitude |
|---|---|---|---|
| Capacitating | 86 | 10.1 Hz | 2.85 um |
| Non-capacitating | 49 | 16.1 Hz | 2.80 um |

A slow oscillation of the head path was tested against AR(1) red-noise surrogates (300 per cell). The positive control, the tail beat band, is significant in 94% of cells. The slow head band is significant in only 13% (17 cells; 5% expected by chance), so for most cells the slow period is just the largest noise peak. Of those 17, 15 show the same period 2 um further down the tail, about 4.4 Hz (2.8x the beat period), as a whole-cell motion should. The records are too short to say more; "not detectable" is not "absent", and the test removes only linear drift so it is conservative. Do not read the 2.8x ratio or a capacitation difference into this.

### Flagellin filament structure (PDB 1UCU)

`examples/pdb_helices.rs` treats each alpha helix as a 3D helix. Ideal geometry is radius 2.3 A, 100 degrees and 1.5 A per residue; with a window equal to the measured (fractional) period, `moving_axis_period_3d` recovers 2.20 A, 100.9 degrees and 1.49 A across five helices. An integer window overestimated radius and rise by 20%. It is a static structure, so it tests geometry, not tracking over time.

### Pea-shoot circumnutation

`examples/real_track.rs` on 20 tracks (Zenodo 10.5281/zenodo.22015401): fixed windows of 100 to 400 samples were inconsistent, because the motion is an irregular ellipse with a period of 340 to 1900 samples. Period-matched windows cut the window-to-window spread to 7 to 31% (radius) and 10 to 40% (turns). That is consistency only.

### Rejected datasets

A bacterial track set was unusable, and the flagellin sequence tables (Nature Communications 2024, figshare) contain no coordinates; see `data/README.md`.

## Spirals around a moving axis

`moma::spiral` (2D) and `moma::helix` (3D) re-express a path relative to a moving origin (the local barycentre), with the 3D phase measured in a parallel-transported frame so it does not twist on its own.

- **Axis estimation:** `moving_axis` / `moving_axis_3d` (fixed window), `adaptive_axis` / `adaptive_axis_3d` (window sized to the local rotation period), `refine_axis` (2D joint fit).
- **Period finding:** `dominant_period` / `longest_period` (and `_3d`), `auto_axis`, `coverage_axis`, `moving_axis_period`, and `helix::shared_period_3d` (one period pooled across simultaneous tracks).
- **Dirty data:** `fill_gaps`, `moving_median_axis`, `reject_outliers`, and `helix::repair_3d`.
- **Ring link:** `quantise_phase` maps phase to `Z_m`; `quantised_steps`, `backward_steps` and `winding_number` give slip and winding signatures.
- **Your data:** `parse_xy` / `helix::parse_xyz`.

### Barycentric consensus

`moma::consensus` merges several answers into one.

- **Ring:** `circular_mean` and `ring_consensus` return the weighted circular mean with a *concentration* (0 to 1) saying how strongly the votes agree; `ring_barycentre` returns just the residue.
- **Space:** `weighted_barycentre` is the mean; `geometric_median` is the robust version, barely moved by a minority of wild points. `plane_barycentre` gives an `OriginShift`.
- **Paths:** `consensus_axis` (and `consensus_axis_2d`) merges several estimates of one moving axis. Each sample is the geometric median of the estimates, and each estimator is weighted by its prior divided by its rms distance to the consensus, so an estimator that is wrong throughout is outvoted rather than averaged in. It returns the weights (which estimators were trusted) and a per-sample spread (where the axis is poorly determined).

`examples/consensus_bench.rs` runs six axis estimators on the synthetic elliptical benchmark (known axis) and merges them. Axis rms error:

| Case | Best single | Median single | Plain mean of all | Consensus |
|---|---|---|---|---|
| steady | 0.047 | 0.17 to 0.21 | 0.66 | **0.047** |
| irregular | 0.29 to 0.30 | 0.46 | 0.65 | **0.29 to 0.30** |
| reversing | 0.47 to 0.48 | 1.62 | 1.43 | 1.62 |

For steady and irregular rotation the consensus finds the best estimator without being told which it is, and the plain mean is 4 to 14 times worse. For **reversing** rotation it does not: four of the six estimators share the same bias (they miss that the orbit does not complete a loop), so they form the majority and the two unbiased ones (midrange, coverage) are outvoted. Consensus protects against a *minority* of estimators failing; it cannot rescue a case where most of them fail alike. If you know which estimators are trustworthy, pass `prior` weights.

**Real data: the travelling wave on a sperm tail** (`examples/consensus_wave.rs`, the 135 cells above). The beat phase at consecutive points along a tail advances by a roughly constant step, so those steps are votes that `circular_mean` can combine. Beat phase was measured at 48 positions between 5% and 60% of the tail, weighted by amplitude; nothing about the wave is given to the tools.

* **Negative control:** with the positions in their true order the phase steps have concentration of about 1.00 (median); shuffling the same positions drops it to about 0.49. Real order beats its own shuffle in all 135 cells.
* **Leave-one-out prediction** of a held-out position's phase from its neighbour (mean absolute error in radians): copying the neighbour 0.071; the single previous step 0.0078; consensus of *all* other steps 0.0147; consensus of the *four nearest* steps 0.0056.
* **What this says:** a global consensus step is *worse* than the single previous step, because the wave speed varies along the tail, so one constant step is the wrong model. A local consensus beats the single step by about 28%, so consensus helps when the voters estimate the same quantity and noise is the problem. Choose the voters to match the model.
* **Calibration:** concentration is almost 1 everywhere (0.9998 to 1.0000 across the thirds), but the ordering is right: the lowest-concentration third has the largest error (0.021 against 0.009 for the global consensus). Treat concentration as a ranking signal, not an absolute scale.
* **Limits:** the published centrelines are smooth and interpolated, which makes this an easy case, and it checks a consensus of one kind (circular) on one dataset. `consensus_axis` has not been tested on real data.
## Author

Neil Crago

## Contributing

Contributions are welcome. If you have an idea for an analysis tool or a dataset to validate against, or find a bug, please feel free to open an issue or submit a pull request.

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
