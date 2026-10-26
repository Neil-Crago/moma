# Changelog

## 1.0.0

### Breaking: removed
The original ring and prime modules, which did not fit the moving-axis direction and had no real-data validation: `analysis`, `biosig`, `codon`, `composite_field`, `core` (`MomaRing`, `OriginStrategy`), `entropy`, `goldbach`, `influence`, `massfield`, `mutation`, `origin_drift`, `primes`, `resonance`, `score` and `strategy`, with their example crates (`prime_gaps`, `goldbach_projector`, `origin_drift`, `mass_field`). Use 0.5.x if you need them. The `rand` dependency is gone; the crate has no dependencies.

### Added
- `consensus`: `circular_mean` and `ring_consensus` (with concentration), `weighted_barycentre`, `geometric_median` (Weiszfeld with the Vardi-Zhang correction), and `consensus_axis` / `consensus_axis_2d`, which merge several axis estimates and weight each estimator by agreement with the rest. `ring_barycentre` and `plane_barycentre` are unchanged.
- `examples/consensus_bench.rs`: consensus of six axis estimators on the synthetic elliptical benchmark. It matches the best estimator for steady and irregular rotation (plain mean 4 to 14 times worse), but not for reversing rotation, where most estimators share a bias and the consensus follows them.
- `helix::shared_period_3d`: one period pooled across simultaneous tracks (one vote per track), which stops a track dominated by a harmonic from setting the period.
- Against AR(1) red-noise surrogates (300 per cell) the tail beat band is significant in 94% of cells (positive control) but the slow head band in only 13% (5% expected by chance). Of those 17 cells, 15 show the same period 2 um further down the tail, about 4.4 Hz (2.8x the beat period).
- `examples/sperm_beat.rs`: real-data check on 135 human sperm flagellum centrelines (Zenodo 10.5281/zenodo.15299846). Positions agreeing with the pooled beat period: 94% (about 75% of cells when each position chose its own); tip amplitude over 1.5x head amplitude in 96% of cells. Volume rate 180 per second (90 Hz piezo, one volume per sweep, from the dataset paper). Median beat: capacitating 10.1 Hz, non-capacitating 16.1 Hz; tip amplitude about 2.8 um in both. Cells share 10 donors, so this is descriptive only.

## 0.5.0

Refocus on moving-axis analysis of spirals and helices.

### Added
- `spiral` (2D): `moving_axis`, `adaptive_axis`, `refine_axis` (joint axis and
  rotation fit with self-reported residual), `to_polar`, `unwrap_phase`,
  `radius_spread`, `winding_number`, `quantise_phase`, `quantised_steps`,
  `backward_steps`, `fill_gaps`, `moving_median_axis`, `reject_outliers`, `parse_xy`.
- `helix` (3D): `moving_axis_3d`, `adaptive_axis_3d`, `tangents`, `to_cylindrical`
  (parallel-transported frame), `repair_3d`, `fill_gaps_3d`,
  `moving_median_axis_3d`, `parse_xyz`.
- `consensus`: `ring_barycentre` (weighted circular mean on `Z_m`) and `plane_barycentre`.
- Fractional-window axis: `moving_axis_period` / `moving_axis_period_3d` (mean over exactly
  `period` samples, fractional end weights). On the alpha helices of the flagellar
  filament PDB entry 1UCU (`examples/pdb_helices.rs`) it recovers radius 2.20 A,
  100.9 deg and 1.49 A per residue (ideal 2.3, 100, 1.5); an integer window of 5 gave
  radius 2.83 A and rise 1.71 A.
- Reversing and irregular rotation: `longest_period` / `longest_period_3d` (robust to
  harmonics; now used by `auto_axis`), `moving_midrange_axis` (+ `_3d`; centre of a
  window's coordinate range, not biased by dwelling), `coverage_axis` (+ `_3d`; window
  grows until one full turn is seen). On the synthetic reversing benchmark axis error
  fell from 1.6 to about 0.5.
- Period-aware axis: `dominant_period`, `window_for_period`, `auto_axis`, `total_turns`,
  `total_turns_smoothed`
  (spiral) and `dominant_period_3d`, `auto_axis_3d` (helix). On 20 real pea tracks this
  cut the window-to-window spread of radius and turn count from about 85% to 7-40%.
- Examples: `elliptical_bench` (synthetic elliptical, irregular and reversing rotation with ground truth), `snaking_spiral`, `helix_stress`, `flagellum`, and `real_track`
  (consistency check on real data; see `data/README.md` for the source).

### Fixed
- `MomaRing::residue` no longer overflows for values near `u64::MAX`; operands are
  reduced before adding. Results can differ from 0.4 for such inputs.

### Removed
- The `cosmology`, `bioinformatics` and `key_derivation_function` example crates.
  The library modules they used (`biosig`, `codon`, `mutation`, ...) remain.
