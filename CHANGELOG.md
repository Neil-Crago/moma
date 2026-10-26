# Changelog

## Unreleased

- `sperm_beat`: of the 17 cells (13%) with a significant slow head oscillation, 15 show the same period (+-20%) 2 um further down the tail, as a whole-cell motion should. Slow period about 4.4 Hz (2.8x the beat period) in both capacitating (9 cells) and non-capacitating (8 cells); too few cells to compare conditions.
- `sperm_beat` adds AR(1) red-noise surrogate tests (300 per cell). Positive control: the tail beat band is significant (p<0.05) in 94% of cells. The head's slow-band peak is significant in only 13% of cells (5% expected from noise), so the slow head-path period is mostly not distinguishable from noise at these record lengths; a minority of cells show a real slow oscillation.
- `sperm_beat` now reports capacitating vs non-capacitating medians and a head-path period. Volume rate is 180 per second (90 Hz piezo, one volume per up or down sweep: 5000-8000 fps / 29-45 planes; from the 3D-SpermFlagella paper's acquisition section). Median beat: capacitating 10.1 Hz, non-capacitating 16.1 Hz; tip amplitude about 2.8 um in both. Cells are not independent (10 donors), so this is descriptive only.
- helix::shared_period_3d: one period pooled across simultaneous tracks (one vote per track), which stops a track dominated by a harmonic from setting the period.
- xamples/sperm_beat.rs: real-data check on 135 human sperm flagellum centrelines (Zenodo 10.5281/zenodo.15299846). Positions agreeing with the pooled beat period: 94% (vs about 75% of cells when each position chose its own); tip amplitude over 1.5x head amplitude in 96% of cells.

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
