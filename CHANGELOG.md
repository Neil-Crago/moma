# Changelog

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
