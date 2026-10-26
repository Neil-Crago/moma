# Test data

`pea_exp3box8AT1.csv` is one file (`3DtrajectorySmooth_P1_6 EXP3box8AT1.csv`, renamed) from:

> Francisco de Carvalho Oliveira, T.; Niemeyer Reissig, G.; Ardini, D.; Nava, D. E.;
> Maia Souza, G.; Castiello, U.; Guerra, S. *Raw kinematical recordings of Pisum
> sativum during circumnutation under support and no-support conditions.* Zenodo
> (2026). DOI [10.5281/zenodo.22015401](https://doi.org/10.5281/zenodo.22015401).
> Licensed CC BY 4.0.

Plain `x,y,z` rows, no header; y is the growth direction. The other 19 recordings
are at the DOI above. This directory is not included in the published crate.

## Datasets tried and rejected

**Magnetotactic bacteria 3D tracks** (Bente et al., eLife 2020, Dryad
[10.5061/dryad.r2nd550](https://doi.org/10.5061/dryad.r2nd550), CC0, 456 MB). Not
usable as a helix test: coordinates are integer pixels of 0.219 µm (xy) and
0.111 µm (z), the speeds in the tracks are about 100-170 µm/s (not the 500+ µm/s
of the paper's headline), and after removing the net direction the residual
(0.5-0.7 µm rms, 3-8 cycles per track) is 7-11x the quantisation floor, so a
helix of 0.2-1 µm radius is not reliably resolvable. Many tracks are also
duplicates of the same cell, and the 400 fps readme disagrees with the stored
frame rate (200 fps). It is a population-tracking dataset, not a flagellar-helix one.
