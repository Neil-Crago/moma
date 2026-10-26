# Test data

`pea_exp3box8AT1.csv` is one file (`3DtrajectorySmooth_P1_6 EXP3box8AT1.csv`, renamed) from:

> Francisco de Carvalho Oliveira, T.; Niemeyer Reissig, G.; Ardini, D.; Nava, D. E.;
> Maia Souza, G.; Castiello, U.; Guerra, S. *Raw kinematical recordings of Pisum
> sativum during circumnutation under support and no-support conditions.* Zenodo
> (2026). DOI [10.5281/zenodo.22015401](https://doi.org/10.5281/zenodo.22015401).
> Licensed CC BY 4.0.

Plain `x,y,z` rows, no header; y is the growth direction. The other 19 recordings
are at the DOI above. This directory is not included in the published crate.

**`1UCU.pdb`** is the PDB entry 1UCU (Yonekura et al., *Nature* 2003, "R-type
straight flagellar filament made of full-length flagellin", 4 Å cryo-EM; one
Salmonella flagellin chain). Public domain (PDB). Used by `examples/pdb_helices.rs`
as a real-coordinate accuracy check: the alpha helices in it have known geometry.

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

**Flagellin sequence clusters** (`41467_2024_53923_MOESM3_ESM.csv`, Nature
Communications 2024 supplementary table 3; 16,948 rows of UniProt ID, cluster,
organism, flagellin length, AlphaFold score, Pfam domains) and
**`Flagellum_71_protein_dataset.tar.gz`** (figshare; 71 FASTA files of flagellar
proteins). Neither contains coordinates or time series, so neither can test the
axis or rotation tools. They are sequence and annotation data. Both are kept
untracked and are not committed.

**3D-SpermFlagella** (`Sperm_flagellum_centerlines_3D.zip`, 117 MB, untracked):
3D+t human sperm flagellum centrelines, Hernandez Herrera, Hernandez,
Bribiesca-Sanchez, Montoya, Diaz-Guerrero, Darszon, Corkidi; Sci Data 2026,
doi:10.1038/s41597-026-06876-2; Zenodo 10.5281/zenodo.15299846; CC BY 4.0.
135 cells, 1-3.5 s each. Each cell folder has `X.csv`/`Y.csv`/`Z.csv` with one
row per point along the tail (head to tip) and one column per time point; NaNs
cluster at the tip. The frame rate is not in the Zenodo record, so periods are
in samples. Used by `examples/sperm_beat.rs`.

Frame rate (from the 3D-SpermVid paper, doi:10.1038/s41597-025-05177-4): objective oscillates at 90 Hz, hyperstacks at 90 per second. Record lengths (65-499 columns) imply 0.7-5.5 s at that rate versus the stated 1-3.5 s, so treat Hz as provisional.
