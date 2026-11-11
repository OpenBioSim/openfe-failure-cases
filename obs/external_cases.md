# pxr and asap cases

This note looks at the free (solvent) leg of seven cases from `failures/2026-08-14_pxr_user_set` and `failures/2026-08-25_asap_pubilc_failures` with SOMD2 and Ghostly. Unlike the janssen bace edges in [hybridisation_changes.md](hybridisation_changes.md), their ligands aren't in the openfe-benchmarks data, so they are built from each case's transformation JSON.

| Case | OpenFE's description |
| --- | --- |
| pxr OCNT-2311573~OADMET-0006255 | solvent failure, triple bond near the alchemical region |
| pxr OADMET-0006142~OCNT-2316690 | one mapped heavy atom, large MBAR error |
| asap eva71_2apro_first_test | solvent and complex failure, hybridisation change |
| asap EV712AproBLRN edge 1 | solvent and complex failure, hybridisation change |
| asap EV712AproBLRN edge 2 | solvent and complex failure, poor alignment from a 2D mapping |
| asap eva71_2apro_prospective edge 1 | solvent failure, overlapping hydrogens |
| asap SHP2XPMN | complex failure only, possible hybridisation change |

## Preparation

`prepare_external.py` rebuilds both ligands from the transformation JSON (atoms, bonds, formal charges, coordinates and OpenFE's partial charges), parameterises them with OpenFF 2.3.0 (unconstrained), and merges and solvates them with the openfe-benchmarks code, once with OpenFE's mapping and once with BioSimSpace's. The cases are listed in `cases_external.json`. The transformation JSONs for EV712AproBLRN and SHP2XPMN contain no partial charges, so these ligands use the charges assigned during parameterisation.

| Case | OpenFE | BioSimSpace |
| --- | --- | --- |
| pxr OCNT-2311573~OADMET-0006255 | 38 atoms | 39 atoms |
| pxr OADMET-0006142~OCNT-2316690 | 4 atoms | 36 atoms |
| asap eva71_2apro_first_test | 37 atoms | 37 atoms, identical |
| asap EV712AproBLRN edge 1 | 37 atoms | 37 atoms, identical |
| asap EV712AproBLRN edge 2 | 40 atoms | 40 atoms |
| asap eva71_2apro_prospective edge 1 | 41 atoms | 42 atoms |
| asap SHP2XPMN | 43 atoms | 44 atoms |

## Anchor dihedrals through a linear nitrile

In four cases, eva71_2apro_first_test, EV712AproBLRN edge 1, pxr OCNT-2311573 and SHP2XPMN (the last two with BioSimSpace's mapping), a group on one ligand is mapped onto a nitrile on the other, so its hydrogens become ghosts attached to the nitrile nitrogen. That nitrogen has a single physical neighbour, the nitrile carbon, so Ghostly treats it as a terminal junction and keeps one anchor dihedral per ghost, ghost–N–C–C. At the nitrile end state the angle C–C≡N is linear, so these dihedrals pass through three collinear atoms, where they are undefined and their forces diverge:

| Case | Kept anchor dihedrals |
| --- | --- |
| eva71_2apro_first_test | three, k = 0.43 kcal/mol, at λ = 1 |
| EV712AproBLRN edge 1 | three, k = 0.32 kcal/mol, at λ = 1 |
| pxr OCNT-2311573 | two, amplitude up to about 2.9 kcal/mol, at λ = 0 |
| SHP2XPMN | two, amplitude up to about 2.9 kcal/mol, at λ = 1 |

This happens at the end state itself, so SOMD2's `torsion_window` can't prevent it. Ghostly now removes dihedrals and impropers involving ghost atoms that pass through an angle that is linear at the end state, and prefers anchors that aren't reached through a linear angle ([ghostly#24](https://github.com/OpenBioSim/ghostly/issues/24), fixed in [ghostly#23](https://github.com/OpenBioSim/ghostly/pull/23)). The ghosts remain held by their angle terms, so they can only rotate about the linear axis.

With that fix, `torsion_window` also applies to five cases, where an angle becomes linear at one end state with non-zero torsions through it at the other: eva71_2apro_first_test and EV712AproBLRN edge 1 (C–S–C and C–O–C becoming linear at λ = 1), pxr OCNT-2311573 (λ = 0) and SHP2XPMN (λ = 1). No angle becomes linear in pxr OADMET-0006142, EV712AproBLRN edge 2 or eva71_2apro_prospective.

## Results

Free leg, BioSimSpace mapping. Each case was equilibrated with the openfe-benchmarks protocol (`obs/equilibrate.py`), using the fixed Ghostly, then run with SOMD2 using replica exchange, 11 evenly spaced λ windows of 1 ns each and a 4 fs timestep. "Old Ghostly" is 0d456c0, before the fix for anchor dihedrals through linear angles. "New Ghostly" is ae91e52, with it. `torsion_window` is `auto` unless stated.

| Case | Old Ghostly | New Ghostly | New Ghostly, `torsion_window false` |
| --- | --- | --- | --- |
| eva71_2apro_first_test | `NaN` at λ = 1.0, no data | completed | completed |
| EV712AproBLRN edge 1 | `NaN` at λ = 1.0, no data | completed | completed, one recovered failure |
| pxr OCNT-2311573 | `NaN` at λ = 0.0, no data | completed | |
| SHP2XPMN | `NaN` at λ = 0.9, no data | completed | |
| EV712AproBLRN edge 2 | | completed | |
| pxr OADMET-0006142 | | completed | |
| eva71_2apro_prospective edge 1 | | failed in equilibration | |

With the old Ghostly, all four nitrile cases crash before saving any data, at or next to the nitrile end state. With the new Ghostly, all of them run, so there's no free energy to compare between the two versions.

eva71_2apro_prospective edge 1 can't be minimised as submitted: each ligand has two hydrogens at identical positions (atoms 36 and 37, and atoms 26 and 27), matching OpenFE's note about overlapping hydrogens. This is a problem with the input rather than the engine, so we didn't change the input to run it.

The MBAR free energies of the free leg, in kcal/mol, with the minimum overlap between neighbouring windows:

| Case | ΔG | Minimum overlap |
| --- | --- | --- |
| eva71_2apro_first_test | 5.09 ± 0.43 | 0.010 |
| eva71_2apro_first_test, `torsion_window false` | 4.81 ± 0.42 | 0.012 |
| EV712AproBLRN edge 1 | 0.04 ± 0.34 | 0.025 |
| EV712AproBLRN edge 1, `torsion_window false` | 1.11 ± 0.33 | 0.025 |
| pxr OCNT-2311573 | 42.99 ± 0.16 | 0.112 |
| SHP2XPMN | 25.58 ± 0.14 | 0.149 |
| EV712AproBLRN edge 2 | 4.37 ± 0.09 | 0.233 |
| pxr OADMET-0006142 | −105.25 ± 0.15 | 0.050 |

A free-leg ΔG isn't meaningful on its own, but these show that the runs give usable data. The window doesn't change the end states, so it can't change the free energy once sampling has converged. For eva71_2apro_first_test the values with and without it agree, while for EV712AproBLRN edge 1 they differ by 1.1 kcal/mol. Both are single 1 ns runs with very poor overlap, so the MBAR error estimates are likely too small, and these two edges would need more λ windows for production runs.

For pxr OADMET-0006142, BioSimSpace's mapping has 36 atoms where OpenFE's has 4, so most of the ligand is mapped rather than a ghost, and the free leg has a small MBAR error. EV712AproBLRN edge 2, which failed in both legs with OpenFE, runs without problems.

## Summary

- Four of the cases have a ghost group attached to a nitrile, where Ghostly kept anchor dihedrals through a linear angle at the end state. They crash with the old Ghostly and run with the fix.
- `torsion_window` applies to five of the cases, and the fix to Ghostly is needed as well as the window for the four nitrile cases.
- EV712AproBLRN edge 2 and pxr OADMET-0006142 run with SOMD2 without any special treatment. For pxr OADMET-0006142 our mapping is much larger than OpenFE's.
- eva71_2apro_prospective edge 1 has overlapping hydrogens in the submitted structures and can't be minimised.
- The comparisons with and without the window are single short runs with poor overlap, so they don't show whether the free energy is unchanged.

## Reproducing

Software: SOMD2 devel at 6dea460, Ghostly `fix_22` at ae91e52 and 0d456c0, Sire and BioSimSpace 2026.3.0.dev, OpenMM 8.6, on an NVIDIA RTX 5070 Ti.

From this repository, with openfe-benchmarks (`feature_obs` branch) checked out next to it:

```
python obs/prepare_external.py --benchmarks ../openfe-benchmarks
```

This writes the free legs to `obs_output/external/<case>/{openfe,biosimspace}/`. Each was equilibrated and run as in [hybridisation_changes.md](hybridisation_changes.md), with `PYTHONPATH` pointing at a copy of Ghostly 0d456c0 for the old-Ghostly runs. The scripts and configs used are in `obs_output/runs_external/`, which isn't tracked.
