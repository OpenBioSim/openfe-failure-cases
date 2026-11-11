# Hybridisation changes: janssen_bace bace_p3_arg368_in

The case `failures/2026-07-29_janssen_bace_bace_p3_arg368_in` reports that three OpenFE RBFE edges, 8f~8g, 17f~17g and 27f~27g, failed with `NaN` in most attempts, in both the solvent and complex legs, sometimes only after about 30% of the sampling time. The submitters suspected the change in hybridisation in the mapping.

This note looks at the free (solvent) leg of the three edges with SOMD2 and Ghostly. It covers whether our default mapping avoids the change in hybridisation, why the simulations fail, and whether SOMD2's `torsion_window` option fixes them.

## The perturbation

All three edges replace a methoxyethoxy group with a butynyloxy group on the same scaffold:

| Ligand | Substituent |
| --- | --- |
| 8f, 17f, 27f | `CH3–O–CH2–CH2–O–Ar` |
| 8g, 17g, 27g | `CH3–C≡C–CH2–O–Ar` |

Both OpenFE's mapping and ours map the ether oxygen and the first CH2 carbon onto the two alkyne carbons. So the angles C–O–C and O–C–C, at about 110–115° in the f ligands, become C–C≡C and C≡C–C at 180° in the g ligands. The two hydrogens on the CH2 carbon become ghost atoms in the g end state.

## Mappings

`make_mappings.py` takes OpenFE's atom mapping from each case's transformation JSON, checks its atom order against the benchmark ligands, and writes two mappings files for openfe-benchmarks' `obs/prepare_rbfe.py`. One uses OpenFE's mapping. The other maps the atoms with BioSimSpace, our default.

| Edge | OpenFE | BioSimSpace |
| --- | --- | --- |
| 8f~8g | 49 atoms | 50 atoms |
| 17f~17g | 53 atoms | 53 atoms, identical |
| 27f~27g | 56 atoms | 56 atoms, identical |

For 8f~8g, the mappings differ only at the terminal methyl. BioSimSpace maps all three hydrogens, while OpenFE leaves one unmapped on each side. In every case the change in hybridisation is in the mapped core, so our mapping doesn't avoid it.

## Why the simulations fail

Torsions and impropers that pass through an angle are interpolated between their values at the two end states, as is the angle itself. At the f end state the torsions through C–O–C and O–C–C are non-zero, and at the g end state, where the angle is linear, they are zero. Near λ = 1 the angle's equilibrium value is close to 180° while the torsions still carry part of their force constant. As the angle's atoms approach collinearity, the torsion's force diverges, and the simulation produces `NaN`. Whether this happens in a given run depends on how close a thermal fluctuation comes to 180°, which explains failures that appear only partway through sampling.

Ghostly removes the torsions that run from ghost atoms through the linear unit. That handles the dummy-atom side, but the torsions between mapped atoms are physical terms of the f ligand, so no dummy-atom scheme can remove them.

## SOMD2's torsion window

SOMD2's `torsion_window` option ([somd2#287](https://github.com/OpenBioSim/somd2/pull/287)) finds angles of a perturbable molecule that are linear at only one end state and have a non-zero torsion or improper through them at the other. When it finds one, the torsions and impropers of that molecule change over only part of each stage of the λ schedule: the first half when the angle is linear at λ = 1, the second half when it's linear at λ = 0. Within that window they follow a smoothstep. So by the time the angle approaches 180°, the torsions through it are already at their final value of zero. The option is `auto` by default.

For all three edges, with both mappings, it detects C–O–C and O–C–C becoming linear at λ = 1, with non-zero torsions through them at λ = 0. So it applies the window [0, 0.5]. The 180° angle between the two ghost hydrogens, which Ghostly sets, has no torsions through it and doesn't trigger it. No torsions are non-zero at the linear end state, where the schedule couldn't help.

## Results

Free leg, BioSimSpace mapping. Each edge was equilibrated with the openfe-benchmarks protocol (`obs/equilibrate.py`), then run with SOMD2 using replica exchange, 11 evenly spaced λ windows of 1 ns each and a 4 fs timestep, once with `--torsion-window false` and once with the default `auto`.

| Edge | `torsion_window false` | `torsion_window auto` |
| --- | --- | --- |
| 8f~8g | `NaN` at λ = 0.9 after 400 ps | completed, no errors |
| 17f~17g | completed, no errors | completed, no errors |
| 27f~27g | `NaN` at λ = 0.9 after 200 ps | completed, no errors |

17f~17g completing without the window is chance. Its alkyne angles reach 179.5° at λ = 0.9 too, but the torsions through them are weak (in 8f~8g the largest is about 0.6 kcal/mol), so reaching near-collinear only sometimes produces `NaN`. This matches OpenFE's observation that the failures didn't always happen immediately.

With the window, the alkyne angles behave as expected on all three edges:

| λ | C–O–C → C–C≡C (mean, max) | O–C–C → C≡C–C (mean, max) |
| --- | --- | --- |
| 0.8 | 166°, 176–178° | 165°, 175–177° |
| 0.9 | 172°, 178–179° | 172°, 177–180° |
| 1.0 | 175°, 179–180° | 176°, 179–180° |

The ranges cover the three edges. The window doesn't change exchange between neighbouring windows: for 17f~17g it was 1–2% for most pairs both with and without it. That overlap is low for all three edges whatever the schedule, because of the size of the perturbation, an O → C change in element together with sp3 → sp. Production runs of these edges would need more λ windows.

The MBAR free energies, in kcal/mol, with the data saved before the crash for the two runs that failed:

| Edge | `torsion_window false` | `torsion_window auto` |
| --- | --- | --- |
| 8f~8g | −16.24 ± 0.32 (400 ps) | −15.49 ± 0.20 |
| 17f~17g | −15.90 ± 0.20 | −15.16 ± 0.20 |
| 27f~27g | −15.88 ± 0.46 (200 ps) | −15.87 ± 0.20 |

The end states are the same with and without the window, so it can't change the free energy once sampling has converged. The differences of up to 0.75 kcal/mol for 8f~8g and 17f~17g come from single 1 ns runs with a minimum overlap between neighbouring windows of 0.06–0.07, where the MBAR error estimates are likely too small. Showing agreement would need repeats, more λ windows, or both.

## Summary

- Our default mapping doesn't avoid the change in hybridisation. It maps the same core as OpenFE.
- The failure mode reproduces in SOMD2: without the window, two of the three edges crash at λ = 0.9, the same edges and the same mechanism that OpenFE reports.
- With SOMD2's default `torsion_window auto`, all three edges run.

## Reproducing

Software: SOMD2 devel at 6dea460 (the merge of #287), Ghostly `fix_22` at 0d456c0, Sire and BioSimSpace 2026.3.0.dev, OpenMM 8.6, on an NVIDIA RTX 5070 Ti. Ligands are parameterised with OpenFF 2.3.0 (unconstrained) and `nagl_openff-gnn-am1bcc-1.0.0` charges, the openfe-benchmarks defaults.

From this repository, with openfe-benchmarks (`feature_obs` branch) checked out next to it:

```
python obs/make_mappings.py --benchmarks ../openfe-benchmarks

cd ../openfe-benchmarks/obs
for m in openfe biosimspace; do
    python prepare_rbfe.py --system-group janssen_bace --system-name bace_p3_arg368_in \
        --mappings ../../openfe-failure-cases/obs/mappings_$m.json --legs free \
        --output-dir ../../openfe-failure-cases/obs_output/$m
done
```

Then for each edge directory `E` under `obs_output/biosimspace/janssen_bace/bace_p3_arg368_in/rbfe/`:

```
python ../openfe-benchmarks/obs/equilibrate.py $E/free.bss --platform cuda
somd2 $E/free_lambda_0.bss --perturbed-system $E/free_lambda_1.bss --replica-exchange \
    --runtime "1 ns" --timestep "4 fs" --num-lambda 11 --torsion-window false
```

and the same with `--torsion-window auto`. The script and config used are in `obs_output/runs/`, which isn't tracked.
