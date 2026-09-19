# Meso-erythritol/threitol MCS analysis

## Inputs

`create_inputs.py` preserves the existing `(2S,3S)` coordinates as threitol
and generates meso-erythritol by reflecting the `O3`, `H5`, and `H9` group
through the `C2-C3-C4` plane. This preserves the local bond geometry while
changing C3 configuration.

RDKit 2024.09.6 assignments after writing and reloading both formats:

- `meso_erythritol.sdf` and `.mol2`: `(2S,3R)`
- `threitol.sdf` and `.mol2`: `(2S,3S)`

Both stereoisomers use the same atom ordering, atom types, charges, and bonded
connectivity. This isolates stereochemistry as the only MCS difference.

## Results

### Meso-erythritol versus threitol only

The stock `msld_mcs_rdecomp.py` MCS matched the entire molecule (18 atoms and
17 bonds). `core_names.csv` therefore contains every atom for both ligands.
R-group decomposition returned no groups and the code terminated with
`KeyError: 'Core'`. No alchemical distinction was constructed.

With `matchChiralTag=True`, the MCS contained 15 atoms. It assigned
`O3 H5 H9` to one site with C3 as its anchor. Charge renormalization then moved
C3 into each fragment, producing two four-atom alternatives (`C3 O3 H5 H9`)
that bridge the shared C2 and C4 atoms. The final common core contains 14 atoms.

This demonstrates that the meso/threo transformation is topologically
constructible without reducing the shared core to almost nothing, but the
stock MCS does not discover it.

### Corrected original four-ligand series

For `meso-erythritol, glycerol, xylitol, sorbitol`, the stock and
chirality-aware calculations produced byte-identical MCS and RTF outputs.
The structural differences among the chain lengths already reduce the core to
the terminal `C1/O1` region. C2 is an anchor and is moved into the alchemical
fragment during charge renormalization, so all erythritol stereocentres are in
the meso-erythritol fragment.

### Five-state series including both diastereomers

For `meso-erythritol, threitol, glycerol, xylitol, sorbitol`, both modes again
produced byte-identical MCS and RTF outputs. Meso-erythritol and threitol are
retained as separate substituents. Each entire C2-and-beyond region is
duplicated, so each lambda end state carries its own coordinates.

This success is contextual rather than evidence that the stock MCS handles
chirality: the shorter/longer ligands force the stereochemical region out of
the common core. The isolated-pair test proves that the stock MCS still
collapses stereoisomers when connectivity and atom types permit it.

## Patched workflow and C3 restraint

The production `msld_mcs_rdecomp.py` now:

- assigns stereochemistry from 3D coordinates and rejects unassigned centres;
- uses `matchChiralTag=True` in the MCS;
- uses the MCS query molecule directly and requires chirality in substructure
  matching and R-group decomposition; and
- fails explicitly when decomposition returns neither groups nor unmatched
  molecules, instead of indexing a missing `Core` result.

Running that patched implementation directly on the pair produced the same
15-atom MCS as the earlier diagnostic monkeypatch. Charge renormalization moved
the C3 anchor into each alternative, leaving a 14-atom core and two
`C3/O3/H5/H9` fragments. The unmodified input force field contained no
impropers; this is preserved in `runs/patched_pair_no_improper` as the control.

A conventional 0-degree CHARMM improper would incorrectly flatten tetrahedral
C3. The regenerated pair inputs therefore use signed, endpoint-specific atom
orders:

- meso: `IMPR C3 C4 C2 O3`
- threo: `IMPR C3 C2 C4 O3`

Both initial structures have an improper of +32.60 degrees under their
respective order. The added restraint is 50.0 kcal/mol/rad^2 with a +32.60
degree target. The generated fragment patches contain translated terms
`IMPR C018 C009 C002 O015` and `IMPR C022 C002 C009 O019`, and the assembled
hybrid PSF contains exactly two impropers.

## Fixed-lambda endpoint smoke test

Both endpoints were tested for 100 ps at 300 K with a 1 fs timestep using
OpenMM 8.2 on CPU. At each endpoint, the inactive fragment's electrostatic and
Lennard-Jones interactions were set to zero while all bonded terms remained
active, matching the LaDyBUGS bonded-force treatment. This was a vacuum
`NoCutoff` smoke test with hydrogen-bond constraints, not a production
solvated/bound validation.

| Endpoint | Lambda | C3 signed-volume range (nm^3) | Improper range | Sign flips |
| --- | --- | --- | --- | --- |
| meso | `[1, 0]` | -0.002862 to -0.001884 | 27.01 to 37.56 deg | 0/101 |
| threo | `[0, 1]` | +0.001981 to +0.002870 | 25.74 to 37.06 deg | 0/101 |

Both endpoints retained their starting C3 configuration. Full statistics are
in `runs/patched_pair_with_improper/endpoint_summary.json`, with sampled values
in `endpoint_meso.csv` and `endpoint_threo.csv`.

## Solvated lambda/Gibbs validation

The patched hybrid was also run in a 2.0 nm periodic box containing 258 TIP3P
waters. After 5 ps of endpoint equilibration, a bias-updated all-state Gibbs
sampler moved over 11 lambda states for 50 ps at 300 K with a 1 fs timestep.
Signed C3 volume and the endpoint-specific improper were recorded after every
10 MD steps (5,000 samples total).

| Metric | Result |
| --- | --- |
| Lambda states visited | 11/11 |
| Lambda-state transitions | 2,450 |
| Alternations between exact endpoints | 16 |
| Meso C3 signed volume | -0.003065 to -0.001663 nm^3 |
| Threo C3 signed volume | +0.001656 to +0.003107 nm^3 |
| Meso/threo sign flips | 0/5,000 and 0/5,000 |
| Separated throughout | yes |

The construction therefore passes this validation: lambda crosses the state
grid repeatedly while the two endpoint geometries retain opposite C3 signs.
This is a construction/stability test with an adapting bias, not evidence of
converged equilibrium populations or a production-quality free-energy result.
The complete statistics are in
`runs/patched_pair_with_improper/solvated_lambda/summary.json`, the per-cycle
trace is `lambda_stereo.csv`, and the coordinates are in `trajectory.dcd`.

## Reproduction

Run with the environment used here:

```bash
/home/boloyede/envs/charmm_2910746/bin/python create_inputs.py
/home/boloyede/envs/charmm_2910746/bin/python prepare_pair.py runs/patched_pair_with_improper
/home/boloyede/envs/charmm_2910746/bin/python build_hybrid.py runs/patched_pair_with_improper
/home/boloyede/envs/charmm_2910746/bin/python endpoint_test.py runs/patched_pair_with_improper --steps 100000
/home/boloyede/envs/charmm_2910746/bin/python solvated_lambda_test.py runs/patched_pair_with_improper --output-name solvated_lambda --box-nm 2.0 --minimization-iterations 500 --equilibration-steps 5000 --cycles 5000 --steps-per-cycle 10 --timestep-fs 1.0 --initial-bias-increment 5.0 --minimum-bias-increment 0.01 --flatness-interval 250 --flatness-fraction 0.6 --progress-interval 100 --trajectory-interval 500 --threads 8
/home/boloyede/envs/charmm_2910746/bin/python plot_solvated_lambda.py runs/patched_pair_with_improper/solvated_lambda
```

The preparation runners deliberately refuse to overwrite an existing result.
