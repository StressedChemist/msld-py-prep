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

## Applying the patch to other small molecules

### Scope and stereoisomer classification

The current patch recognizes assigned, atom-centred tetrahedral
stereochemistry. It does not yet provide equivalent safeguards for E/Z double
bonds, atropisomers, or other axial/planar stereogenic elements.

For two molecules with the same constitution, compare corresponding
stereocentres after establishing a chemically correct atom mapping:

| Relationship | Pattern across mapped tetrahedral stereocentres |
| --- | --- |
| Same stereoisomer | All assignments describe the same configuration |
| Enantiomer candidates | Every stereocentre is inverted |
| Diastereomers | At least one stereocentre is retained and at least one is inverted |

“Every centre is inverted” is only a screening rule. Enantiomers are
non-superimposable mirror images, so molecular symmetry must also be checked.
Meso compounds, equivalent atoms, pseudoasymmetric centres, and a change in
CIP priority can make a simple comparison of `R`/`S` text labels misleading.
Use a mirror-image/isomorphism check when symmetry is possible. In an achiral
environment, enantiomers should have identical equilibrium thermodynamics;
they can differ in a chiral protein or other chiral environment.

### 1. Prepare stereochemically explicit inputs

Create one input stem per lambda endpoint, for example:

```text
ligand_RR
ligand_RS
```

List those stems, without extensions, in `mol_list.txt`. For every stem,
provide at least the SDF structure and matching CGenFF stream/parameter input
required by the preparation workflow. Retaining matching MOL2, RTF, and PRM
files is recommended for traceability and downstream construction.

The following input rules are essential:

- Use a 3D SDF with the intended configuration, or a 2D SDF with explicit
  wedge/hash bonds. Do not use an unannotated 2D structure.
- Retain explicit hydrogens consistently across the series.
- Keep SDF atoms and force-field atoms in exactly the same order. The current
  workflow associates atom types, charges, and names by index.
- Give atoms unique, stable names and use those names consistently in the SDF,
  MOL2, and CHARMM files.
- Parameterize every endpoint before the MCS step. A stereochemical match does
  not compensate for missing or inconsistent atom types and charges.

The patched loader assigns chiral tags from 3D coordinates when available,
assigns CIP labels, and stops if RDKit finds a potentially chiral but
unassigned tetrahedral centre. Treat that failure as an input error; do not
disable it merely to make preparation continue.

### 2. Run the chirality-aware MCS

Run from a working directory containing `mol_list.txt` and the ligand files.
Set these paths for the local installation rather than copying the absolute
paths used for this example:

The default `msld_py_prep.py` entry point on this branch now calls
`MCSS_RDecomp`, so its normal two-pass workflow is chirality-aware. Configure
the system variables at the top of that file before using it. To run only the
MCS/decomposition stage explicitly, use:

```bash
export MSLD_PREP_ROOT=/path/to/msld-py-prep
export MSLD_PYTHON=python
export PYTHONPATH="$MSLD_PREP_ROOT${PYTHONPATH:+:$PYTHONPATH}"

"$MSLD_PYTHON" -c \
  'from msld_mcs_rdecomp import MCSS_RDecomp; MCSS_RDecomp("mol_list.txt", "MCS_for_MSLD.txt")'
```

The command writes `MCS_for_MSLD.txt`, `core_names.csv`, `core_charges.csv`,
`stereoisomer_pairs.csv`, and `stereochemistry_audit.csv`. The example
`prepare_pair.py` shows how to connect this step to `MsldCRN` and `MsldPRM`;
change its hard-coded `stems` list and input directory for another molecule
set.

For constitutionally identical pairs, `stereoisomer_pairs.csv` reports the
mapped centres that are retained and inverted. It labels a mixture of retained
and inverted centres as `diastereomers`; inversion of every mapped centre is
reported conservatively as `enantiomer_candidate`. Classification requires an
unambiguous graph mapping selected by unique, stable atom names. The workflow
stops rather than guessing when symmetry leaves that mapping ambiguous.

`stereochemistry_audit.csv` contains one row per assigned tetrahedral centre.
For every centre whose configuration changes, it reports whether the centre
and all of its neighbours remained in the MCS and whether the endpoint RTF has
a center-first `IMPR` record. It resolves that record's atom types against the
endpoint PRM and records the force constant and target angle. Preparation fails
by default if the full local stereochemical environment remains in the core,
an improper or its parameter is missing, the parameter is ambiguous, the force
constant is not positive, or its target lies within 5 degrees of a planar
0/180-degree geometry. The reports are written before the exception so they
can be inspected to correct the inputs.

For diagnostic comparisons only, the improper check can be relaxed with
`MCSS_RDecomp(..., require_stereo_impropers=False)`. This does not make the
result safe for dynamics and must not be used for production preparation.

### 3. Inspect the core and endpoint fragments

Do not accept the output solely because the program completed. Check
`stereoisomer_pairs.csv`, `stereochemistry_audit.csv`, `core_names.csv`,
`MCS_for_MSLD.txt`, and the generated
`site*_sub*_pres.rtf` files.

For every stereocentre that differs between endpoints:

- the stereochemistry-defining atoms must not all remain in one shared core;
- each endpoint must have its own coordinates for enough of the local
  tetrahedral environment to define its configuration; and
- the alternatives must be separate substituents at the same lambda site.

If connectivity-identical stereoisomers produce a full shared core with no
distinguishing fragments, stereochemistry has been ignored. The patched code
fails closed in this case. If nearly the entire ligand becomes substituent,
the construction may still be valid, but it approaches separate absolute
calculations and loses much of the efficiency of a relative transformation.
Highly symmetric molecules require additional scrutiny because several graph
mappings can be chemically equivalent while giving different atom-name maps.

### 4. Add a signed restraint for each duplicated stereocentre

MCS separation alone does not stop a tetrahedral centre from inverting during
dynamics. Each endpoint-specific centre therefore needs a signed improper or
an equivalent chiral restraint.

For a centre `C` and three ordered neighbours `A`, `B`, and `D`, monitor the
signed volume

```text
V = ((r_A - r_C) x (r_B - r_C)) . (r_D - r_C)
```

The sign depends on the neighbour order, so define that order once and record
it. With the same chemically corresponding order, inverted configurations
must have opposite signs. The magnitude should remain comfortably separated
from zero; crossing zero is an inversion or a nearly planar configuration.

Add an endpoint-specific CHARMM improper using the native tetrahedral angle,
not a conventional zero-degree target that would flatten the centre. Two
equivalent conventions are possible:

1. retain the same neighbour order and use opposite signed target angles; or
2. swap two neighbours for one endpoint so both impropers use the same signed
   target, as done for erythritol here.

For example, this test uses:

```text
IMPR C3 C4 C2 O3   ! meso endpoint
IMPR C3 C2 C4 O3   ! threo endpoint; two neighbours reversed
```

Measure the target from each minimized input geometry and choose the force
constant deliberately. Do not copy the erythritol target or force constant to
an unrelated chemistry without testing it. After building the hybrid PSF,
confirm that one intended improper exists for every endpoint-specific
stereocentre and that it refers to the translated hybrid atom names.

### 5. Validate before production sampling

Use the following acceptance sequence for every differing stereocentre:

1. Run a short fixed-lambda simulation at each endpoint.
2. Track the signed volume using one consistent chemical neighbour order.
3. Require the expected sign at every saved sample and inspect the minimum
   absolute volume; merely comparing the first and last frames is insufficient.
4. Run a solvated lambda/Gibbs test and verify that lambda visits all intended
   states and makes repeated transitions while every endpoint-specific signed
   volume remains on its expected side of zero.
5. Inspect the trajectory for distorted bonds, flattened centres, solvent
   overlap, and inactive-fragment instability.

A passing construction test establishes that the dual topology can move in
lambda without racemizing. It does not establish converged bias weights,
equilibrium populations, or a production-quality free-energy difference.
Those require longer sampling, fixed or converged biases, independent
replicates, and the usual free-energy convergence checks.

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
- classifies connectivity-identical pairs using a stable atom-name mapping,
  writes stereochemical pair/centre audit CSV files, and requires a
  center-first endpoint improper with an active, non-planar parameter at each
  changing centre; and
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

## Conventional CHARMM continuous-MSλD validation

The same hybrid was independently tested with CHARMM Developmental 50a2's
native BLOCK/QLDM implementation rather than discrete LaDyBUGS Gibbs updates.
The test used an F2SI two-state continuous lambda coordinate, `EXCL 2 3`, and
`RMLA BOND THET IMPR`, so the signed endpoint impropers remained fully active
throughout lambda propagation. The molecular system stayed at 300 K in the
same 2.0 nm box of 258 TIP3P waters.

For this construction stress test, the fictitious lambda thermostat was raised
to 3000 K and a +2.2 kcal/mol fixed meso bias was applied. These choices force
repeated endpoint traversal; they are not proposed production bias or
temperature settings and cannot be used to estimate a free-energy difference.

| Metric | Full 20 ps | First 10 ps | Second 10 ps |
| --- | ---: | ---: | ---: |
| Saved coordinate/lambda frames | 2,000 | 1,000 | 1,000 |
| Lambda range | 0.0–1.0 | 0.0–1.0 | 0.0–1.0 |
| Crossings of lambda = 0.5 | 83 | 33 | 50 |
| Alternations between lambda ≤ 0.1 and ≥ 0.9 | 47 | 15 | 32 |
| Meso endpoint frames, lambda ≥ 0.9 | 906 | 461 | 445 |
| Threo endpoint frames, lambda ≤ 0.1 | 173 | 77 | 96 |

Across all 2,000 frames, meso C3 remained between -0.002996 and
-0.001838 nm³ and threo C3 remained between +0.001869 and +0.002947 nm³.
Neither endpoint had a sign flip. Lambda motion and endpoint visits persisted
in both halves of the trajectory, so this native CHARMM run also passes the
construction criterion. The compact trace, plot, and machine-readable summary
are in `runs/patched_pair_with_improper/charmm_msld_validation`.

## Reproduction

From this `stereo_analysis` directory, select the Python interpreter that has
RDKit and OpenMM installed. Choose new output-directory names if the recorded
results already exist, because the runners refuse to overwrite them.

```bash
export MSLD_PYTHON=python

"$MSLD_PYTHON" create_inputs.py
"$MSLD_PYTHON" prepare_pair.py runs/patched_pair_with_improper
"$MSLD_PYTHON" build_hybrid.py runs/patched_pair_with_improper
"$MSLD_PYTHON" endpoint_test.py runs/patched_pair_with_improper --steps 100000
"$MSLD_PYTHON" solvated_lambda_test.py runs/patched_pair_with_improper --output-name solvated_lambda --box-nm 2.0 --minimization-iterations 500 --equilibration-steps 5000 --cycles 5000 --steps-per-cycle 10 --timestep-fs 1.0 --initial-bias-increment 5.0 --minimum-bias-increment 0.01 --flatness-interval 250 --flatness-fraction 0.6 --progress-interval 100 --trajectory-interval 500 --threads 8
"$MSLD_PYTHON" plot_solvated_lambda.py runs/patched_pair_with_improper/solvated_lambda
```

For the native CHARMM test, first make a CHARMM executable available on
`PATH`. On the test cluster that was done with `module load
charmm/charmm/c50a2`; other installations should use their local environment.
Some MPI builds may also require site-specific runtime settings. Then run:

```bash
"$MSLD_PYTHON" charmm_msld_test.py \
  runs/patched_pair_with_improper \
  --output-name charmm_msld_validation_reproduction \
  --charmm "$(command -v charmm)" \
  --steps 20000 --save-frequency 10 --threads 8 \
  --meso-bias 2.2 --lambda-temperature 3000 --require-pass

"$MSLD_PYTHON" plot_solvated_lambda.py \
  runs/patched_pair_with_improper/charmm_msld_validation_reproduction
```

`charmm_msld_test.py` refuses to overwrite an existing result directory. It
also requires equal coordinate and lambda frame counts, reports motion in each
half of the trajectory, and sets `working_construction` only when lambda
alternations persist and both signed volumes retain their expected signs with
an absolute-volume margin of at least 0.001 nm³. `--require-pass` makes a failed
acceptance test return a nonzero status.

The preparation runners deliberately refuse to overwrite an existing result.
