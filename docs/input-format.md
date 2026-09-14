# Importing independent calculations

The boundary between upstream physics and this package is an explicit set of
site-aligned arrays. The package does not claim to launch Twister, kp_tblg,
pytwist, CP2K, or Julia. Export their outputs into this common format, recording
the exact code commit, parameters, convergence checks, and normalization.

## CSV route

`sites.csv` has one row per sampled local environment:

```csv
site_id,registry,x_nm,y_nm,lambda_ev,coupling_ev
AB-center,AB,0,0,1.0,0.001
AA-center,AA,5,5,0.85,0.0012
SP-center,SP,2.5,2.5,0.95,0.0011
```

These numbers are **format examples only**, not calculated values.

`ldos.csv` is long-form, with a full common energy grid for every site:

```csv
site_id,energy_ev,dos_states_per_ev_nm2
AB-center,-2.0,1.0
AA-center,-2.0,1.1
SP-center,-2.0,1.05
```

Continue with all remaining energy samples for every site. The fragment above
is intentionally not a complete dataset. At least three energy points are
required; realistic calculations need many more. Rows may be unsorted. The
importer joins by site ID and sorts energy, rejecting duplicate energies, missing
site spectra, inconsistent grids, unknown IDs, nonfinite values, and negative DOS.

```bash
tbg-et import-csv --sites sites.csv --ldos ldos.csv \
  --metadata metadata.json --output results/imported
```

## Native NPZ route

Use `numpy.savez_compressed` with these exact names. Object arrays/pickle are
not accepted.

| Array | Shape | Meaning |
| --- | --- | --- |
| `energy_ev` | (M,) | Finite, strictly increasing energies |
| `dos` | (N, M) | Nonnegative states/eV/nm²; nonzero spectrum per site |
| `xy_nm` | (N, 2) | Physical sampled coordinates |
| `registry` | (N,) Unicode | AA, AB (including BA), or SP |
| `lambda_ev` | (N,) | Positive reorganization energy |
| `coupling_ev` | (N,) | Positive coupling magnitude |
| `site_ids` | (N,) Unicode | Unique nonempty IDs |

For atomic projected DOS, explicitly convert per-atom or per-cell units to the
areal normalization, and record the projection weight, degeneracies and reference
area. Do not normalize each local spectrum independently to unity: that removes
part of the physical DOS variation under study. Use one compatible energy zero.

## Metadata

`metadata.json` is required for both routes:

```json
{
  "schema_version": 1,
  "kind": "external",
  "source": "Describe input files, solver versions, commits and calculation IDs",
  "dos_units": "states/eV/nm^2",
  "energy_reference": "Describe the common energy zero and Fermi level",
  "structure_model": "Describe relaxation, force field and geometry provenance",
  "dos_model": "Describe Hamiltonian, mesh, projection, broadening and normalization",
  "lambda_model": "Describe screening calculation or sampled energy-gap ensembles",
  "coupling_model": "Describe diabatic states, molecular pose, basis and extraction",
  "limitations": ["List convergence limitations and missing physical effects"]
}
```

Replace every descriptive placeholder before interpreting a result. `external`
means supplied independently, not verified by this package. Additional JSON
fields are preserved. Use `synthetic` for constructed demonstration inputs.
No automatic importer can determine whether provenance statements are true.

## Upstream microscopic helpers

```python
from tbg_et.structure import classify_registry
from tbg_et.microscopic import two_state_coupling, gap_reorganization

# Local relative-displacement field from an independently relaxed structure:
registry = classify_registry(displacement_fractional)

# Matrix elements from the same normalized real two-state diabatic basis:
h_ev = two_state_coupling(h_aa_ev, h_bb_ev, h_ab_ev, s_ab)

# Same vertical-gap convention in both equilibrated charge-state ensembles:
estimates = gap_reorganization(gaps_oxidized_ev, gaps_reduced_ev)
lambda_ev = estimates["lambda_mean_ev"]
```

The quantities in this snippet must come from calculations you run independently.
These helpers do not generate or simulate those data.
