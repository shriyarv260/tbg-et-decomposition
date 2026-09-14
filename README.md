# tbg-et-decomposition

**Separate local electron-transfer enhancement into density of states,
reorganization energy, and electronic coupling.**

A reproducible Python research workbench for comparing AA, AB/BA, and SP sites
in twisted bilayer graphene. It integrates local DOS with nonadiabatic
Marcus–Hush–Chidsey kinetics, evaluates all eight reference/target combinations,
and attributes the log-rate enhancement over all six replacement orders.

**Status: working analysis software with a synthetic demonstration.** This
repository does not yet contain relaxed TBG structures, calculated TBG LDOS,
RuHex CDFT/POD calculations, or experimental validation. The default twist-angle
dependence is prescribed. It must not be used as evidence for the AA anomaly,
a microscopic mechanism, or research novelty.

## Run it

Python 3.10 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
tbg-et demo --output results/demo
```

The default run covers **0.8°, 1.1°, 1.5°, and 2.0°**, with 576 spatial samples
and 1,601 energy samples per angle. No credentials, solver installation, or
network access is needed after package installation. On Windows, activate
with `.venv\Scripts\activate`.

Each angle produces:

| Artifact | Contents |
| --- | --- |
| `overview.png`, `overview.pdf` | Registry and rate maps, lambda, local DOS, potential sweep, attribution |
| `rates.csv` | Every site, inputs, full rate, and DOS-only comparison |
| `summary.json` | All eight hybrids, six paths, Shapley contributions, numerical checks |
| `report.md` | Readable interpretation and limitations |
| `input.npz`, `metadata.json` | Exact reusable inputs and provenance |
| `manifest.json` | SHA-256 hashes and runtime versions |

The parent directory includes `twist_sweep.csv` and `twist_sweep.png`.
Outputs require a new or empty directory to avoid overwriting previous runs.

## Analyze your own calculations

Use the [input specification](docs/input-format.md) to export compatible
site and LDOS arrays from your upstream workflow:

```bash
tbg-et import-csv --sites sites.csv --ldos ldos.csv \
  --metadata metadata.json --output results/imported

tbg-et analyze --input results/imported/input.npz \
  --metadata results/imported/metadata.json \
  --reference AB-center --target AA-center \
  --temperature 298.15 --eta 0 --mu 0 --area 1 \
  --output results/analysis
```

Site IDs are explicit: the tool never silently chooses an AA/AB reference or
averages dissimilar local DOS spectra. To try this command immediately, use
the `input.npz` and `metadata.json` generated in `results/demo/theta_1.1/`.
Synthetic provenance remains synthetic when reimported.

## What the calculation means

```text
independent structure / LDOS / lambda / coupling inputs
                         │
                         ▼
             validated site-aligned dataset
                         │
                         ▼
         MHC-DOS integral for each local site
                         │
                         ▼
      eight hybrids → six paths → log-rate attribution
                         │
                         ▼
               maps, tables, audit trail
```

DOS and lambda generally cannot be separated into unique independent factors
because lambda changes the energy weighting inside the DOS integral. We report
an explicit convention: the average contribution over all six orders in
**log-rate space**. Contributions sum to `ln(k_target / k_reference)`; their
exponentials multiply to the rate ratio. Negative contributions are allowed.
Lambda is held fixed in a direct-DOS intervention even if it was derived from DOS.
See [methods and units](docs/methods.md).

| Component | Implemented | Research work still needed |
| --- | --- | --- |
| Local structure | Periodic registry classifier; synthetic rigid displacement grid | Relaxed geometry and local displacement from an upstream solver |
| Electronic structure | Validated local DOS imports | Converged Hamiltonian/wavefunction/LDOS calculations |
| Kinetics | Stable reduction/oxidation integrals with physical invariants tested | Experimental calibration and appropriate electrochemical boundary conditions |
| Screening/lambda | Imported lambda; explicitly unfitted sensitivity surrogate; energy-gap estimator | Published nonlocal screening implementation or constant-potential simulations |
| Coupling | Imported magnitudes; real two-state orthogonalization utility | Molecule/electrode CDFT or POD calculations and basis consistency |
| Decomposition | Eight hybrids, all orders, exact log-rate closure | Interpretation and uncertainty on real microscopic inputs |
| SECCM | No current model | Transport, meniscus averaging, and concentration/standard-state conversion |

## Scientific context

The project is motivated by spatial MHC-DOS modeling that uses fitted local
prefactors ([Babar & Viswanathan, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11284846/))
and evidence that electrode screening modifies reorganization energy
([Maroo et al., 2026](https://www.nature.com/articles/s41586-026-10311-2)).
This implementation does **not** reproduce either paper's results.

The [research plan](docs/research-plan.md) describes the remaining calculations
and evidence required before interpreting a result as a microscopic explanation.
Upstream code is not copied or vendored here.
