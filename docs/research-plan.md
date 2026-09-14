# From a runnable pipeline to a research result

## Question

For a fixed redox species, potential, temperature, and molecular sampling protocol,
how does the predicted AA/AB rate contrast change when the direct DOS, local
reorganization energy, and coupling are replaced independently?

The current repository implements the analysis and tests its mathematical
invariants. It does not establish that the proposed combination is novel.

## Stage 1: independently calculated structure and LDOS

Prepare relaxed cells near 0.8°, 1.1°, 1.5°, and 2.0°. Record the actual
commensurate angles rather than relabeling them as exact requested angles.
Archive the cell, atom positions, boundary conditions, force field, and relaxation
residuals. Derive the local displacement field and inspect AA/AB/SP assignments.

Compute LDOS from a chosen Hamiltonian on the same sites, converging the basis,
k-point grid, energy grid, and broadening. Record spin/valley/layer projections.
Test whether spatially integrated LDOS reproduces the total DOS. A nearest-registry
label does not substitute for atomic structural information.

Candidate upstream projects from the supplied brief:

- [Twister](https://github.com/qtm-iisc/Twister)
- [kp_tblg](https://github.com/stcarr/kp_tblg)
- [pytwist](https://github.com/sturk111/pytwist)

These are integration targets, not installed or validated dependencies.

## Stage 2: baseline and experimental comparison

Import local DOS and use fixed reference lambda and coupling to establish an
MHC-DOS baseline. Compare a matching normalization and potential convention with
[ElectrochemicalKinetics.jl](https://github.com/BattModels/ElectrochemicalKinetics.jl).
Preserve the original measured data, uncertainty, digitization procedure if any,
and area weighting before comparing to experiment. Do not encode a desired
order-of-magnitude discrepancy as a synthetic observation.

[Babar & Viswanathan, 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11284846/)
combines local DOS-derived kinetics with SECCM transport and fitted spatial
prefactors. This repo has no transport solver, so its effective rates should not
be compared directly to nanopipette currents.

## Stage 3: screening and reorganization

Implement and validate the published nonlocal dielectric screening treatment, or
obtain local reorganization energies from constant-potential sampling. Compare
the metallic and low-DOS limits and check dependence on ion distance, dielectric
parameters and chemical potential. The current surrogate is only a sensitivity
exercise and must be replaced for a claim about microscopic screening.

[Maroo et al., 2026](https://www.nature.com/articles/s41586-026-10311-2)
connects electrode electronic structure to reorganization through screening. Its
methods are a starting point for this validation; the local extension to a
spatially inhomogeneous moire cell must be justified independently.

## Stage 4: microscopic coupling

For representative AA, AB, and SP environments, calculate diabatic couplings for
the same redox species with a controlled distribution of molecular positions,
orientations, charge states, solvent configurations, and distances. Record basis,
functional, charge localization constraints, and finite-size corrections.

Candidate references from the brief include the
[CP2K CDFT tutorial](https://github.com/nholmber/cp2k-cdft-tutorial) and
[graphene charge-transfer workflow](https://github.com/Cheng75913/MP11-Graphene_Charge_Transfer).
This repository does not include validated input decks for either.

Check whether energy-independent H is adequate and whether couplings remain in
the weak-coupling regime. Average |H|² over appropriate configurations rather
than squaring an averaged signed matrix element.

## Stage 5: inference and uncertainty

Run the imported datasets through the eight-hybrid analysis. Compare all ordered
attributions and the Shapley summary. Preserve dependence between sampled DOS,
lambda, and coupling; independent parameter resampling can break their physical
correlations. Repeat across independent geometries and equilibrated trajectory
blocks. Report uncertainty on the rate ratio and attribution, not only input means.

Account for area weights before comparing domain averages; the rate evaluated at
an averaged DOS/lambda/H does not generally equal the average local rate.

## Evidence required before claiming completion of the scientific project

- Converged relaxed structures and computed local spectra with traceable provenance.
- Independently determined lambda and coupling, with uncertainty and physical limits checked.
- A quantitatively justified experimental comparison at matched conditions.
- Sensitivity to attribution convention, sampling, numerical resolution, and omitted physics.
- A literature review supporting any novelty claim.

None of those evidence requirements is satisfied merely by running the synthetic demo.
