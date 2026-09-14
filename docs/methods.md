# Methods, conventions, and limitations

## Rate integral

For an energy-independent coupling magnitude H in the nonadiabatic limit:

\[
k_{\rm red}=\frac{2\pi}{\hbar}|H|^2 A
\int D(E) f(E-\mu)
\frac{\exp[-(\lambda+\eta-(E-\mu))^2/(4\lambda k_BT)]}
{\sqrt{4\pi\lambda k_BT}}\,dE.
\]

E, mu, lambda, H, and kBT use eV. Eta is supplied in volts; for a one-electron
reaction its numerical value is the energy shift in eV. **Positive eta favors
oxidation.** The oxidation expression replaces f by (1-f) and the squared
gap by `(lambda - eta + E - mu)^2`.

For the same inputs the integrands obey `k_red/k_ox = exp(-eta/kBT)` even for
asymmetric DOS. This identity is a regression test. At eta=0 the forward and
reverse rate coefficients agree. These are rate coefficients, not net current;
redox activities have not been included.

The DOS is in states/eV/nm², including whatever spin/valley degeneracy is recorded
by the upstream calculation. The effective coupling area A is explicitly in nm².
The DOS/area product is a local channel density, giving effective per-molecule
rates in s⁻¹. A is a model convention, **not a pipette area**. H and DOS must use
compatible normalizations; a molecule-to-aggregate coupling must not be combined
with an aggregate DOS that counts the same states a second time. Comparing to
heterogeneous rate constants in cm/s requires a separately justified conversion.

E and mu share a fixed reference. The code does not solve quantum capacitance,
electrostatic potential division, self-consistent doping, or a Frumkin correction.
At each eta/mu it uses the supplied DOS, lambda, and H unchanged. For potential-
dependent inputs, rerun with independently calculated datasets for each potential.

## Numerics

Trapezoidal quadrature supports nonuniform strictly increasing energy grids.
The summation uses log-space accumulation and stable Fermi factors. No DOS
normalization, smoothing, interpolation, extrapolation, or negative-DOS clipping
is silently performed. A zero DOS/coupling yields zero rate at the core API level;
the attribution pipeline requires positive rates for all counterfactual hybrids.

The report compares the full and coarsened grids. This is a sensitivity check,
not proof of convergence: peaks missing on both grids remain missing. Repeat
upstream calculations on denser k meshes, with varied spectral broadening and
wider energy windows. The supplied finite energy window is never automatically
extended. The reference and target must use the same full energy grid.

## Attribution

Let the three factors be D, lambda, and H. Define v(S) as the natural logarithm
of the rate with target values for factors in S and reference values otherwise.
Compute every subset S. For each of the six permutations, replace one factor at
a time; average each factor's marginal log-rate change over those permutations.

These are exact three-variable Shapley values. They sum to v(all)-v(empty).
Because coupling enters only as |H|², its contribution is exactly
`2 ln(H_target/H_reference)`; tests check this identity. The other two shares
depend on their interaction in the integral. All ordered paths are exported,
so a reader can inspect that dependence instead of relying on a single ordering.

If lambda is a function of DOS, this operation **freezes** lambda while replacing
DOS. It separates direct phase-space and mediated reorganization effects by
choice, not by an independently realizable intervention. Shapley values do not
prove a causal mechanism, and there is no universal ordering-independent physical
product of three unique enhancement factors.

## Screening and microscopic extraction

The synthetic demonstration uses
`lambda = lambda_metal + penalty/(1 + D_thermal/D_scale)`, where D_thermal is the
integral of DOS against minus the Fermi derivative. The 0.82 eV metal baseline
is motivated by the RuHex discussion in
[Maroo et al., 2026](https://www.nature.com/articles/s41586-026-10311-2).
The penalty, scale, and spatial dependence are illustrative and unfitted.
**This rational interpolation is not that paper's nonlocal dielectric theory.**

For scientific analysis, supply a lambda field computed independently. The
`gap_reorganization` utility estimates lambda as half the difference in mean
vertical gaps from oxidized and reduced ensembles, using `E_reduced-E_oxidized`
in both. It also reports variance/(2kBT) in each ensemble as a linear-response
diagnostic. It does not correct autocorrelation, assess equilibration, or provide
confidence intervals.

`two_state_coupling` uses `(H_ab - S_ab*(H_aa+H_bb)/2)/(1-S_ab^2)` for two real,
normalized diabatic states, returning its magnitude. A test compares it against
matrix Löwdin orthogonalization. Near-linearly-dependent states are ill-conditioned
and must be rejected by the upstream physical analysis even when |S|<1.

## Structure

The registry classifier uses fractional displacement in a 60° graphene basis,
with periodic distances to AA, AB/BA, and three SP representatives. It returns the
nearest ideal registry, not a domain-wall width, local strain, or twist estimate.
The ideal moire period is `0.246 nm / (2 sin(theta/2))`. Synthetic plots use a rigid
displacement grid; no reconstruction/relaxation solver has been invoked.

## Physical assumptions

- One-electron outer-sphere transfer and weak, energy-independent coupling.
- Classical Gaussian nuclear energy-gap statistics with positive lambda.
- Independently supplied local parameters evaluated at a common temperature and potential.
- No solvent dynamics, adiabatic crossover, many-body corrections, transport, or fitted experimental prefactors.

The computational pipeline is usable now. Testing the proposed scientific
hypothesis requires the upstream calculations and validation described in the
[research plan](research-plan.md).
