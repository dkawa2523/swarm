# Multi-term Boltzmann Entry-Point Design

This document describes the integrated `multiterm_boltzmann` entry point.
It absorbs the useful design from `swarm_bolsig_beyond_solver.zip` into the
existing `electron_swarm` package instead of keeping a separate `swarm_beyond`
namespace.

## Scope

The integrated entry point keeps the shared solver scope narrow:

- uniform DC electric field
- B = 0
- axisymmetric m = 0 Legendre expansion
- flux transport in the shared result schema
- lmax>1 coefficients as a reference-anchored integral-cross-section closure
- shared `electron_swarm` YAML input and CSV output schema

The default `moment_closure` path remains an estimate and records
`meta_physical_validity=moment_closure_estimate_not_multiterm_operator`.
`method=operator, lmax=1` delegates to the two-term reference path, while
`method=operator, lmax>1` uses a reference-anchored closure only after explicit
`allow_experimental_operator: true` opt-in. `method=hybrid` still fails fast.

## Transport Semantics

Flat legacy fields keep their historical meaning as flux coefficients:

- `drift_velocity_m_s`
- `mobility_m2_V_s`
- `diffusion_L_m2_s`
- `diffusion_T_m2_s`

The typed result model carries a `TransportSet`:

- `transport.flux`: spatial average of particle velocity and flux diffusion.
- `transport.bulk`: swarm center-of-mass drift and bulk diffusion only when a
  validated hydrodynamic or MC bulk mode is available.
- `transport.source`: source-gradient corrections, including
  `gradient_velocity_m_s = W_bulk - W_flux`, only when bulk exists.

Moment-closure does not populate standard bulk/source-gradient transport. Its
reference values remain isolated as estimates:

- `meta_estimated_bulk_drift_velocity_m_s`
- `meta_estimated_bulk_reduced_diffusion_L_m-1_s-1`
- `meta_estimated_source_gradient_velocity_m_s`

## Input and Data Flow

`multiterm_boltzmann` uses the existing unified cross-section loader:

1. `electron_swarm.core.config.load_config()` reads YAML.
2. `electron_swarm.core.cross_sections.load_cross_sections()` loads canonical
   CSV or legacy LXCat/BOLSIG-like text.
3. `cross_sections.high_energy_extrapolation` controls behavior above the
   tabulated energy range. The default examples use `zero`; `hold` is available
   for legacy comparisons, and `error` is useful for strict benchmark runs.
4. `build_multiterm_case()` converts the shared config and `CrossSectionSet`
   into the internal energy grid, number density, and electric field.
5. Cross sections, velocities, collision frequencies, and signed energy losses
   are projected once onto the energy grid.
6. `MultiTermBoltzmannSolver.solve_case()` returns `SwarmCaseResult`, so the
   existing writers, plots, and COMSOL export path remain reusable.

The standalone parser from the zip package is not part of the official input
surface for this repository.

Internal modules keep extension points small:

- `grid.py`: energy grids, quadrature widths, electron speeds.
- `models.py`: solver-internal case, solution, and rate containers.
- `projection.py`: cross-section normalization and grid projection.
- `closure.py`: moment-closure estimate backend.
- `operator.py`: operator backend surface, `lmax=1` two-term compatibility
  path, lmax>1 reference-anchored closure, and comparison helpers.

## Numerical Model

Moment-closure computes a normalized EEDF from a power-balance closure using
either a Maxwellian or Druyvesteyn shape. Reaction rates are evaluated by
convolution:

```text
k_r = integral sigma_r(epsilon) v(epsilon) F(epsilon) d epsilon
```

The flux drift uses an effective momentum relaxation frequency. Superelastic
processes are treated as energy gains in the power balance; attachment has no
energy loss unless explicit process metadata is added in the future. The power
balance must bracket or minimize to a small residual; otherwise the solver
fails instead of returning a boundary estimate.

The current lmax>1 path does not emit standard bulk/source-gradient
coefficients. If `hydrodynamic: true` is requested, the result records that the
request was not computed.

Process normalization is applied before projection. If a species has an
`EFFECTIVE` cross section, same-species `ELASTIC`/`MOMENTUM` cross sections are
not added again to the transport momentum frequency. They remain available for
rate tables, but the normalization warning is recorded so users can inspect
input ambiguity.

`method=operator, lmax=1` delegates to the existing two-term
Scharfetter-Gummel implementation and adapts the result to the multi-term
result schema with `meta_multiterm_method_used=operator_lmax1_two_term`.
`method=operator, lmax>1` is executable as
`meta_multiterm_method_used=operator_reference_anchored_lmax_gt1`. It requires
`allow_experimental_operator: true`. This path anchors f0, rates, and flux
transport to the native two-term Scharfetter-Gummel reference, then adds bounded
higher-Legendre coefficients using the available integral momentum relaxation.
Runs are marked
`meta_physical_validity=reference_anchored_lmax_gt1_integral_cross_section_closure`.
This is intentionally not a full multi-term differential-scattering solution;
it prevents the previous cold, nonphysical lmax>1 eigenmode from entering
transport tables while keeping a clean extension point for future physics.

A true `lmax>1` implementation needs information that an LXCat-style integral
cross-section table usually does not contain:

- Legendre moments of the differential elastic scattering kernel, not only an
  elastic or momentum-transfer integral.
- Consistent anisotropic source/sink blocks for excitation, ionization,
  attachment, and superelastic channels.
- Conservative energy-space discretization that reduces exactly to the native
  two-term Scharfetter-Gummel block for `lmax=1`.
- A finite-k or equivalent hydrodynamic benchmark before bulk diffusion and
  source-gradient coefficients are published as standard transport data.

The shared boundary is split into two public native helpers:
`BoltzmannTwoTermSolver.solve_native_distribution()` returns the native energy
grid, quadrature widths, normalized EEDF, diagnostics, and metadata before
public CSV/result postprocessing. `assemble_native_operator_block()` returns
the Scharfetter-Gummel sparse matrix, projected collision data, electric field,
and density on the same grid. Future multi-term operator assembly should reuse
that block as the `lmax=1` regression anchor. A future true `l>1` solver should
be added only after the input model carries the missing higher-order collision
moments; until then there is no public direct sparse Legendre block assembly.

Operator runs also emit compact quality indicators rather than a separate
benchmark contract: `meta_operator_tail_probability`,
`meta_operator_tail_rate_fraction`, `meta_operator_highest_l_relative_l1`, and
`meta_operator_lmax_convergence_ok`. The last value uses the existing
`multiterm_boltzmann.lmax_convergence_tolerance` setting and is intended as a
screening flag for whether a higher `lmax` check is worth running.
The collision model is reported with
`meta_operator_ionization_source_model` and
`meta_operator_l_gt_0_elastic_model` and
`meta_operator_l_gt_0_inelastic_model`. In the current integral-cross-section
scope, `l=0` ionization follows the native two-term
`ionization_energy_sharing` setting. Higher-order elastic damping is reported
as `reference_anchored_integral_momentum_closure`, and higher-order inelastic
terms use a sink-only model with isotropic `l=0` source closure.

## Extension Path

The next implementation stages should be:

1. Keep `operator.py`'s `lmax=1` two-term compatibility path as the regression
   baseline.
2. Keep lmax>1 public results anchored to the native two-term reference while
   only integral cross sections are available.
3. Use `assemble_native_operator_block()` as the reusable Scharfetter-Gummel
   block boundary for future direct block assembly.
4. Add differential sparse Legendre collision/source blocks only after the input
   model includes higher-order scattering moments.
5. Benchmark finite-k diffusion extraction over Ar, Ar/N2, attachment, and
   superelastic cases before exposing bulk/source-gradient output.
6. Use `tools/validate_multiterm_operator.py --quick` only as a compact
   executable smoke sweep for Ar, Ar/N2, attachment, superelastic, and lmax
   trends. It is not a physics validation gate.
7. Use `tools/benchmark_operator_gate.py --quick` as the small reference gate.
   BOLOS is optional unless `--require-bolos` is passed; MC is included only
   with `--with-mc`.
