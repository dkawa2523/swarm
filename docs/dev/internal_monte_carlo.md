# Internal Monte Carlo Backend

The internal `monte_carlo` backend is a lightweight product backend for
electron swarm comparisons. It is not an embedded MCIG clone and does not emit
validated bulk transport.

Its implementation is contained in `electron_swarm/solvers/monte_carlo`.
`solver.py` is the public adapter, `batch.py` owns independent multi-case execution,
`setup.py` validates and prepares static run state, and `case.py` coordinates
one-case execution. `case_state.py` owns mutable state and initialization,
`case_phases.py` owns fixed-population, weighted, and tail state transitions,
and `case_result.py` owns observable and result assembly. A direct batch
requires an explicit base seed and derives a separate stable uint32 seed from
that base and each E/N value. Case streams therefore do not depend on list
order or on whether neighboring anchors are present. Workflow replica seeds
are independently derived from the workflow base seed, mixture, E/N, replica,
and solver identity before the solver derives its case stream. The combined
two-stage seed schema and both seed values are recorded with the result.
Collision, orbit, population, histogram, rate,
weighted-runtime, and compiled-kernel modules own their corresponding numerical
responsibilities. Fixed-population transport is implemented in
`direct_transport.py`; synchronized weighted-growth transport is implemented in
`weighted_transport.py`; only their result value type and weighted-mean formula
are shared through `transport_common.py`.

The dependency direction is `solver -> batch -> case -> state/phases/result`,
with case internals and `setup` depending only on lower-level numerical owners.
Those owners do not import the public adapter, and moved private helpers are not
re-exported. The
package exports only `MonteCarloSolver`; workflow code may additionally consume
the stable public evidence contract in `evidence.py` without importing private
execution modules.

## Angular Scattering

Supported internal samplers:

- `isotropic`
- `maxent_p1`

`maxent_p1` requires both total/effective elastic and momentum-transfer cross
sections for every active scattering species. If either side is missing, the
backend fails explicitly instead of falling back to isotropic scattering.

## EEDF Estimator

The EEDF is a normalized density in `1/eV`:

```text
sum(eedf * energy_width_eV) == 1
```

Its finite-volume grid is part of the estimator contract.  From zero to
20 eV it is uniform in electron speed, giving quadratic energy edges and
sub-eV cells much narrower than the high-energy cells.  Above 20 eV the cell
width grows smoothly with a bounded relative width.  Every physical reaction
threshold plus the support transitions and extrema of each active
`sigma(E) sqrt(E)` kernel is retained as a cell boundary.  This construction
resolves low-field EEDF gradients and narrow rate kernels without importing a
deterministic distribution or making the grid depend on arbitrary LXCat table
point density.  Workflow provenance records the EEDF-estimator schema and
rejects reuse of a database produced by a different or unrecorded grid.

For a DC, zero-magnetic-field flight, internal MC partitions EEDF residence at
exact energy-bin crossings and integrates every piecewise-linear
`sigma(E) v(E)` kernel after splitting at cross-section knots and thresholds.
Magnetic flights retain an endpoint-energy rule that preserves magnetic no-work
behavior. The canonical EEDF CSV includes bin widths plus nullable
`sample_count`, `effective_sample_count`, and `relative_standard_error`.

When a reaction rate remains unresolved, the independent tail ensemble keeps
exact statistical weight in every occupied energy stratum and allocates its
exploration particles from the independently normalized reaction kernels as
well as the advancing threshold front. No deterministic solver distribution is
imported or blended into the MC estimate.

Tail sampling is enabled by an explicit Monte Carlo `tail_max_collisions`
budget. It is independent of `physics.energy_grid_policy.threshold_refinement`,
which belongs to deterministic energy-grid construction. A configured MC tail
budget can therefore never be silently disabled by a grid setting from another
solver mode. Omitting the budget (or setting it to YAML `null`) disables the
tail phase and records a configured tail budget of zero in execution
provenance.
The solve plan records the request as `configured`; after production, each case
records exactly one actual state: `disabled`, `configured_not_triggered`, or
`executed`.

## Transport Meaning

Internal MC reports flux-like particle-tracking estimates:

- drift: total weighted displacement over total weighted time
- longitudinal diffusion: residual axial variance over total weighted time
- transverse diffusion: transverse variance over total weighted time

These are not bulk/source-gradient coefficients. In ionizing swarms,
fixed-particle and weighted-branching population models represent different
observables and should not be mixed without checking `transport_definition`.

For low-field current response, `transport_estimator: paired_field_parity`
runs positive- and negative-field legs from identical bit-generator states.
It averages only the field-aligned odd responses (drift, mobility, and energy
mobility). The positive-field leg remains the sole owner of EEDF, reaction
rates, mean energy, and diffusion, and its tail budget is unchanged. Loss of
path correlation after a stochastic branch can reduce the variance benefit.
Under field-inversion symmetry, the aligned average adds no algebraic bias;
finite-warmup bias in either leg remains subject to the stationarity gate. This mode is
defined only for homogeneous DC calculations with `B=0`; it is never selected
from an internal E/N threshold.

Configured warmup and production budgets apply to each field leg. Thus the
mode preserves the complete positive-field sampling evidence and adds one
transport-only mirror leg, approximately doubling the orbit work for the
explicitly selected case. The case RNG advances through the positive leg
only, preserving the primary case result relative to `single_field`.

## Development Audit

Detailed MC diagnostics are available through `tools/benchmarks/benchmark_internal_mc_audit.py`
and `tools/benchmarks/benchmark_internal_mc_physics_triage.py`. Product summary outputs keep
only interpretation metadata.
