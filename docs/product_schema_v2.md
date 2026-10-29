# Product Schema v2

Product YAML files must set `schema_version: 2` and select solvers with
`run.solvers`. Each solver entry is a mapping with `id` and optional `enabled`;
string entries and labels are not part of the product schema.

Canonical solver ids:

- `two_term`
- `multi_term`
- `monte_carlo`
- `propagator`

Product solver configuration is intentionally small. `two_term` accepts only
`backend: native_sg` and `min_momentum_cross_section_m2`. Ionization sharing
belongs to `physics.ionization.energy_sharing`; `loss_only` explicitly selects
the single-daughter source. `multi_term` accepts only `method` and `lmax`. `monte_carlo` accepts
only sampling controls listed below. `propagator` accepts its
bounded grid and steady-solve controls. Experimental backend selectors, PN
formulation knobs, and convergence heuristics are development details, not
public YAML fields.

Duplicate keys, merge keys, numeric strings, and booleans in numeric fields are
rejected. Mixture fractions must sum to one; workflow loading does not normalize
arbitrary weights. Workflow `mc` settings require an enabled MC solver, and
deterministic execution settings require an enabled deterministic solver.
`conditions.length_scale_m` has no homogeneous-swarm meaning and is not accepted.
DC requests cannot specify a time-dependent frequency.

Gas conditions use exactly one of pressure in Pa or number density in `m^-3`.
Pressure and temperature determine density through the ideal-gas relation;
an explicitly supplied number density is used directly. `gas_mixture` defines
component fractions and explicit masses in amu; every component requires
`mass_amu`. The old single-gas `conditions.species`/`mass_amu` shorthand is
removed. Cross-section masses must agree with the configured species within
relative rounding tolerance 5e-4; every operator then receives the exact
configured mass. A material mismatch fails before solving. The resolved
physical context records density, temperature, mixture, collision-data digest,
fields, population/source physics, and e-e model independently of solver choice.


The common entry is `swarm-workflow run <config.yaml>` for both direct runs and
campaigns. [Input and execution guidance](swarm_inputs.md) defines the optional
`campaign` section and aligned mixture-fraction arrays. The old external base
reference is removed. E/N belongs to `run.e_over_n_Td`, composition to
`conditions.gas_mixture`, and the master seed to `solvers.monte_carlo.seed`.
Direct single-mixture calculation YAML remains accepted by the low-level
`electron-swarm` CLI and core Python API. Campaign comparison output is
unsupported and is rejected during planning rather than silently discarded.

## Multi-Term Methods

- `pn_closure_direct`: a B=0/DC/axisymmetric direct PN solve for every
  `lmax >= 1`. All adjacent moments are coupled in both directions. Ordinary
  integral cross sections supply an explicit angular closure, not DCS data.
- `pn_dcs`: runnable with `physics.angular_scattering.model: moment_table`.
  It uses normalized Legendre moments from the table in the direct PN angular
  damping. Ordinary XS closure alone cannot run `pn_dcs`.

Ordinary-integral-XS `multi_term` results record
`ordinary_integral_xs_closure=true` and `exact_dcs_based=false`.
Multi-term drift is the solved F1 velocity moment; diffusion is an F0-gradient
approximation. The corresponding transport definition is
`pn_f1_flux_drift_f0_gradient_diffusion`.

## Angular Scattering

Valid model/closure pairs:

- `isotropic` / `zero`
- `momentum_power` / `power`
- `maxent_p1` / `maxent`
- `moment_table` / `table`

Moment-table input is a wide CSV with `energy_eV,m0,m1,...`. Energy must be
nonnegative and strictly increasing, `m0` must be 1, and higher moments must be
finite values in `[-1, 1]`. The complete supplied sequence must be realizable
by a nonnegative angular measure on `[-1, 1]`; range checks alone are insufficient.
No extrapolation is performed.
`exact_dcs_based=true` applies only when `provenance: dcs_derived`.

Solver-plan support describes implementation behavior, not physical-model
fidelity. A supported selected closure therefore does not make a run degraded.
The plan and case results separately record `angular_scattering_fidelity` and
`angular_scattering_assumption`; ordinary cross sections are always identified
as `integral_xs_closure`, including a fully supported isotropic kernel. Only a
moment table whose provenance is `dcs_derived` is identified as DCS-derived.

```yaml
physics:
  angular_scattering:
    model: moment_table
    higher_moment_closure: table
    moment_table:
      path: moment_tables/argon_demo_moments.csv
      format: normalized_legendre_moments
      provenance: model_derived
      extrapolation: error
```

## Monte Carlo

Product `monte_carlo` always means the limited internal particle backend.
External BOLSIG+, MCIG, and external MC ingest are benchmark/reference tooling
concerns and are not accepted by product YAML.

The internal MC sampler uses and reports the angular model requested under
`physics.angular_scattering`. It supports `isotropic` and `maxent_p1` angular
samplers plus DC magnetic dynamics through a Boris pusher. Unsupported angular
samplers are handled by `feature_policy` and are never silently approximated.
Attachment trajectory loss is not implemented for either population model;
an active attachment cross section therefore fails or skips that solver through
`feature_policy` before particle evolution. Propagator superelastic requests
are handled in the same way. Planning uses only active mixture species.
`maxent_p1` requires both explicit total elastic and momentum-transfer cross
sections for each active scattering species. Magnetic trajectories are
substepped only inside a sampled null-collision trial interval; collision
acceptance still occurs once at the trial event.

Public `solvers.monte_carlo` fields are limited to `seed`,
`particles`, `production_trial_periods`, `warmup_trial_periods`, `tail_trial_periods`,
`tail_rate_rse_trigger`, `transport_correlation_lag_barriers`,
`transport_estimator`, and `numeric_kernel`. `numeric_kernel: auto` uses the
compiled full-barrier kernel for supported DC/B=0/isotropic/global-clock
weighted branching and records a reference-kernel fallback for other models.
An explicit unsupported `numeric_kernel: numba` request fails. The default
collision clock uses a global marked majorant: stationary-target `v sigma`
groups for inelastic channels and Maxwellian relative-speed `g sigma` groups
for elastic species. Elastic events use exact center-of-mass binary kinematics.
Population physics is selected by `physics.ionization.population_model`.
The default `weighted_branching` uses a bounded two-daughter temporal-growth
model shared with the deterministic solvers. `fixed_particle_single_daughter`
selects fixed-population MC; the deterministic solvers explicitly reject or
skip that request under `feature_policy`. The old `solvers.monte_carlo.population_model`
key is removed with a migration error. Version 0.4 MC configs that relied on
its implicit fixed-population default must now request
`physics.ionization.population_model: fixed_particle_single_daughter` explicitly
to retain that physical experiment; otherwise the 0.5 temporal-growth default
applies.
`energy_sharing: loss_only` explicitly retains one daughter and gives no
ionization population growth; the reported daughter count is one. The
weighted MC implementation rejects that combination as unsupported, while
supported deterministic treatments retain their explicit plan status.

The two daughter energies follow `physics.ionization.energy_sharing`; the
ionization threshold is always recorded as an energy loss. When the particle
ensemble grows past twice `particles`, an emergency cap resamples it back to
`particles`. Weighted production additionally synchronizes and normalizes the
population at every common-time barrier, resampling there when needed. Both
operations preserve total weight.

The three `*_trial_periods` budgets count common-time null-collision trial
periods, not accepted physical collisions. Defaults are 256 particles, zero
warmup, 100 production periods, and no tail phase. These defaults and the
sampling unit are visible in `--plan` and the run manifest. The old names
`max_collisions`, `warmup_collisions`, and `tail_max_collisions` are rejected
with a migration error. Budgeting paired parity includes two production legs
and, when configured, one conditional tail leg.

Direct product execution requires an explicit nonnegative `seed`. It is a base
identity: each E/N case receives a deterministic SHA-256-derived uint32 stream,
so reordering or running an anchor alone does not change its samples. Canonical
summary rows record both the configured base seed and the effective case seed.
Campaign replicas first derive an independent seed from `solvers.monte_carlo.seed`, mixture,
E/N, replica index, and solver id; the stored seed-schema version identifies
that workflow stage together with the solver's case-seed stage.

For zero-magnetic-field flights, internal MC partitions EEDF residence at exact
energy-bin crossings and integrates piecewise-linear reaction kernels after
splitting at cross-section knots and thresholds. Magnetic flights use the
endpoint-energy rule of the current Boris trajectory implementation. Product metadata reports the observable definition with
`transport_definition` plus solver-comparison interpretation fields.
Population-model details such as
`swarm_population_treatment`, branching counters, and estimator diagnostics are
available only when benchmark tools request diagnostic collection; they are not
part of normal product summary metadata. Weighted branching is still a
flux-like particle-tracking estimator; validated bulk transport remains out of
scope.

For `fixed_particle_single_daughter`, MC particle/energy flux transport is
measured on common physical-time planes with multiple finite-lag time origins.
The product reports direct energy mobility plus longitudinal and transverse
energy diffusivity. The scalar energy-diffusion field remains unset when the
two components differ; consumers must select an explicit tensor projection.
For `weighted_branching`, the v6 configured-lag estimator reports
direct residence energy mobility and a restricted density-packet energy-current
diffusion tensor. A single `transport_correlation_lag_barriers: L` derives
samples at `L/4,L/2,L,2L`, blocks of `2L`, lineage ESS gates at `L/2,L`,
and production-equivalent early/late residence estimates. Every replica must
retain at least `max(128, 0.05*particles)` effective lineages at both hard lags;
`2L` is supplemental and is included only when it meets the same threshold.
This is suitable for the declared standard-local-
energy projection after ensemble qualification; it does not identify the
independent full two-gradient energy response matrix.

`transport_estimator` defaults to `single_field`. For homogeneous DC, B=0
cases, `paired_field_parity` starts positive- and negative-field legs from the
same random-generator state and averages their field-aligned drift, mobility,
energy mobility, and corresponding stationarity windows. EEDF, reaction rates,
mean energy, and diffusion remain owned by the complete positive-field leg.
Configured warmup and production counts apply to each field leg, so this mode
has two transport legs; the negative-field leg does not run rare-tail sampling.

The external sweep workflow may assign more sampling to statistically difficult
E/N anchors without changing the product solver schema.  `campaign.mc.sampling_plan` is
a list of complete per-anchor rows:

```yaml
campaign:
  mc:
    e_over_n_Td: [100.0, 750.0, 3000.0]
    replicas: 12
    sampling_plan:
      - e_over_n_Td: 3000.0
        particles: 2048
        warmup_trial_periods: 2000
        production_trial_periods: 32768
        tail_trial_periods: 16384
        replicas: 16
        transport_correlation_lag_barriers: 64
        transport_estimator: single_field
```

Every row must name an existing MC anchor. `particles`, `production_trial_periods`, and
`replicas` are positive integers; `warmup_trial_periods` is nonnegative. An
explicit `tail_trial_periods` value is positive. When no tail budget is set in
the base solver or a row explicitly uses YAML `null`, the resolved canonical
plan records zero and no tail phase is run. The
correlation lag is optional in a workflow row and uses
the base solver value when omitted; when present it must be a power of two at
least 4. `transport_estimator` is likewise optional per row and inherits the
base solver value; it accepts `single_field` or `paired_field_parity`. An anchor
without a row uses `solvers.monte_carlo` sampling values
from the same calculation document and `campaign.mc.replicas`.  The workflow resolves all
anchors into one canonical plan, stores it in database provenance, and refuses
resume when that plan changes. The lag is also a public
`solvers.monte_carlo.transport_correlation_lag_barriers` field for direct
product runs. An enabled Monte Carlo workflow also
requires the same explicit nonnegative `solvers.monte_carlo.seed`; independent replicas are never
qualified from implicit or duplicated random streams.  COMSOL-ready table
generation requires the canonical plan, verifies it against every stored
replica's run settings, and records the plan plus its SHA-256.  A passing
`aggregate_quality` row alone is not MC transport qualification: the
solver-specific convergence evidence in `build-tables` must also pass.
When `campaign.mc.convergence` is enabled,
`maximum_particle_barriers_per_replica` and
`maximum_total_particle_barriers` bound the nominal warmup, production, and
conditional-tail plan before workers start. Follow-up time or replica
extensions must fit the same bounds.
`build-tables --mc-qualification-profile full_transport` applies that complete
transport gate. The explicit `function_eedf_restricted_lmea` profile instead
gates weighted-growth v6 lineage/population evidence, EEDF/rate-tail
consistency, direct elastic loss, and the active mobility. Distribution
stationarity is measured by the paired 95% bound of the mean energy. Mobility
uses the paired early/late mean log ratio for systematic drift and the
Student-t 95% half-width across independent production replicas for coefficient
precision. This avoids counting random early/late window noise twice while
retaining a 10% bound on both drift and reported mobility. Particle/energy
diffusion and energy-mobility failures remain in the same quality artifact and
full-transport metadata; they are not relabeled, repaired, or used to remove
an otherwise valid active anchor.
For weighted-growth qualification, warmup must in fact be positive and
`production_trial_periods` must be an integer multiple of `2L` containing at least 128
complete blocks. In this population model the name `production_trial_periods` denotes
common-time null-collision trial periods/barriers, not the realized number of
accepted physical collisions.

Every workflow records the selected solvers' runtime source hashes and numerical
library versions. Resume and table generation require matching implementations,
base physics and collision/angular input hashes. Databases predating this
identity must be recomputed; changing metadata cannot qualify old results.

For an extended sweep, `mc.reuse_database` may name an existing workflow
SQLite database. The sweep imports only raw MC cases with the same base
physics, cross sections, estimator version, mixture, sampling controls, and
derived replica seed; aggregate results are always rebuilt in the destination.
This lets a new E/N range compute only genuinely new anchors without combining
incompatible Monte Carlo evidence.

Deterministic workflow sweeps remain single-process by default. Parallel
execution is an external workflow concern, not a product solver option:

```yaml
execution:
  deterministic:
    workers: 8
    global_memory_budget_mb: 4096
```

The workflow caps the request by available CPUs and a conservative per-worker
memory reservation. The propagator reservation is its configured
`solvers.propagator.max_memory_mb` plus process overhead; solvers without a
public memory ceiling use
a conservative workflow reservation. E/N values are split into contiguous
same-solver chunks so continuation remains available inside each chunk.
Workers return results only; the parent process writes SQLite checkpoints in
workflow order. Changing this execution contract changes workflow provenance,
and worker failures are propagated without publishing a completed aggregate.

A deterministic workflow with one mixture and one enabled solver can request
bounded mean-energy support without inventing an EEDF row:

```yaml
mean_energy_support:
  required_max_mean_energy_eV: 35.45
  relative_guard: 0.10
  maximum_steps: 2
  maximum_e_over_n_Td: 5000.0
```

After the configured E/N sweep, the workflow inverts the calculated upper
`(E/N, mean energy)` envelope, runs the same solver at each proposed field, and
writes those complete cases before aggregation. Existing continuation anchors
are reused on resume. The sweep fails if the guarded target is still unmet at
the step or field bound. Monte Carlo is excluded because its support must be
extended through the replica campaign and sampling budget.

## Propagator

`propagator` is the P1 deterministic velocity-cell solver. Its public scope is
spatially homogeneous DC, `B=0`, positive E/N, no electron-electron
collisions, and no finite-k request. It holds a normalized population on an
energy x polar-angle grid. Acceleration and local angular collisions are solved
together as positive stationary shell responses; the velocity origin uses a
characteristic-aligned cut-cell response rather than a frozen spherical shell.

Public `solvers.propagator` fields are:

| Field | Default | Validation |
| --- | ---: | --- |
| `method` | `stationary_response` | this value only |
| `energy_cells` | 600 | integer in `[64, 4096]` |
| `polar_cells` | 72 | even integer in `[8, 360]` |
| `max_iterations` | 2000 | integer in `[10, 20000]` |
| `convergence_tolerance` | `1e-8` | finite value in `[1e-12, 1e-4]` |
| `max_memory_mb` | 1024 | integer at least 128 |

The energy grid uses a nested sinh-stretched speed coordinate over the
configured core support, inserts active reaction thresholds, and adds a
stretched tail to the configured ceiling in one construction without rescaling
the core. The map is linear in speed at the origin while allocating more cells
to the thermal and sub-eV region than a uniform-speed grid. The solver measures physical
outgoing characteristic escape at that boundary; it does not append, resolve,
or remap a converged distribution. Response applications, the energy ceiling,
and memory are hard bounds. A bounded failure raises an error and does not
switch scheme or substitute another solver.

The supported angular kernels are `isotropic` and `maxent_p1`. The latter
requires explicit total and momentum-transfer elastic cross sections for every
active species. A KL projection of the cell-integrated maxent phase law builds
a positive reciprocal source--destination joint measure whose discrete P1
moment is exact. The total elastic integral controls angular event loss/gain;
the momentum-transfer integral independently controls P1 relaxation and the
finite-temperature recoil operator. Elastic energy exchange uses a reversible
SG generator with exact Maxwell cell masses. The origin integrates both elastic
rates over radial bands, while regular shells analytically average their `1/v`
geometry coefficient and refine the response to a fixed internal error bound.
Inelastic transfers and reported rates share the same donor-cell quadrature.
P1 returns EEDF, energy-angle population, reaction rates, growth, flux
drift/mobility, and a same-operator elastic energy-loss moment. Bulk drift,
particle diffusion, and electron-energy mobility/diffusion remain unavailable.

The Propagator-only artifact plus the exact-P1/separated-rate operator
invariants pass the maxent-Ar, inelastic-Ar, mixed-gas, low-field, energy-boundary,
and medium--fine numerical-core gates. This qualifies `propagator` as one bounded
calculation model in its declared scope; acceptance does not depend on agreement
with another solver. Raw-DCS support, additional gases, strongly optically thick
synthetic origin limits, and P2 transport are scope extensions rather than P1
release blockers. The repository COMSOL mapping is a `physical_target` with
the same explicit restricted local-mean-energy closure boundary used by the
other solver routes. The GEC-CCP workflow separately binds its 300 x 48 table
grid to 600 x 72 checks over 3000--4000 Td before COMSOL preflight. See the fixed
[Propagator-only qualification artifact](dev/results/propagator_p1_core_qualification.json)
and the generic-schema
[GEC-CCP target artifact](dev/results/propagator_gec_ccp_target_qualification.json).

## Physics Features

- e-e: `none`, `relaxation_postprocess`, or `fp_energy`.
- ionization source: `equal`, `primary_secondary`, or `loss_only`.
- tail refinement: `physics.energy_grid_policy.max_eV_limit` bounds the
  complete reported cells for both Boltzmann solvers. With `adaptive: true`,
  two-term and multi-term extend the energy domain and re-solve up to this
  ceiling and a finite cycle budget. Multi-term re-solves all retained PN
  moments. Measured tail probability, target, and refinement cycles live in
  `case.diagnostics`; multi-term also records whether the target was met.
  Reaching a bound does not imply that the target was met. Propagator
  constructs its bounded tail before solving.
- finite-k: request is validated and policy-handled, but no solver emits
  validated finite-k coefficients yet.

Unsupported physics is handled by:

```yaml
feature_policy:
  unsupported: fail        # fail | skip_solver
  degraded: record         # fail | record
```

Unsupported physics is never approximated or ignored as a fallback. Degraded
support is either recorded in the solver plan or rejected with `degraded: fail`.

Workflow qualification additionally rejects two-term tables that reach the
configured energy ceiling, even when their measured tail is small. A valid
direct solver result therefore need not qualify as a COMSOL input.

Elastic energy relaxation requires an elastic momentum-transfer integral.
An explicit total elastic integral is sufficient only with an isotropic
closure or a supplied first angular moment. Effective momentum integrals
cannot identify elastic recoil. PN isotropic scattering requires matching
total and momentum integrals when both are supplied; use `maxent_p1` to
declare anisotropy. Two-term retains its momentum integral and reports the
resulting first-moment assumption, so comparisons can identify differing
angular physics even when both requests name `isotropic`.

## Outputs

Canonical product outputs:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_energy_angle_distribution.csv` when a solver returns one
- `<base>_rf_phase.csv` and `<base>_rf_phase_eedf.csv` for phase-resolved RF
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled
- `<base>_manifest.json`

The EEDF CSV always includes `energy_width_eV` plus nullable `sample_count`
and `weight_effective_count` columns. Boltzmann solvers leave these empty.
Internal Monte Carlo reports raw histogram counts and weight concentration;
these correlated samples do not establish a standard error. Independent
workflow replicas provide statistical uncertainty. The EEDF is a
normalized density in `1/eV`, so `sum(eedf * energy_width_eV) == 1` up to
floating-point tolerance.

The energy-angle CSV stores one row per energy/polar cell with representative
`energy_eV`, `mu`, cell widths `energy_width_eV`, `mu_width`, and
`energy_angle_density_eV_inv`. Multiplying density by both widths and summing
gives one. It is emitted by `propagator`; solvers that
do not return an angle-resolved distribution do not produce this file.

Summary metadata includes only product-facing interpretation fields: solver
method, angular model/source, ionization/e-e/magnetic/tail/transport treatment,
MC base/effective case seeds, and DCS/closure flags. For MC, tail treatment is
`disabled`, `configured_not_triggered`, or `executed`; the solve plan uses
`disabled` or `configured` because activation is known only after production.
Detailed grid, convolution, tail, e-e operator,
direct-PN, and internal-MC audit values are diagnostics or development-tool
outputs, not summary columns.

Solver-comparison rows include `angular_model_status`, `comparison_status`,
each solver's `numerical_status` and `transport_definition`, per-observable
availability status, scalar relative differences, and optional `eedf_l1_error`.
Stale transport has no relative difference. Direct MC coefficients are marked
as requiring independent ensemble qualification. A required comparison can still
require known, matching angular metadata, but external reference agreement is
benchmark tooling rather than a product schema contract.

Publication stages a complete generation, backs up previous managed files,
and commits the manifest last under an OS lock. A failed replacement restores
the previous generation. The manifest contains the executed configuration,
resolved solver settings, implementation/library identities, input hashes,
attribution status, quantity definitions, and each CSV's hash and typed schema.
`electron_swarm.io.read_product_run(path)` verifies those hashes and schemas
under the same lock, returning `.metadata` and `.tables` together. Text
identifiers and integer seeds retain their exact values; missing numeric values
remain missing. Applications needing a consistent generation should use this
reader; unrelated CSV reads cannot establish that guarantee.

For `solver_execution`, `outcome.status` describes execution coverage: `complete`, `partial`, or
`skipped`. When every solver is disabled or skipped, no simulation CSVs are
created, previous managed simulation CSVs are removed, and CLI exit code is 3.
The plan and manifest explain why no calculation ran. Partially executed requests return CLI exit code 4. Disabled entries have
plan state `disabled`, `requested=false`, and `skipped=false`; they do not
make an otherwise complete run partial. Comparison reports `available`, `partial`, `unavailable`, `incompatible`,
`unverified`, or `disabled`, with separate `coverage_status` and
`qualification_status` fields. Optional
missing pairs remain visible; required missing pairs fail the run. Completion
does not establish convergence, model agreement, or physical validation.

Input attribution is reported as `declared`, `unattributed`, `not_provided`,
or `hash_mismatch`. A matching sidecar with a source URL is a supplied citation,
not independent verification of the physical data. Numerical source-table
qualification and external experimental/reference validation are distinct.

Version 0.5 requires explicit masses and moves population physics out of
solver settings. Model identity (`swarm.numerical_model_sha256.v2`) tracks
numerical sources. Execution identity (`swarm.execution_pipeline_sha256.v1`)
tracks input interpretation, planning, and runtime processing. Output format
has its own `swarm.product_run.v3` schema. CLI/output-only edits leave numerical
qualification valid; changed model or execution identity requires new workflow
calculations. Model changes also require fresh numerical qualification. Product YAML remains schema v2 with the same four solver ids.
To recompute, set the workflow's `database` to a new SQLite filename and run
the sweep again. The old database remains evidence for its recorded source;
its cases cannot be merged into the new implementation's campaign.

The `output` block configures `directory`, `base_name`, optional
`figures: [eedf, transport, reaction_rates]`, and optional
`export: {solver: multi_term, quantities: [reaction_rates]}`. The selected
canonical solver must be enabled and runnable. This first export supports
DC reaction-coefficient points only; RF and other quantities fail explicitly.
See [the coefficient output workflow](swarm_inputs.md#結果を他のモデルへ渡す).
Binary64 values are always written losslessly; lossy formatting is not a product option.
Product YAML does not accept external reference inputs.
BOLSIG+ / MCIG files are benchmark-tool inputs, not product schema fields.

In 0.10, figure requests are validated against runnable solvers before execution.
`energy_angle` requires Propagator data; `rf_phase` requires time-dependent
two-term data. Reaction-coefficient figures currently require DC. PNG/SVG,
curve data and native grids share the canonical manifest and atomic publication.
Saved-data plotting preserves numerical source identities and does not grant
new execution or qualification evidence. See [figure usage](swarm_inputs.md#出力と図).

Executed results carry an immutable execution record. A changed case, array,
physical context, observable definition, or run metadata is published with
`evidence_scope: derived_results`, the frozen `source_execution`, and derived
numerical status. Old comparison rows are discarded. Caller-created results
use `provided_results` and cannot acquire execution evidence merely by supplying
a configuration. Their outcome statuses are `derived` and `provided`; the
original execution coverage stays inside the immutable parent record. Unchanged results use `solver_execution` and its content ID.
The writer does not mutate caller data. Publication location is separate.

Each manifest case records typed `transport_observables` with quantity, method,
unit, and approximation flag. Two-term/PN `D_L` is the isotropic F0-gradient
approximation: its equality to `D_T` does not establish a solved longitudinal
response. MC weighted density-packet coefficients have distinct quantities.
Comparisons suppress different observable meanings (`definition_mismatch`)
and different physical conditions (`physics_mismatch`, with differing fields).
Unknown contexts cannot satisfy a required comparison. Matching approximate
observables remain labeled `available_model_approximation`. Workflow databases
and exported table/bundle source policies preserve the same context and meanings.

Observable definitions use actual coefficients. Every `reduced_*` transport
column equals the actual coefficient multiplied by the reported gas number
density N in m^-3: mobility times N has units m^-1 V^-1 s^-1, and diffusion
times N has units m^-1 s^-1. Energy transport is normalized by mean electron
energy, so its mobility/diffusion units match their particle counterparts.
Unavailable coefficients stay missing rather than becoming zero.

Version 0.6 rejects type-coercing or ragged product tables. The CSV schema is
ordered by the actual header; numeric, integer, boolean, and text cells follow
their declared types. Whole-workflow physical context is schema v2 and shares
its normalized physics definitions with direct results. Angular table content
and named e-e parameters are included. Different angular labels or table
contents cannot produce an available comparison; error metrics are suppressed.
`coverage_status: complete` only describes pair presence, while comparison
`status: available` describes compatible input coverage. Numerical convergence,
stale transport, and independent MC qualification remain separate.

Output directories must be nonempty paths. Explicit null or empty strings do
not select defaults or the configuration directory. The Python output API uses
the same validation before creating files. Workflow databases bind a separate
execution identity, including seed/configuration transformations, so version
0.5 databases must be recomputed under new filenames for the 0.6 pipeline.
Numerical model identities are unchanged in 0.6; the existing 23-case core
qualification remains valid without editing or relabeling its source hashes.

Publication snapshots the complete result before checking its execution record.
Mutations of caller-owned arrays, scalars, or metadata during file generation
do not change the published snapshot. This protects both numeric consistency
and the association with immutable execution evidence.

## 0.7 campaign identity and artifacts

Workers receive the prepared calculation snapshot. Their executed input hashes
must match the campaign before cases are saved. `calculation_sha256` identifies
resolved calculation settings and mixtures, excluding publication and campaign
quality controls. The complete YAML is checked for changes during execution.
Old databases cannot be relabeled with this identity. Canonical table and COMSOL
bundle manifests use format version 2; Propagator target evidence uses
`swarm.propagator_target_qualification.v2`. Product manifests remain v3.
