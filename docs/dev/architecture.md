# Product Architecture

The product accepts schema v2 only. Canonical solver ids are `two_term`,
`multi_term`, `monte_carlo`, and `propagator`, selected with `run.solvers`.
Physics requests live under `physics.*` and remain separate from solver-mode
selection.

## Product execution path

Every product run follows one direction:

1. [`core/config_parser.py`](../../electron_swarm/core/config_parser.py)
   loads schema-v2 YAML and composes the typed section parsers under
   [`core/config_sections`](../../electron_swarm/core/config_sections).
   [`core/config_validation.py`](../../electron_swarm/core/config_validation.py)
   is the single owner of semantic, range, and cross-field invariants.
2. [`orchestration/plan.py`](../../electron_swarm/orchestration/plan.py)
   resolves requested physics against solver capabilities and feature policy.
3. [`core/solver_registry.py`](../../electron_swarm/core/solver_registry.py)
   constructs the canonical solver selected by each runnable plan item.
4. [`orchestration/solver_product_contracts.py`](../../electron_swarm/orchestration/solver_product_contracts.py)
   validates each result against its typed solver product contract. Metadata
   records what the solver actually executed; the contract does not invent
   missing evidence.
5. [`io/writers.py`](../../electron_swarm/io/writers.py) writes canonical
   per-solver and comparison outputs. Legacy aliases are not product outputs.

[`runner.py`](../../electron_swarm/runner.py) and the orchestration executor are
composition layers. Capability support answers whether a model can run;
fidelity metadata separately records the physical closure and assumptions.

## Solver and Monte Carlo evidence ownership

`electron_swarm/core` owns product-domain types, configuration, capability and
feature resolution, transport, and results. `electron_swarm/solvers` owns
numerical implementations and returns core result types. Canonical solver ids
live in the neutral `core/solver_ids.py` owner, so configuration does not
depend back on the executable solver registry.

Each canonical numerical model owns one package under `solvers`: `two_term`,
`multi_term`, `monte_carlo`, and `propagator`. A package `__init__.py` exposes
only its solver class; numerical helpers remain model-local. Shared gas-density
and rate-convolution primitives live in `physics/kinetics.py`. Energy-grid
geometry, neutral-gas collision assembly, finite-volume operators, and EEDF
observables shared by the Boltzmann solvers live under
`solvers/boltzmann_common/`. Two-term growth closure and transport remain under
`solvers/two_term/`; propagator and Monte Carlo do not import the Boltzmann
common package.

The common solver contract is batch execution through `solve_all()`, matching
the orchestration boundary. `two_term`, `multi_term`, and `propagator` use the
small `IndependentCaseSolver` adapter for per-E/N evaluation. Monte Carlo owns
its batch directly because its run setup and random stream span that boundary.

Within `solvers/two_term`, `solver.py` is the public adapter, `models.py` owns
typed numerical results, `grid.py` owns grid construction and initialization,
`steady.py` owns the stationary solve and adaptive refinement,
`time_periodic.py` owns RF-cycle propagation, and `observables.py` owns
transport, rate, and canonical case assembly. The final collision/operator
block is assembled once and reused by result and transport evaluation.
Internal modules never import the public adapter.

Within `solvers/multi_term`, `case.py` and `grid.py` construct the numerical
case, `operator.py` owns the staggered full-PN collision/field matrix,
`steady.py` owns the fail-closed stationary eigenmode, `observables.py` owns
F1 drift, F0 diffusion, and rates, and `result.py` owns product projection.
`solver.py` only composes those responsibilities.

Within `solvers/propagator`, `solver.py` owns operator caching, warm starts,
and steady-case composition; `memory.py` owns bounded-memory estimates and
`case_result.py` owns observable evaluation, diagnostics, and canonical result
assembly. Grid, collision, shell/origin response, angular kernel, steady solve,
and P1 input qualification remain separate numerical owners. Result assembly
does not import the solver adapter.

[`solvers/monte_carlo/evidence.py`](../../electron_swarm/solvers/monte_carlo/evidence.py)
is the public, stable evidence-schema surface consumed by workflow code. The
rest of `solvers/monte_carlo` is private implementation. `solver.py` is the
public solver adapter; `batch.py` derives an independent deterministic stream
for each E/N case from the explicit base seed; `setup.py` validates and prepares shared,
random-free run state; and `case.py` coordinates one case. Mutable case state
and initialization live in `case_state.py`, phase transitions in
`case_phases.py`, and observable/evidence/result assembly in `case_result.py`.
The dependency direction is `solver -> batch -> case -> state/phases/result`;
none of those private owners imports the adapter or re-exports a moved helper.
`audits.py` owns
energy/run audits, `population.py` owns particle-population and resampling
state, and `result_evidence.py` assembles reaction-rate and population
evidence. `weighted_runtime.py` owns the weighted-branching phase and compiled
barrier dispatch. `direct_transport.py` owns fixed-population flux/Helfand
observation, while `weighted_transport.py` owns synchronized weighted-growth
moments and lag evidence. Their shared value contract and weighted-mean formula
live in the small `transport_common.py` owner. Collision, orbit, histogram, and
compiled kernels remain internal numerical details; `kinematics.py` owns the
pure particle velocity, energy, direction, and Boris operations shared by
those details.

## Workflow ownership

- [`selection`](../../swarm_workflow/selection) owns solver-independent physical
  context, immutable whole-closure selection evidence, and bounded low-E/N
  anchor-fallback planning. Campaign, quality, and table layers depend on this
  neutral policy boundary; COMSOL bundle export only consumes its decisions.
- [`commands`](../../swarm_workflow/commands) is the CLI adapter. `parser.py`
  defines the command surface, `handlers.py` translates parsed arguments into
  domain calls, and [`cli.py`](../../swarm_workflow/cli.py) is only the stable
  console entry point.
- [`campaign`](../../swarm_workflow/campaign) owns campaign configuration,
  persistence, aggregation, and execution. `config.py` owns typed workflow
  configuration; `provenance.py` owns the canonical workflow identity and the
  only quality-only resume decision; `repository.py` and `store.py` own SQLite
  reads and writes; `quality.py` binds the pure quality-policy codec to
  immutable SQLite provenance; `statistics.py` and `aggregate.py` own
  replicate statistics and materialization; `sweep.py` owns execution
  planning, scheduling, and sweep summaries. Workflow hashes include
  solver-specific controls only when that solver is enabled, so changes to an
  inactive MC sampling schema cannot invalidate two-term or propagator data.
  Obsolete workflow hashes are rejected; stored metadata is never rewritten to
  manufacture a match.
- [`quality/policy.py`](../../swarm_workflow/quality/policy.py) owns thresholds,
  canonical policy encoding, and pure provenance payloads. It has no campaign
  or SQLite dependency; stored-policy resolution belongs to `campaign/quality.py`.
  [`quality/table.py`](../../swarm_workflow/quality/table.py) owns the exported
  quality-table schema.
- [`quality/monte_carlo/contracts.py`](../../swarm_workflow/quality/monte_carlo/contracts.py)
  owns common pure MC qualification rules.
  [`quality/monte_carlo/policy.py`](../../swarm_workflow/quality/monte_carlo/policy.py)
  owns bounded campaign decisions and sampling-budget policy.
  [`quality/monte_carlo/direct_transport.py`](../../swarm_workflow/quality/monte_carlo/direct_transport.py)
  and
  [`quality/monte_carlo/weighted_transport.py`](../../swarm_workflow/quality/monte_carlo/weighted_transport.py)
  own direct fixed-population and weighted-growth transport evaluation,
  respectively. Each normalizes raw replica input into typed ensemble evidence,
  evaluates stationarity and lag convergence, evaluates coefficient precision,
  and then assembles the canonical quality summary. Their public summarize
  functions only orchestrate those stages. These modules do not own SQLite or
  output writing.
- [`tables/contracts.py`](../../swarm_workflow/tables/contracts.py) owns table
  schemas, constants, errors, and typed datasets.
  [`tables/repository.py`](../../swarm_workflow/tables/repository.py) owns SQLite
  reads, aggregate freshness, and sampling-provenance validation.
  [`tables/monte_carlo.py`](../../swarm_workflow/tables/monte_carlo.py) owns MC
  table qualification and statistical evidence materialization. Its top-level
  builder composes explicit exclusion, support, data-loading, transport-policy,
  and metadata stages without changing the underlying MC estimates, while
  [`tables/energy_loss.py`](../../swarm_workflow/tables/energy_loss.py) owns
  solver-specific elastic energy-loss materialization.
  [`tables/builder.py`](../../swarm_workflow/tables/builder.py) only composes qualification,
  dataset assembly, and canonical table writing.
- [`comsol/input`](../../swarm_workflow/comsol/input) owns provenance-bound
  bundle export, portable selection binding validation, Function-EEDF
  construction and audit, and the mean-energy interpolation coordinate. Its package API
  exposes only `export_comsol_bundle` and its public summary/error contracts.
  Input construction never imports model adapters or the COMSOL execution
  runtime; model adapters may depend on input contracts in that direction.
  Within `function_eedf`, `contracts.py` owns grids and import contracts,
  `io.py` reads and validates artifacts, `moments.py` owns quadrature and tilt,
  `source.py` normalizes source cells, `kernels.py` owns provenance-bound
  collision-rate kernels, `c1.py` builds the continuous family and coordinates
  projection, and `adaptive_grid.py` owns the independently refined physical
  energy and mean-energy axes. Within `export`, `bundle.py`
  orchestrates export, `manifest.py` validates provenance and builds manifests,
  `writers.py` owns common serialization, and the Function-EEDF source, C1,
  and physical-table materializers have separate owners. Within `eedf_audit`,
  `contracts.py` owns the immutable plan and artifact columns, `planning.py`
  builds provenance-bound probes, `sampling.py` selects deterministic native-grid
  samples, `java.py` renders and extracts the saved-model audit, `analysis.py`
  owns acceptance calculations, and `io.py` owns validated artifact reads.
  No legacy monolithic module is retained beside these packages.

## COMSOL model-adapter boundary

[`comsol/runtime`](../../swarm_workflow/comsol/runtime) is the other shared
COMSOL layer. It owns executable and license discovery, Java compilation and
batch execution, process logs, freshness checks, and artifact provenance. It
does not know any plasma model. Each package under `comsol/models/` composes
the public input/runtime services with one MPH topology, mapping schema, solve
sequence, and acceptance contract. CLI registration remains a thin explicit
parser/handler entry for each implemented adapter; there is no universal map
or implicit fallback model. See
[`docs/comsol_model_adapters.md`](../comsol_model_adapters.md).

## GEC-CCP ownership

- [`comsol/models/gec_ccp`](../../swarm_workflow/comsol/models/gec_ccp) is the
  complete model adapter. Its package root exposes only the public run
  contracts, prepare/execute entry points, and user summaries.
- `contracts.py` and `data.py` own stable types and tabular parsing;
  `mapping.py` composes YAML loading, while `config_values.py` and
  `mapping_sections.py` validate primitive values and independent mapping
  sections. `closure.py` owns physical closure decisions and `mph.py` inspects
  the saved-model contract.
- [`validation`](../../swarm_workflow/comsol/models/gec_ccp/validation)
  composes typed manifest/table evidence through provenance, active-table,
  elastic-loss, EEDF, Monte Carlo quality, and final quality stages. Guard,
  joint-consistency, and source-specific checks have concrete owners in that
  package. Validation depends on model contracts and never on execution or
  post-solve audits.
- `planning/evidence.py` assembles qualified closure evidence. `prepare.py`
  only orders validation, layout, Java materialization, and manifest writing;
  `prepare_stages.py` owns those materialization stages. `java.py` composes the
  public GEC COMSOL Java generators, while `java_apply_sections.py` and
  `java_export_sections.py` own application/export sections. Shared Java
  primitives live in `java_support.py`; the audited mean-energy argument
  contract lives independently in `closure_arguments.py`. The
  [`execution`](../../swarm_workflow/comsol/models/gec_ccp/execution) runs the
  plan. Its public package boundary exposes only `execute_gec_ccp_run`;
  `context.py` owns stage results, `inputs.py` freezes execution provenance,
  and `preflight.py`, `solver.py`, `runtime.py`, `postsolve.py`, and
  `status.py` own the typed execution stages. `pipeline.py` only orders them.
- [`audits`](../../swarm_workflow/comsol/models/gec_ccp/audits) owns post-solve
  conservation, Function-EEDF, saved-model, and transport acceptance checks.
- [`plots`](../../swarm_workflow/comsol/models/gec_ccp/plots) is a downstream
  consumer of accepted artifacts. `workflow.py` is the plotting entry point;
  it only coordinates provenance, manifests, and renderer results. Contracts,
  readers/metrics, spatial rendering, and closure/operating-EEDF rendering
  remain independent modules inside that package.

## Positive-column ownership

- [`comsol/models/positive_column`](../../swarm_workflow/comsol/models/positive_column)
  is the independent 1D positive-column adapter. `config.py` owns mapping and
  model-input validation, `java.py` owns generated COMSOL source, `verify.py`
  owns saved-model/function verification, and `comparison.py` owns downstream
  profile comparison and plotting.
- `workflow.py` is the model-specific composition boundary for prepare,
  execute, result validation, and canonical summaries. It depends on the
  shared COMSOL runtime and input contracts; neither shared package imports
  this model adapter.

A new plasma model adds an independent package and thin CLI entry point; it
does not extend either existing mapping with unrelated optional fields or add
model dispatch to the shared input/runtime layers.

## Dependency and qualification rules

- Dependencies point from workflow composition toward public core/solver
  contracts. Solver and core modules must never import `swarm_workflow`.
- Workflow code may consume MC schema constants from public
  `electron_swarm.solvers.monte_carlo.evidence`; it must not import other
  `electron_swarm.solvers.monte_carlo` implementation modules.
- An extracted owner must not import its façade. Façades may compose owners and
  expose documented public entry points, but moved private helpers must not be
  re-exported or restored through compatibility wrappers.
- Propagator qualification artifacts bind the listed source files and inputs
  through the canonical implementation fingerprint in
  [`quality/propagator_source.py`](../../swarm_workflow/quality/propagator_source.py),
  validated by
  [`quality/solver.py`](../../swarm_workflow/quality/solver.py).
  The core runner calls `PropagatorSolver` directly, so the evidence qualifies
  the numerical solver boundary rather than unrelated product orchestration.
  Full qualification and target refinement both compare source/input identity
  before and after their long calculations; a changed source, config,
  cross-section file, workflow evidence row, or qualification algorithm aborts
  publication. Quick smoke output is kept outside the formal artifact path.
- Propagator workflow databases record `propagator_solver_source_sha256`.
  Target refinement requires it to match the qualified core and carries the key
  through table and COMSOL-bundle provenance. Databases created before this
  identity existed must be recomputed rather than relabelled.
  A solver-specific table exports only the selected solver's implementation
  metadata even when its workflow database contains several solvers; aggregate
  database evidence may retain every solver key. Generic COMSOL export
  revalidates the copied core evidence and its source hash before publishing a
  Propagator bundle.
  Editing any fingerprinted numerical-core file invalidates that qualification;
  never bypass the check or hand-edit evidence. Re-run the qualification
  workflow and issue new evidence after an intentional core change.
- The Propagator fingerprint is a source-and-input identity, not a promise of
  bitwise reproducibility across Python, NumPy, SciPy, operating systems, or
  hardware. The qualification artifact records that runtime context for review;
  a stricter environment lock is intentionally not imposed by the product
  contract.
