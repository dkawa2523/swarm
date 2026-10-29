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
   [`orchestration/preparation.py`](../../electron_swarm/orchestration/preparation.py)
   snapshots configuration, loaded inputs, concrete solver settings, and source
   identities before execution. The runner verifies that snapshot again before
   publishing; callers keep ownership of their unchanged configuration.
3. [`core/solver_registry.py`](../../electron_swarm/core/solver_registry.py)
   constructs the canonical solver selected by each runnable plan item.
4. [`orchestration/solver_product_contracts.py`](../../electron_swarm/orchestration/solver_product_contracts.py)
   validates each result against its typed solver product contract. Metadata
   records what the solver actually executed; the contract does not invent
   missing evidence.
5. [`io/writers.py`](../../electron_swarm/io/writers.py) writes canonical
   per-solver and comparison outputs. Legacy aliases are not product outputs.
   [`io/artifacts.py`](../../electron_swarm/io/artifacts.py) owns the typed,
   lossless scalar-table encoding and the result object that retains manifest
   context alongside DataFrames. Comparisons are computed once and published
   from the result, so output publication cannot silently disable them.

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
only its solver class; numerical helpers remain model-local. Shared gas conditions and density live in `core/gas.py`; rate-convolution
primitives live in `physics/kinetics.py`. Shared
energy-grid construction and authoritative reported-cell geometry live in
`grids/energy.py`. Neutral-gas collision assembly, finite-volume operators, and
EEDF observables shared by the Boltzmann solvers live under
`solvers/boltzmann_common/`. Two-term growth closure and transport remain under
`solvers/two_term/`; propagator and Monte Carlo do not import the Boltzmann
common package.

Canonical EEDF cells cover a contiguous energy domain from zero, with reported
widths defining their support. Solver comparison, workflow aggregation, plots,
and Function-EEDF construction consume that same geometry; they never infer
new widths from the evaluation nodes or assume that nodes are cell midpoints.

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
`steady.py` owns the fail-closed stationary eigenmode, `refinement.py` owns
bounded whole-PN energy-domain re-solves, `observables.py` owns
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

## Single-document execution boundary

`campaign/document.py` owns the schema-v2 calculation document plus optional
campaign section. It resolves aligned composition columns and delegates scalar
calculation validation to the product parser. `campaign/config.py` validates
persistence, replicas, resources, and bounded continuation. It retains the
resolved calculation snapshot rather than a reference to another YAML file.
The common `swarm-workflow run` CLI accepts both direct and campaign documents;
`electron_swarm` stays independent of workflow and SQLite.

Campaign preparation and execution share one plan. Job payloads carry the
prepared calculation and expected input hashes; workers do not reload editable
YAML. The canonical calculation identity excludes publication and campaign
quality controls, while execution also checks the complete source document.
Removing a required field or an unimplemented campaign comparison fails before
creating the database. The old external base-file loader path is removed.

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
  immutable SQLite provenance; `aggregate.py` owns materialization and shares
  independent-replica estimates and conservative probability rebinning from
  [`core/statistics.py`](../../electron_swarm/core/statistics.py);
  `sweep.py` owns execution
  planning, scheduling, and sweep summaries. Workflow hashes include
  solver-specific controls only when that solver is enabled, so changes to an
  inactive MC sampling schema cannot invalidate two-term or propagator data.
  Obsolete workflow hashes are rejected; stored metadata is never rewritten to
  manufacture a match.
- [`reaction_rates.py`](../../swarm_workflow/reaction_rates.py) reads product
  or campaign snapshots into neutral DC coefficient points. It has no COMSOL
  or full-transport qualification dependency. The core
  [`io/reaction_rates.py`](../../electron_swarm/io/reaction_rates.py) owns that
  table schema, fixed physical series, uncertainty, and the read/write API.
  [`io/publication.py`](../../electron_swarm/io/publication.py) owns the shared
  locked generation replacement. Numerical solvers do not import these layers.
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
  dataset assembly, and native format-v3 table writing on E/N. Reaction rows
  retain observed mean energy; elastic-loss points retain their sign and
  statistical evidence. Mean-energy monotonicity is not a native-data gate.
- [`comsol/input`](../../swarm_workflow/comsol/input) owns provenance-bound
  bundle export, portable selection binding validation, Function-EEDF
  construction and audit, and the mean-energy interpolation coordinate. Its package API
  exposes only `export_comsol_bundle` and its public summary/error contracts.
  Input construction never imports model adapters or the COMSOL execution
  runtime; model adapters may depend on input contracts in that direction.
  `validation.py` owns contained artifact paths, immutable table hashes,
  and exact CSV schemas. Export verifies every source table before copying
  and verifies the copied bytes against that same identity. Model adapters
  add source quality and physical-context checks without redefining byte
  integrity. Within `function_eedf`, `contracts.py` owns grids and import contracts,
  `io.py` reads and validates artifacts, `moments.py` owns quadrature and tilt,
  `source.py` normalizes source cells, `kernels.py` owns provenance-bound
  collision-rate kernels, `c1.py` builds the continuous family and coordinates
  projection, and `adaptive_grid.py` owns the independently refined physical
  energy and mean-energy axes. Within `export`, `bundle.py`
  orchestrates export, `manifest.py` validates provenance and builds manifests,
  `writers.py` owns common serialization, and the Function-EEDF source, C1,
  and physical-table materializers have separate owners.
  `coefficient_tables.py` reindexes native points to mean-energy tables for
  COMSOL, validates each quantity's own support, and records source-table
  hashes. It does not join rates to transport nodes, fit, interpolate or repair
  values. Native table and COMSOL bundle versions are independent; bundle v2
  retains its model-facing bindings. Within `eedf_audit`,
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
  is the independent 1D positive-column adapter. `contracts.py` owns public
  run/result types, `config.py` owns mapping validation, and `bundle.py` owns
  qualified source, physical-context, and table validation. `lookup.py` owns
  reaction identity and unique interpolation coordinates; Java generation
  and saved-model verification consume the same validated reaction selection.
  `java.py` owns generated COMSOL source, `verify.py` owns saved-model/function
  verification, `results.py` owns one strict native CSV/summary reader, and
  `result_units.py` owns explicit quantity meanings and supported conversion
  scales. Audit consumes one normalized snapshot; it does not reread raw CSV
  or repeat numeric/selection validation. `comparison.py` selects saved
  solutions and integrates piecewise-linear spatial profiles, publishing
  curve CSV, metrics and optional PNG/SVG in one new directory. Its current
  uniformity report is not a charge-conservation residual.
- `workflow.py` is the model-specific composition boundary for prepare,
  execute, and canonical summaries. It depends on the
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
- All solver workflows bind installed runtime sources and numerical-library
  versions through
  [`core/source_identity.py`](../../electron_swarm/core/source_identity.py).
  The same identity is checked before/after execution, resume, and table export.
  Tests, CLI qualification tools, and local repository directories are outside
  the installed runtime boundary. Propagator qualification uses this identity,
  validated by
  [`quality/solver.py`](../../swarm_workflow/quality/solver.py).
  The core runner calls `PropagatorSolver` directly. Its immutable artifact
  includes input text and hashes so validation also works from an installed
  wheel. Target comparison and validation code are installed workflow modules.
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

## Physical meaning and publication

`core/physical_context.py` describes one resolved physical problem independently
of numerical settings. `core/observables.py` owns typed quantity/method/unit
meanings. The runner attaches both to each case; workflow storage and table
policies preserve them. `core/cross_sections.py` reconciles source mass rounding
once and passes one authoritative configured mass to all numerical engines.

`core/provenance.py` seals an execution's settings and result content. Publication
uses that immutable record to distinguish solver execution, derived data, and
provided data. Modified/provided results do not inherit convergence or comparison
success. The original execution remains inspectable for a derived result.

`core/source_identity.py` separates the numerical model from the execution
pipeline. CLI and output formatting have their own adapter/schema; changing
them does not invalidate numerical model qualification. `cli.py` owns argument
parsing and exit codes; `runner.py` owns execution composition. Run/feature-policy
parsing and output/comparison parsing have separate section owners. Output and
comparison dataclasses live in `core/output.py`; input/numerical default types
remain covered by execution identity in `core/config.py`.

The physical-context owner also defines typed angular/e-e views and the common
workflow context. Angular table bytes are identified independently of paths;
inactive controls are normalized. Selection and model adapters consume that
shared schema rather than reconstructing physics locally.

`io/artifacts.py` uses a strict record reader and schema-driven scalar conversion.
All output producers preserve Python/NumPy scalar types before encoding; pandas
inference cannot erase large integer counts. Output config validation is owned
by its dataclass and rechecked on the publication snapshot.

`campaign/identity.py` owns a portable workflow execution identity distinct from
solver models and core execution. It binds seed/config projection, persistence,
and aggregation. The same guard protects execution, resume, and table reuse.
Numerical source qualification is not invalidated by these I/O-only changes.

Publication deep-copies results before checking the execution record and before
building any table. This closes the check/use gap between result identity and
serialization; the caller may retain mutable result objects without mutating
the generation being published.

## Consumer configuration in 0.7

COMSOL model settings remain separate from Swarm calculation settings. All
hand-written mapping paths resolve from the mapping directory, while bundle
references resolve from the bundle. Effective positive-column maps freeze
resolved paths when relocated; `results.output_directory` owns their output
location. A mapping's bundle reference is sufficient without a duplicate CLI
argument. Canonical table/bundle format v2 carries `calculation_sha256`.
The format version has one owner in `tables/contracts.py`.

Simulator-independent DC reaction coefficients are implemented separately from
COMSOL table requirements. A non-COMSOL coefficient consumer example remains
I4 work; no empty Global-model package or chemistry integrator is introduced.

## Standard figures in 0.10

`io/tables.py` owns native result-table projection. `io/publication.py` owns the
single staged generation, typed CSVs, image hashes and rollback for product,
coefficient and figure outputs. The previous duplicated CSV write paths are
removed. `io/figure_data.py` groups fixed physical conditions and independent MC
replicas; `io/figures.py` owns rendering from those numerical tables. Matplotlib
is optional and loaded only for requested figures.

`swarm_workflow/figures.py` reads a verified manifest or a read-only SQLite
snapshot. It does not call a solver or current-source qualification for analysis.
It retains historical source identities, while coefficient export still applies
its existing current-implementation and planned-replica requirements. Campaign
workers do no plotting; their complete saved snapshot is published once with
both figures and any selected coefficient export.

The public replot API takes a saved source path, rather than accepting mutable
artifact objects that could retain an execution claim after numerical edits.
Source numeric CSV hashes remain mandatory when regenerating missing images.
