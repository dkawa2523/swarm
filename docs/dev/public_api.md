# Public Product API

Public entry points (0.10):

- Python package: `electron_swarm`
- CLI module: `python -m electron_swarm <config.yaml>`
- Preflight without calculation: `python -m electron_swarm <config.yaml> --plan`
- Console script: `electron-swarm <config.yaml>`
- Python helpers: `load_config`, `run`, `run_from_config`, `SwarmConfig`

The common calculation/campaign CLI is
`swarm-workflow run <config.yaml> [--plan [--json]]`. A single Swarm YAML adds
`campaign` only when persistence, independent replicas, or multiple mixtures
are needed. `swarm_workflow.campaign.document.read_calculation` returns a
resolved document; its `template` is explicitly the first case template, while
its `mixtures` describe all cases. It is not an implicit all-cases runner.
`load_workflow` validates campaign controls and retains their calculation
snapshot. `plan_sweep` and `run_sweep` share the same plan. Workers consume
snapshots and verify executed input hashes before storage.

The lower-level `electron-swarm` CLI and the core Python API continue to serve
direct single-mixture calculations. They do not silently ignore a campaign.
See [Swarm inputs](../swarm_inputs.md) for the maintained input surface.

Stable schema surface:

- required `schema_version: 2`
- solver ids: `two_term`, `multi_term`, `monte_carlo`, `propagator`
- `multi_term` methods: `pn_closure_direct`, `pn_dcs`
- solver selection under `run.solvers`
- physics requests under `physics.*`
- policy under `feature_policy`

Stable output files:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_energy_angle_distribution.csv` when returned by a solver
- `<base>_rf_phase.csv` and `<base>_rf_phase_eedf.csv` for phase-resolved RF
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled
- `<base>_reaction_rates.csv` for a selected neutral export
- `<base>_figure_data.csv` and `<base>_figure_<kind>_<NNNN>.png/svg` for requested figures
- `<base>_manifest.json`

`electron_swarm.io.read_product_run(manifest_path)` returns a
`ProductRunArtifact` after verifying the generation, content hashes, and table
schemas. Its `.metadata` includes the outcome, comparison coverage, resolved
settings, input attribution, configuration, input hashes, and implementation
identities. Its `.tables` maps names such as `summary_csv` to DataFrames.
CSV floats round-trip binary64 values; strings such as `NA` and `001` and large
integer seeds keep their original meanings. Missing numeric values remain
missing. Formatting for presentation belongs in a consumer of these tables.

`write_outputs(result, output)` publishes the already computed comparison in
the result; it has no separate comparison configuration. The publication
location is recorded independently of the original execution configuration.
With no executed solver, only the plan and manifest are emitted. For an executed result, manifest
`outcome.status` is `complete`, `partial`, or `skipped`; CLI `skipped` returns 3 and `partial` returns 4. Disabled solvers are not requested.
Completion describes execution coverage, not numerical or physical validation.

Preflight `--plan` returns an object with status, `solver_plan`,
`resolved_solvers`, `physical_cases`, and `expected_case_count`. Resolved settings include MC
defaults, trial-period units, transport-leg count, and Propagator grid semantics.

Comparison summary rows include `angular_model_status`, `comparison_status`,
numerical status and transport definition for both solvers, per-observable
availability, scalar relative differences, and optional `eedf_l1_error`.
Optional missing results produce explicit rows and comparison coverage rather
than disappearing. A required comparison rejects incomplete coverage and incompatible or unknown physical/angular inputs. The
default reference is the requested MC solver when present, otherwise the first
requested canonical solver. A nonzero value against a zero reference has an
undefined relative difference, recorded as missing with `reference_zero`.

Version 0.4 removes the context-dropping `read_product_outputs` API,
`output.float_format`, unused `conditions.length_scale_m`, and the public
`two_term.nonconservative_model` switch. The physical ionization choice is
`physics.ionization.energy_sharing`; old MC collision-count names must migrate
to `production_trial_periods`, `warmup_trial_periods`, and `tail_trial_periods`.
Duplicate YAML keys and unused workflow MC/execution blocks fail at loading.

Stable summary metadata:

- solver method and angular model/source
- ordinary-XS closure / exact-DCS flags
- `multi_term` `lmax`
- ionization source treatment
- e-e treatment and transport-stale flag
- transport definition
- magnetic treatment
- tail refinement treatment
- active collision-process treatment and numerical status
- propagator velocity-space, scattering-role, recoil, and transport treatment

Direct-PN residuals, negative-mass diagnostics, MC audit counters, benchmark
failure categories, and solver-plan capability details are internal or
benchmark outputs unless documented in the product schema.

Version 0.5 requires explicit gas masses and moves population selection to
`physics.ionization.population_model` (default `weighted_branching`). The CLI
adapter is `electron_swarm.cli`; orchestration stays in `runner`.

Every executed `SwarmCaseResult` carries a frozen `physical_context` and a
snapshot of typed `transport_observables`. Comparison requires a common physical
problem and selects only observables with matching quantities/units. F0-gradient
longitudinal diffusion remains explicitly approximate.

`write_outputs` checks the immutable execution record before publication.
Unchanged results have `evidence_scope=solver_execution`; changed results have
`derived_results`, their original `source_execution`, and no old comparison;
caller-created results have `provided_results`. Their `outcome.status` is
`derived` or `provided`, respectively. Derived/provided cases cannot
retain a numerically-converged execution claim. The manifest output schema is
`swarm.product_run.v3`. Numerical model identity, execution identity, and output
schema are versioned separately so CLI/output edits do not invalidate numerical
qualification. Workflow storage and table policies preserve case context.

## Typed boundary in 0.6

The configuration remains schema v2. Product manifests use
`swarm.product_run.v3`; the reader rejects older product manifests rather than
silently assigning the new meanings. Whole-workflow physical context uses
`swarm_physical_context.v2`. Its normalized physics view is owned by
`core/physical_context.py` and exposed as `workflow_physical_context` there;
the duplicate `selection.physical_context` helper is removed.

Each case includes named angular and e-e physics. The angular view binds the
moment-table content hash, provenance, and extrapolation policy. E-e targets
use named parameters; an unused target temperature does not create a different
physical problem when mean energy is conserved.

Comparison `status` is `available`, `partial`, `unavailable`, `incompatible`,
`unverified`, or `disabled`. `coverage_status` separately says whether all
requested pairs exist. `qualification_status` records convergence, stale
transport, or the need for an independent MC ensemble. `available` is not a
claim of independent physical validation. Required and optional comparisons
share input eligibility; incompatible or unknown inputs have no relative error
or EEDF norm.

Typed CSV reading validates every record width before conversion, handles
quoted multiline text, and never infers a dataframe index. Writing follows
actual column order and rejects type/precision loss before replacing a file.
Boolean/text values cannot masquerade as numbers. Large integers retain integer
storage through all producer dataframes and the CSV codec. Exact integral
floating counts below 2**53 are accepted; larger counts require integer storage.

`OutputConfig` validates paths and filename components for both YAML and Python
calls. Explicit null/empty output directories are invalid. `write_outputs`
validates a snapshot even if a caller mutated the settings after construction.

Workflow execution identity (`swarm.workflow_execution_sha256.v2`) covers
configuration projection, seed derivation, storage, and aggregation plus the
pandas version and the shared scalar/probability statistics. Resume and table building reject missing or changed identity.
CLI/publication changes remain outside numerical model and workflow identities.

Publication owns an independent snapshot of the complete result (including
arrays and metadata) and output settings. Changes to caller-owned objects
during serialization cannot alter that generation or borrow its execution
identity. The comparison count `available_pairs` counts pairs with both input
cases present; input eligibility is reported separately by `status` and the
physical/angular statuses. A missing text cell is canonicalized to empty text;
missing numeric/boolean cells remain nullable.

## Campaign input in 0.7

The external base-file reference and the CLI `sweep` alias are removed. Use one
schema-v2 YAML and `swarm-workflow run`. E/N is under `run.e_over_n_Td`, gas
fractions under `conditions.gas_mixture`, and the master seed under
`solvers.monte_carlo.seed`. Equal-length fraction arrays represent aligned
mixtures. The campaign owns persistence and sampling controls, not a second
copy of the physical inputs. Comparison output for campaigns is explicitly
unsupported and rejected while planning.

`calculation_sha256` replaces the old file-based calculation identity in
workflow/table/bundle provenance. Table and bundle manifests use format v2;
Propagator target evidence uses v2. Old identities and formats are rejected,
without aliases or edited historical hashes. Product scalar output remains v3.
Positive-column mapping paths are file-relative; effective maps retain resolved
paths. Its `results.output_directory` selects the run artifacts. All three
COMSOL run CLIs use the mapping's bundle unless an explicit override is given.

## Result persistence and selection in 0.8

`WorkflowStore.read_case(mixture_id=..., solver=..., e_over_n_Td=...,
replicate=...)` restores one typed `SwarmCaseResult`, including its physical
context, observable definitions, rates, EEDF/counts, angular distribution, and
RF phase data when present. This restores data, not an execution record;
publishing it as a caller-created run has `provided_results` evidence scope.
RF phases and their EEDFs are stored atomically with the case. Resume excludes
an RF case with missing phase rows or cells. Old stores without the required
RF schema are rejected and must be recalculated into a new database.
The unused per-bin `1/sqrt(ESS)` diagnostic is removed; independent-replica
uncertainty remains in aggregate outputs.

Positive-column results use typed `ResultOutput` and `ResultField` settings:
`results.output_directory`, `dataset`, required `solutions: all|last`, and
`fields: [{name, expression, unit}]`. The old `probes` key and unknown
coordinate/time settings fail with explicit errors. Export reads the saved
output MPH without running a study. CSV identity columns are `dataset`,
`solution_index`, `point_index`, `time_s`, `parameter_name`, `parameter_unit`,
and `parameter_value`, followed by the requested fields. Stationary solutions
have an empty `time_s`. Only one outer solution and at most one inner solution
parameter are supported. Spatial profile comparison requires one selected
solution per result and matching operating conditions, including time.

## Neutral coefficient points in 0.9

`OutputConfig.export` is an optional frozen `RateExportConfig(solver=...,
quantities=("reaction_rates",))`. YAML uses the matching `output.export`
mapping with a quantities list. The selected solver must be enabled and
runnable for every requested mixture. DC E/N is the only supported axis.
Plans report the selection and canonical CSV path; they create no files.

Direct publication adds `<base>_reaction_rates.csv` to the same product v3
manifest. Campaigns without figures and the saved-data export command use
`<base>_manifest.json` with `swarm.reaction_rates.v1`.
`electron_swarm.io.read_reaction_rates(path)` reads either manifest and returns
`ReactionRatesArtifact.metadata` and `.table`. It preserves the typed CSV,
hashes, fixed-context series, native points/range, and source identities.
`write_reaction_rates(artifact, output)` republishes a neutral artifact;
edited data lose their numerical claim and retain source lineage.
All publishers atomically replace one generation per directory/base name,
including removal of obsolete canonical outputs when switching to campaign
coefficient output. Saved-data export refuses to overwrite its input generation.

`swarm_workflow.reaction_rates.export_reaction_rates(source_path, output,
solver=...)` accepts a product manifest or a current workflow SQLite database.
The CLI is `swarm-workflow export-reaction-rates <source> --solver <id>
--output <directory> [--base-name swarm]`. Source databases are read only.
Rate-only exports do not invoke full-table aggregation, inspect EEDF grids,
require transport, convert Townsend coefficients, or reorder mean energy.
They never substitute another solver or fill missing channels with zero.

Rows identify `series_id`, solver, E/N, species/process/type/threshold,
target fraction, per-target `k`, `x k`, independent replica count, SE/CI95,
estimate status, and quantity-relevant numerical status. Fixed physical
context and original per-case solver/feature metadata live in each series.
`postprocessed_rates` preserves updated e-e rates while retaining the
source's stale-transport flag. MC precision remains explicitly unassessed;
independent-replica intervals are not a plasma-model input qualification.
One-replica and all-zero MC intervals remain missing; all zero is censored.
The format contains coefficient points, with no interpolation or extrapolation.

## Standard figures in 0.10

`OutputConfig.figures` is an optional tuple of `eedf`, `transport`,
`reaction_rates`, `energy_angle`, and `rf_phase`; YAML uses a unique list.
It defaults to empty. Plans check source availability, DC reaction-rate axes
and RF field type before creating files. Install the optional `[plot]` extra.
Neither `run` nor the plot command adds a chemical-species integrator.

`swarm_workflow.figures.plot_results(source_path, output)` and
`swarm-workflow plot <manifest-or-db> --figures eedf transport
--output <directory> [--base-name swarm]` plot saved numerical data. An existing
manifest can be regenerated in place, including missing images, after all CSV
hashes and schemas have been verified. Numerical execution identities remain
unchanged. Compatible historical DBs can be analyzed without relabeling their
source implementation or granting a new qualification.

`electron_swarm.io.read_figures(path)` verifies PNG/SVG and typed tables,
returning `ProductRunArtifact`. Direct plots share product v3; DB snapshots
use `swarm.figures.v1`, role `swarm_figures`, evidence `saved_results_analysis`.
Neutral-coefficient plots retain `swarm.reaction_rates.v1`. Every role uses
`<base>_manifest.json`. Selected campaign coefficient exports and figures are
published from one DB snapshot in the same generation; `read_reaction_rates`
also accepts that figure manifest. Analysis does not relax coefficient export
requirements. A campaign requesting neither export nor figures writes its DB
only, with no new analysis publication.

`<base>_figure_<kind>_<NNNN>.png/svg` and `<base>_figure_data.csv` are owned
outputs, atomically replaced or removed with the generation. Manifest
`figures` records panels, units, curve source cases, numerical status,
display limits, generator version and image hashes. Curve points retain
cell bounds, replica counts and pointwise t95 intervals; 2D panels refer to
the unchanged native angular/RF tables. No duplicate grid is introduced.

EEDF cell densities are drawn as steps without renormalization. Signed
native tails remain signed; log omissions and cropped absolute mass are
declared. Independent MC replicas are averaged equally, conservatively
rebinned to common energy support when necessary. No missing interval becomes
zero. Different transport quantities have separate panels; stale or approximate
transport stays labeled. Reaction curves use per-target k with native E/N
markers; connecting lines are visual guides. RF heatmaps use periodic nearest
saved phase samples, preserving the source phase coordinates.

Readers verify image hashes by default. Numeric consumers may call
`read_product_run(path, verify_figures=False)`; CSV verification is mandatory.
Saved reaction-rate export also tolerates missing images without weakening
its numeric checks. Plotting accepts saved paths, rather than mutable artifact
objects that could borrow execution evidence after edits.

## Positive-column results in 0.11

`swarm_workflow.comsol.models.positive_column.results.read_positive_column_results`
accepts `positive_column_run_summary.json` (`swarm.positive_column_result.v1`)
and returns `PositiveColumnResults(rows, units, metadata)`. It verifies the
relative CSV and hash, retains dataset/solution/point/time identities, and
normalizes explicitly supported field units. Custom fields keep declared
units. Nonfinite fields, duplicate identities, missing units and negative
density or excitation/ionization sources fail. Signed currents are allowed.

`compare_comsol_profiles(external_results, reference_results, output_dir=...,
external_solution=None, reference_solution=None, create_plot=False)` now takes
two current completed result summaries. Multiple saved solutions require the
original COMSOL indices. Gas context, voltage, pressure and physical time must
match. The comparison exports exact curve values, spatially integrated
metrics, optional PNG/SVG, and source summary/CSV hashes. Its output directory
must be new; publication is atomic. Bare CSV/old-summary compatibility and
unweighted nodal diagnostics are removed.

The single map still selects fields and export units; no plotting YAML is
added. Readback and comparison require no new solve. A standalone command to
re-export an arbitrary saved MPH, including built-in reference results, is
not implemented. See [the adapter workflow](../comsol_workflow.md#profile-comparison)
for supported units and CLI selection.

## Native tables and COMSOL conversion in 0.12

`build_tables` keeps its existing API and three qualified sources. Its source
manifest is now format v3, independent of the COMSOL bundle's format v2.
Native tables use E/N; `rates_vs_en.csv` also contains `mean_energy_eV`.
`elastic_energy_loss.csv` is solver-native data with the same sign and
independent-replica uncertainty. Native source manifests no longer claim a
secondary interpolation axis or a COMSOL input role.

`export_comsol_bundle` derives `transport_vs_mean_energy.csv`,
`rates_vs_mean_energy.csv`, and the available elastic-loss mean-energy table.
It checks declared units and the single-valued mean-energy axis per quantity,
retains different native support ranges, and copies source tables unchanged.
Each derived table records `source_table.file`, `source_table.sha256`, and the
projection method. Nonmonotonic native points are retained at the source;
conversion fails instead of repairing them. Rebuild old source tables rather
than relabeling their version. The neutral four-solver reaction-rate export,
schema-v2 configuration and model-specific closure requirements are unchanged.
