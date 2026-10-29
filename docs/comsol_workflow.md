# COMSOL external Swarm tables

This workflow applies externally calculated electron-swarm data to the COMSOL
6.4 one-dimensional argon positive-column model. It is an auditable
electron-data interface, not a claim that a spatial plasma model is reduced to
a homogeneous swarm calculation.

The reusable execution layer lives in `swarm_workflow/comsol/runtime/` and
owns executable discovery, command construction, process logs, and artifact
provenance. Positive-column mapping, Java generation, verification, workflow,
and comparison code live together in
`swarm_workflow/comsol/models/positive_column/`. A new COMSOL model should add
its own package under `models/` and pass concrete paths to the shared runtime;
the runtime must not import a model package. The cross-model extension
contract is documented in [COMSOL plasma-model adapters](comsol_model_adapters.md).

This adapter accepts a qualified, unprocessed `two_term` bundle with the full
transport closure required by this MPH. `multi_term`, `monte_carlo`, and
`propagator` are explicitly unsupported here. The common bundle exporter and
the other model adapters retain their own solver-source contracts.

`bundle.source: two_term` and `physical_context` are required in the schema-v2
map. The latter declares the model's gas temperature and pure-Ar mixture
(`fraction: 1`, `mass_amu: 39.948`). Preparation requires matching bundle
physical-context and mixture evidence. A `--gas-temperature-K` override is
recorded in the effective mapping and must match the selected bundle; it is
then applied and read back in COMSOL.

Each reaction lookup identifies `species`, `process`, `process_type`, and
`threshold_eV`. The canonical bindings are `e+Ar=>e+Ars` at 11.5 eV for `eir2`
and `e+Ar=>2e+Ar+` at 15.8 eV for `eir4`. A process type alone cannot select a
reaction. Apply, property readback, and spatial audits use the same validated
lookup; duplicate coordinates in the selected channel fail.

COMSOL export validates every manifest-listed source CSV's SHA256 and exact
column inventory before replacing bundle artifacts. Unlisted files are ignored;
unsupported listed tables, unsafe paths, altered coefficients, or changed
quality evidence fail. Copied tables retain the verified source hash, so export
cannot rebind edited tables to the original workflow provenance. Positive-column
preparation additionally requires the current export manifest, source hashes,
physical context, and passing quality CSV covering every table point. Execution
pins that bundle generation across apply, verify, solve, and export.

Since 0.12, `build-tables` writes native format-v3 E/N points. It no longer
writes `transport_vs_mean_energy.csv` or `rates_vs_mean_energy.csv`.
`export-comsol` creates these tables in the format-v2 bundle. Each derived
table records its source filename and hash; the native table is copied too.
Native elastic loss is `elastic_energy_loss.csv`; its COMSOL projection is
`elastic_energy_loss_vs_mean_energy.csv`, retaining the sign and uncertainty.

Every reaction row carries its own observed mean energy. Conversion preserves
the support of each quantity, including MC rates whose E/N range exceeds the
qualified transport range. It never joins or fills values from another grid.
Source units must match the declared canonical units; missing or different
units fail instead of receiving a new SI label. Nonmonotonic native E/N points
remain readable, but a non-single-valued mean-energy function fails COMSOL
conversion. No smoothing, clipping or interpolation repairs are applied.
Rebuild old format-v2 source tables from the recorded database. The Swarm YAML,
map schema v2 and CLI sequence remain the same.

A failed export of a multi-mixture generation clears the old root success
manifest; partially processed children cannot represent a completed root.

The adapter separates typed run contracts, bundle validation, channel lookup
validation, and post-solve audits into `contracts.py`, `bundle.py`, `lookup.py`,
and `audit.py`; `workflow.py` composes the stages. Old maps without the explicit
source, physical context, or reaction identity fail with a migration error.

## Formulation contract

Direct inspection of the archived MPH files found
`MeanElectronEnergyModel = LocalEnergyApproximationE` in the base model, the
external-table model, and the built-in Boltzmann reference. These models use
the local energy approximation (LEA): COMSOL retains and solves the mean
electron energy balance. The archived external path is therefore not a local
field approximation (LFA).

The bundle writes seven arrays:

| Written array | Independent variable | COMSOL property | LEA status |
|---|---|---|---|
| mean electron energy | E/N | `pes1.enrgXdata/enrgYdata` | inactive |
| reduced electron mobility | mean energy | `pes1.muNXdata/muNYdata` | candidate active |
| longitudinal reduced electron diffusion | mean energy | `pes1.deNXdata/deNYdata` | candidate active |
| reduced electron-energy mobility | mean energy | `pes1.mueNXdata/mueNYdata` | candidate active |
| reduced electron-energy diffusion | mean energy | `pes1.denNXdata/denNYdata` | candidate active |
| `eir2` direct-excitation Townsend coefficient | mean energy | `eir2.xtownratedata/ytownratedata` | candidate active |
| `eir4` direct-ionization Townsend coefficient | mean energy | `eir4.xtownratedata/ytownratedata` | candidate active |

Thus, seven written arrays do not mean seven simultaneously active closures.
LEA has six candidate-active arrays; its E/N-to-mean-energy lookup is inactive
because mean energy is a solved field. If the model were deliberately changed
to LFA, the E/N-to-mean-energy lookup would become active, while
electron-energy mobility and diffusion would be inactive because LFA does not
solve a mean-energy transport equation. LFA would therefore have five
candidate-active arrays.

“Candidate active” describes the equation/formulation contract. It is not
proof that an archived executable class actually used the array. Formal
activation requires a clean compile, formulation readback, dependent-variable
inspection, property readback, and a response test in which each candidate
table is perturbed.

COMSOL remains responsible for electron and heavy-species balance equations,
electrostatics, walls and secondary emission, the external circuit, geometry,
mesh, and the nonlinear/time-dependent solver. The mapping replaces selected
electron transport and direct-reaction data only. `eir2` and `eir4` are
individual direct reaction channels, not total excitation and ionization
sources.

The canonical mapping is
`Model/maps/positive_column_external.yaml`. It starts from
`Model/positive_column_1d.mph`; a hand-edited MPH is not a canonical input.
The mapping declares `mean_energy_formulation: local_energy`, and the Java
adapter fails if COMSOL does not read back
`LocalEnergyApproximationE`.

## End-to-end workflow

Use one schema-v2 Swarm YAML and keep solver identity separate from requested physics
features. Add `campaign.database` to that same file for persistence/resume;
there is no external base configuration. See [Swarm inputs](swarm_inputs.md).
The model map uses file-relative paths and its own bundle reference; `--bundle`
is an optional explicit override. Positive-column result output is configured
with `results.output_directory`. A calculation-settings excerpt:

```yaml
schema_version: 2
run:
  solvers:
    - id: two_term
conditions:
  gas_temperature_K: 293.15
  pressure_Pa: 13.3322
feature_policy:
  unsupported: fail
  degraded: record
```

Generate the Swarm result, canonical tables, COMSOL bundle, and spatial result:

```powershell
swarm-workflow run examples\workflow_argon_positive_column_two_term.yaml
swarm-workflow build-tables outputs\positive_column\two_term.sqlite --output outputs\positive_column\two_term_tables --source two_term
swarm-workflow export-comsol outputs\positive_column\two_term_tables\mixture_0000 --output outputs\positive_column\two_term_bundle\mixture_0000
swarm-workflow run-positive-column Model\maps\positive_column_external.yaml --bundle outputs\positive_column\two_term_bundle\mixture_0000
```

`run-positive-column` uses a new class-output directory for every invocation. It treats
compiler error text as failure even when `comsolcompile` returns exit code
zero, and it refuses to execute a stale class. The run records the Java and
class SHA256 values, COMSOL build, source MPH, mapping, bundle, configuration,
and executed voltage sequence.

The run first attempts the configured 200 V target from the source model's
stored state. If needed, it falls back to 20, 50, 100, and 200 V continuation.
`--dry-run` validates the execution plan without requiring a COMSOL license.

### Selecting saved results

The result map explicitly identifies the saved dataset, solution selection,
and field expressions. For example:

```yaml
results:
  output_directory: ../../outputs/positive_column/comsol_external
  dataset: dset1
  solutions: last
  fields:
    - {name: x, expression: x, unit: m}
    - {name: electron_density, expression: plas.ne, unit: 1/m^3}
```

This excerpt illustrates selection; a workflow run still requires all fields
in the maintained map for its numerical acceptance checks. `fields` replaces
the old `probes` expression list. Unknown settings, including point coordinates
and interpolated times, are rejected. `solutions` is required: `all` exports
every saved inner solution, while `last` exports only the final one.
`output_directory` defaults to `outputs/comsol_positive_column` relative to the
map; the maintained map sets it explicitly.

The archived source model has a transient study and dataset `dset1` bound to
`sol1`. `run.voltages_V` is a continuation sequence to its final voltage,
not a parameter sweep producing one result per voltage. The maintained map
therefore explicitly selects `last`. Export loads the saved output MPH and
does not run a study. CSV rows retain `dataset`, `solution_index`, `point_index`,
`time_s`, `parameter_name`, `parameter_unit`, and `parameter_value` before the
requested fields. Stationary results have an empty `time_s`; transient times
are evaluated in seconds. More than one outer solution or inner parameter is
unsupported and fails explicitly. The summary binds the selection, field
units, and CSV hash. Export must create a new file; an earlier CSV or stdout
cannot satisfy it.

The exporter reads back the evaluation's `unit` property after `getData()`
and rejects a change from the requested unit (ignoring whitespace). Complex
fields fail instead of silently exporting their real part. The CSV retains
the requested units; audit and comparison share one reader that converts
known quantities to canonical units. Density and excitation/ionization
sources must be nonnegative; signed currents remain signed.

The portable result entry point is `positive_column_run_summary.json`
(`swarm.positive_column_result.v1`), alongside `positive_column_results.csv`.
It binds the actual CSV hash, gas context, dataset, saved-solution selection
and field units. Move both files together. Python users can call
`swarm_workflow.comsol.models.positive_column.results.read_positive_column_results`
to obtain canonical `.rows`, `.units`, and original `.metadata` without COMSOL.

The API shape and solution metadata follow the official
[Eval reference](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_api_results.52.047.html)
and [SolutionInfo API](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/api/com/comsol/model/SolutionInfo.html).
Software tests execute a two-solution/three-point fixture, and generated Java
compiles against the installed 6.4 API. Licensed model execution remains a
separate acceptance step.

Each new execution replaces any prior completed summary with `running` before
preflight. A preflight rejection writes `blocked`; a later stage failure writes
`failed` with its stage and reason. Both clear acceptance, so an older successful
run cannot be mistaken for the latest request.

### Verification levels

The verifier intentionally separates three questions:

1. **Written-property audit:** read back all seven property arrays and compare
   them with the source bundle.
2. **Formulation audit:** confirm LEA/LFA and the corresponding active and
   inactive sets.
3. **Spatial numerical-consistency audit:** for LEA, compare the stable
   exported spatial variables for reduced mobility, longitudinal reduced
   diffusion, and the `eir2`/`eir4` Townsend coefficients with their table
   evaluations at every spatial point.

Mean energy is not compared with the E/N lookup in the LEA spatial contract;
it is a solved dependent variable. Electron-energy mobility and diffusion are
property-verified but are not claimed to have independent stable spatial
export variables in this model.

Table range is an interpolation-domain check, not a physical-validity
certificate. The production bundle must cover at least 0.05--5000 Td. Values
below the lowest computed point may use only the explicit, manifested
low-field floor policy; negative E/N and upper-range escape fail. Extending a
table to 5000 Td prevents extrapolation but does not, by itself, establish
local-equilibrium or drift-diffusion validity at high field.

## Profile comparison

```powershell
swarm-workflow compare-positive-column external\positive_column_run_summary.json reference\positive_column_run_summary.json --output outputs\comparison_01 --plot
```

The inputs are completed results from the current adapter. Bare CSVs and old
summaries are rejected: do a new current-adapter execution, rather than adding
unit or acceptance labels to historical data. This interface compares, for
example, two meshes or formulations. A standalone export of a built-in
reference MPH is not yet connected to this result interface.

If a result contains multiple saved solutions, pass `--external-solution N`
and `--reference-solution M` using the original COMSOL solution indices.
Selected profiles require the same gas context and temperature, constant
matching voltage and pressure, and matching physical time when present.
Times and operating conditions are never averaged to make profiles match.
Unknown inner parameters, repeated spatial coordinates and disjoint domains
fail. Comparison reconstructs both profiles on the union of native knots and
the endpoints of their common interval, without extrapolation.

Known quantities convert automatically from explicitly supported units:
position m/cm/mm/um; density `1/m^3` or `1/cm^3` (also inverse-power spelling);
mean energy eV/J; potential/voltage V/kV/mV; pressure Pa/kPa/hPa; E/N
`V*m^2`/Td; current density A/m2, A/cm2, mA/m2, mA/cm2 (use `^2` in YAML);
sources `1/(m^3*s)` or `1/(cm^3*s)`; reduced mobility `1/(V*m*s)` or
`1/(V*cm*s)`; reduced diffusion `1/(m*s)` or `1/(cm*s)`; Townsend m2/cm2
(use `^2`). Unsupported units fail during map loading. Custom field names
retain their declared units and are not inferred as additional physical
quantities. Plots label canonical units; reduced field is displayed in Td.

The output directory must be new. Metrics, canonical curves in
`comparison_profiles.csv`, a source-bound `comparison_summary.json`, and
optional PNG/SVG are published together. A failed plot does not leave a
completed comparison. Saved run-stage times are used when available;
`--external-runtime` and `--reference-runtime` can override both solve times.
Partial supplied timing fails; missing recorded times leave speedup missing.

For nonuniform spatial grids, the primary relative L2 metric is based on exact
integration of the piecewise-linear reconstructions:

\[
\epsilon_2 =
\left[
\frac{\int (q_{\rm ext}-q_{\rm ref})^2\,dx}
     {\int q_{\rm ref}^2\,dx}
\right]^{1/2}.
\]

The comparison also reports signed spatial-integral ratio and spatially
weighted correlation. Zero-reference relative errors and zero-variance
correlations remain missing. The unweighted nodal diagnostic is removed.

The CSV column `total_current_density` is specifically
\(J_{\rm e}+J_{\rm Ar^+}\), the electron-plus-argon-ion conductive current.
It is not terminal total current and excludes displacement current and
circuit contributions. `conductive_current_uniformity` reports spatially
integrated mean and standard deviation. Current relative standard deviation (RSD) measures
spatial constancy only; it is not a charge-conservation residual. The
“central 80%” statistic is a geometric window that excludes 10% of the domain
at each end. It must not be relabeled as the plasma bulk, nor may the excluded
regions be identified as sheaths without additional diagnostics.
Current statistics use each result's native full domain; curve errors use
the shared comparison interval stated in the summary.

## Acceptance criteria

- clean compile and class provenance, with no stale-class fallback;
- matching formulation and candidate-active table contract;
- all seven written properties reproduced at readback;
- the four LEA spatial numerical-consistency checks passed;
- identical gas temperature, pressure, voltage, reaction inventory, boundary
  conditions, mesh/order, initial state, and solver tolerances where a
  closure-focused comparison is intended;
- finite outputs, nonnegative densities and physical source terms, and no
  unapproved table-range escape;
- response-based functional activation tests for each candidate-active table;
- raw 200- and 400-element profiles before any grid-sensitivity claim;
- repeated independent processes for timing, with matched timer scopes and
  apply/verify/run/export stages retained separately.

Historical run products are not version-controlled. A benchmark is accepted
only from a fresh run whose generated provenance satisfies the criteria above.
