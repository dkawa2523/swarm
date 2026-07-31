# COMSOL external Swarm tables

This workflow applies externally calculated electron-swarm data to the COMSOL
6.4 one-dimensional argon positive-column model. It is an auditable
electron-data interface, not a claim that a spatial plasma model is reduced to
a homogeneous swarm calculation.

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

Use schema version 2 and keep solver identity separate from requested physics
features:

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
swarm-workflow sweep examples\workflow_argon_comsol.yaml
swarm-workflow build-tables outputs\argon_comsol\swarm.sqlite --output outputs\argon_comsol\tables --source two_term
swarm-workflow export-comsol outputs\argon_comsol\tables\mixture_0000 --output outputs\comsol_bundle\mixture_0000
swarm-workflow run-comsol Model\maps\positive_column_external.yaml --bundle outputs\comsol_bundle\mixture_0000
```

`run-comsol` uses a new class-output directory for every invocation. It treats
compiler error text as failure even when `comsolcompile` returns exit code
zero, and it refuses to execute a stale class. The run records the Java and
class SHA256 values, COMSOL build, source MPH, mapping, bundle, configuration,
and executed voltage sequence.

The run first attempts the configured 200 V target from the source model's
stored state. If needed, it falls back to 20, 50, 100, and 200 V continuation.
`--dry-run` validates the execution plan without requiring a COMSOL license.

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
swarm-workflow compare-comsol external.csv builtin.csv --output outputs\comsol_comparison --external-runtime 10 --reference-runtime 40 --plot
```

Both CSV files must contain constant, matching `applied_voltage` and
`gas_pressure` columns. The built-in result is a reference, not experimental
ground truth.

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
weighted correlation. The former unweighted nodal L2 is retained only as a
legacy diagnostic.

The CSV column `total_current_density` is specifically
\(J_{\rm e}+J_{\rm Ar^+}\), the electron-plus-argon-ion conductive current.
It is not terminal total current and excludes displacement current and
circuit contributions. Current relative standard deviation (RSD) measures
spatial constancy only; it is not a charge-conservation residual. The
“central 80%” statistic is a geometric window that excludes 10% of the domain
at each end. It must not be relabeled as the plasma bulk, nor may the excluded
regions be identified as sheaths without additional diagnostics.

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

The present archive does not satisfy all of these criteria: its historical
class provenance is stale, and the clean rerun is blocked by COMSOL license
error `-10` (`Product has expired`). See `docs/comsol_benchmark.md` for the
provisional numerical results and their interpretation limits.
