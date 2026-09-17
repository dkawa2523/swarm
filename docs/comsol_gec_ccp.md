# Argon GEC-CCP external-Swarm contract

This document defines the current two-term, Monte Carlo, and P1 propagator
physical inputs to the time-periodic argon GEC-CCP COMSOL model. This is a
model-specific adapter contract. Shared COMSOL services and the
extension path for other plasma models are defined in
[COMSOL plasma-model adapters](comsol_model_adapters.md).

The public entry points are:

| Purpose | Two term | Monte Carlo | Propagator P1 |
| --- | --- | --- | --- |
| Base configuration | `examples/argon_gec_ccp_base.yaml` | `examples/argon_gec_ccp_monte_carlo_weighted_branching.yaml` | `examples/argon_gec_ccp_propagator.yaml` |
| Sweep workflow | `examples/workflow_argon_gec_ccp_two_term.yaml` | `examples/workflow_argon_gec_ccp_monte_carlo.yaml` | `examples/workflow_argon_gec_ccp_propagator.yaml` |
| COMSOL mapping | `comsol_modes/maps/argon_gec_ccp_two_term_function_eedf.yaml` | `comsol_modes/maps/argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml` | `comsol_modes/maps/argon_gec_ccp_propagator_function_eedf.yaml` |
| Result role | physical target | physical target after MC qualification | physical target after core and target-grid qualification |

These product configurations use `schema_version: 2`, `run.solvers`, and
their canonical ids `two_term`, `monte_carlo`, and `propagator`.

## Physics ownership

The Swarm calculations are steady, spatially homogeneous DC calculations.
COMSOL consumes their results through local mean electron energy during its RF
time-periodic plasma solve. This is a restricted local-mean-energy closure; it
is not an RF kinetic solve.

| Quantity | Two-term target | MC target | Propagator target | COMSOL responsibility |
| --- | --- | --- | --- | --- |
| EEDF | Adaptive moment- and rate-controlled spreadsheet Function EEDF | Adaptive moment- and rate-controlled spreadsheet Function EEDF | Adaptive moment- and rate-controlled spreadsheet Function EEDF | Evaluate the active local-energy closure |
| Reduced electron mobility | External two-term value | External weighted-MC value | External flux-moment value | Apply the imported value |
| Particle diffusion | Not imported | Not imported | Unavailable | Einstein relation from imported mobility |
| Energy mobility and diffusion | Not imported | Not imported | Unavailable | Restricted local-energy transport |
| Elastic electron-energy loss | External solver-native moment | External direct trajectory moment | External same-operator moment | Apply as the electron-energy sink |
| Excitation and ionization | COMSOL integrates cross sections with the Function EEDF | Rates are preintegrated from the same MC EEDF | COMSOL integrates cross sections with the Function EEDF | Apply the source-specific rate closure |
| RF field, Ar/Ar+, walls, secondary emission | Not a Swarm output | Not a Swarm output | Not a Swarm output | Native COMSOL model |

All three solvers use one COMSOL Spreadsheet serialization. Both tensor axes
are adaptive: every actual solver anchor is retained on the mean-energy axis,
and additional C1 rows are inserted only where COMSOL's linear interpolation
needs them. Every materialized row preserves normalization and mean energy.
Refinement continues until the EEDF total-variation error and every available
active cross-section-rate integral meet the declared projection budget. Rate
importance floors belong only to the error norms; they never change an EEDF
or rate value, and exact zeros remain zero.

The qualified MC production route imports inelastic rate tables computed from
the same EEDF and argon cross sections. It does not substitute two-term
coefficients or add an artificial rate floor. The P1 propagator is projected
from its finite-volume energy cells through the same source-independent
adaptive representation. It otherwise uses the same restricted
local-mean-energy closure as the two-term target,
with `swarm_mobility_einstein` and `external_solver_native` elastic loss. It
does not supply particle or energy diffusion. Its deterministic P0/P1 core
qualification passes, including low-field, medium--fine, mixture, energy-
ceiling, and weak-transfer gates. A separate target artifact binds the
300 x 48 production grid to independent 600 x 72 checks at 3000, 3500, and
4000 Td; this evidence is mandatory for the GEC preflight. The mapping
declares a `physical_target`, but a run is accepted only after its
Function-EEDF, operating-support, transport, conservation, and saved-model
audits pass. That target is limited to the explicit restricted
local-mean-energy closure; it does not imply P2, nonlocal-RF, or unavailable
diffusion claims.

Ordinary integral cross sections do not define differential scattering. The
current MC qualification is explicitly limited to DC, zero magnetic field,
isotropic scattering, no electron-electron collisions, equal ionization energy
sharing, and zero high-energy cross-section extrapolation.

## Monte Carlo work budget

`workflow_argon_gec_ccp_monte_carlo.yaml` contains the complete 15-anchor plan
from 1 to 4000 Td. It uses four independent replicas and 16 workers. Low-field
anchors use small ensembles and long relaxation windows; high-field anchors
use larger ensembles and shorter correlation lags. A separate bounded tail
phase activates only where excitation or ionization statistics require it.

The workflow enforces both per-replica and total particle-barrier ceilings. A
failed statistical gate may produce one bounded follow-up decision through
`decide-mc` and `advance-mc`; the code cannot extend work indefinitely. The
canonical workflow does not name a `reuse_database`, so its result does not
silently depend on a prior campaign.

## Execution sequence

For two term:

```powershell
swarm-workflow sweep examples\workflow_argon_gec_ccp_two_term.yaml
swarm-workflow build-tables outputs\argon_gec_ccp\two_term_function_eedf_restricted_lmea.sqlite --output outputs\argon_gec_ccp\two_term_function_eedf_restricted_lmea_tables --source two_term
swarm-workflow export-comsol outputs\argon_gec_ccp\two_term_function_eedf_restricted_lmea_tables --output outputs\argon_gec_ccp\two_term_function_eedf_restricted_lmea_bundle
swarm-workflow run-gec-ccp comsol_modes\maps\argon_gec_ccp_two_term_function_eedf.yaml --bundle outputs\argon_gec_ccp\two_term_function_eedf_restricted_lmea_bundle\mixture_0000
```

For Monte Carlo:

```powershell
swarm-workflow sweep examples\workflow_argon_gec_ccp_monte_carlo.yaml
swarm-workflow build-tables outputs\argon_gec_ccp\monte_carlo_function_eedf_restricted_lmea.sqlite --output outputs\argon_gec_ccp\monte_carlo_function_eedf_restricted_lmea_tables --source monte_carlo --mc-qualification-profile function_eedf_restricted_lmea
swarm-workflow decide-mc --mc-tables outputs\argon_gec_ccp\monte_carlo_function_eedf_restricted_lmea_tables\mixture_0000 --two-term-tables outputs\argon_gec_ccp\two_term_function_eedf_restricted_lmea_tables\mixture_0000 --attempt 1 --output outputs\argon_gec_ccp\mc_solver_selection.json
swarm-workflow export-comsol outputs\argon_gec_ccp\monte_carlo_function_eedf_restricted_lmea_tables --selection outputs\argon_gec_ccp\mc_solver_selection.json --output outputs\argon_gec_ccp\monte_carlo_function_eedf_restricted_lmea_bundle
swarm-workflow run-gec-ccp comsol_modes\maps\argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml --bundle outputs\argon_gec_ccp\monte_carlo_function_eedf_restricted_lmea_bundle\mixture_0000
```

For the bounded propagator target:

```powershell
swarm-workflow sweep examples\workflow_argon_gec_ccp_propagator.yaml
py -3 tools\qualify_propagator_target.py `
  --target argon_gec_ccp_restricted_local_mean_energy `
  --fields-Td 3000 3500 4000 `
  --operating-bracket-Td 3000 3500 `
  --table-support-cap-Td 4000 `
  --medium-grid 300 48 `
  --fine-grid 600 72 `
  --mixture-id 0 `
  --fine-max-memory-mb 1024 `
  --database outputs\argon_gec_ccp\propagator_function_eedf_restricted_lmea.sqlite `
  --config examples\argon_gec_ccp_propagator.yaml `
  --core-qualification docs\dev\results\propagator_p1_deterministic_qualification_20260908.json `
  --workers 3 `
  --memory-budget-mb 3072 `
  --output docs\dev\results\propagator_gec_ccp_target_qualification_20260915.json
swarm-workflow build-tables outputs\argon_gec_ccp\propagator_function_eedf_restricted_lmea.sqlite --output outputs\argon_gec_ccp\propagator_function_eedf_restricted_lmea_tables --source propagator --solver-qualification docs\dev\results\propagator_p1_deterministic_qualification_20260908.json --target-qualification docs\dev\results\propagator_gec_ccp_target_qualification_20260915.json
swarm-workflow export-comsol outputs\argon_gec_ccp\propagator_function_eedf_restricted_lmea_tables --output outputs\argon_gec_ccp\propagator_function_eedf_restricted_lmea_bundle
swarm-workflow run-gec-ccp comsol_modes\maps\argon_gec_ccp_propagator_function_eedf.yaml --bundle outputs\argon_gec_ccp\propagator_function_eedf_restricted_lmea_bundle\mixture_0000
```

The generic qualification CLI writes
`swarm.propagator_target_qualification.v1`. The GEC-CCP adapter then validates
that artifact against its own explicit fields, bracket, support cap, mixture,
grid pair, and refinement limits during preflight.

Add `--dry-run` to `run-gec-ccp` to inspect the MPH contract, validate the
bundle, and materialize Java without starting COMSOL.

## COMSOL run behavior

A normal physical-target run performs these operations:

1. Validate the schema, bundle hashes, solver-selection record, MPH contract,
   closure support, and source-specific statistical evidence.
2. Apply the external closure to a separate MPH copy.
3. Audit the active Function-EEDF binding when the closure uses one.
4. Run the external time-periodic study.
5. Export the eleven canonical result tables under `swarm_tables/`.
6. Audit conservation, transport, Function-EEDF behavior, support, MC rate
   censoring, and saved-model provenance.

The COMSOL built-in Druyvesteyn calculation is not required to accept an
external target and is skipped by default. This removes one complete periodic
solve and eleven duplicate CSV exports from every production execution. To run
it as an explicit comparison, copy a canonical map and add:

```yaml
model:
  baseline_output_mph: ../work/reference.mph
run:
  include_builtin_reference: true
  period_dataset: dset1
  baseline_waveform_dataset: dset2
```

That opt-in run adds `gec_baseline_run`, `gec_baseline_export`, and the
`builtin_druyvesteyn/` result directory. Cross-solver two-term versus MC plots
compare their accepted external results directly and do not recompute this
reference.

## Canonical outputs and acceptance

The run directory contains:

- `gec_ccp_plan.json`, with immutable input and generated-Java hashes;
- `swarm_tables/`, with the eleven COMSOL CSV exports;
- source-specific audit JSON/CSV files;
- `conservation_audit.json` and any applicable transport/support audits;
- `gec_ccp_run_status.json`, the final acceptance record.

The canonical retained run consists of that result directory, its output MPH,
and the four matching apply, Function-EEDF audit, external-run, and export log
directories. Timestamped failed attempts and superseded preflight material are
diagnostic scratch data; after their cause is represented by a regression test
or maintained contract, they are not part of the product result.

The declared closure is accepted only when the solve completed and every
applicable numerical and physics gate passed. Its terminal record then contains
`quality_accepted_for_declared_closure: true`. A result is promoted separately
to `physical_target_accepted: true` only when the mapping declares a physical
target and supplies every required physical closure. A successful COMSOL
process exit alone is insufficient for either decision.

A `diagnostic_control` run has a separate terminal status and is never
promoted to `physical_target_accepted: true`. Its purpose is to verify the
input mapping and inspect solver sensitivity without presenting a restricted
local-mean-energy closure as full physical validation.

## Scope limits

The current GEC-CCP target does not provide a full kinetic gradient-response
closure, nonlocal RF kinetics, arbitrary crossed electric and magnetic fields,
or state-resolved surface evolution. The wall model remains the shared native
COMSOL drift-diffusion boundary closure. See `docs/physics_limitations.md` for
the wider solver capability limits.
