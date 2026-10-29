# electron-swarm

Electron-swarm solver comparison tools, simulator-neutral reaction coefficients,
and numerically qualified input workflows for COMSOL plasma models.

[Download COMSOL models through Git LFS](comsol_models/README.md). The inventory
separates original MPH files, historical Swarm-configured ICP models, and
current inputs that have not yet been applied to an MPH.

## Run a comparison

Product YAML uses `schema_version: 2` and selects `two_term`, `multi_term`,
`monte_carlo`, or `propagator` through `run.solvers`. Physics features are
requested separately; unsupported or degraded treatment is always recorded.
See [the schema](docs/product_schema_v2.md) and [physics
limits](docs/physics_limitations.md).

Minimal run:

```powershell
py -3 -m swarm_workflow.cli run examples\two_term.yaml
```

Use one Swarm YAML for direct runs and campaigns. Add `campaign.database` for
persistence/resume and `campaign.mc.replicas` for independent MC replicas;
E/N, mixtures, and `solvers.monte_carlo.seed` have one input owner.
Add `--plan` for a concise preview or `--plan --json` for effective settings,
physical cases, sampling units, and resource bounds. See
[configuration and execution](docs/swarm_inputs.md).
Optional `output.export: {solver: multi_term, quantities: [reaction_rates]}`
adds simulator-neutral DC coefficient points with physical context and MC
replica uncertainty, independently of transport or COMSOL requirements.
Use `swarm-workflow export-reaction-rates` to produce the same data from saved results.
With `pip install 'electron-swarm[plot]'`, add
`output.figures: [eedf, transport, reaction_rates]` to the same YAML for
PNG/SVG and their numerical data. Optional `energy_angle` and `rf_phase`
figures use native solver results. `swarm-workflow plot <manifest-or-db>
--figures eedf transport --output <directory>` redraws saved data without
running a solver. Chemical-species time evolution belongs to the plasma model.
The run writes summary, EEDF, rates, solver-plan CSVs and a manifest under
`outputs/examples/`. The manifest records settings, input hashes, implementation
identity, input attribution, and output hashes. `electron_swarm.io.read_product_run`
returns the verified `.metadata` and `.tables` together, preserving numeric
precision and text identifiers. A run with no executed solver reports `skipped`,
returns CLI exit code 3, and writes only the plan and manifest.
A partially executed request returns exit code 4; disabled solvers are excluded
from requested coverage. Changed results are published as derived data with
their immutable source execution, rather than carrying forward its success.
Use `examples/compare_three_solvers.yaml` when a direct multi-model comparison
is wanted; Monte Carlo examples intentionally declare their finite sampling
budgets separately. Comparison availability, pair coverage, and numerical qualification are
reported separately. Rows include physical-condition matching, numerical status, and
per-observable meanings. Gas masses are explicit. Ionization population physics
is selected under `physics.ionization.population_model` (default: two-daughter
weighted growth). A direct MC run needs independent replicas for uncertainty;
`swarm-workflow` supplies that campaign and qualification workflow.

`propagator` is a deterministic energy--polar-angle P1 solver for homogeneous
positive-E/N DC cases with `B=0`. Its declared P0/P1 core is qualified; it
returns EEDF, energy-angle population, rates, growth, flux drift/mobility, and
same-operator elastic energy loss. It does not return bulk drift, particle
diffusion, or electron-energy transport. The core qualification evidence is
[machine-readable](docs/dev/results/propagator_p1_core_qualification.json);
each consuming plasma-model adapter owns an explicit target requirement. The
GEC-CCP high-field range is covered by a separate generic-schema
[target artifact](docs/dev/results/propagator_gec_ccp_target_qualification.json).

## GEC-CCP routes

| Solver | Workflow | COMSOL map | Role |
| --- | --- | --- | --- |
| Two term | `examples/workflow_argon_gec_ccp_two_term.yaml` | `comsol_modes/maps/argon_gec_ccp_two_term_function_eedf.yaml` | two-term input |
| Monte Carlo | `examples/workflow_argon_gec_ccp_monte_carlo.yaml` | `comsol_modes/maps/argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml` | MC ensemble-qualified input |
| Propagator | `examples/workflow_argon_gec_ccp_propagator.yaml` | `comsol_modes/maps/argon_gec_ccp_propagator_function_eedf.yaml` | P1 grid-qualified input |

Two term and Propagator supply an active Function EEDF. The MC route retains
its Function-EEDF grid as source evidence and supplies inelastic rates
preintegrated from that same distribution. All three supply reduced mobility
and solver-native elastic energy loss. COMSOL owns the RF solve, heavy species,
boundaries, Einstein particle diffusion, and restricted local-energy transport.
Source qualification, independent reference validation, and licensed COMSOL
model acceptance are separate steps.

The exact sweep, table, bundle, dry-run, and execution commands are maintained
once in [the GEC-CCP contract](docs/comsol_gec_ccp.md). A built-in Druyvesteyn
reference is opt-in and is not repeated during a normal external-Swarm run.

COMSOL integration is not limited to GEC-CCP. Model-independent input and
execution services are reused by independent model adapters; the supported
adapters and the contract for adding another plasma model are documented in
[COMSOL plasma-model adapters](docs/comsol_model_adapters.md).

## Install and test

```powershell
py -3 -m pip install -e ".[dev]"
py -3 -m pytest -q
```
