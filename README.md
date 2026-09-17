# electron-swarm

Electron-swarm solver comparison tools with qualified external-Swarm inputs
for COMSOL plasma models.

## Product contract

Product YAML uses `schema_version: 2` and selects `two_term`, `multi_term`,
`monte_carlo`, or `propagator` through `run.solvers`. Physics features are
requested separately; unsupported or degraded treatment is always recorded.
See [the schema](docs/product_schema_v2.md) and [physics
limits](docs/physics_limitations.md).

Minimal run:

```powershell
py -3 -m electron_swarm.runner examples\two_term.yaml
```

This writes `two_term_summary.csv`, `two_term_eedf.csv`,
`two_term_rates.csv`, and the solver-plan CSV under `outputs/examples/`.
Use `examples/compare_three_solvers.yaml` when a direct multi-model comparison
is wanted; Monte Carlo examples intentionally declare their finite sampling
budgets separately.

`propagator` is a deterministic energy--polar-angle P1 solver for homogeneous
positive-E/N DC cases with `B=0`. Its declared P0/P1 core is qualified; it
returns EEDF, energy-angle population, rates, growth, flux drift/mobility, and
same-operator elastic energy loss. It does not return bulk drift, particle
diffusion, or electron-energy transport. The core qualification evidence is
[machine-readable](docs/dev/results/propagator_p1_deterministic_qualification_20260908.json);
each consuming plasma-model adapter owns an explicit target requirement. The
GEC-CCP high-field range is covered by a separate generic-schema
[target artifact](docs/dev/results/propagator_gec_ccp_target_qualification_20260915.json).

## GEC-CCP routes

| Solver | Workflow | COMSOL map | Role |
| --- | --- | --- | --- |
| Two term | `examples/workflow_argon_gec_ccp_two_term.yaml` | `comsol_modes/maps/argon_gec_ccp_two_term_function_eedf.yaml` | physical target |
| Monte Carlo | `examples/workflow_argon_gec_ccp_monte_carlo.yaml` | `comsol_modes/maps/argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml` | qualified physical target |
| Propagator | `examples/workflow_argon_gec_ccp_propagator.yaml` | `comsol_modes/maps/argon_gec_ccp_propagator_function_eedf.yaml` | deterministic P1 physical target |

Two term and Propagator supply an active Function EEDF. The MC route retains
its Function-EEDF grid as source evidence and supplies inelastic rates
preintegrated from that same distribution. All three supply reduced mobility
and solver-native elastic energy loss. COMSOL owns the RF solve, heavy species,
boundaries, Einstein particle diffusion, and restricted local-energy transport.

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
