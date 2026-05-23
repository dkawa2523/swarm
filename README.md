# swarm

Product electron swarm solver comparison tools.

Schema v2 is the only public YAML schema. Select solvers with canonical ids:
`two_term`, `multi_term`, and `monte_carlo`. Legacy `run.mode`, `both`, `all`,
`boltzmann_two_term`, `multiterm_boltzmann`, and legacy output aliases are not
supported.

## Run

```powershell
py -3 -m electron_swarm configs\swarm_product_v2_minimal.yaml
py -3 -m electron_swarm configs\swarm_product_v2_example.yaml --no-write
```

## Minimal Schema

```yaml
schema_version: 2
run:
  solvers:
    - id: two_term
    - id: multi_term
  e_over_n_Td: [50, 100]
physics:
  angular_scattering:
    model: momentum_power
    higher_moment_closure: power
  electron_electron:
    enabled: false
    model: none
feature_policy:
  unsupported: fail
  degraded: warn
solvers:
  multi_term:
    method: pn_closure_surrogate
    lmax: 4
output:
  directory: outputs
  base_name: swarm
```

## Outputs

Canonical files:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

Summary metadata is intentionally compact: solver method, angular model,
ordinary-XS closure status, e-e treatment, magnetic treatment, and tail metrics.

## Physics Scope

`multi_term.method: pn_closure_surrogate` is the runnable ordinary-cross-section
multi-term path. It is an angular-closure/surrogate path, not an exact DCS-based
PN solver. Results record `physics_level=surrogate`,
`exact_dcs_based=false`, and `ordinary_integral_xs_closure=true`.

`pn_closure_direct`, `pn_dcs`, magnetic dynamics, RF/time-dependent fields,
finite-k transport, full e-e Fokker-Planck, Coulomb MC, and YAML-side
state-resolved generation are unsupported or `NotImplemented` until implemented
end to end.

## Tests

```powershell
py -3 -m pytest -q
```
