# swarm

Electron swarm simulation tools for low-pressure plasma modeling. The repository
contains the legacy particle Monte Carlo workflow and a newer unified
`electron_swarm` workflow that can run multiple solvers through one YAML schema.

## Solvers

- `monte_carlo`: existing particle Monte Carlo implementation in `swarm_mc`.
- `boltzmann_two_term`: native BOLSIG-like two-term Boltzmann backend.
- `multiterm_boltzmann`: axisymmetric DC, B=0, m=0 entry point. The default
  path is a moment-closure EEDF estimate. `method: operator, lmax: 1` is the
  two-term reference path. `method: operator, lmax > 1` is experimental and
  requires `allow_experimental_operator: true` before it can run.

`run.mode: both` keeps the historical pairing of `boltzmann_two_term` and
`monte_carlo`. Use `run.mode: all` to run all three solvers.

## Setup

With `uv`:

```powershell
uv venv
uv sync
```

Without `uv`, use the Python launcher directly:

```powershell
py -3 -m pip install numpy scipy pandas matplotlib pyyaml json5 molmass lxcat_data_parser pytest
```

## Unified Workflow

Boltzmann two-term only:

```powershell
py -3 -m electron_swarm configs\unified\boltzmann_only.yaml
```

Multi-term Boltzmann only:

```powershell
py -3 -m electron_swarm configs\unified\multiterm_boltzmann.yaml
```

Operator `lmax=1` compatibility path:

```powershell
py -3 -m electron_swarm configs\unified\multiterm_operator_lmax1.yaml
```

Reference-anchored lmax>1 validation path:

```powershell
py -3 -m electron_swarm configs\unified\multiterm_operator_experimental.yaml
```

This YAML opts in with `allow_experimental_operator: true`. Treat it as a
development and benchmark target, not the default production path.

All unified solvers:

```powershell
py -3 -m electron_swarm configs\unified\all_template.yaml
```

Use `--no-write` to run without writing CSV or plot outputs.

```powershell
py -3 -m electron_swarm configs\unified\multiterm_boltzmann.yaml --no-write
```

## Input Schema

Unified YAML files are organized into:

- `run`: solver mode, E/N values in Td, case prefix.
- `conditions`: gas temperature, pressure or number density, gas mixture, and
  optional `length_scale_m` for locality warnings.
- `cross_sections`: canonical CSV or legacy `txt|lxcat|bolsig|dat` files,
  plus an explicit `high_energy_extrapolation` rule (`zero`, `hold`, or
  `error`) for energies above the tabulated data.
- solver-specific sections: `boltzmann_two_term`, `multiterm_boltzmann`,
  `monte_carlo`.
- `output`: output directory, base file name, plot/rate/EEDF switches, legacy
  compatibility aliases.

The multi-term solver uses the existing unified cross-section loader. It does
not use the standalone parser from `swarm_bolsig_beyond_solver.zip` as the
official repository input path.

For multi-term runs, cross sections are normalized before projection. If an
`EFFECTIVE` cross section exists for a species, same-species `ELASTIC` and
`MOMENTUM` entries are not double-counted in the transport momentum frequency.
The `operator` backend has compact assembly diagnostics and an `lmax=1`
two-term comparison harness. `method: operator, lmax: 1` delegates to the validated
two-term Scharfetter-Gummel backend and returns the result through the
`multiterm_boltzmann` schema. `lmax > 1` currently uses a
reference-anchored integral-cross-section closure: f0, rates, and flux
transport come from the validated two-term Scharfetter-Gummel solve, while
higher Legendre coefficients are bounded diagnostics based on an integral
momentum relaxation closure. This avoids presenting an underdetermined
integral-cross-section problem as a full multi-term differential-scattering
solve. `hydrodynamic: true` is recorded as requested, but finite-k
bulk/source-gradient transport is not emitted for this anchored lmax>1 path.
The shared native two-term boundary for this path is
`BoltzmannTwoTermSolver.solve_native_distribution()` for the solved EEDF and
`BoltzmannTwoTermSolver.assemble_native_operator_block()` for the reusable
Scharfetter-Gummel energy-space matrix. These APIs expose the grid, quadrature
widths, collision projection, diagnostics, and metadata before public result
postprocessing, so future operator work has a stable regression anchor.
There is no public direct sparse `l>1` Legendre block solver in this phase:
integral cross sections do not determine the higher-order differential
scattering moments needed for one. A future implementation should first prove
that the general block system reduces to the native two-term block for
`lmax=1`, then add anisotropic collision/source blocks behind the benchmark
gate.
Operator summaries also include lightweight quality indicators such as
`meta_operator_tail_rate_fraction` and
`meta_operator_highest_l_relative_l1`; these help decide when to extend the
energy grid or rerun with a larger `lmax` without adding a separate benchmark
workflow.
For a compact executable smoke sweep of the experimental `lmax > 1` closure
across Ar, Ar/N2, high-E/N, attachment, and superelastic cases, run:

```powershell
py -3 tools\validate_multiterm_operator.py --quick
```

This confirms numerical execution and trend reporting only; it is not a physics
validation gate. For the compact reference gate against native two-term, optional
BOLOS, and optional MC, run:

```powershell
py -3 tools\benchmark_operator_gate.py --quick
```

Use `--require-bolos` only on environments where BOLOS is installed and should
be mandatory. Use `--with-mc` for the slower stochastic MC comparison.

## Optional Collision Extensions

The unified runner supports small optional collision extensions without adding
solver-specific branches to `runner.py`.

- `state_resolved` can generate superelastic processes from existing excitation
  cross sections and `species_states` populations. Generated superelastic
  cross sections use a small low-energy floor so detailed-balance singularities
  do not dominate grid-center quadrature.
- `molecular_states.vibrational.levels` is a shorthand for state populations
  that later state-resolved transitions can reference.
- `electron_electron` currently provides a postprocess relaxation model for
  Boltzmann-family EEDFs only. It is not a full Coulomb/Fokker-Planck operator;
  rates and mean energy are recomputed after relaxation, while transport
  coefficients remain marked as stale.
- `energy_grid.refine.enabled: true` optionally adds points around process
  thresholds. It is off by default, so existing grids and results are unchanged.

Example:

```powershell
py -3 -m electron_swarm examples\state_resolved_electron_electron.yaml --no-write
```

## Output Files

Unified outputs:

- `<base>_summary.csv`
- `<base>_eedf.csv`
- `<base>_rates.csv`
- comparison plots when `output.write_plots: true`

Legacy compatibility outputs:

- `summary_mc.csv`, `summary_boltzmann.csv`, `summary_multiterm.csv`
- `eedf_table_mc.csv`, `eedf_table_boltzmann.csv`, `eedf_table_multiterm.csv`
- `energy_table_mc.csv`, `energy_table_boltzmann.csv`,
  `energy_table_multiterm.csv`
- aliases `summary.csv`, `eedf_table.csv`, `energy_table.csv` for
  `output.compatibility.primary_solver`

Flat columns such as `drift_velocity_m_s` and `diffusion_L_m2_s` contain flux
coefficients. The lmax>1 anchored operator does not emit standard
bulk/source-gradient transport. Moment-closure reference bulk values remain
isolated in `meta_estimated_*` columns.

## COMSOL Export

Existing COMSOL export remains compatible with legacy aliases and now also
recognizes multi-term fallback tables. When multiple solver-specific tables
are present and no `summary.csv` alias exists, set
`comsol_export.input.primary_solver` explicitly.

```powershell
py -3 -m electron_swarm configs\unified\multiterm_boltzmann.yaml
```

Then run the exporter against the output directory using the existing
`swarm_comsol_exporter` workflow.

## Physical Scope

The `multiterm_boltzmann` operator scope is:

- uniform DC electric field
- B = 0
- axisymmetric m = 0 Legendre expansion
- elastic, momentum/effective, excitation, ionization, attachment, and
  superelastic process categories
- flux transport output from the moment-closure estimate or the
  reference-anchored lmax>1 closure
- no validated finite-k bulk/source-gradient output for lmax>1 until a
  differential-collision multi-term operator is implemented and benchmarked

Not included in this MVP:

- arbitrary E-B angle
- m != 0 full multi-harmonic solver
- RF/time-dependent fields
- full differential cross-section anisotropic scattering
- state-resolved plasma chemistry ontology
- finite-k transport outside the local hydrodynamic validity range

## Tests

```powershell
py -3 -m pytest tests/test_boltzmann_two_term.py tests/test_unified_integration.py tests/test_multiterm_boltzmann.py -q
```

Optional BOLOS comparison remains available through:

```powershell
py -3 tools\benchmark_operator_gate.py --quick --require-bolos
```
