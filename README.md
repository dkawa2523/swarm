# swarm

Electron swarm simulation tools for low-pressure plasma modeling. The repository
contains the legacy particle Monte Carlo workflow and a newer unified
`electron_swarm` workflow that can run multiple solvers through one YAML schema.

## Solvers

- `monte_carlo`: existing particle Monte Carlo implementation in `swarm_mc`.
- `boltzmann_two_term`: native BOLSIG-like two-term Boltzmann backend.
- `multiterm_boltzmann`: axisymmetric DC, B=0, m=0 entry point. The default
  path is a moment-closure EEDF estimate; `method: operator, lmax > 1`
  enables the production sparse operator for integral cross-section input in
  this B=0/DC/m=0 scope. Set `hydrodynamic: true` to request finite-k
  bulk/source-gradient transport.

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

Sparse operator path (`lmax > 1`, flux by default):

```powershell
py -3 -m electron_swarm configs\unified\multiterm_operator.yaml
```

`configs\unified\multiterm_operator_experimental.yaml` is kept as a
compatibility example for earlier workflows.

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
`multiterm_boltzmann` schema. `lmax > 1` is executable as a production flux
sparse operator for this B=0/DC/m=0 integral cross-section scope with
`meta_physical_validity=operator_flux_b0_dc_m0_integral_cross_sections`.
When `hydrodynamic: true`, finite-k extraction must succeed before production
bulk/source-gradient transport is emitted.
The shared native two-term boundary for this path is
`BoltzmannTwoTermSolver.solve_native_distribution()` for the solved EEDF and
`BoltzmannTwoTermSolver.assemble_native_operator_block()` for the reusable
Scharfetter-Gummel energy-space matrix. These APIs expose the grid, quadrature
widths, collision projection, diagnostics, and metadata before public result
postprocessing, so future operator work has a stable regression anchor.
The `l=0..lmax` sparse system uses `LegendreBlockLayout` for explicit
term/energy indexing. The executable operator system keeps an E=0 `l=0`
collision block, momentum/effective plus inelastic sink relaxation for `l>0`,
and a conservative upwind finite-volume electric-field coupling in energy
space. Anisotropic inelastic source terms remain an integral-cross-section
approximation because differential scattering data are outside the current
input model.
`build_operator_assembly_diagnostics()` fixes the developer contract for the
flattened coefficient order and the `integral f0 dE = 1` normalization row.
Metadata names containing `operator_assembly_scaffold_*` are legacy CSV names;
they now describe the executable operator matrix path, not a separate scaffold.
Operator summaries also include lightweight quality indicators such as
`meta_operator_tail_rate_fraction` and
`meta_operator_highest_l_relative_l1`; these help decide when to extend the
energy grid or rerun with a larger `lmax` without adding a separate benchmark
workflow.
For compact multi-term validation across Ar, Ar/N2, high-E/N, attachment, and
superelastic smoke cases, run:

```powershell
py -3 tools\validate_multiterm_operator.py --quick
```

Add `--with-bolos` on environments where the optional BOLOS package is
installed.

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
coefficients. Operator runs emit standard bulk/source-gradient transport only
when `hydrodynamic: true` and the finite-k fit succeeds. Moment-closure
reference bulk values remain isolated in `meta_estimated_*` columns.

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
- flux transport output from the moment-closure estimate or sparse operator
- finite-k bulk/source-gradient output for `method: operator` when
  `hydrodynamic: true`

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
py -3 tools\validate_against_bolos.py configs\unified\boltzmann_only.yaml
```
