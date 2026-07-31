# electron-swarm

Electron swarm solver comparison tools for low-pressure plasma modeling, with
a small external-table path for COMSOL.

The product schema requires `schema_version: 2` and canonical solver ids
`two_term`, `multi_term`, and `monte_carlo`. Solver selection is independent of
requested physics features, and every result records how those features were
handled.

```powershell
py -3 -m pip install -e ".[dev]"
electron-swarm examples\two_term.yaml --no-write
```

COMSOL workflow:

```powershell
swarm-workflow sweep examples\workflow_argon_comsol.yaml
swarm-workflow build-tables outputs\argon_comsol\swarm.sqlite --output outputs\argon_comsol\tables --source two_term
swarm-workflow export-comsol outputs\argon_comsol\tables\mixture_0000 --output outputs\comsol_bundle\mixture_0000
swarm-workflow run-comsol Model\maps\positive_column_external.yaml --bundle outputs\comsol_bundle\mixture_0000 --dry-run
```

The 1D positive-column COMSOL path injects seven quantities: mean electron
energy, reduced mobility, longitudinal reduced diffusion, reduced
electron-energy mobility (`reduced_electron_energy_mobility_m2_V_s_m3`,
`1/(V m s)`), reduced electron-energy diffusion
(`reduced_electron_energy_diffusion_m2_s_m3`, `1/(m s)`), and excitation and
ionization Townsend coefficients. COMSOL retains the spatial balance
equations, electrostatics, walls, heavy-species chemistry, geometry, mesh, and
nonlinear solver. Seven written arrays do not imply seven active closures:
the canonical mapping uses COMSOL's local-energy formulation, where the
E/N-to-mean-energy table is inactive and the other six arrays are only
candidate active until a clean COMSOL run verifies activation. See
`docs/comsol_workflow.md`,
`docs/comsol_benchmark.md`, and
`reports/comsol_swarm_benchmark_2026/`.

The time-periodic argon GEC CCP path is separate because COMSOL must keep mean
electron energy as an RF-dependent solved variable:

```powershell
swarm-workflow sweep examples\workflow_argon_gec_ccp.yaml
swarm-workflow build-tables outputs\argon_gec_ccp\swarm_schema_v2.sqlite --output outputs\argon_gec_ccp\tables_schema_v2 --source two_term
swarm-workflow export-comsol outputs\argon_gec_ccp\tables_schema_v2 --output outputs\argon_gec_ccp\comsol_bundle_schema_v2
swarm-workflow run-gec-ccp comsol_modes\maps\argon_gec_ccp_swarm.yaml --bundle outputs\argon_gec_ccp\comsol_bundle_schema_v2\mixture_0000 --dry-run
```

The checked-in model uses Druyvesteyn, not Maxwellian. See
`docs/comsol_gec_ccp.md` for the exact mapping and comparison procedure.

```powershell
py -3 -m pytest -q
```
