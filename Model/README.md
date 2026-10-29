# Local COMSOL models

The three original `.mph` inputs in this directory are distributed through Git
LFS. Generated work files stay local and are ignored by git. See the
[model inventory and download instructions](../comsol_models/README.md) for the
distinction between original models, historical configured models, and current
inputs whose MPH application has not completed. The product mapping is
`maps/positive_column_external.yaml`.

Required local files:

```text
Model/positive_column_1d.mph
```

Keep `positive_column_1d.mph` immutable. The mapping writes seven x-y lookup quantities
(14 COMSOL property arrays), applies the formulation-relevant subset, then
saves its result under `Model/work/`.

The map requires `bundle.source: two_term`, a gas temperature and pure-Ar
composition matching the bundle, and explicit species/process/type/threshold
identity for each direct reaction. This adapter currently supports the full
unprocessed two-term closure only. Other solver closures require their own
qualified model contract. Dry-run and execution reject altered table hashes,
incomplete source provenance, and failed or incomplete quality CSV evidence.

Run the complete external-table path with:

```powershell
swarm-workflow run-positive-column Model\maps\positive_column_external.yaml
```

The mapping declares the bundle path. Use `--bundle` only for an explicit
replacement; the effective mapping records the resolved input.

`results` selects `output_directory` and requires the saved `dataset`, explicit
`solutions: all|last`, and `fields` with names, expressions, and units.
The maintained map uses `dset1` and `last` for the final transient state.
Voltage continuation is not a parameter sweep. Result CSVs retain solution,
point, and time identities; arbitrary coordinate/time requests are rejected.

The external bundle supplies electron mobility, longitudinal diffusion,
electron-energy mobility, electron-energy diffusion, and selected direct
reaction Townsend coefficients. COMSOL evaluates these coefficients inside
the electron and mean-energy flux equations and remains responsible for the
spatial PDEs, residual chemistry, heavy species, wall and surface reactions,
electrostatics, circuit, geometry, mesh, and solver settings. In the mapped
local-energy formulation, the written E/N-to-mean-energy lookup is inactive.
