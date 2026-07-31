# Local COMSOL models

Licensed `.mph` inputs and generated work files stay local and are ignored by
git. The product mapping is `maps/positive_column_external.yaml`.

Required local files:

```text
Model/positive_column_1d.mph
Model/positive_column_1d_boltzmann.mph
Model/work/positive_column_1d_swarm_template.mph
```

Create the template as a copy of `positive_column_1d.mph`; do not overwrite the
Application Library source. The mapping writes seven x-y lookup quantities
(14 COMSOL property arrays), applies the formulation-relevant subset, then
saves its result under `Model/work/`.

Run the complete external-table path with:

```powershell
swarm-workflow run-comsol Model\maps\positive_column_external.yaml --bundle outputs\comsol_bundle\mixture_0000
```

The external bundle supplies electron mobility, longitudinal diffusion,
electron-energy mobility, electron-energy diffusion, and selected direct
reaction Townsend coefficients. COMSOL evaluates these coefficients inside
the electron and mean-energy flux equations and remains responsible for the
spatial PDEs, residual chemistry, heavy species, wall and surface reactions,
electrostatics, circuit, geometry, mesh, and solver settings. In the mapped
local-energy formulation, the written E/N-to-mean-energy lookup is inactive.
