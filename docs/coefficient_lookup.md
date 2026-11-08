# Use Swarm coefficients in another simulator

`build-tables` produces qualified DC coefficients on the original E/N grid.
Use one `mixture_0000` directory with its manifest and CSVs. Temperature,
number density, and mole fractions must match the external model's fixed
background gas. The tables also record the solver and its physical assumptions.

This example reads the distributed pure-Ar GEC-ICP Swarm result at 300 K and
2.66644 Pa. It works with the `electron-swarm` package installed, from the
repository root after `git lfs pull`.

```powershell
New-Item -ItemType Directory -Force outputs\global_argon | Out-Null
Copy-Item comsol_models\gec_icp\swarm_results.sqlite outputs\global_argon\swarm.sqlite
py -3 -m swarm_workflow.cli build-tables outputs\global_argon\swarm.sqlite --source two_term --output outputs\global_argon\tables
py -3 -m examples.lookup_argon_coefficients outputs\global_argon\tables\mixture_0000 --e-over-n-Td 45
```

Use `swarm-workflow run comsol_models/gec_icp/swarm.yaml` to regenerate the
source with the current implementation when the saved database requires migration.
For another gas, build from that case's schema-v2 YAML and workflow database.
Use `swarm-workflow export-reaction-rates` for reaction coefficient points
from any of the four solver modes, including `multi_term`.

```python
from swarm_workflow.tables import Reaction, load_coefficients

tables = load_coefficients(
    "outputs/global_argon/tables/mixture_0000",
    gas_temperature_K=300.0,
    gas_number_density_m3=2.66644 / (1.380649e-23 * 300.0),
    mole_fractions={"Ar": 1.0},
)
ionization = tables.series(
    "rate_coefficient_m3_s",
    reaction=Reaction("Ar", "e+Ar=>2e+Ar+", "ionization", 15.8),
)
rate_m3_s = ionization.at(45.0)
```

- `tables.reactions` lists the exact reaction keys. The two-body source per
  volume is `n_e * n_target * k`, using the individual rate in `m^3/s`.
- `tables.series(...)` also reads canonical mean-energy and transport columns
  and `elastic_energy_loss_rate_coefficient_eV_m3_s` where available. Consult
  `series.unit` and `tables.source_policy` for each quantity's units and
  transport or energy-loss definition. The source policy also records the
  physical assumptions of the DC sweep. The elastic energy exchange retains its sign.
- An exact E/N anchor preserves the CSV value. Between adjacent qualified
  anchors the reader interpolates linearly in Td. Out-of-range requests,
  missing values, excluded anchors and censored MC rates raise
  `CoefficientLookupError`. Restricted MC tables expose mean energy,
  reduced mobility, and signed elastic loss through this reader.

Keep the `manifest.json` with the CSVs when moving the tables. Swarm supplies
coefficients at the specified gas state. The consuming plasma model advances
its densities, species composition, and energy balances.
