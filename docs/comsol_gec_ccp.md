# Argon GEC CCP: Swarm table comparison workflow

This path targets `comsol_modes/argon_gec_ccp.mph` without modifying the
source model. It produces two solved copies when a valid COMSOL Plasma Module
license is available:

1. the model's built-in Druyvesteyn closure;
2. the same model with Swarm transport and reaction-rate lookup tables.

The local MPH contract is inspected before Java is generated. The workflow
requires the time-periodic plasma interface `ptp`, plasma feature `pes1`,
electron-impact reactions `eir1` through `eir3`, studies `std1` and `std2`,
and the existing line datasets. A mismatch fails before COMSOL is launched.

## What the source model actually assumes

The checked-in MPH contains `eedf = Druyvesteyn`, not Maxwellian. Its three
electron-impact reactions use cross-section data and inherit that interface
EEDF. The plasma model has a prescribed mobility path
(`SpecifyMueOnly`), while the model solves electron density and mean electron
energy in the time-periodic RF formulation.

For the external case, COMSOL must continue solving mean energy. Setting
`SpecifyMeanElectronEnergy = MeanEnergyTable` would collapse the RF energy
dynamics into a local-field prescription and is intentionally not used.
Instead the workflow sets `SpecifyElectronDensityAndEnergy = UseLookupTables`
and imports these quantities versus mean energy:

- reduced electron mobility `muN`;
- reduced particle diffusion `DN`;
- reduced electron-energy mobility `mueN`
  (`reduced_electron_energy_mobility_m2_V_s_m3`, `1/(V m s)`);
- reduced electron-energy diffusion `DeN`
  (`reduced_electron_energy_diffusion_m2_s_m3`, `1/(m s)`);
- elastic, excitation, and ionization rate coefficients.

Each reaction is switched from `UseCrossSectionData` to `UseLookupTable` with
`RateConstantForm = UseRate`. This is essential: leaving a reaction on
cross-section mode would still let COMSOL evaluate it from the built-in EEDF
and would mix the two closures.

`eedf.csv` and `eedf_f0.csv` are exported for diagnostics. The latter uses
the COMSOL EEPF convention `f0 = EEDF/sqrt(energy)`. Neither file is imposed
as a spatially uniform RF EEDF. A steady DC swarm EEDF is not, in general,
the phase-resolved EEDF of a 13.56 MHz CCP.

## Reproducible commands

From the repository root:

```powershell
swarm-workflow sweep examples\workflow_argon_gec_ccp.yaml
swarm-workflow build-tables outputs\argon_gec_ccp\swarm_schema_v2.sqlite `
  --output outputs\argon_gec_ccp\tables_schema_v2 --source two_term
swarm-workflow export-comsol outputs\argon_gec_ccp\tables_schema_v2 `
  --output outputs\argon_gec_ccp\comsol_bundle_schema_v2
swarm-workflow run-gec-ccp comsol_modes\maps\argon_gec_ccp_swarm.yaml `
  --bundle outputs\argon_gec_ccp\comsol_bundle_schema_v2\mixture_0000 --dry-run
```

The configured sweep contains 28 points from 0.05 to 10000 Td. The
time-periodic solve uses power continuation at 0.1, 0.25, 0.5, and 1 W,
then runs the Time Periodic to Time Dependent conversion study.

To execute and plot:

```powershell
swarm-workflow run-gec-ccp comsol_modes\maps\argon_gec_ccp_swarm.yaml `
  --bundle outputs\argon_gec_ccp\comsol_bundle_schema_v2\mixture_0000 `
  --comsol "C:\Program Files\COMSOL\COMSOL64\Multiphysics\bin\win64\comsol.exe"
swarm-workflow plot-gec-ccp `
  --bundle outputs\argon_gec_ccp\comsol_bundle_schema_v2\mixture_0000 `
  --comsol-results outputs\argon_gec_ccp\comsol `
  --output outputs\argon_gec_ccp\plots
```

If a coarse transport-continuation step fails after saving a converged external
MPH, rerun `run-gec-ccp` with `--resume-external`. This skips the baseline,
loads the last converged external MPH, skips already saved coefficient stages,
advances electron mobility, energy mobility, energy diffusion, and electron
diffusion one at a time, and then activates the elastic, excitation, and
ionization rate tables sequentially. The external case enables COMSOL's source and
reaction-source stabilization for the logarithmic plasma formulation. Its
automatic damped-Newton solver permits damping down to `1e-8` and disables the
otherwise destabilizing minimum-step recovery jump. The built-in Druyvesteyn
baseline remains unchanged.

The plot command always produces a Swarm closure overview. It produces axial
profile and electrode-waveform comparison plots only when both COMSOL result
sets exist; missing COMSOL outputs are recorded as skipped and never replaced
with synthetic data.

## Comparison and acceptance checks

Keep geometry, mesh, wall coefficients, secondary emission, RF frequency,
power-control/self-bias settings, tolerances, and initial conditions identical.
Compare at minimum:

- period-averaged axial and radial electron density;
- mean electron energy, potential, ionization source, and absorbed power;
- powered-electrode voltage and current over an RF period;
- convergence history, periodicity residual, and runtime;
- whether the solved mean energy remains inside the 0.576 to 273.5 eV table
  range produced by the current sweep.

Do not treat the built-in Druyvesteyn result as experimental truth. It is the
controlled reference for isolating closure changes. A stronger validation
needs an external Boltzmann solver such as BOLSIG+ for the same cross sections
and, for the RF nonlocal regime, PIC/MCC or measured GEC data. Ordinary
integral cross sections also do not establish an exact angular-scattering
model, so two-term versus multi-term/Monte Carlo comparisons must state their
angular closure.

COMSOL background:

- [GEC CCP application model](https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.argon_gec_ccp/argon_gec_ccp.html)
- [Electron energy distribution functions](https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_boltzmann.06.09.html)
- [Time-periodic EEDF formulation](https://doc.comsol.com/6.3/doc/com.comsol.help.plasma/plasma_ug_plasma.09.08.html)
