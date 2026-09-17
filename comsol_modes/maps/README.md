# COMSOL model maps

## GEC-CCP

This directory contains three physical-target mappings:

- `argon_gec_ccp_two_term_function_eedf.yaml`
- `argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml`
- `argon_gec_ccp_propagator_function_eedf.yaml`

All three import an external Function EEDF, reduced mobility, and solver-native
elastic energy loss. COMSOL retains Einstein particle diffusion, electron
energy transport, the RF solve, heavy species, and boundary physics. The
two-term map lets COMSOL integrate inelastic cross sections. The Monte Carlo
map imports inelastic rates preintegrated from the same MC EEDF because its
earlier direct Function-EEDF target was not a robustly converged production
closure.
The propagator map also uses COMSOL Function-EEDF rate integration. Its
finite-volume EEDF is imported as a compact structured spreadsheet. The map
is a physical target only when both maintained Propagator qualification
artifacts and all declared downstream closure audits pass.

A normal run solves and audits only `swarm_tables`. To request an additional
COMSOL built-in Druyvesteyn reference, add all four fields below to a copy of
a map:

```yaml
model:
  baseline_output_mph: ../work/reference.mph
run:
  include_builtin_reference: true
  period_dataset: dset1
  baseline_waveform_dataset: dset2
```

See `docs/comsol_gec_ccp.md` for the full input ownership and command sequence.

## GEC-ICP

`argon_gec_icp_model.yaml` records the immutable topology of the repository
2D axisymmetric ICP MPH. The following source maps layer a qualified Swarm
bundle, one cold-start operating point, and source-specific outputs over it:

- `argon_gec_icp_two_term_function_eedf.yaml`
- `argon_gec_icp_monte_carlo_function_eedf_restricted_lmea.yaml`
- `argon_gec_icp_propagator_function_eedf.yaml`

All three use a Function-EEDF plus reduced-mobility restricted LMEA. COMSOL
owns diffusion, energy transport, cross-section integrals, argon/metastable
chemistry, walls, and the frequency-transient inductive field. The maps run no
COMSOL coefficient sweep. See `docs/comsol_gec_icp.md` for the exact ownership,
qualification, execution, and convergence contract.
