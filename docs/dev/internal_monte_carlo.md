# Internal Monte Carlo Backend

The internal `monte_carlo` backend is a lightweight product backend for
electron swarm comparisons. It is not an embedded MCIG clone and does not emit
validated bulk transport.

## Angular Scattering

Supported internal samplers:

- `isotropic`
- `maxent_p1`

`maxent_p1` requires both total/effective elastic and momentum-transfer cross
sections for every active scattering species. If either side is missing, the
backend fails explicitly instead of falling back to isotropic scattering.

## EEDF Estimator

The EEDF is a normalized density in `1/eV`:

```text
sum(eedf * energy_width_eV) == 1
```

Internal MC records time-residence probability with midpoint energy for each
orbit substep. The canonical EEDF CSV includes bin widths plus nullable
`sample_count`, `effective_sample_count`, and `relative_standard_error`.

## Transport Meaning

Internal MC reports flux-like particle-tracking estimates:

- drift: total weighted displacement over total weighted time
- longitudinal diffusion: residual axial variance over total weighted time
- transverse diffusion: transverse variance over total weighted time

These are not bulk/source-gradient coefficients. In ionizing swarms,
fixed-particle and weighted-branching population models represent different
observables and should not be mixed without checking `transport_definition`.

## Development Audit

Detailed MC diagnostics are available through `tools/benchmark_internal_mc_audit.py`
and `tools/benchmark_internal_mc_physics_triage.py`. Product summary outputs keep
only interpretation metadata.
