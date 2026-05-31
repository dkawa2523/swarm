# Ar EEDF Consistency Benchmark

This benchmark checks the Ar / BOLSIG-style control case where `two_term` is the
native Scharfetter-Gummel reference and `multi_term pn_closure_direct lmax: 1`
must reduce to the same EEDF convention.

Run the lightweight benchmark:

```powershell
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency.yaml
```

Outputs are written under `outputs/benchmarks`:

- `ar_eedf_consistency_eedf_metrics.csv`
- `ar_eedf_consistency_failure_analysis.csv`
- `ar_eedf_consistency_report.md`

The report status is the direct `lmax: 1` gate status. Diagnostic rows for
higher-l stability or light Monte Carlo statistics do not make the direct gate
fail.

## Gate

The direct gate compares normalized EEDF `F(E)` in `1/eV`, not EEPF:

- `eedf_relative_l1 < 0.01`
- `mean_energy_relative_difference < 0.005`
- `drift_velocity_relative_difference < 0.01`
- `major_rate_relative_difference < 0.02`
- normalization error `< 1e-8`
- `negative_mass_fraction < 1e-8`

`pn_closure_direct lmax: 1` uses the shared kinetic projection data and solves a
coupled sparse SG-reduction block. It does not copy the `two_term` EEDF or
reuse `two_term` rates.

## Higher lmax

`pn_closure_direct` values above `lmax: 1` run in the limited B=0, DC,
axisymmetric ordinary-XS angular-closure scope. The benchmark records
lmax-to-lmax EEDF differences, residuals, negative mass, and any instability in
the failure-analysis CSV.

## Monte Carlo

The default benchmark keeps MC intentionally light and can report
`MC_statistical_uncertainty`. A heavier manual config is available:

```powershell
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency_mc_manual.yaml
```

This manual run is not part of default pytest or CI.
