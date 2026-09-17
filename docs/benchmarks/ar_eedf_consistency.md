# Ar Boltzmann Consistency Benchmark

This benchmark compares the independent `two_term` and `multi_term` models on
an Ar control case and records the convergence of the PN truncation.

```powershell
py -3 tools\benchmarks\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency.yaml
```

The normalized EEDF, mean energy, drift, and rates are compared across models.
Those differences are diagnostic: the multi-term result is not copied from or
repaired against the two-term result. The multi-term numerical gate instead
requires a normalized nonnegative F0, a converged complete-PN residual, and
stable changes as `lmax` is raised.

`lmax: 1` and higher orders use the same staggered PN operator. The reported
drift is always the solved F1 velocity moment. The benchmark records
lmax-to-lmax EEDF changes, full residuals, and negative mass.

A light Monte Carlo row can be included for diagnostic comparison. A heavier
manual configuration is available at
`configs/benchmarks/ar_bolsig_eedf_consistency_mc_manual.yaml`; it is not part
of default tests.
