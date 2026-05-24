# Testing and Benchmarks

Phase 15 keeps the default developer path lightweight while adding explicit
markers for optional checks.

## Test Commands

Default local run:

```powershell
py -3 -m pytest -q
```

CI-equivalent fast run:

```powershell
py -3 -m pytest -q -m "not slow and not mc"
```

Focused numerical regression:

```powershell
py -3 -m pytest -q -m regression
```

Optional Monte Carlo checks:

```powershell
py -3 -m pytest -q -m mc
```

`slow`, `mc`, and `regression` are pytest markers. Lightweight regression tests
stay in CI; stochastic or heavier MC checks are marked `mc`.

## Regression Baselines

Regression baselines are small JSON files under `tests/regression_baselines`.
They capture accepted product behavior for representative scalar outputs with
deliberately loose tolerances. Update a baseline only when a solver or physics
change intentionally changes the accepted numerical behavior.

## Benchmarks

Run one product benchmark config:

```powershell
py -3 benchmarks\run_product_benchmarks.py --config configs\benchmarks\ar_simple.yaml
```

Run all benchmark configs:

```powershell
py -3 benchmarks\run_product_benchmarks.py
```

The benchmark runner writes only `outputs/benchmarks/benchmark_summary.csv`.
Canonical solver CSV outputs are not written during benchmarks.

The focused Ar EEDF consistency benchmark is:

```powershell
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency.yaml
```

It writes `ar_eedf_consistency_eedf_metrics.csv`,
`ar_eedf_consistency_failure_analysis.csv`, and
`ar_eedf_consistency_report.md` under `outputs/benchmarks`.
`pn_closure_direct` is intentionally excluded until a true coupled f0/f1 block
PN implementation exists.
