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
`pn_closure_direct` `lmax: 1` is included as the strict direct reduction gate.
Higher direct PN `lmax` rows are included as limited-scope angular-closure
smoke/regression diagnostics.
The report status is the direct `lmax: 1` gate status; higher-l instability and
light-MC statistical uncertainty remain in the failure CSV as diagnostic rows.
See `docs/benchmarks/ar_eedf_consistency.md` for the gate definition and output
interpretation.

The Ar BOLSIG+ equivalence benchmark compares `two_term` and
`multi_term pn_closure_direct lmax=1` against ingested BOLSIG+ output:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig
```

It writes `ar_bolsig_plus_equivalence_summary.csv`,
`ar_bolsig_plus_equivalence_eedf_metrics.csv`,
`ar_bolsig_plus_equivalence_failure_analysis.csv`, and
`ar_bolsig_plus_equivalence_report.md`.  See
`docs/benchmarks/ar_bolsig_plus_equivalence.md`.

The MCIG reference benchmark uses the same external-reference CLI:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig
```

It writes `ar_mcig_reference_summary.csv`,
`ar_mcig_reference_eedf_metrics.csv`,
`ar_mcig_reference_failure_analysis.csv`, and
`ar_mcig_reference_report.md`.  See
`docs/benchmarks/ar_mcig_reference.md`.

The combined BOLSIG+ / MCIG triage benchmark runs the product solvers once and
classifies disagreements across BOLSIG+, MCIG, `two_term`, direct `multi_term`,
and optional internal MC:

```powershell
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml
```

It writes `ar_triage_matrix.csv`, `ar_triage_eedf_metrics.csv`,
`ar_triage_transport_metrics.csv`, `ar_triage_rate_metrics.csv`,
`ar_triage_failure_analysis.csv`, and `ar_triage_report.md`. Use `--plot` for
an optional EEDF figure when external reference files are present. Use
`--fail-on-code-regression` for CI-style gating of product implementation
regressions, and reserve `--fail-on-physics-mismatch` for stricter manual
reference studies. Use `--run-bolsig` / `--run-mcig` with explicit command
templates when local external binaries should be executed before ingest. See
`docs/benchmarks/ar_bolsig_mcig_triage.md`.

External BOLSIG+ / MCIG output files can also be ingested by adding
`references.external` entries to the benchmark config.  The benchmark then
writes `ar_reference_comparison_summary.csv`,
`ar_reference_eedf_metrics.csv`, `ar_reference_failure_analysis.csv`, and
`ar_reference_report.md`.  See
`docs/benchmarks/ar_external_references.md` for the canonical CSV format.
The limited implementation scope for direct PN higher-order moments is documented in
`docs/dev/direct_pn_lmax_gt1_gate.md`.

The default config keeps Monte Carlo intentionally light. For manual MC
comparison, use:

```powershell
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency_mc_manual.yaml
```

That run is not part of default pytest or CI.
