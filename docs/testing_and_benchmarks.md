# Testing and Benchmarks

The default developer path stays lightweight. Heavier Monte Carlo and manual
external-reference checks are opt-in.

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

Benchmarks are development aids for solver regressions, physics triage, and
manual comparison against external references. BOLSIG+ and MCIG files or local
binaries are not product solver modes. Reference parsing, external command
execution, and EEDF comparison helpers live under `tools`, outside the runtime
`electron_swarm` package.

Run product benchmark configs:

```powershell
py -3 benchmarks\run_product_benchmarks.py --config configs\benchmarks\ar_simple.yaml
py -3 benchmarks\run_product_benchmarks.py
```

The benchmark runner writes only `outputs/benchmarks/benchmark_summary.csv`.
Canonical solver CSV outputs are not written during benchmarks.

Focused benchmark CLIs:

```powershell
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency.yaml
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml
```

The focused tools write CSV/report files under `outputs/benchmarks`. See the
matching page under `docs/benchmarks/` for command-specific inputs and outputs.

External BOLSIG+ / MCIG output files can also be ingested by benchmark-only
`references.external` entries. Product YAML loaded by
`electron_swarm.load_config` rejects those entries as unknown top-level fields.
See `docs/benchmarks/ar_external_references.md` for the canonical CSV format.

The default config keeps Monte Carlo intentionally light. For manual MC
comparison, use:

```powershell
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency_mc_manual.yaml
```

That run is not part of default pytest or CI.

For internal MC audit-only runs, `tools\benchmark_internal_mc_audit.py` writes
the raw histogram bin width, probability mass, cumulative probability, survival
probability, and sample-quality columns. It does not smooth or fit the EEDF to
external references.

For internal-MC physics triage, use
`tools\benchmark_internal_mc_physics_triage.py`. It compares raw bin mass,
quantiles, tail survival, and model-definition metadata without declaring any
external result to be the truth.
