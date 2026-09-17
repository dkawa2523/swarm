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
py -3 tools\benchmarks\run_product_benchmarks.py --config configs\benchmarks\ar_simple.yaml
```

Pass `--config` more than once to run several product benchmark files. Focused
external-reference YAML contains additional tool-owned sections and is run only
through the matching focused command below.

The benchmark runner writes only `outputs/benchmarks/benchmark_summary.csv`.
Canonical solver CSV outputs are not written during benchmarks.

Focused benchmark CLIs:

```powershell
py -3 tools\benchmarks\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency.yaml
py -3 tools\benchmarks\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig
py -3 tools\benchmarks\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig
py -3 tools\benchmarks\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml
```

The focused tools write CSV/report files under `outputs/benchmarks`. See the
matching page under `docs/benchmarks/` for command-specific inputs and outputs.

The P1 propagator has a bounded, reproducible **Propagator-only** deterministic
qualification runner:

```powershell
py -3 tools\qualify_propagator_p1_deterministic.py
py -3 tools\qualify_propagator_p1_deterministic.py --quick --output outputs\benchmarks\propagator_p1_deterministic_quick.json
```

It records the code revision and dirty state, Python/NumPy/SciPy and hardware,
cross-section SHA-256, per-case grid/work/memory/timing diagnostics, maxent-P1
and inelastic medium--fine differences, mixed-gas angular refinement,
energy-ceiling independence, and weak-transfer conservation in
`docs/dev/results/propagator_p1_deterministic_qualification_20260908.json`.
It never instantiates another solver and explicitly excludes COMSOL and P2.
The current artifact records `p1_deterministic_core_qualified: true` and no
blocking gates. `--quick` is a smoke run and is not a replacement for the full
matrix.

External BOLSIG+ / MCIG output files can also be ingested by benchmark-only
`references.external` entries. Product YAML loaded by
`electron_swarm.load_config` rejects those entries as unknown top-level fields.
See `docs/benchmarks/ar_external_references.md` for the canonical CSV format.

The default config keeps Monte Carlo intentionally light. For manual MC
comparison, use:

```powershell
py -3 tools\benchmarks\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency_mc_manual.yaml
```

That run is not part of default pytest or CI.

For internal MC audit-only runs, `tools\benchmarks\benchmark_internal_mc_audit.py` writes
the raw histogram bin width, probability mass, cumulative probability, survival
probability, and sample-quality columns. It does not smooth or fit the EEDF to
external references.

For internal-MC physics triage, use
`tools\benchmarks\benchmark_internal_mc_physics_triage.py`. It compares raw bin mass,
quantiles, tail survival, and model-definition metadata without declaring any
external result to be the truth.
