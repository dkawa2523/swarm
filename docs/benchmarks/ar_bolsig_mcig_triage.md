# Ar BOLSIG+ / MCIG Triage

This benchmark compares product Ar swarm results against external BOLSIG+ and
MCIG references in one matrix.  BOLSIG+ and MCIG are reference sources, not
product solver ids.

Run:

```powershell
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml
```

Provide external files explicitly:

```powershell
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml --bolsig-output data\references\ar_bolsig_plus_equivalence.csv --mcig-output data\references\ar_mcig_reference.csv
```

Or run installed external tools first and then ingest their output:

```powershell
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml --run-bolsig '"C:\path\to\bolsig.exe" --input "{input}" --output "{output}"' --bolsig-input data\bolsig.in --bolsig-output data\references\ar_bolsig_plus_equivalence.csv --run-mcig '"C:\path\to\mcig.exe" --input "{input}" --output "{output}"' --mcig-input data\mcig.in --mcig-output data\references\ar_mcig_reference.csv --plot
```

The external commands are templates. Supported placeholders are `{output}`,
`{input}`, `{config}`, and `{reference_id}`. The commands must produce files in
the configured ingest format; missing binaries or failed commands stop the run
with a clear error.

Use `--require-bolsig` or `--require-mcig` when missing external references
should fail the benchmark.  Otherwise missing files are reported and no fake
reference data is generated.

Useful options:

- `--plot`: write `ar_triage_eedf.png` when external reference files are
  available.
- `--fail-on-code-regression`: fail the command when the triage detects a
  product implementation regression such as the direct lmax=1 gate breaking.
- `--fail-on-physics-mismatch`: fail the command for physics-model mismatch
  rows as well; this is stricter and is intended for manual reference studies.

Outputs:

- `ar_triage_matrix.csv`
- `ar_triage_eedf_metrics.csv`
- `ar_triage_transport_metrics.csv`
- `ar_triage_rate_metrics.csv`
- `ar_triage_failure_analysis.csv`
- `ar_triage_report.md`
- optional `ar_triage_eedf.png` with `--plot`

The metric CSVs use `angular_model_status` for angular compatibility. Detailed
evidence appears in `ar_triage_failure_analysis.csv` and the markdown report.
The triage is a development aid: MCIG is not treated as an error-free reference
when angular metadata or statistical uncertainty are missing.

For an internal-MC-only audit table, run:

```powershell
py -3 tools\benchmark_internal_mc_audit.py --config configs\benchmarks\ar_bolsig_eedf_consistency_mc_manual.yaml
```
