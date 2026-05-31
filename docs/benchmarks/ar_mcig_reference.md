# Ar MCIG Reference Benchmark

This benchmark ingests external MCIG output as a Monte Carlo swarm reference for
Ar single gas.  MCIG is a reference source, not a product solver id.

Run without requiring an external file:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig
```

Provide a generated MCIG file explicitly:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig --mcig-output data\references\ar_mcig_reference.csv --angular-model isotropic
```

Run a local MCIG binary first by passing an explicit command template:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig --run-mcig '"C:\path\to\mcig.exe" --input "{input}" --output "{output}"' --mcig-input data\mcig.in --mcig-output data\references\ar_mcig_reference.csv --angular-model isotropic --plot
```

The command must create the declared output file in a supported ingest format.

Use `--require-mcig` to fail when the MCIG file is missing.  Without it, missing
output is recorded as a reference-ingest failure and no fake MCIG data is
generated.

The supported canonical CSV columns are:

```text
case_id,E_over_N_Td,energy_eV,eedf_eV_inv,mean_energy_eV,drift_velocity_m_s
```

Use `eepf_eV_m32` instead of `eedf_eV_inv` when the MCIG export is EEPF.  The
benchmark converts all curves to normalized EEDF `F(E)` in `1/eV`.

Optional uncertainty columns use a `_ci95` suffix, for example
`mean_energy_eV_ci95` or `drift_velocity_m_s_ci95`.  If no uncertainty is
available, comparison confidence is marked `unknown`; differences are not
automatically treated as implementation bugs.

Angular scattering matters.  Set `angular_model` in `references.external`, or
pass `--angular-model isotropic|mcig_default|unknown`.  Rows with unknown or
mismatched angular metadata are marked degraded.

Outputs:

- `ar_mcig_reference_summary.csv`
- `ar_mcig_reference_eedf_metrics.csv`
- `ar_mcig_reference_failure_analysis.csv`
- `ar_mcig_reference_report.md`
- optional `ar_mcig_reference_eedf.png`

MCIG vs BOLSIG+ rows are included when both reference files are available.  Use
those rows to separate two-term approximation differences from product code
mismatches.
