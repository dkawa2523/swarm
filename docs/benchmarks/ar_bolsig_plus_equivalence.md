# Ar BOLSIG+ Equivalence

This benchmark compares product results against external BOLSIG+ output for Ar
single gas.  BOLSIG+ is treated as a two-term weakly ionized swarm reference,
not as a product solver mode.

Run:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig
```

Use an explicit output file when the config does not point to an existing
reference:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig --bolsig-output data\references\ar_bolsig_plus_equivalence.csv
```

Run a local BOLSIG+ binary first by passing an explicit command template:

```powershell
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig --run-bolsig '"C:\path\to\bolsig.exe" --input "{input}" --output "{output}"' --bolsig-input data\bolsig.in --bolsig-output data\references\ar_bolsig_plus_equivalence.csv --plot
```

The command must create the declared output file in a supported ingest format.

Add `--require-bolsig` to fail when the external file is missing.  Without it,
missing BOLSIG+ output is reported in the failure CSV and no fake reference is
generated.

The external reference file should use `electron_swarm_reference_csv`:

```text
case_id,E_over_N_Td,energy_eV,eedf_eV_inv,mean_energy_eV,drift_velocity_m_s
```

Use `eepf_eV_m32` instead of `eedf_eV_inv` when the BOLSIG+ export is EEPF.
The benchmark converts all curves to normalized EEDF `F(E)` in `1/eV`.

Outputs:

- `ar_bolsig_plus_equivalence_summary.csv`
- `ar_bolsig_plus_equivalence_eedf_metrics.csv`
- `ar_bolsig_plus_equivalence_failure_analysis.csv`
- `ar_bolsig_plus_equivalence_report.md`

Pass thresholds are mean energy `<1%`, major rates `<3%`, and EEDF relative L1
`<5%` against BOLSIG+.  The existing strict `multi_term pn_closure_direct
lmax=1` versus `two_term` gate is also checked.
