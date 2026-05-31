# Ar External References

External BOLSIG+ and MCIG results are benchmark reference inputs, not solver
modes.  Product solver ids remain only `two_term`, `multi_term`, and
`monte_carlo`.

Use schema v2 `references.external` to ingest existing output files:

```yaml
references:
  external:
    - id: bolsig_plus
      path: data/references/ar_bolsig_50Td.csv
      format: electron_swarm_reference_csv
      eedf_convention: eedf
    - id: mcig
      path: data/references/ar_mcig_50Td.csv
      format: electron_swarm_reference_csv
      eedf_convention: eedf
      uncertainty: unavailable
```

The supported canonical CSV format is one row per energy bin with required
columns:

```text
case_id,E_over_N_Td,energy_eV,eedf_eV_inv
```

Use `eepf_eV_m32` instead of `eedf_eV_inv` when the external output is EEPF.
The ingest path converts EEPF with `F(E)=EEPF*sqrt(E)` and normalizes all
curves to EEDF `F(E)` in `1/eV` with `integral F(E)dE = 1`.

Optional scalar columns include `mean_energy_eV`, `drift_velocity_m_s`,
`mobility_m2_V_s`, `diffusion_L_m2_s`, `diffusion_T_m2_s`, and
`net_ionization_frequency_s`.  Optional rate columns use `rate__<name>`.

The Ar benchmark writes:

- `ar_reference_comparison_summary.csv`
- `ar_reference_eedf_metrics.csv`
- `ar_reference_failure_analysis.csv`
- `ar_reference_report.md`

If BOLSIG+ or MCIG binaries are not installed, provide their output files in
the canonical CSV format.  Missing files are reported as external-reference
failures; no fake reference data is generated.

Optional external execution is available through explicit command templates.
The benchmark runs the command, verifies that the declared output file exists,
then ingests that file exactly like a pre-generated reference:

```powershell
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml --run-bolsig '"C:\path\to\bolsig.exe" --input "{input}" --output "{output}"' --bolsig-input data\bolsig.in --bolsig-output data\references\ar_bolsig_plus_equivalence.csv --run-mcig '"C:\path\to\mcig.exe" --input "{input}" --output "{output}"' --mcig-input data\mcig.in --mcig-output data\references\ar_mcig_reference.csv --plot
```

Supported placeholders are `{output}`, `{input}`, `{config}`, and
`{reference_id}`.  The command must create a file in a supported ingest format;
the repository does not fabricate reference curves or assume a vendor-specific
CLI.

BOLSIG+ is treated as a two-term reference.  MCIG is treated as a Monte Carlo
swarm reference.  Angular scattering assumptions must still be checked: a
reference comparison is not a same-angular validation unless the external
metadata proves the angular model is equivalent.
