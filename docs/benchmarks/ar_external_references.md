# Ar External References

External BOLSIG+ and MCIG results are benchmark reference inputs, not solver
modes. Product solver ids are `two_term`, `multi_term`, `monte_carlo`, and
`propagator`.

Benchmark configs may include a benchmark-only `references.external` block to
ingest existing output files. This block is consumed by the benchmark tools and
is not part of product YAML accepted by `electron_swarm.load_config`:

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
The configured `eedf_convention` must agree with the distribution column;
contradictions fail. Missing identifiers/E/N, negative distributions, invalid
numbers, and conflicting repeated scalar values fail rather than dropping rows.
The input normalization and whether mean energy was reported or derived remain
in reference metadata.

Supply `energy_width_eV` when cell geometry is known. Otherwise midpoint cells
are inferred explicitly and that assumption is recorded. Comparisons integrate
the union of the reported supports; they do not truncate and renormalize one
curve onto the other's grid. Missing transport coefficients remain unavailable;
reference data is not converted into a solver result with zero transport.

Optional scalar columns include `mean_energy_eV`, `drift_velocity_m_s`,
`mobility_m2_V_s`, `diffusion_L_m2_s`, `diffusion_T_m2_s`, and
`net_ionization_frequency_s`.  Optional rate columns use `rate__<name>`.
`rate__ionization`, `rate__excitation`, and `rate__attachment` mean the sums of
mixture-weighted coefficients of that process type. Other names must match one
candidate process exactly. Unmatched or ambiguous rates never become zero.
Named reaction errors are calculated independently; the largest rate error is
not copied into every reaction column.
The major-rate gate uses reference rates at least 1% of the largest reported
reference rate and requires every selected rate to be comparable. Source file
hashes, geometry, input normalization, and the mean-energy source are retained
in summary rows. Declared unavailable uncertainty is not used for acceptance.

The Ar benchmark writes:

- `ar_reference_comparison_summary.csv`
- `ar_reference_eedf_metrics.csv`
- `ar_reference_failure_analysis.csv`
- `ar_reference_report.md`

Summary and metric CSVs use the compact `angular_model_status` field with
`match`, `mismatch`, or `unknown`. Detailed angular reasons belong in
failure-analysis evidence, not in separate summary columns.
`PARTIAL` and `unavailable_metrics` expose missing threshold evidence. A reported
MC scalar confidence interval may qualify the mean-energy difference, but it
cannot override failed EEDF or reaction-rate gates. Unknown angular assumptions
remain degraded. These decisions describe the listed numerical comparisons,
not universal validation of the solver or collision data.

If BOLSIG+ or MCIG binaries are not installed, provide their output files in
the canonical CSV format.  Missing files are reported as external-reference
failures; no fake reference data is generated.

Optional external execution is available through explicit command templates.
The benchmark verifies the declared output file before ingesting it:

```powershell
py -3 tools\benchmarks\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml --run-bolsig '"C:\path\to\bolsig.exe" --input "{input}" --output "{output}"' --bolsig-input data\bolsig.in --bolsig-output data\references\ar_bolsig_plus_equivalence.csv --run-mcig '"C:\path\to\mcig.exe" --input "{input}" --output "{output}"' --mcig-input data\mcig.in --mcig-output data\references\ar_mcig_reference.csv --plot
```

Supported placeholders are `{output}`, `{input}`, `{config}`, and
`{reference_id}`. The repository never fabricates reference curves. BOLSIG+ is
treated as a two-term reference, and MCIG as a Monte Carlo swarm reference.
