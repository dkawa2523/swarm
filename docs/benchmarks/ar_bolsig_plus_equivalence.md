# Ar BOLSIG+ Equivalence

## Distributed pure-Ar positive-column Swarm input

On 2026-11-08, BOLSIG+ 07/2024 and the distributed two-term Swarm input were
compared at 12 fields from
1 to 5000 Td. Both used the same three ground-state Ar cross sections,
293.15 K, number density from 13.3322 Pa, equal ionization energy sharing,
zero ionization degree for Coulomb collisions, and temporal growth. The
reference used the 1600-point automatic grid and disabled cross-section
extrapolation. [Release downloads](https://www.bolsig.laplace.univ-tlse.fr/download.html)
and [reference manual](https://www.bolsig.laplace.univ-tlse.fr/manual.html).

| Comparison | Maximum difference |
| --- | ---: |
| Mean electron energy | 0.0146% |
| Reduced mobility and scalar particle/energy transport | 0.1025% |
| Major reaction coefficients | 0.0257% |
| Normalized EEDF L1 | 0.8846% |

For the reference run, use the distributed CSV's three processes as BOLSIG+
elastic momentum-transfer, excitation, and ionization data. Set the elastic
electron/Ar mass ratio from 39.948 amu and use thresholds 11.5 and 15.8 eV.
In the BOLSIG+ input use `READCOLLISIONS`, the conditions in the [case data](../dev/results/argon_bolsig_two_term_20261108.json),
one `RUN` per listed E/N, and `SAVERESULTS` with rates and the energy distribution.
Convert the exported EEPF to `F(E) = sqrt(E) * EEPF` and compare the normalized
densities by cell-overlap integration. A major reaction here is at least 1%
of the largest reference coefficient at the same field.

The scalar Swarm F0-gradient diffusion agrees with the ordinary BOLSIG+
diffusion output in this check. A separate BOLSIG+ gradient-expansion run
gives `D_L / D_T` from 0.154 to 0.590 at these fields. It demonstrates the
distinct longitudinal response. The comparison qualifies the matched
two-term calculation and its stated transport definitions.

## Compare a supplied reference export

This benchmark compares product results against external BOLSIG+ output for Ar
single gas.  BOLSIG+ is treated as a two-term weakly ionized swarm reference,
not as a product solver mode.

Run:

```powershell
py -3 tools\benchmarks\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig
```

Use an explicit output file when the config does not point to an existing
reference:

```powershell
py -3 tools\benchmarks\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig --bolsig-output data\references\ar_bolsig_plus_equivalence.csv
```

To automate a local reference export, pass a wrapper command to `--run-bolsig`
with `{input}` and `{output}` for the specified input and output paths. The
wrapper writes the reference in the ingest format below.

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

Pass thresholds for the two-term BOLSIG+ control are mean energy `<1%`, major
rates `<3%`, and EEDF relative L1 `<5%`. Multi-term rows are independent PN
comparisons; their numerical acceptance comes from full-system residual and
lmax-convergence evidence rather than forced equality to the two-term model.
