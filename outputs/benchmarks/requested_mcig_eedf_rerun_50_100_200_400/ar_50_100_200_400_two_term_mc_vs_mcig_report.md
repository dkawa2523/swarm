# Product two_term/internal MC vs MCIG GUI EEDF sweep

Fresh product run with schema_version=2.

## Product run

- Config: `outputs\benchmarks\requested_mcig_eedf_rerun_50_100_200_400\ar_50_100_200_400_two_term_mc_vs_mcig.yaml`
- Product solvers: `two_term`, `monte_carlo`
- E/N: 50, 100, 200, 400 Td
- Ar, 300 K, pressure 100 Pa, no magnetic field, dc field, no e-e collisions
- Cross sections: `examples/cross_sections/argon_minimal.csv`
- High-energy handling: `high_energy_extrapolation: hold`, `max_eV_limit: 1000.0`
- Internal MC: 1024 particles, 300 warmup collisions, 3000 production collisions, seed 20260624

## MCIG availability

- 100 Td: actual MCIG GUI EEDF archived in `mcig_reference_inputs/`
- 200/400 Td: actual MCIG GUI EEDF archived in `mcig_reference_inputs/`
- 50 Td: no saved MCIG GUI EEDF found; product curves are plotted without MCIG comparison metrics

## Metrics vs MCIG GUI actual

| E/N Td | comparison | rel L1 | candidate mean eV | MCIG mean eV | mean rel diff | tail prob diff >15.76 eV |
|---:|---|---:|---:|---:|---:|---:|
| 100 | product two_term vs MCIG GUI actual | 0.0172106 | 6.49665 | 6.42756 | 0.0107499 | 0.00103416 |
| 100 | product internal MC vs MCIG GUI actual | 0.0266762 | 6.46462 | 6.42756 | 0.00576681 | -0.000245281 |
| 200 | product two_term vs MCIG GUI actual | 0.0421962 | 7.95411 | 7.75772 | 0.0253149 | 0.0066373 |
| 200 | product internal MC vs MCIG GUI actual | 0.0457581 | 7.98326 | 7.75772 | 0.0290734 | 0.00274012 |
| 400 | product two_term vs MCIG GUI actual | 0.0826695 | 10.5695 | 10.1091 | 0.0455392 | 0.0156706 |
| 400 | product internal MC vs MCIG GUI actual | 0.069099 | 10.6842 | 10.1091 | 0.0568911 | 0.00879689 |

## Evaluation

The product curves track the MCIG bulk shape well at 100 Td and degrade as E/N increases. The increase in L1 and mean-energy difference at 200/400 Td is consistent with the known non-equivalent physics assumptions rather than a clean solver-only numerical discrepancy.

The most important caveats are unchanged across the available MCIG cases: MCIG angular scattering model `11.000` is not mapped to the product `isotropic` sampler/closure, MCIG energy sharing model/parameter is not equivalent to product `energy_sharing: equal`, and the product run uses hold extrapolation above the 100 eV cross-section table to keep the 1000 eV EEDF range available.

50 Td should not be included in MCIG error conclusions until an actual MCIG GUI EEDF export for 50 Td is produced with the same collision file/settings.

## Outputs

- Linear overview: `outputs\benchmarks\requested_mcig_eedf_rerun_50_100_200_400\ar_50_100_200_400_two_term_mc_vs_mcig_eedf_overview_linear.png`
- Tail overview: `outputs\benchmarks\requested_mcig_eedf_rerun_50_100_200_400\ar_50_100_200_400_two_term_mc_vs_mcig_eedf_overview_tail.png`
- Curves CSV: `outputs\benchmarks\requested_mcig_eedf_rerun_50_100_200_400\ar_50_100_200_400_two_term_mc_vs_mcig_curves.csv`
- Metrics CSV: `outputs\benchmarks\requested_mcig_eedf_rerun_50_100_200_400\ar_50_100_200_400_two_term_mc_vs_mcig_metrics.csv`
- Settings CSV: `outputs\benchmarks\requested_mcig_eedf_rerun_50_100_200_400\ar_50_100_200_400_two_term_mc_vs_mcig_mcig_settings.csv`
