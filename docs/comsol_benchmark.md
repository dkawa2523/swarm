# COMSOL benchmark status

## Evidence status

As of 2026-07-31, the adapter uses a clean class directory for every run,
detects compiler error text even when `comsolcompile` returns exit code zero,
rejects stale classes, reads back the mean-energy formulation, and records
Java/class/MPH/bundle/configuration hashes. A clean external run was attempted,
but COMSOL stopped before model application with license error `-10`
(`Product has expired`). The planned independent three-repeat benchmark is
therefore unavailable.

All numerical comparisons below are a **provisional historical technical
benchmark** based on profiles archived on 2026-07-14. That run's compiler log
contains `Failed to compile java file` and `Compilation failed` despite a
recorded return code of zero. The executed class cannot be proven to correspond
to the archived Java source. Numerical agreement between exported properties
and intended tables is useful for diagnosis, but it is not functional
activation proof.

Direct MPH inspection also corrects the model description: the base model,
archived external-applied model, and built-in Boltzmann reference all use
`LocalEnergyApproximationE`. The external case retained COMSOL's mean-energy
equation. Seven arrays were written, but only six are candidate active in LEA;
the E/N-to-mean-energy table is inactive. An LFA variant would instead have
five candidate-active arrays because electron-energy mobility and diffusion
are inactive when no mean-energy equation is solved.

## Historical model conditions and confounders

- Applied voltage: 200 V.
- COMSOL gas state: 293.15 K and 13.3322 Pa.
- Historical Swarm table state: 300 K and approximately 13.3 Pa.
- Main archived spatial profiles: 200 elements, 801 exported points.
- The built-in saved solver includes both the mean-energy degree of freedom and
  an EEDF degree of freedom. The comparison is not a single-MPH,
  closure-only A/B test.
- The archived 400-element result survives only as aggregate statistics; its
  raw profile CSV was not retained. It cannot support a recomputed mesh
  sensitivity, GCI, or convergence-order claim.

The built-in Boltzmann result is a reference calculation, not experimental
truth.

## Exact piecewise-linear profile comparison

The primary metric integrates the piecewise-linear profile reconstructions
exactly on the common nonuniform grid. “Integral ratio” is
external/reference. For signed quantities it uses signed integrals.

| Quantity | Relative L2 | Spatial correlation | Integral ratio |
|---|---:|---:|---:|
| electron density | 0.530 | 0.995 | 0.484 |
| mean electron energy | 0.662 | 0.778 | 1.457 |
| electric potential | 0.153 | 0.994 | 1.152 |
| E/N | 0.156 | 0.993 | 1.188 |
| electron current density | 0.366 | 0.592 | 0.867 |
| electron + Ar+ conductive current | 0.351 | -0.116 | 0.854 |
| `eir2` direct-excitation source | 0.808 | 0.646 | 0.939 |
| `eir4` direct-ionization source | 0.802 | 0.785 | 0.718 |

The density correlation is high while its integral is only 48.4% of the
reference: the normalized shape is similar, but the amplitude/global particle
balance differs substantially. Potential and E/N are also not independent
validation targets because E/N is derived from the potential gradient.
Potential's integral ratio is gauge-dependent and is retained only as a
descriptive statistic.

The weak conductive-current correlation should not be overinterpreted:
near-constant signed profiles have little variance for a correlation
denominator. Relative L2, mean level, and spatial RSD are more informative for
that quantity.

## Spatial localization and high-field regime

Exact piecewise-linear error attribution over geometric
cathode-side-10% / central-80% / anode-side-10% windows gives:

| Quantity | cathode 10% | central 80% | anode 10% |
|---|---:|---:|---:|
| electron density | 0.139% | 99.591% | 0.271% |
| mean electron energy | 44.082% | 55.187% | 0.731% |
| electric potential | 5.659% | 86.217% | 8.124% |
| E/N | 98.287% | 0.932% | 0.780% |
| electron current | 83.473% | 14.645% | 1.882% |
| conductive current | 84.548% | 13.180% | 2.272% |
| `eir2` direct excitation | 54.084% | 45.743% | 0.172% |
| `eir4` direct ionization | 97.187% | 2.811% | 0.002% |

These are geometric partitions, not a diagnosed bulk/sheath decomposition.
The result localizes much of the E/N, current, and direct-ionization
difference near the cathode-side boundary, while the density-amplitude error
is distributed through the central window.

Using COMSOL's typical 500 Td drift-diffusion guidance as a contextual warning,
not a hard cutoff, E/N above 500 Td occupies 3.178% of external-profile length
and 3.075% of reference-profile length. Yet that region contributes 26.50%
and 22.34% of the respective `eir2` absolute integrals, and 98.97% and 50.92%
of the respective `eir4` absolute integrals. The direct-ionization comparison
is therefore dominated by a short high-field region where locality and
drift-diffusion assumptions require special scrutiny. Providing tables through
5000 Td prevents numerical extrapolation; it does not validate the model form.

## Conductive-current constancy

The archived CSV field named `total_current_density` is
\(J_{\rm e}+J_{\rm Ar^+}\), not terminal total current. It excludes displacement
current and any circuit contribution.

For the 200-element profiles, exact piecewise-linear population RSD is:

| Path | full domain | geometric central 80% |
|---|---:|---:|
| external tables | 37.370% | 0.366% |
| built-in reference | 0.234% | 0.0763% |

RSD quantifies spatial constancy; it is not a discrete continuity residual or
proof of current conservation. The large external full-domain value is
edge-localized, but the archive lacks charge-density, displacement-current,
terminal-current, boundary-flux, and elementwise-residual diagnostics needed
to classify it physically or numerically.

A legacy summary reports point-count-weighted RSD values of 76.96%/0.5137%
for the 200-element external run and 35.37%/0.2036% for a 400-element run
(full domain/geometric central 80%). Because the 400-element raw profile is
missing and the weighting differs from the primary exact spatial metric,
these numbers are historical context only, not a mesh-sensitivity result.

## Reaction and activation audit

The base and external MPH archives contain the same feature and reaction
inventory. The external tables target only:

| Feature | Reaction role | External input |
|---|---|---|
| `eir2` | direct excitation, e + Ar -> e + Ars | Townsend lookup |
| `eir4` | direct ionization, e + Ar -> 2e + Ar+ | Townsend lookup |

Elastic scattering, superelastic de-excitation, stepwise ionization, Penning
ionization, metastable quenching, and surface reactions remain separate.
Consequently, the plotted `eir2`/`eir4` values are direct-channel source terms,
not total excitation or total ionization.

Against the intended bundle, the historical profile's reduced mobility and
`eir2`/`eir4` Townsend values agreed at all 801 points within 1%. The
E/N-to-mean-energy lookup did not: median relative difference was 57.9%, the
95th percentile was 97.6%, the maximum was 959%, and 792/801 points exceeded
1%. That mismatch is expected under LEA because mean energy is solved rather
than imposed by the E/N lookup. It supports the formulation diagnosis but,
because of stale class provenance, cannot prove archived executable
activation. The current LEA verification contract therefore checks four
stable spatial quantities: reduced mobility, longitudinal reduced diffusion,
and `eir2`/`eir4` Townsend coefficients.

## Timing

| Scope | external path | built-in reference | descriptive ratio |
|---|---:|---:|---:|
| run-stage wall time | 42.914 s | 370 s | 8.62 |
| end-to-end wall time | 80.120 s | 379 s | 4.73 |

The 42.914 s observation is a run-stage wall time, not a PDE-kernel or
solve-only time. The external end-to-end archive contains apply 13.780 s,
verify 10.838 s, run-stage 42.914 s, and export 12.588 s. Each path has
`n = 1`; internal timer scopes are not proven identical, and the external
class provenance is defective. The ratios are descriptive historical
observations, not speedup estimates or evidence of solver superiority.

## Conditions for a formal result

A formal rerun must:

1. restore a valid COMSOL license and start from
   `Model/positive_column_1d.mph`;
2. compile into an empty output directory and verify Java/class hashes;
3. read back `LocalEnergyApproximationE`, dependent variables, all seven
   written properties, and the six-candidate-active LEA contract;
4. pass the four LEA spatial numerical-consistency checks and perturbation
   tests for functional activation;
5. match Swarm and COMSOL at 293.15 K and 13.3322 Pa;
6. preserve raw profiles at 200 and 400 elements, adding a finer mesh if
   convergence is claimed;
7. export terminal/displacement/conductive currents, charge density, boundary
   fluxes, and global particle/charge balances;
8. run external and reference paths as independent processes at least three
   times with identical timer contracts;
9. report median and range for apply, verify, run, export, and end-to-end time;
10. retain full logs and Java/class/MPH/bundle/configuration/cross-section
    hashes.

Until those gates pass, the 2026-07-14 values remain provisional and must not
be used to claim experimental accuracy, equation-level equivalence, or a
general performance advantage.
