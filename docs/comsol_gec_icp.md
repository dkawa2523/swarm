# Argon GEC-ICP external-Swarm closure

The GEC-ICP adapter applies one independently calculated `two_term`,
`monte_carlo`, or `propagator` closure to
`comsol_modes/argon_gec_icp.mph`. These are three canonical solver modes, not
variants of one composite solver. COMSOL performs one cold-start
frequency-transient solve; it does not sweep E/N, pressure, coefficients, or
coil power.

## Inputs and ownership

| Source | Current workflow | COMSOL map |
| --- | --- | --- |
| Two term | `examples/workflow_argon_gec_icp_two_term.yaml` | `comsol_modes/maps/argon_gec_icp_two_term_function_eedf.yaml` |
| Monte Carlo | `examples/workflow_argon_gec_icp_monte_carlo.yaml` | `comsol_modes/maps/argon_gec_icp_monte_carlo_function_eedf_restricted_lmea.yaml` |
| Propagator P1 | `examples/workflow_argon_gec_icp_propagator.yaml` | `comsol_modes/maps/argon_gec_icp_propagator_function_eedf.yaml` |

The shared model topology is declared in
`comsol_modes/maps/argon_gec_icp_model.yaml`. Every map uses 1500 W,
13.56 MHz, 300 K, 2.66644 Pa, a 0.01 s final time, and four logarithmic output
intervals per decade.

Each external solver owns its steady-DC EEDF and reduced mobility. COMSOL
retains Einstein particle diffusion, electron-energy transport, all
cross-section rate and elastic-loss integrals, argon/metastable chemistry,
walls, heavy-species transport, the inductive RF field, conductivity coupling,
and the electron heat source. This is a restricted local-mean-energy closure,
not a nonlocal RF kinetic solve. No coefficient is repaired after either the
Swarm or COMSOL calculation.

The model-independent bundle records `source_policy.qualification_outputs`:
solver quantities covered by its qualification evidence. The GEC-ICP map and
plan separately record the subset actually bound here. Thus MC's direct elastic
loss remains qualified solver evidence, while this ICP calculation deliberately
uses COMSOL's cross-section integral for elastic energy exchange.

The 13 comparison anchors span 1–2500 Td. Solver-native support endpoints are:

- two term: 3000 Td;
- Monte Carlo: 2726.6505889 Td;
- propagator: 2744.2936035 Td.

They are real solver cases used to cover the transient mean-energy guard, not
synthetic extrapolation points.

## P1 improvements

The Function-EEDF input now starts from a continuous moment-conserving C1
family. Its physical electron-energy and mean-energy axes are both refined
adaptively. Every solver anchor is retained, and every materialized row is
nonnegative, normalized, and has the requested first moment.

Refinement stops only when the actual tensor-product table meets 0.5% bounds
for both the sqrt(E)-weighted total-variation error and every active
cross-section-rate kernel. Rates above 0.1% of a kernel's family peak use a
relative bound; smaller rates use a peak-scaled absolute bound. This avoids
forcing dense rows in the vanishing-rate limit without altering any rate.
Out-of-support values still require a new core-solver anchor; clamping and
synthetic tail extension are not used.

Historical pre-refactor measurements used these adaptive tables:

| Source | Solver anchors | Mean nodes | Energy nodes | Rows | Max mean-axis TV/rate error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Two term | 14 | 93 | 605 | 56,265 | 0.495% / 0.475% |
| Monte Carlo | 14 | 114 | 797 | 90,858 | 0.497% / 0.493% |
| Propagator | 14 | 89 | 543 | 48,327 | 0.455% / 0.499% |

The former P1 input used only the 14 solver anchors on the mean-energy axis.
Although its anchor rows were valid, COMSOL's piecewise-linear interpolation
between them produced rate discontinuity error large enough to make the
transient solve extremely slow. The adaptive 89-row P1 table removes that
core representation error; that historical COMSOL solve stage completed in about
60 s.

Monte Carlo was also corrected in the solver itself: reaction-kernel weighted
ensembles target rare rates, residence/rate integration uses quadrature, the
tail cadence is fixed, and the 1 Td mobility uses common-random-number
opposite-field parity. The 5 Td trajectory window was extended. The canonical
workflow computes all 14 anchors directly with one bounded sampling plan and
does not import cases from an earlier attempt database. All final MC anchors
must pass the declared restricted-LMEA gate, so the active ICP result remains
pure MC without a two-term fallback.

The canonical `monte_carlo_function_eedf.sqlite` has not yet been generated.
The older attempt and supported databases are retained only as historical
numerical evidence; their workflow provenance does not qualify the canonical
input. A new MC/COMSOL release therefore requires running the sequence below.

Propagator no longer imports the two-term implementation. Shared gas-density
and rate-convolution primitives live in `electron_swarm/physics/kinetics.py`;
the MC, two-term, and propagator numerical packages remain independently
runnable.

## Reproduction

Two term:

```powershell
py -3 -m swarm_workflow.cli sweep examples\workflow_argon_gec_icp_two_term.yaml
py -3 -m swarm_workflow.cli build-tables outputs\argon_gec_icp\two_term_function_eedf_extended_support.sqlite --output outputs\argon_gec_icp\two_term_tables_tail_converged --source two_term
py -3 -m swarm_workflow.cli export-comsol outputs\argon_gec_icp\two_term_tables_tail_converged --output outputs\argon_gec_icp\two_term_bundle
py -3 -m swarm_workflow.cli run-gec-icp comsol_modes\maps\argon_gec_icp_two_term_function_eedf.yaml --bundle outputs\argon_gec_icp\two_term_bundle\mixture_0000
```

The MC workflow computes all 14 anchors from an empty canonical database:

```powershell
py -3 -m swarm_workflow.cli sweep examples\workflow_argon_gec_icp_monte_carlo.yaml
py -3 -m swarm_workflow.cli build-tables outputs\argon_gec_icp\monte_carlo_function_eedf.sqlite --output outputs\argon_gec_icp\monte_carlo_tables --source monte_carlo --mc-qualification-profile function_eedf_restricted_lmea
py -3 -m swarm_workflow.cli export-comsol outputs\argon_gec_icp\monte_carlo_tables --output outputs\argon_gec_icp\monte_carlo_bundle
py -3 -m swarm_workflow.cli run-gec-icp comsol_modes\maps\argon_gec_icp_monte_carlo_function_eedf_restricted_lmea.yaml --bundle outputs\argon_gec_icp\monte_carlo_bundle\mixture_0000
```

Propagator qualification is bound to the current implementation before target
refinement and export:

```powershell
py -3 tools\qualify_propagator_p1_deterministic.py --output docs\dev\results\propagator_p1_deterministic_qualification_20260908.json
py -3 -m swarm_workflow.cli sweep examples\workflow_argon_gec_icp_propagator.yaml
py -3 tools\qualify_propagator_target.py --target argon_gec_icp_restricted_local_mean_energy --fields-Td 1500 2500 2744.2936035 --operating-bracket-Td 1500 2500 --table-support-cap-Td 2744.2936035 --medium-grid 300 48 --fine-grid 600 72 --database outputs\argon_gec_icp\propagator_function_eedf_mean_energy_guarded.sqlite --config examples\argon_gec_icp_propagator.yaml --core-qualification docs\dev\results\propagator_p1_deterministic_qualification_20260908.json --workers 2 --memory-budget-mb 2048 --output docs\dev\results\propagator_gec_icp_mean_energy_guarded_qualification_20260913.json
py -3 -m swarm_workflow.cli build-tables outputs\argon_gec_icp\propagator_function_eedf_mean_energy_guarded.sqlite --output outputs\argon_gec_icp\propagator_tables --source propagator --solver-qualification docs\dev\results\propagator_p1_deterministic_qualification_20260908.json --target-qualification docs\dev\results\propagator_gec_icp_mean_energy_guarded_qualification_20260913.json
py -3 -m swarm_workflow.cli export-comsol outputs\argon_gec_icp\propagator_tables --output outputs\argon_gec_icp\propagator_bundle
py -3 -m swarm_workflow.cli run-gec-icp comsol_modes\maps\argon_gec_icp_propagator_function_eedf.yaml --bundle outputs\argon_gec_icp\propagator_bundle\mixture_0000
```

Use `--dry-run` on `run-gec-icp` to validate topology, hashes, physical
context, support, and qualifications without starting a licensed solve.

## Acceptance and historical measurements

Apply plus saved-closure readback and native Function-EEDF evaluation run in
one COMSOL process; the physical solve is a second process. Acceptance requires
all 26 planned times through 0.01 s, exact saved-input provenance, native
Function-EEDF moments/rates, the 10% transient mean-energy guard, finite
positive states, passive absorbed power, and the final-four-point stationarity
limits. A successful COMSOL exit alone is insufficient.

The last completed pre-refactor runs passed their then-current acceptance
checks. Their numerical values are recorded here only as historical context:

| Source | Closure/audit | Solve | Electron inventory | Mean energy | Ar* inventory | Absorbed power |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Two term | 33.0 s | 58.3 s | 1.90171e14 | 7.51879 eV | 3.16378e15 | 1392.653 W |
| Monte Carlo | 36.9 s | 68.9 s | 1.91435e14 | 7.48730 eV | 3.13100e15 | 1393.035 W |
| Propagator | 32.1 s | 60.3 s | 1.91905e14 | 7.47528 eV | 3.13262e15 | 1393.180 W |

These numbers are not current release evidence until the canonical workflows,
qualification artifacts, bundles, and COMSOL runs are regenerated. In the
historical comparison, the pure-MC EEDF differs from two term by at most 3.81%
total variation; P1 differs from two term by at most 3.94%, both at 2500 Td.
MC and P1 differ by at most 1.86%, at 2 Td.

The like-for-like model comparison uses the exact common 1 ms state. Relative
to the original Maxwellian/COMSOL-mobility closure, two term, MC, and P1 change
electron inventory by -35.39%, -34.96%, and -34.80%, and mean energy by
+46.73%, +46.12%, and +45.88%, respectively. These are whole-closure
differences because EEDF and mobility change together. The original saved
1 ms state fails the current Ar* stationarity limit (11.543% versus 3%), so it
is not a converged 10 ms reference.

Current inspectable Propagator input evidence:

- `docs/dev/results/propagator_p1_deterministic_qualification_20260908.json`
- `docs/dev/results/propagator_gec_icp_mean_energy_guarded_qualification_20260913.json`
- `outputs/argon_gec_icp/propagator_tables/mixture_0000/manifest.json`
- `outputs/argon_gec_icp/propagator_bundle/mixture_0000/manifest.json`
- `outputs/argon_gec_icp/comsol_propagator_function_eedf/gec_icp_plan.json`

The obsolete pre-refactor Propagator comparison and COMSOL-result directories
were removed. A licensed run of the final command above must create a fresh
`run_status.json` before COMSOL convergence or physical-target acceptance is
claimed; `--dry-run` proves only input and model preflight readiness.
