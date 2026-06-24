# MCIG Benchmark Method Note

MCIG is a development benchmark target, not a product solver mode and not a
ground-truth fitting target. Use it to compare well-defined swarm observables
after the model assumptions are made explicit.

## What To Check

- EEDF convention: compare normalized `F(E)` in `1/eV` with
  `integral F(E) dE = 1`.
- Cross-section input: use the same species, thresholds, interpolation range,
  and high-energy extrapolation policy.
- Angular model: compare only rows whose `angular_model_status` is `match`
  before assigning a solver implementation issue.
- Ionization and growth model: record whether the run is fixed-population,
  weighted-growth, temporal growth, spatial growth, or loss-only.
- Statistics: inspect MC uncertainty, bin counts, effective sample counts,
  and power or energy-balance flags before interpreting tail differences.

## Product Interpretation

The product comparison CSVs intentionally expose only compact interpretation
columns. Benchmark tools may write extra metrics, but they should keep angular
compatibility to the single `angular_model_status` field and put detailed
reasons in failure-analysis evidence.

Internal MC uses a null-collision trial clock. Magnetic motion is substepped
only inside the sampled trial interval, and collision acceptance is evaluated
once at the trial event. Its fixed-population and weighted-growth outputs are
different observables; do not compare them as if they were the same bulk
transport quantity.

## Non-Goals

This benchmark note does not define a new product contract. It does not add
PN coefficient-derived transport, full DCS transport, Coulomb MC, or arbitrary
crossed-field PN support. Those remain future physics work.
