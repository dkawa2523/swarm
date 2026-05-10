# Swarm Solver Project Notes

This repository is a compact basis for comparing electron swarm calculations
from the same YAML conditions and cross-section data.

The maintained entry point is `electron_swarm`. It can dispatch:

- `monte_carlo`: particle Monte Carlo through `swarm_mc`.
- `boltzmann_two_term`: native BOLSIG-like two-term solver.
- `multiterm_boltzmann`: moment-closure estimate by default, with lmax=1 as a
  two-term reference adapter and lmax>1 as a reference-anchored closure for
  integral-cross-section validation work.

Current development priorities:

- Keep solver dispatch and result writing simple enough for third-party review.
- Keep optional collision extensions outside individual solver internals.
- Mark approximate or experimental physics clearly in result metadata.
- Use focused tests for shared YAML behavior, result metadata, and regression
  safety instead of broad contracts around every internal helper.

For user-facing commands and schema notes, see `README.md`.
