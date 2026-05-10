# Production Boltzmann two-term backend design

This overlay adds a production-oriented electron Boltzmann two-term backend to the
existing particle Monte Carlo swarm workflow.  The goal is not to replace BOLSIG+
as a reference code, but to provide a BOLSIG+/Hagelaar-Pitchford-class native
solver that shares input, output, plotting, and regression infrastructure with
the particle Monte Carlo implementation.

## Solver choices

`boltzmann_two_term.backend` accepts the following values:

- `native_bolsig`: built-in production backend.  This is the default.
- `bolos`: optional independent reference backend using the third-party BOLOS
  package when installed.
- `auto`: use BOLOS when importable, otherwise use `native_bolsig`.
- `internal`: backward-compatible alias for `native_bolsig`.

## Governing approximation

The electron distribution is expanded in the direction cosine relative to the
uniform electric field,

```text
f(v, mu) = f0(v) + mu f1(v),     mu = cos(theta),
```

and all outputs are derived from the normalized energy distribution function
`F(epsilon)` satisfying

```text
integral_0^inf F(epsilon) d epsilon = 1.
```

The native backend solves the steady temporal-growth eigenproblem

```text
L[F] = lambda F
```

where `L` is the energy-space collision/heating operator and `lambda` is the net
population growth frequency caused by ionization and attachment.  Conservative
runs have `lambda = 0` up to residual tolerance.

## Energy-space finite-volume operator

The energy flux is written in fitted convection-diffusion form,

```text
J = A(epsilon) F - D(epsilon) dF/depsilon,
```

with zero numerical flux at the energy-domain boundaries.  The electric-field
heating term contains both diffusion and the geometric drift term required when
solving for the EEDF rather than for `f0` directly:

```text
D_field = (2/3) (e E)^2 / m_e * epsilon_J / nu_m / e^2
A_field = D_field / (2 epsilon_eV)
```

Elastic energy exchange with a gas at temperature `T_g` is represented as a
Fokker-Planck contribution that preserves the Maxwell energy PDF when `E = 0`:

```text
D_elastic = 2 (m_e/M) nu_m epsilon kT_g
A_elastic = 2 (m_e/M) nu_m (0.5 kT_g - epsilon)
```

The total nearest-neighbour flux is discretized using the Scharfetter-Gummel
exponential scheme.  This avoids the oscillations and negative EEDF tails common
with central differencing when the local Peclet number is large.

## Inelastic and nonconservative collisions

Excitation, superelastic/de-excitation, ionization, and attachment processes are
assembled as sparse source/sink operators.  Excitation shifts electrons from
`epsilon` to `epsilon - threshold`; superelastic collisions use a negative shift;
attachment is a pure sink.  Ionization supports three models:

- `equal`: two outgoing electrons share the excess energy equally.
- `primary_secondary`: one cold secondary at `secondary_electron_energy_eV` and
  one primary carrying the remaining energy.
- `loss_only`: threshold energy loss without net multiplication, useful for
  sensitivity checks.

## Transport coefficients

The native backend computes flux transport coefficients from the two-term
anisotropic correction, not from a single averaged collision frequency.  With
`nu_m/N = sum_i x_i sigma_m,i v`,

```text
mu N = -e/(3 m_e) integral [2 epsilon/(nu_m/N) dF/depsilon - F/(nu_m/N)] d epsilon
D N  =  1/3       integral [v^2/(nu_m/N) F] d epsilon
W    = (mu N) (E/N)
```

The scalar diffusion output is used for both `diffusion_L` and `diffusion_T` in
the shared schema.  It corresponds to the standard two-term flux diffusion
coefficient.  If density-gradient/bulk longitudinal diffusion is required, it
should be added as a separate higher-order transport module rather than hidden
inside this local-field coefficient.

## Adaptive energy grid

The native solver can regrid automatically.  After each solution it checks

- the integral probability in the top `tail_cells_fraction` of the grid,
- the final-cell EEDF value relative to the peak,
- `mean_energy_multiplier * mean_energy`.

If the tail is not negligible, the upper energy bound is increased and the EEDF
is interpolated onto a new grid.  This mirrors the common BOLSIG/BOLOS workflow
of using a domain several tens of mean energies wide while keeping strong
low-energy resolution.

## Validation workflow

1. Run `configs/boltzmann_only.yaml` and check that every row has
   `meta_converged = True`, small `meta_residual_L1`, and negligible
   `meta_tail_probability`.
2. When BOLOS is installed, run
   `tools/benchmark_operator_gate.py --quick --require-bolos` and compare mean
   energy, reduced mobility, reduced diffusion, and rate coefficients.
3. Compare MC and Boltzmann results with `run.mode: both`; differences outside
   expected two-term limitations indicate anisotropy, nonlocal effects,
   insufficient MC statistics, or cross-section/model inconsistencies.

## Known physical limitations

The two-term approximation is accurate for many weakly ionized, spatially
uniform swarm/fluid-table calculations.  Accuracy degrades for strongly
anisotropic EEDFs, strong magnetic fields, strong electron-electron collisions,
spatially nonlocal sheaths, RF phase-resolved kinetics, and cases requiring
higher Legendre terms or a full multi-term Boltzmann treatment.  Those cases
should be routed to MC or a future multi-term module while retaining the same
shared I/O schema.
