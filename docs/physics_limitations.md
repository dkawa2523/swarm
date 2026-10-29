# Physics Limitations

This product is a solver comparison platform, not a complete plasma-chemistry
or differential-scattering package.

## Solver Meaning

- `two_term` is a local-field two-term Boltzmann solver. It reports flux
  transport from the solved EEDF.
- Ordinary-cross-section `multi_term` is a direct PN angular-closure solver.
  Integral cross sections do not determine full differential scattering, so it
  must not be read as an exact DCS multi-term solver. Its full PN operator uses
  bidirectional adjacent-moment field coupling. Drift comes from solved `F1`;
  diffusion remains an `F0`-gradient approximation and is labeled as such.
- `pn_dcs` consumes normalized Legendre moment tables. It is marked
  DCS-based only when the table provenance is `dcs_derived`.
- Internal `monte_carlo` is a particle backend. For the fixed-population model,
  particle and energy flux transport are sampled directly from synchronized
  trajectory phase-space moments. It does not report validated bulk transport.
- `propagator` P1 is a deterministic energy--polar-angle finite-volume solver
  for homogeneous DC, positive E/N, and `B=0`. It resolves angular population
  without a two-term truncation, but its answer still depends on the selected
  integral-XS angular closure and velocity-grid discretization. It reports flux
  drift/mobility; bulk drift and all diffusion coefficients are unavailable.

## Angular Scattering

Supported angular models are `isotropic`, `momentum_power`, `maxent_p1`, and
`moment_table`.

`momentum_power` and `maxent_p1` are closure assumptions based on ordinary
cross sections. Internal MC can sample `isotropic` directly. It can sample
`maxent_p1` only when each active scattering species has both explicit total
elastic and momentum-transfer cross sections. `propagator` uses the same
active-species requirement for `maxent_p1`. `momentum_power` and
`moment_table` do not define unique internal-MC samplers in this release.

Implementation support and physical fidelity are independent. Executing an
explicitly selected isotropic or maximum-entropy closure without a fallback is
supported behavior, while its result still carries the
`integral_xs_closure` fidelity label and the selected closure assumption. This
prevents a declared approximation from being confused with an implementation
degradation, without promoting ordinary integral cross sections to DCS data.
The maintained argon application table records this interpretation and its
unrecorded original external citation in a
[provenance sidecar](../examples/cross_sections/argon_application_library.provenance.json);
no DCS claim is inferred from that legacy asset.

## Monte Carlo Estimators

Internal MC uses a null-collision event clock. For `B=0`, EEDF residence is
partitioned at exact energy-bin crossings and piecewise-linear reaction kernels
are integrated after splitting at cross-section knots and thresholds. Magnetic
motion uses endpoint-energy residence with a Boris push between trial events.
Collision acceptance is evaluated at the trial event, not at every orbit
substep.

Attachment particle loss is not yet implemented. If an active attachment
cross section is present, both fixed-particle and weighted-branching requests
fail or skip according to `feature_policy` before evolving trajectories;
attachment is never included only in reported
rates while omitted from the EEDF and transport dynamics.

Fixed-particle ionization tracks one daughter and records fixed-population flux
transport. Common physical-time planes and multiple finite-lag time origins
provide the position-velocity, energy-flux, and position-energy correlations.
Energy mobility and longitudinal/transverse energy diffusion therefore come
from MC trajectories, not an EEDF integral or two-term reconstruction. Raw
coefficients must pass positivity, within-run plateau, and replica-RSE gates
before table qualification.

`weighted_branching` tracks both daughters with weights and synchronized
resampling. Its v6 configured-lag estimator reports direct residence energy mobility
and restricted density-packet energy diffusion without an EEDF or two-term
reconstruction. Qualification requires the minimum effective lineage count at
both `L/2` and the configured production lag `L` to be at least 128 and at least 5%
of the particle population in every replica. Lag `2L` is supplemental evidence
only and is used only when it meets the same threshold. This remains a flux
estimator, not validated bulk transport or a complete two-gradient response.

## Elastic Energy Loss

`two_term` reports `K_epsilon_el` as the energy moment of the same discrete
elastic drift/diffusion operator used in its EEDF solve. This includes the gas-
temperature terms. The elastic operator is isolated by zero-field reassembly
from the same `elastic_A`/`elastic_D` coefficients and SG discretizer; it is not
a cold-gas formula evaluated after the solve.

The Boltzmann solvers drive neutral elastic energy relaxation with the elastic
momentum-transfer integral, independently of the total collision count and
the effective momentum integral used for transport. Inelastic cross sections
never enter the elastic recoil coefficients. A total-only isotropic input
uses `sigma_m = sigma_0`; a total-only PN moment-table input uses
`sigma_m = sigma_0 (1-m1)`. When both integrals and a moment table are supplied,
they must agree within `rtol=1e-6, atol=1e-30 m2`. Effective-only scattering
inputs cannot identify elastic energy relaxation and are rejected.

Isotropic PN requires matching explicit total and momentum-transfer integrals;
unequal integrals require an explicit anisotropic closure such as `maxent_p1`.
Two-term can consume the momentum integral alone, but if a supplied total
reveals a nonzero first angular moment, its assumption metadata records that
difference from an isotropic particle model. Missing elastic recoil inputs
and inconsistent isotropic PN requests are resolved through `feature_policy`
before solver execution.

Internal MC samples each elastic target velocity from the configured Maxwellian,
accepts collisions with the relative-speed rate `g sigma(E_rel)`, and applies
exact two-body center-of-mass kinematics. Held cross-section tails use a
speed-weighted Maxwell proposal with a global analytic envelope; neutral
velocities are not truncated. Its elastic energy-loss coefficient is therefore
the direct signed event-energy change divided by target density and sampled
trajectory residence time, including gas-temperature energy return. Independent
replicas provide the uncertainty exported in
`elastic_energy_loss_vs_mean_energy.csv`. Missing or mixed solver artifacts are
not repaired from another solver.

`propagator` P1 uses a finite-temperature, first-mass-ratio elastic
Fokker--Planck generator. Its Scharfetter--Gummel discretization preserves the
cell-integrated Maxwell equilibrium and detailed balance. Its drift and
diffusion coefficients use `nu_m = N v sigma_m`, not the total event frequency
`nu_0 = N v sigma_0`; this distinction is essential for anisotropic scattering.
The angular loss/gain operator separately uses `nu_0` and realizes
`nu_0 (1-g) = nu_m` through its exact discrete P1 moment. The reported
`K_epsilon_el` is the energy moment of that same generator, so neutral thermal
energy return is part of both the solved distribution and the reported loss;
it is not a cold-gas expression applied after the solve.

## Ionization Energy Sharing

The Boltzmann, particle, and Propagator transfer operators use the same daughter
energy rules. Equal sharing divides the incident energy above threshold in
half. Primary/secondary sharing caps the configured secondary energy at the
available incident energy above threshold, leaving the remainder to the
primary. Loss-only sharing tracks one electron with all remaining energy.
The daughter energy sum is therefore always the incident energy minus the
threshold before projection onto the energy grid. Finite-volume Boltzmann
projection preserves this identity when the daughter energies lie inside the
represented energy centers; daughters below the first center enter the first
cell and retain the corresponding finite-grid energy error.

Two-term temporal growth is evaluated from the particle-number source assembled
with those same transfers: two-daughter ionization creates one net electron,
attachment removes one, and conservative or ignored processes create none.
Conservative field and elastic terms do not enter this update; their large
discrete fluxes must not generate artificial growth through floating-point
cancellation. The full eigenproblem residual remains an independent gate.

## Propagator Qualification

The P1 core uses stationary collision-coupled shell responses, a
characteristic-aligned origin ball, a positive Perron solve, and cell-integrated
collision observables. The core grid is a nested sinh-stretched speed map,
linear at the origin and concentrated in the thermal/sub-eV range. Regular
shells analytically average `1/v`; three-level observed contraction controls
coefficient error, while a separate estimate controls the rational response
error. Inelastic rates and daughter maps use one positive quadrature partitioned
at cross-section knots, thresholds, and ionization-sharing kinks.

The deterministic qualification is Propagator-only. Explicit Magboltz 11.17
GAS2 total and momentum-transfer integrals exercise the `maxent_p1` closure at
`0.05, 0.1, 1 Td`; isotropic inelastic Ar covers `10, 100, 500 Td`.
The `300 x 36` to `600 x 72` comparison gates major scalars at 1% and
cell-density EEDF L1 at 2%. A 90% Ar / 10% O2 case passes 36-to-72-angle refinement,
and independent 7.5/15/30 keV solves establish energy-ceiling independence.

The qualification records `p1_deterministic_core_qualified: true` with no
blocking gates. Its scope excludes every other solver, COMSOL, and P2; it makes
no claim about bulk drift or diffusion. Like the other solver modes,
`propagator` is accepted as a calculation model with explicit closure,
discretization, and observable limits. Universal accuracy for untested gases,
raw-DCS models, or strongly optically thick synthetic origin balls is not a
completion requirement. Those cases require target-specific validation only
when they are placed in scope. The fixed environment, source fingerprint,
provenance hashes, timings, and gates are in the
[Propagator-only qualification artifact](dev/results/propagator_p1_core_qualification.json).

The GEC-CCP range has a separate target gate rather than enlarging the general
solver claim. A split-grid check identified angular, not energy, resolution as
the high-field limiter. The `300 x 48` production grid is checked against
independent `600 x 72` solves with the same 1% scalar and 2% EEDF limits.
Measured differences, the complete 3000--4000 Td decision, and the database
evidence digest are recorded once in the generic-schema
[GEC-CCP target qualification artifact](dev/results/propagator_gec_ccp_target_qualification.json).

## Unsupported Physics

- Stationary high-frequency RF; time-dependent fields outside the two-term
  sinusoidal, instantaneous-F1 approximation
- finite-k hydrodynamic transport
- full Landau electron-electron collisions
- Coulomb particle-particle MC
- raw angle-resolved DCS parsing
- arbitrary crossed-field PN dynamics
- YAML-side state-resolved or generated-superelastic chemistry

Requested unsupported physics is handled by `feature_policy`; it is not silently
ignored.

## Electron-Electron Treatment

`relaxation_postprocess` and `fp_energy` update Boltzmann-family EEDF and rates
after the solver run. Transport is not recomputed and is marked stale.
The rate-derived effective Townsend proxy uses that retained drift and inherits
the transport-stale limitation. It is not a newly measured spatial-growth result.
Solver residuals and elastic-loss diagnostics from the original distribution
are discarded; stale transport cannot qualify for downstream tables.
Time-dependent e-e postprocessing is unsupported because phase-resolved
distributions would also require an update.

`relaxation_postprocess` distinguishes a conserved bounded-grid mean from a
specified temperature. A requested mean outside the positive-temperature
Maxwell family's attainable range raises an error. With mean conservation
disabled, the target uses the specified temperature directly; its finite-domain
mean need not equal the infinite-domain value of 1.5 times that temperature.

`fp_energy` uses a positive conservative finite-volume step. When requested,
the discrete inverse temperature is solved so the mean energy is conserved
along with probability. A domain that cannot represent the required positive
temperature fails rather than clipping or correcting the distribution afterward.
PN may contain tiny signed roundoff within its native acceptance bound. Before
e-e relaxation, negative probability and energy fractions below `1e-8` may be
projected onto nonnegative input while preserving probability and mean energy;
the input L1 change is bounded by `1e-7` and recorded when applied. Larger or
unrepresentable corrections fail. This changes only the relaxation input and
does not repair or relabel the native PN solution.

Angular closures are evaluated independently for each gas species before
combining collision operators. This preserves mixture linearity and the same
closure assumptions across PN, MC, and Propagator. Ordinary integral cross
sections still do not identify the true DCS.

Direct MC histogram `weight_effective_count` measures weight concentration.
It is not an independent sample count or a standard-error estimate: residence
samples are correlated and branching shares ancestry. Table uncertainty is
estimated across independently seeded workflow replicas.
