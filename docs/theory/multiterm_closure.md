# Multi-Term Closure Notes

The runnable ordinary-cross-section multi-term method is
`pn_closure_direct`. It uses angular closure assumptions because integral
cross sections do not determine a full
differential scattering distribution. The direct path uses the `lmax: 1`
SG-reduction gate and a limited higher-l angular-closure block solve for B=0,
DC, axisymmetric m=0 cases.

`momentum_power` derives the first Legendre moment from total and
momentum-transfer cross sections and closes higher moments with a power law.
`maxent_p1` uses the same first moment and computes higher moments from a
maximum-entropy P1 distribution.

`pn_dcs` consumes normalized Legendre moments from a table and therefore does
not use ordinary-XS closure metadata. It runs the direct PN block with
`angular_moment_source=moment_table`; `exact_dcs_based=true` is reserved for
tables whose provenance is `dcs_derived`. Raw angle-resolved DCS parsing remains
roadmap work.
