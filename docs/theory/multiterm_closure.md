# Multi-Term Closure Notes

The runnable ordinary-cross-section multi-term method is
`pn_closure_surrogate`. It uses angular closure assumptions because integral
cross sections do not determine a full differential scattering distribution.

`momentum_power` derives the first Legendre moment from total and
momentum-transfer cross sections and closes higher moments with a power law.
`maxent_p1` uses the same first moment and computes higher moments from a
maximum-entropy P1 distribution.

`pn_dcs` consumes precomputed Legendre moments from a table and therefore does
not use ordinary-XS closure metadata. It is still a table-moment closure path,
not the direct PN block operator.
