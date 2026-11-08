# Function EEDF projection benchmark

The finite-volume solver owns the tabulated reaction coefficients. The
exporter owns the continuous representation passed to COMSOL. The adapter
checks the imported Function and the plasma solution with separate decisions.
A matched experiment would test the model and its closures together.

## Reproducible numerical case

The [fixture](../../tests/fixtures/references/argon_p1_projection/) retains
1188 finite-volume EEDF bins and nine reaction coefficients at 5, 20, and
3500 Td. These values come from the 54-point pure-Ar CCP Swarm calculation
at 300 K and 13.3322 Pa. The inputs are the maintained
[P1 workflow](../../examples/workflow_argon_gec_ccp_propagator.yaml) and the
ground-state [argon cross sections](../../examples/cross_sections/argon_application_library.provenance.json).
The source table manifest SHA-256 was
`cf7a395ee6ddeeca00f3b54614631be6d5a816e18189c54a95076a8d525937ad`.

| File | SHA-256 |
| --- | --- |
| `eedf.csv` | `e927397c3657cc3f834f4aa7da1e8dd21147019e95d2e4956f6ebe6eaab5d969` |
| `rates_vs_en.csv` | `15f7e7b1e831cdb2e7364af62009176a40dbad00eef1eec4bfa376a41ffdf719` |

```powershell
py -3 -m pytest -q tests/test_comsol_export.py -k argon_p1_projection
```

The test runs the normal bundle export and uses independent quadrature for
both representations. The raw P1 reference uses constant probability density
in each cell. The import uses a continuous piecewise-linear `f0` and
`F(E) = sqrt(E) * f0(E)`. Integrate `sigma(E) * v(E) * F(E)` separately for
each reaction, splitting intervals at the section knots and threshold.
The exported shape is nonnegative, normalized, and has the requested first
energy moment. The projection matches source rates significant at `1e-5` of
the process peak with a scaled residual of at most `1e-7`. Significance
selects comparisons and leaves the values intact.

At 5 Td, the raw cells reproduce `7.108620556133124e-19 m³/s` for excitation.
The previous shape-only projection gave `8.977346426987355e-19 m³/s`.
The fixture covers this difference and the higher-energy reactions.

## Between-anchor sensitivity

The exporter retains the shape-based mean-energy interpolation. It adjusts
its rates by a factor interpolated between the stored anchors, then fits the
nonnegative projected rows. COMSOL interpolates the resulting table linearly.
Each row below is a separate holdout from the 14-point P1 ICP sweep. Both
predictions use the remaining anchors and the omitted point's mean energy.

| E/N (Td) | Coefficient | Held-out source (m³/s) | Linear in mean energy | Shape-based correction |
| ---: | --- | ---: | ---: | ---: |
| 3 | Excitation | 1.092e-20 | 3.452e-19 | 1.948e-20 |
| 30 | Ionization | 5.874e-19 | 2.804e-18 | 8.225e-19 |
| 300 | Ionization | 1.369e-15 | 3.169e-15 | 1.065e-15 |

These differences motivate additional Swarm points where reactions turn on
or change rapidly. The fit residual describes the stored anchors. The
between-anchor prediction needs its own evidence.

## Other references

The [12-point BOLSIG+ comparison](ar_bolsig_plus_equivalence.md) provides an
independent solver control for the two-term Ar input and its transport
definitions.

NIST reports Langmuir-probe and microwave-interferometer results for the Ar
GEC cell at 13.33 Pa, including a comparison at 250 V RF amplitude. Match
the excitation, geometry, measurement position and diagnostic definition
before applying those data to a voltage-driven run. The present adapter
acceptance uses 1 W. [NIST GEC measurements](https://doi.org/10.6028/jres.100.031)

The Dutton and scanning-drift-tube references on LXCat serve a different
purpose. Compare bulk or flux transport only when the reported definition
matches the solver observable. Retrieve those datasets from the source with
their contributor headers. [LXCat contributors](https://us.lxcat.net/contributors/)
[LXCat redistribution terms](https://master.lxcat.net/instructions/redistribution.php)
