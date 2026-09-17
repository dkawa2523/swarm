# COMSOL plasma-model adapters

The repository separates electron-swarm data from each spatial COMSOL plasma
model. A canonical Swarm bundle is model-independent; a model adapter owns how
that bundle is applied to one concrete MPH topology and how the resulting solve
is accepted. This supports additional plasma models without turning GEC-CCP or
positive-column assumptions into global behavior.

## Supported adapters

| Adapter | Package | Run command | Map | Target contract |
| --- | --- | --- | --- | --- |
| 1D positive column | `swarm_workflow/comsol/models/positive_column/` | `run-positive-column` | `Model/maps/positive_column_external.yaml` | `docs/comsol_workflow.md` |
| Argon GEC-CCP | `swarm_workflow/comsol/models/gec_ccp/` | `run-gec-ccp` | `comsol_modes/maps/argon_gec_ccp_*.yaml` | `docs/comsol_gec_ccp.md` |
| Argon GEC-ICP | `swarm_workflow/comsol/models/gec_icp/` | `run-gec-icp` | `comsol_modes/maps/argon_gec_icp_*.yaml` | `docs/comsol_gec_icp.md` |

No adapter is a generic fallback. A different ICP, CCP, positive-column,
global, fluid, or hybrid model needs its own adapter when its MPH tags,
formulation, closure, run sequence, or acceptance evidence differs.

## Dependency direction

The dependency direction is fixed:

```text
campaign/quality/tables -> selection <- comsol/input
workflow commands -> model adapter -> comsol/input
                                `-> comsol/runtime
```

- `selection` owns physical context, immutable whole-closure decisions, and
  bounded anchor-fallback plans without knowing COMSOL bundles or MPH models.
- `comsol/input` builds provenance-bound tables and bundles. It knows no MPH
  tags, geometry, studies, datasets, or solver sequence.
- `comsol/runtime` discovers COMSOL, compiles and executes generated Java, and
  records logs and artifact provenance. It knows no plasma model.
- `comsol/models/<model_id>` owns the mapping schema, MPH inspection, closure
  validation, Java generation, run stages, result exports, and acceptance for
  one model family.
- `commands/parser.py` and `commands/handlers.py` expose thin model-specific CLI
  entry points. They contain no model physics.

Shared layers must not import `comsol.models`. Model adapters may depend on the
public shared layers in that direction. Do not introduce a universal mapping
whose optional fields encode several unrelated MPH topologies.

A bundle's `source_policy.qualification_outputs` lists solver quantities covered
by its qualification evidence; it does not claim that every downstream COMSOL
model consumes them. The model adapter's mapping and `closure_ownership` record
the subset actually bound into that MPH. This keeps solver qualification
independent of target-specific closure ownership.

## Adapter contract

A new model adapter must define and test these items together:

1. A typed schema-v2 mapping with the actual model, bundle, closure, run,
   results, and log contract. Preserve the adapter's documented path rule:
   positive-column paths resolve from the repository root, while GEC-CCP paths
   resolve from the mapping directory. A new adapter must choose one explicit
   rule and test it for both its source map and any materialized effective map.
2. An immutable source MPH, a different output MPH, and explicit non-colliding
   generated-code, result, and log locations. An adapter may intentionally
   colocate its generated Java and run artifacts.
3. Explicit ownership of EEDF, transport, reaction rates, elastic energy loss,
   spatial equations, fields, chemistry, walls, circuit, and solver settings.
4. Bundle validation against the exact source solver, hashes, physical context,
   interpolation support, table inventory, and required qualification evidence.
5. Fail-closed validation of bindings whose mismatch could silently select the
   wrong physics, formulation, reaction, study, solution, or dataset. Perform
   offline checks where supported and licensed COMSOL readback where required;
   do not add broad diagnostics for values whose absence already fails clearly.
6. A deterministic dry-run that validates offline-checkable inputs, records
   remaining licensed readbacks, and materializes an auditable plan and Java
   without compiling or launching COMSOL.
7. A staged execution path using `comsol/runtime`, with no stale-class fallback
   and no in-place modification of the source MPH.
8. Model-specific post-solve checks and one canonical terminal summary/status.
   Process success alone is not physical acceptance.
9. A `run-<model-id>` CLI command (snake-case package names become kebab-case
   commands), a maintained target document, and behavior tests for both
   success and fail-closed paths.

The package root should expose only typed errors, plans/summaries, prepare and
execute functions, and user-facing formatters. Mapping, Java, validation,
execution, audit, and plotting concerns remain separate when their size
justifies separate modules; one-use wrappers and compatibility aliases are not
kept.

## Swarm-input compatibility

The current COMSOL table sources are `two_term`, `monte_carlo`, and
`propagator`. `multi_term` is a canonical calculation model but has no COMSOL
table-export contract yet. An adapter must consume only quantities present and
qualified in its selected bundle. Missing diffusion, energy transport, rates,
or EEDF evidence must be closed explicitly by COMSOL or rejected; they are not
borrowed from another solver or manufactured after the solve.

Target-specific refinement or statistical gates belong to the bundle/map
contract. GEC qualification thresholds, Ar reaction tags, RF study sequencing,
positive-column voltage continuation, and LEA/LFA activation are model-local
examples, not shared defaults.

The GEC-ICP adapter is likewise model-local. It binds a steady-DC Function
EEDF and reduced mobility into the 13.56 MHz frequency-transient ICP model;
COMSOL retains diffusion, energy transport, embedded ground/metastable
chemistry, inductive fields, and target-specific convergence acceptance.
