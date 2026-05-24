# Product Architecture

The product runner follows one path:

1. Load schema v2 YAML.
2. Build a solver plan from requested solvers, requested physics, capabilities,
   and feature policy.
3. Execute runnable plan rows.
4. Attach compact result metadata.
5. Write canonical CSV outputs.

Public solver ids are `two_term`, `multi_term`, and `monte_carlo`. Physics
features are requested under `physics.*` and are not solver ids.

Runner orchestration stays thin. Capability and policy decisions live in the
plan layer, solver execution lives in orchestration/executor, and output shape
lives in the writer.
