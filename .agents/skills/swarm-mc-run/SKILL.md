---
name: swarm-mc-run
description: "Execute, resume, assess, and boundedly advance internal Monte Carlo electron-swarm campaigns in this repository. Use for MC execution, timing, statistical qualification, and authorized finite follow-up attempts; leave COMSOL bundle export and solve execution to COMSOL skills."
---

# Run Monte Carlo Campaigns

Apply this skill only in the electron-swarm repository. Read the root AGENTS.md. Load the workflow using swarm-mc-config guidance and print its resolved jobs, seeds, workers, and nominal particle-barrier ceiling before starting.

## Avoid duplicate or unbounded work

- Check live Python processes and the resolved SQLite path. If the same workflow/database is active, observe or report that run instead of launching another.
- Run only the plan encoded in the workflow. Do not increase particles, barriers, replicas, workers, or tail work ad hoc.
- Preserve an existing compatible SQLite database so completed raw cases resume without duplication. Never delete it as a retry mechanism.
- A changed physical context or sampling plan needs a new database, except for an exact advance-mc workflow authorized by the prior immutable decision.
- Stop on structural evidence failures. Do not spend additional sampling on missing provenance, unsupported physics, estimator mismatches, or invalid values.

## Execute one encoded attempt

~~~powershell
$elapsed = Measure-Command {
  py -3 -m swarm_workflow.cli sweep <workflow.yaml>
}
$elapsed.TotalSeconds
~~~

Capture the command, start/end time, wall time, completed and reused cases, worker count, database, and any failing anchor. Distinguish wall time from summed replica CPU time and from the nominal particle-barrier budget.

## Assess statistical qualification

Build tables with the profile required by the downstream observable:

~~~powershell
py -3 -m swarm_workflow.cli build-tables <mc.sqlite> --output <mc-tables> --source monte_carlo --mc-qualification-profile <full_transport-or-function_eedf_restricted_lmea>
~~~

Use full_transport for a full transport claim. Use function_eedf_restricted_lmea only for that declared restricted COMSOL closure. If qualification fails and diagnostic artifacts are needed, repeat table generation with --allow-unqualified-mc; label the result evidence-only and never export it as production input.

Read the per-anchor quality evidence and manifest. An aggregate pass alone is insufficient: require the solver-specific convergence evidence, canonical sampling-plan provenance, independent seeds, and the selected profile's active gates.

## Make a bounded follow-up decision

When the task requires a terminal MC-or-fallback closure decision and qualified two-term tables exist:

~~~powershell
py -3 -m swarm_workflow.cli decide-mc --mc-tables <mc-mixture-tables> --two-term-tables <two-term-mixture-tables> --attempt <n> --output <selection.json>
~~~

Handle the recorded action exactly:

- accept_monte_carlo: stop; the selected source is MC.
- select_two_term: stop; retain MC as validation evidence.
- blocked: stop and report the structural or fallback problem.
- extend_time or add_replicas: generate the exact next workflow with advance-mc; do not hand-edit it.

~~~powershell
py -3 -m swarm_workflow.cli advance-mc <current-workflow.yaml> --decision <selection.json> --output <next-workflow.yaml> --database <new.sqlite>
~~~

Before launching a follow-up, report its incremental and cumulative budget. Execute it when the user's request covers the bounded campaign through terminal qualification; otherwise leave the exact workflow ready. Never exceed maximum_attempts or either particle-barrier ceiling.

Finish with the observed wall time, source database, qualification profile, anchor-level failures, policy action, selected solver if terminal, and artifact paths. A COMSOL bundle is produced later by swarm-comsol-input.
