# Developer tools

- `benchmarks/`: reproducible solver and external-reference comparisons.
- `validation/`: evidence-only audits and plots that leave solver results unchanged;
  `evaluate_mc_refinement.py` compares raw MC replicas and particle populations.
- `qualify_propagator_p1_deterministic.py`: formal Propagator core
  qualification.
- `qualify_propagator_target.py`: model-independent, explicit target-range
  refinement qualification.

Runtime solver code belongs under `electron_swarm/`; COMSOL workflow code
belongs under `swarm_workflow/comsol/`. Generated outputs are written below
the ignored `outputs/` directory and are not source artifacts.
