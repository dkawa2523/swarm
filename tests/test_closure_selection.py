from __future__ import annotations

import csv
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import shutil

import pytest
import yaml

from electron_swarm import load_config
from swarm_workflow.comsol.input import ComsolExportError, export_comsol_bundle
from swarm_workflow.comsol.input.bundle_selection import validate_bundle_selection
from swarm_workflow.selection import ClosureSelectionError, context_from_metadata
from swarm_workflow.quality.monte_carlo.policy import (
    MonteCarloPolicyError,
    decide_monte_carlo_closure,
    parse_convergence_policy,
)
from swarm_workflow.campaign.config import (
    _combined_hash,
    _cross_section_file_hashes,
    load_workflow,
)
from swarm_workflow.campaign.sweep import advance_monte_carlo_workflow, run_sweep
from swarm_workflow.tables import build_tables
from product_helpers import base_product_config, write_config
from test_mc_policy import _decision, _tables


def test_missing_physical_context_is_not_reconstructed_from_legacy_paths(
    tmp_path: Path,
) -> None:
    base = write_config(tmp_path, base_product_config(tmp_path, ["two_term"]))

    assert context_from_metadata(
        {
            "base_config_path": str(base),
            "base_config_sha256": sha256(base.read_bytes()).hexdigest(),
        }
    ) is None


@pytest.mark.parametrize("value", [True, 4.5, "4", 0, -1])
def test_convergence_policy_requires_positive_integer_controls(value: object) -> None:
    with pytest.raises(MonteCarloPolicyError, match="positive integers"):
        parse_convergence_policy({"initial_replicas": value})


def test_convergence_policy_parses_bounded_population_and_tail_factors() -> None:
    policy = parse_convergence_policy(
        {"particle_increase_factor": 3, "tail_extension_factor": 5}
    )

    assert policy.particle_increase_factor == 3
    assert policy.tail_extension_factor == 5
    with pytest.raises(MonteCarloPolicyError, match="particle_increase_factor"):
        parse_convergence_policy({"particle_increase_factor": 1})
    with pytest.raises(MonteCarloPolicyError, match="tail_extension_factor"):
        parse_convergence_policy({"tail_extension_factor": 1})


@pytest.mark.parametrize(
    ("failure", "expected_action", "field", "expected_value"),
    (
        ("mean_energy_stationarity", "extend_time", "max_collisions", 128),
        (
            "mc_weighted_growth_lineage_degenerate",
            "increase_particles",
            "particles",
            256,
        ),
        (
            "required_rate_censored_all_zero:Ar:ionization:ionization",
            "extend_tail",
            "tail_max_collisions",
            128,
        ),
    ),
)
def test_advance_writes_and_enforces_only_the_recorded_followup(
    tmp_path: Path,
    failure: str,
    expected_action: str,
    field: str,
    expected_value: int,
) -> None:
    mc, fallback = _tables(tmp_path, failures={5.0: [failure]})
    base = (
        Path(__file__).resolve().parents[1]
        / "examples/argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    xs_hash = _combined_hash(_cross_section_file_hashes(load_config(base)))
    for directory in (mc, fallback):
        path = directory / "manifest.json"
        manifest = json.loads(path.read_text())
        manifest["hashes"]["cross_sections_sha256"] = xs_hash
        if directory == mc:
            manifest["hashes"]["base_config_sha256"] = sha256(
                base.read_bytes()
            ).hexdigest()
        path.write_text(json.dumps(manifest))
    manifest = json.loads((mc / "manifest.json").read_text())
    original = tmp_path / "first.yaml"
    original.write_text(
        yaml.safe_dump(
            {
                "base_config": str(base),
                "database": "first.sqlite",
                "e_over_n_Td": [5.0, 10.0],
                "mixtures": [{"Ar": 1.0}],
                "mc": {
                    "replicas": 4,
                    "base_seed": 1729,
                    "convergence": {},
                    "sampling_plan": manifest["mc_sampling_plan"]["entries"],
                },
            }
        )
    )
    decision_path = tmp_path / "first-decision.json"
    decide_monte_carlo_closure(mc, fallback, attempt=1, output=decision_path)
    assert json.loads(decision_path.read_text())["action"] == expected_action
    output = advance_monte_carlo_workflow(
        original,
        decision_path,
        output_path=tmp_path / "next/second.yaml",
        database_path=tmp_path / "second.sqlite",
    )
    next_workflow = load_workflow(output)
    assert [asdict(entry) for entry in next_workflow.mc_sampling_plan] == json.loads(
        decision_path.read_text()
    )["next_sampling_plan"]
    assert next_workflow.mc_reuse_database == (tmp_path / "first.sqlite").resolve()
    assert getattr(next_workflow.mc_sampling_plan[0], field) == expected_value
    assert next_workflow.mc_sampling_plan[1].max_collisions == 32
    assert all(entry.replicas == 4 for entry in next_workflow.mc_sampling_plan)

    tampered = yaml.safe_load(output.read_text())
    tampered["mc"]["sampling_plan"][0][field] *= 2
    output.write_text(yaml.safe_dump(tampered))
    with pytest.raises(MonteCarloPolicyError, match="next_sampling_plan"):
        load_workflow(output)
    tampered["mc"]["sampling_plan"][0][field] //= 2
    output.write_text(yaml.safe_dump(tampered))
    decision_path.write_bytes(decision_path.read_bytes() + b"\n")
    with pytest.raises(
        MonteCarloPolicyError, match="changed after workflow generation"
    ):
        load_workflow(output)


def _fallback_bundle(tmp_path: Path) -> tuple[Path, Path, Path]:
    mc, _ = _tables(
        tmp_path / "mc-evidence", failures={5.0: ["mean_energy_stationarity"]}
    )
    database = tmp_path / "source.sqlite"
    config = base_product_config(tmp_path, ["two_term"])
    base = write_config(tmp_path, config)
    workflow = tmp_path / "two_term_workflow.yaml"
    workflow.write_text(
        yaml.safe_dump(
            {
                "base_config": str(base),
                "database": str(database),
                "e_over_n_Td": [5.0, 10.0],
                "mixtures": [{"Ar": 1.0}],
            }
        )
    )
    run_sweep(workflow)
    build_tables(database, tmp_path / "tables", source="two_term")
    source = tmp_path / "tables/mixture_0000"
    fallback = json.loads((source / "manifest.json").read_text())
    candidate = json.loads((mc / "manifest.json").read_text())
    candidate["hashes"]["cross_sections_sha256"] = fallback["hashes"][
        "cross_sections_sha256"
    ]
    candidate["mixture"] = fallback["mixture"]
    candidate["physical_context"] = fallback["physical_context"]
    (mc / "manifest.json").write_text(json.dumps(candidate))
    selection = _decision(tmp_path, mc, source, attempt=2)
    assert selection["selected_solver"] == "two_term"
    return source, tmp_path / "decision.json", tmp_path / "bundle"


def test_selected_closure_is_portable_and_coefficient_changes_are_rejected(
    tmp_path: Path,
) -> None:
    source, decision, bundle = _fallback_bundle(tmp_path)
    export_comsol_bundle(source, bundle, selection_path=decision)
    binding = validate_bundle_selection(bundle, required=True)
    assert binding["selected_solver"] == "two_term"
    assert binding["attempt"] == 2
    for name in json.loads((source / "manifest.json").read_text())["tables"]:
        assert (bundle / name).read_bytes() == (source / name).read_bytes()
    relocated = tmp_path / "relocated"
    shutil.copytree(bundle, relocated)
    (source / "quality.csv").write_text("old source no longer available")
    assert validate_bundle_selection(relocated, required=True) == binding
    table = relocated / "transport_vs_mean_energy.csv"
    with table.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows[0]["reduced_mobility_m2_V_s_m3"] = "999"
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ClosureSelectionError, match="differs from solver selection"):
        validate_bundle_selection(relocated, required=True)


def test_export_rejects_rebuilt_source_after_decision_before_touching_output(
    tmp_path: Path,
) -> None:
    source, decision, bundle = _fallback_bundle(tmp_path)
    bundle.mkdir()
    sentinel = bundle / "manifest.json"
    sentinel.write_text("existing bundle")
    path = source / "manifest.json"
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ComsolExportError, match="generation differs"):
        export_comsol_bundle(source, bundle, selection_path=decision)
    assert sentinel.read_text() == "existing bundle"


def test_canonical_mc_workflow_has_one_bounded_first_attempt() -> None:
    root = Path(__file__).resolve().parents[1]
    workflow = load_workflow(
        root / "examples" / "workflow_argon_gec_ccp_monte_carlo.yaml"
    )
    assert len(workflow.mc_sampling_plan) == 15
    assert sum(row.replicas for row in workflow.mc_sampling_plan) == 60
    assert workflow.mc_convergence.maximum_attempts == 2
    assert workflow.mc_convergence.maximum_replicas == 8
    assert workflow.mc_reuse_database is None
