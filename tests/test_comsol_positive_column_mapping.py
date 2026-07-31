from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from swarm_workflow.comsol_adapter import format_apply_plan
from swarm_workflow.comsol_mapping import load_comsol_mapping


ROOT = Path(__file__).resolve().parents[1]
MAPPING = ROOT / "Model" / "maps" / "positive_column_external.yaml"


def test_canonical_positive_column_mapping_is_minimal() -> None:
    mapping = load_comsol_mapping(MAPPING)

    assert [function.tag for function in mapping.functions] == ["sw_meanE", "sw_muN"]
    assert mapping.model.input_mph.name == "positive_column_1d.mph"
    assert mapping.closure.mean_energy.table == "mean_energy_vs_en.csv"
    assert mapping.run.voltages_V == (20.0, 50.0, 100.0, 200.0)
    assert [(item.feature, item.process_type) for item in mapping.reaction_lookups] == [
        ("eir2", "excitation"),
        ("eir4", "ionization"),
    ]


def test_canonical_mapping_dry_run_describes_only_active_external_inputs() -> None:
    mapping = load_comsol_mapping(MAPPING)
    text = format_apply_plan(SimpleNamespace(mapping=mapping, java_path=Path("apply.java")))

    assert "sw_meanE" in text
    assert "sw_muN" in text
    assert "sw_DLN" not in text
    assert "energy_loss" not in text
    assert "reaction_lookups:" in text
