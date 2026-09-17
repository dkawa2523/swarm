from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path

import pytest

from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SCHEMA,
)
from swarm_workflow.comsol.input.export.bundle import _copy_composite_evidence
from swarm_workflow.comsol.input.export.contracts import ComsolExportError
from swarm_workflow.comsol.input.export.function_source import (
    _read_function_eedf_source,
)


def test_composite_function_eedf_retains_every_selected_anchor_bin(
    tmp_path: Path,
) -> None:
    path = tmp_path / "eedf.csv"
    columns = (
        "electron_energy_eV",
        "energy_width_eV",
        "E_over_N_Td",
        "mean_energy_eV",
        "eedf",
        "pooled_effective_sample_count",
    )
    rows = (
        (0.5, 1.0, 1.0, 1.0, 0.5, ""),
        (1.5, 1.0, 1.0, 1.0, 1.0e-7, ""),
        (0.5, 1.0, 2.0, 2.0, 0.5, 100.0),
        (1.5, 1.0, 2.0, 2.0, 0.1, 2.0),
        (2.5, 1.0, 2.0, 2.0, 1.0e-7, 0.5),
    )
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)
    composition = {
        "schema": ANCHOR_FALLBACK_SCHEMA,
        "primary_solver": "monte_carlo",
        "fallback_solver": "two_term",
        "scope": "low_e_over_n_anchor_fallback",
        "component_mixing_within_anchor": False,
        "whole_closure_replacement": False,
        "high_e_over_n_monte_carlo_preserved": True,
        "postprocess_repair": False,
        "anchors": [
            {"E_over_N_Td": 1.0, "effective_source": "two_term"},
            {"E_over_N_Td": 2.0, "effective_source": "monte_carlo"},
        ],
    }

    prepared = _read_function_eedf_source(
        path,
        source_kind="composite",
        source_composition=composition,
    )

    assert prepared is not None
    grouped, audit = prepared
    assert grouped[1.0][-1][2] == pytest.approx(1.0e-7)
    assert grouped[2.0][-1][2] == pytest.approx(1.0e-7)
    assert audit == {
        "policy": "full_source_support",
        "source_values_modified": False,
    }


def test_composite_evidence_is_copied_only_when_inventory_and_hashes_agree(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    output = tmp_path / "output"
    source.mkdir()
    output.mkdir()
    artifact = source / "plan.json"
    artifact.write_text('{"status":"selected"}\n', encoding="utf-8")
    digest = sha256(artifact.read_bytes()).hexdigest()
    evidence = {"plan.json": {"role": "anchor_fallback_plan", "sha256": digest}}
    manifest = {
        "source": "composite",
        "evidence": evidence,
        "source_composition": {"evidence": evidence},
    }

    copied = _copy_composite_evidence(source, output, manifest)

    assert copied == evidence
    assert (output / "plan.json").read_bytes() == artifact.read_bytes()
    artifact.write_text('{"status":"changed"}\n', encoding="utf-8")
    with pytest.raises(ComsolExportError, match="hash mismatch"):
        _copy_composite_evidence(source, output, manifest)


def test_composite_function_eedf_rejects_unlisted_anchor(tmp_path: Path) -> None:
    path = tmp_path / "eedf.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "electron_energy_eV",
                "energy_width_eV",
                "E_over_N_Td",
                "mean_energy_eV",
                "eedf",
                "pooled_effective_sample_count",
            )
        )
        writer.writerow((0.5, 1.0, 3.0, 1.0, 1.0, 10.0))
    composition = {
        "schema": ANCHOR_FALLBACK_SCHEMA,
        "primary_solver": "monte_carlo",
        "fallback_solver": "two_term",
        "scope": "low_e_over_n_anchor_fallback",
        "component_mixing_within_anchor": False,
        "whole_closure_replacement": False,
        "high_e_over_n_monte_carlo_preserved": True,
        "postprocess_repair": False,
        "anchors": [{"E_over_N_Td": 2.0, "effective_source": "monte_carlo"}],
    }

    with pytest.raises(ComsolExportError, match="selected E/N anchor"):
        _read_function_eedf_source(
            path,
            source_kind="composite",
            source_composition=composition,
        )
