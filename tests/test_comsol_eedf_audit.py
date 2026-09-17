from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

from swarm_workflow.comsol.input.eedf_audit import (
    ComsolEedfAuditError,
    ComsolEedfAuditPlan,
    analyze_comsol_eedf_audit,
    extract_comsol_eedf_audit_log,
    prepare_comsol_eedf_audit,
    read_native_eedf_grid,
)
from swarm_workflow.comsol.input.eedf_audit.contracts import (
    _QUERY_COLUMNS,
    _VALUE_COLUMNS,
)
from swarm_workflow.comsol.input.function_eedf.moments import (
    _tilt_linear_row_to_mean,
)


def _write_eedf(path: Path, *, structure: str = "spreadsheet") -> None:
    assert structure == "spreadsheet"
    energies = np.concatenate(([0.0], np.geomspace(1.0e-5, 40.0, 128)))
    means = np.asarray([1.0, 3.0, 7.0])
    rows = []
    for mean in means:
        seed = np.exp(-energies / max(mean, 1.0e-12))
        rows.append(_tilt_linear_row_to_mean(energies, seed, float(mean)))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("electron_energy_eV", "mean_energy_eV", "eepf_eV_m32"))
        for mean, row in zip(means, rows, strict=True):
            writer.writerows(zip(energies, [mean] * len(energies), row, strict=True))


def _write_contract(plan: ComsolEedfAuditPlan, *, fununit: str = "1") -> None:
    values = {
        "function_names": plan.function_tag,
        "source": "file",
        "struct": plan.import_contract.structure,
        "scaledata": "auto",
        "interp": "linear",
        "extrap": "const",
        "fununit": fununit,
        "argunit": "eV|eV",
        "column_type": "col1|arg|col2|arg|col3|value",
        "imported_name": "eedf.csv",
        "filename": "eedf.csv",
        "problem_count": "0",
        "problems": "",
        "eedf_selection": plan.function_tag,
        "audit_table_sha256": plan.table_sha256,
        "audit_queries_sha256": plan.query_sha256,
    }
    plan.contract_path.write_text(
        "".join(f"{key}\t{value}\n" for key, value in values.items()),
        encoding="utf-8",
    )


def _write_perfect_values(query_path: Path, values_path: Path) -> None:
    with query_path.open("r", encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source))
    with values_path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=_VALUE_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({**row, "comsol_value": row["expected_a"]})


def _scale_values(
    values_path: Path,
    *,
    factor: float,
    predicate: Callable[[dict[str, str]], bool],
    only_largest: bool = False,
) -> None:
    with values_path.open("r", encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source))
    selected = [row for row in rows if predicate(row)]
    if only_largest:
        selected = [max(selected, key=lambda row: abs(float(row["expected_a"])))]
    selected_ids = {row["point_id"] for row in selected}
    with values_path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=_VALUE_COLUMNS)
        writer.writeheader()
        for row in rows:
            if row["point_id"] in selected_ids:
                row["comsol_value"] = str(float(row["comsol_value"]) * factor)
            writer.writerow(row)


def test_prepare_audit_covers_native_interpolation_contract(tmp_path: Path) -> None:
    table = tmp_path / "eedf.csv"
    model = tmp_path / "model.mph"
    _write_eedf(table)
    model.write_bytes(b"test fixture")

    plan = prepare_comsol_eedf_audit(
        model_path=model,
        table_path=table,
        output_directory=tmp_path / "audit",
        component="comp1",
        physics="ptp",
        function_tag="sw_eedf",
        moment_mean_count=3,
    )

    with plan.query_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    assert tuple(reader.fieldnames or ()) == _QUERY_COLUMNS
    assert {row["kind"] for row in rows} == {
        "anchor",
        "cell_center",
        "derivative",
        "energy_midpoint",
        "moment",
        "outside",
    }
    assert plan.point_count == len(rows)
    assert plan.moment_mean_count == 3
    java = plan.java_path.read_text(encoding="utf-8")
    assert 'model.param().evaluate("sw_eedf_audit_value")' in java
    assert 'function.getStringArray("fununit")' in java
    assert 'function.getString("scaledata")' not in java
    assert "function.problem().tags()" in java
    assert ".study(" not in java
    assert ".save(" not in java
    assert ".importData(" not in java


def test_sparse_local_hat_queries_sample_active_support(tmp_path: Path) -> None:
    table = tmp_path / "sparse-hats.csv"
    energies = np.linspace(0.0, 13.0, 27)
    means = (1.0, 2.0, 3.0)
    with table.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("electron_energy_eV", "mean_energy_eV", "eepf_eV_m32"))
        for mean in means:
            values = np.zeros_like(energies)
            values[1::2] = np.exp(-energies[1::2] / mean)
            writer.writerows(zip(energies, [mean] * len(energies), values, strict=True))
    model = tmp_path / "model.mph"
    model.write_bytes(b"test fixture")

    plan = prepare_comsol_eedf_audit(
        model_path=model,
        table_path=table,
        output_directory=tmp_path / "audit",
        component="comp1",
        physics="ptp",
        function_tag="sw_eedf",
        moment_mean_count=3,
    )
    with plan.query_path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    sampled = [row for row in rows if row["kind"] in {"anchor", "derivative"}]
    assert sampled
    assert all(float(row["expected_a"]) > 0.0 for row in sampled)


def test_analyze_audit_checks_values_moments_derivatives_and_rates(
    tmp_path: Path,
) -> None:
    table = tmp_path / "eedf.csv"
    model = tmp_path / "model.mph"
    _write_eedf(table)
    model.write_bytes(b"test fixture")
    plan = prepare_comsol_eedf_audit(
        model_path=model,
        table_path=table,
        output_directory=tmp_path / "audit",
        component="comp1",
        physics="ptp",
        function_tag="sw_eedf",
        moment_mean_count=3,
    )
    _write_contract(plan)
    _write_perfect_values(plan.query_path, plan.values_path)
    cross_sections = tmp_path / "cross_sections.csv"
    cross_sections.write_text(
        "type,energy_eV,cross_section_m2\n"
        "excitation,0,0\n"
        "excitation,2,1e-20\n"
        "excitation,40,2e-20\n",
        encoding="utf-8",
    )

    result = analyze_comsol_eedf_audit(
        plan,
        cross_section_csv=cross_sections,
    )

    assert result["passed"] is True
    assert result["saved_function_contract"]["passed"] is True
    assert result["moments"]["maximum_error"] < 1.0e-7
    assert result["derivatives"]["continuity_contract"].startswith("C0")
    assert result["rates"]["processes"]["excitation"]["normalized_rmse"] == 0.0

    _write_contract(plan, fununit="1/eV^(3/2)")
    failed = analyze_comsol_eedf_audit(plan, cross_section_csv=cross_sections)
    assert failed["passed"] is False
    assert "fununit" in failed["saved_function_contract"]["mismatches"]


def test_native_projection_small_interior_deviation_is_accepted(tmp_path: Path) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    _scale_values(
        plan.values_path,
        factor=1.05,
        predicate=lambda row: row["kind"] == "cell_center",
        only_largest=True,
    )

    result = analyze_comsol_eedf_audit(plan)

    assert result["passed"] is True
    assert result["active_table_binding"]["passed"] is True
    assert result["native_projection_fidelity"]["passed"] is True
    assert (
        result["native_projection_fidelity"]["maximum_relative_error_significant"]
        < 0.10
    )


def test_native_projection_significant_error_above_ten_percent_fails(
    tmp_path: Path,
) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    _scale_values(
        plan.values_path,
        factor=1.20,
        predicate=lambda row: row["kind"] == "cell_center",
        only_largest=True,
    )

    result = analyze_comsol_eedf_audit(plan)

    assert result["passed"] is False
    assert result["active_table_binding"]["passed"] is True
    assert result["native_projection_fidelity"]["passed"] is False
    assert (
        result["native_projection_fidelity"]["maximum_relative_error_significant"]
        > 0.10
    )


def test_native_projection_nrmse_above_two_percent_fails(tmp_path: Path) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    _scale_values(
        plan.values_path,
        factor=1.05,
        predicate=lambda row: (
            row["kind"] == "cell_center"
            or (row["kind"] == "moment" and row["group"].startswith("moment_mid_"))
        ),
    )

    result = analyze_comsol_eedf_audit(plan)

    fidelity = result["native_projection_fidelity"]
    assert result["passed"] is False
    assert result["active_table_binding"]["passed"] is True
    assert fidelity["maximum_relative_error_significant"] < 0.10
    assert fidelity["all_point_normalized_rmse"] > 0.02
    assert fidelity["passed"] is False


def test_strict_anchor_deviation_remains_a_binding_failure(tmp_path: Path) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    _scale_values(
        plan.values_path,
        factor=1.01,
        predicate=lambda row: row["kind"] == "anchor",
        only_largest=True,
    )

    result = analyze_comsol_eedf_audit(plan)

    assert result["passed"] is False
    assert result["active_table_binding"]["passed"] is False
    assert result["native_projection_fidelity"]["passed"] is True


def test_read_native_grid_is_rectangular(tmp_path: Path) -> None:
    table = tmp_path / "eedf.csv"
    _write_eedf(table)
    grid = read_native_eedf_grid(table)
    assert grid.values_eV_m32.shape == (3, 129)
    assert np.all(np.diff(grid.electron_energies_eV) > 0.0)


def _prepared_perfect_audit(
    tmp_path: Path,
) -> ComsolEedfAuditPlan:
    table = tmp_path / "eedf.csv"
    model = tmp_path / "model.mph"
    _write_eedf(table)
    model.write_bytes(b"test fixture")
    plan = prepare_comsol_eedf_audit(
        model_path=model,
        table_path=table,
        output_directory=tmp_path / "audit",
        component="comp1",
        physics="ptp",
        function_tag="sw_eedf",
        moment_mean_count=3,
    )
    _write_contract(plan)
    _write_perfect_values(plan.query_path, plan.values_path)
    return plan


def test_saved_function_format_is_bound_to_actual_file(tmp_path: Path) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    result = analyze_comsol_eedf_audit(plan)
    assert result["passed"] is True
    assert result["saved_function_contract"]["expected"]["struct"] == "spreadsheet"
    assert result["provenance"]["query_sha256"] == plan.query_sha256
    plan.contract_path.write_text(
        plan.contract_path.read_text().replace("struct\tspreadsheet", "struct\tgrid"),
        encoding="utf-8",
    )
    result = analyze_comsol_eedf_audit(plan)
    assert result["passed"] is False
    assert "struct" in result["saved_function_contract"]["mismatches"]


@pytest.mark.parametrize("artifact", ["table_path", "query_path", "model_path"])
def test_audit_rejects_inputs_changed_after_preparation(
    tmp_path: Path, artifact: str
) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    path = getattr(plan, artifact)
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ComsolEedfAuditError, match="changed after preparation"):
        analyze_comsol_eedf_audit(plan)


def test_audit_rejects_values_from_other_probe_coordinates(tmp_path: Path) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    with plan.values_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows[0]["mean_energy_eV"] = "999"
    with plan.values_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=_VALUE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ComsolEedfAuditError, match="mismatched probe coordinates"):
        analyze_comsol_eedf_audit(plan)


def test_old_audit_log_cannot_be_attached_to_new_input(tmp_path: Path) -> None:
    plan = _prepared_perfect_audit(tmp_path)
    contract_lines = plan.contract_path.read_text().splitlines()
    with plan.query_path.open(newline="") as stream:
        queries = list(csv.DictReader(stream))
    text = "\n".join(
        ["SWARM_EEDF_CONTRACT\t" + line for line in contract_lines]
        + [
            f"SWARM_EEDF_VALUE\t{i}\t{row['expected_a']}"
            for i, row in enumerate(queries)
        ]
    )
    log = tmp_path / "audit.log"
    log.write_text(text, encoding="utf-8")
    extract_comsol_eedf_audit_log(plan, log)
    assert analyze_comsol_eedf_audit(plan)["passed"] is True
    log.write_text(text.replace(plan.table_sha256, "0" * 64), encoding="utf-8")
    with pytest.raises(ComsolEedfAuditError, match="different input or query plan"):
        extract_comsol_eedf_audit_log(plan, log)


def test_spreadsheet_cannot_hide_a_missing_grid_point(tmp_path: Path) -> None:
    table = tmp_path / "arbitrary_name.txt"
    _write_eedf(table, structure="spreadsheet")
    lines = table.read_text().splitlines()
    del lines[5]
    table.write_text("\n".join(lines), encoding="utf-8")
    with pytest.raises(ComsolEedfAuditError, match="not rectangular"):
        read_native_eedf_grid(table)
