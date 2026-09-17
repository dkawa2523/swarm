"""Materialize the LXCat mixture validation evidence and concise Japanese report."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from electron_swarm.core.config_parser import load_config  # noqa: E402
from electron_swarm.runner import run  # noqa: E402
from swarm_workflow.campaign.config import load_workflow  # noqa: E402
from swarm_workflow.campaign.sweep import _clone_config  # noqa: E402
from swarm_workflow.plots.eedf_contracts import EedfCase  # noqa: E402
from swarm_workflow.plots.eedf_data import (  # noqa: E402
    load_mc_eedf_database,
    load_table_eedf,
    validated_case,
)
from swarm_workflow.plots.eedf_metrics import compare_eedf_cases  # noqa: E402


VALIDATION_ROOT = REPOSITORY_ROOT / "outputs" / "lxcat_swarm_validation"
OUTPUT = VALIDATION_ROOT / "evaluation"
CONFIG_ROOT = REPOSITORY_ROOT / "configs" / "validation" / "lxcat_mixtures"
TOOLS_ROOT = REPOSITORY_ROOT / "tools" / "validation"
FIELDS = (1.0, 10.0, 100.0, 1000.0)
MIXTURES = (
    {
        "key": "ar_o2",
        "label": "Ar/O₂ (90/10)",
        "deterministic_id": 0,
        "attachment_id": 0,
        "tail_thresholds_eV": (4.5, 11.548, 12.072),
        "target_qualification": "propagator_ar_o2_90_10_20260914.json",
        "medium_database": VALIDATION_ROOT / "propagator.sqlite",
        "medium_mixture_id": 0,
    },
    {
        "key": "ar_n2",
        "label": "Ar/N₂ (90/10)",
        "deterministic_id": 1,
        "attachment_id": None,
        "tail_thresholds_eV": (6.725, 11.188, 15.581),
        "target_qualification": "propagator_ar_n2_90_10_20260914.json",
        "medium_database": VALIDATION_ROOT
        / "qualifications"
        / "propagator_ar_n2_medium.sqlite",
        "medium_mixture_id": 0,
    },
    {
        "key": "ar_cl2",
        "label": "Ar/Cl₂ (90/10)",
        "deterministic_id": 2,
        "attachment_id": 1,
        "tail_thresholds_eV": (3.252, 10.54, 11.49),
        "target_qualification": "propagator_ar_cl2_90_10_20260914.json",
        "medium_database": VALIDATION_ROOT
        / "qualifications"
        / "propagator_ar_cl2_medium.sqlite",
        "medium_mixture_id": 0,
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_qualification(path: Path) -> dict[float, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return {
            float(row["E_over_N_Td"]): row for row in csv.DictReader(stream)
        }


def scalar_rows(database: Path, mixture_id: int) -> dict[float, dict[str, float]]:
    connection = sqlite3.connect(database)
    try:
        rows = connection.execute(
            """
            SELECT e_over_n_Td, mean, standard_error, relative_standard_error
            FROM aggregate_scalars
            WHERE mixture_id=? AND scalar_group='case'
              AND scalar_name='mean_energy_eV'
            ORDER BY e_over_n_Td
            """,
            (mixture_id,),
        ).fetchall()
        return {
            float(field): {
                "mean": float(mean),
                "standard_error": (
                    math.nan if standard_error is None else float(standard_error)
                ),
                "relative_standard_error": (
                    math.nan if relative_standard_error is None else float(relative_standard_error)
                ),
            }
            for field, mean, standard_error, relative_standard_error in rows
        }
    finally:
        connection.close()


def net_growth_rows(database: Path, mixture_id: int) -> dict[float, float]:
    connection = sqlite3.connect(database)
    try:
        return {
            float(field): float(value)
            for field, value in connection.execute(
                """
                SELECT e_over_n_Td, net_ionization_frequency_s
                FROM cases
                WHERE mixture_id=? AND replicate=0
                ORDER BY e_over_n_Td
                """,
                (mixture_id,),
            )
        }
    finally:
        connection.close()


def compact_metrics(reference: EedfCase, candidate: EedfCase) -> dict[str, float]:
    values = compare_eedf_cases(reference, candidate)
    return {
        "total_variation": values["total_variation"],
        "shape_total_variation": values["shape_total_variation"],
        "hellinger_distance": values["hellinger_distance"],
        "jensen_shannon_divergence_nats": values[
            "jensen_shannon_divergence_nats"
        ],
        "wasserstein_1_eV": values["wasserstein_1_eV"],
        "mean_energy_relative_difference": values[
            "reported_mean_energy_relative_difference"
        ],
    }


def load_primary_datasets() -> dict[str, dict[str, Any]]:
    datasets: dict[str, dict[str, Any]] = {}
    for mixture in MIXTURES:
        key = str(mixture["key"])
        deterministic_id = int(mixture["deterministic_id"])
        datasets[key] = {
            "two_term": load_table_eedf(
                VALIDATION_ROOT
                / "tables"
                / "two_term"
                / f"mixture_{deterministic_id:04d}",
                "two_term",
            ),
            "propagator": load_table_eedf(
                VALIDATION_ROOT
                / "tables"
                / "propagator"
                / f"mixture_{deterministic_id:04d}",
                "propagator",
            ),
            "monte_carlo": load_mc_eedf_database(
                VALIDATION_ROOT / key / "monte_carlo_analysis.sqlite",
                mixture_id=0,
            ),
        }
    return datasets


def primary_comparison(
    datasets: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    metrics_rows: list[dict[str, Any]] = []
    mean_rows: list[dict[str, Any]] = []
    qualification_rows: list[dict[str, Any]] = []
    for mixture in MIXTURES:
        key = str(mixture["key"])
        data = datasets[key]
        qualification_path = (
            VALIDATION_ROOT
            / key
            / "tables"
            / "monte_carlo"
            / "mixture_0000"
            / "mc_qualification.csv"
        )
        qualification = read_qualification(qualification_path)
        mc_scalars = scalar_rows(
            VALIDATION_ROOT / key / "monte_carlo_analysis.sqlite", 0
        )
        for field in FIELDS:
            cases = {
                solver: dataset.cases[field] for solver, dataset in data.items()
            }
            active_qualified = int(
                qualification[field]["active_closure_quality_passed"]
            )
            overall_qualified = int(qualification[field]["passed"])
            mean_rows.append(
                {
                    "mixture": key,
                    "E_over_N_Td": field,
                    "two_term_mean_energy_eV": cases["two_term"].reported_mean_energy_eV,
                    "propagator_mean_energy_eV": cases[
                        "propagator"
                    ].reported_mean_energy_eV,
                    "monte_carlo_mean_energy_eV": cases[
                        "monte_carlo"
                    ].reported_mean_energy_eV,
                    "monte_carlo_mean_energy_standard_error_eV": mc_scalars[field][
                        "standard_error"
                    ],
                    "monte_carlo_mean_energy_rse": mc_scalars[field][
                        "relative_standard_error"
                    ],
                    "monte_carlo_active_closure_qualified": active_qualified,
                    "monte_carlo_full_qualification_passed": overall_qualified,
                }
            )
            pairs = (
                ("propagator_vs_two_term", cases["two_term"], cases["propagator"]),
                ("raw_mc_vs_two_term", cases["two_term"], cases["monte_carlo"]),
                ("raw_mc_vs_propagator", cases["propagator"], cases["monte_carlo"]),
            )
            for pair, reference, candidate in pairs:
                metrics_rows.append(
                    {
                        "mixture": key,
                        "E_over_N_Td": field,
                        "pair": pair,
                        **compact_metrics(reference, candidate),
                    }
                )
            qualification_rows.append(
                {
                    "mixture": key,
                    "E_over_N_Td": field,
                    "aggregate_quality_passed": int(
                        qualification[field]["aggregate_quality_passed"]
                    ),
                    "full_qualification_passed": overall_qualified,
                    "active_closure_quality_passed": active_qualified,
                    "failure_reasons_json": qualification[field][
                        "failure_reasons_json"
                    ],
                    "active_closure_failure_reasons_json": qualification[field][
                        "active_closure_failure_reasons_json"
                    ],
                    "mobility_rse": qualification[field]["mobility_rse"],
                    "max_major_rate_rse": qualification[field]["max_major_rate_rse"],
                    "eedf_normalization_error": qualification[field][
                        "eedf_normalization_error"
                    ],
                }
            )
    return metrics_rows, mean_rows, qualification_rows


def propagator_refinement() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for mixture in MIXTURES:
        artifact_path = (
            VALIDATION_ROOT
            / "qualifications"
            / str(mixture["target_qualification"])
        )
        artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        fine = {float(row["E_over_N_Td"]): row for row in artifact["fine_runs"]}
        medium_growth = net_growth_rows(
            Path(mixture["medium_database"]), int(mixture["medium_mixture_id"])
        )
        for refinement in artifact["medium_fine_refinement"]:
            field = float(refinement["E_over_N_Td"])
            differences = refinement["relative_differences"]
            task_eedf_refinement_passed = bool(
                fine[field]["quality_gates_passed"]
                and refinement["eedf_weighted_L1"] <= refinement["eedf_limit"]
                and differences["mean_energy"] <= refinement["scalar_limit"]
                and differences["drift_velocity"] <= refinement["scalar_limit"]
            )
            rows.append(
                {
                    "mixture": mixture["key"],
                    "E_over_N_Td": field,
                    "fine_core_quality_passed": int(
                        fine[field]["quality_gates_passed"]
                    ),
                    "task_eedf_refinement_passed": int(task_eedf_refinement_passed),
                    "formal_all_scalar_refinement_passed": int(refinement["passed"]),
                    "eedf_weighted_L1": refinement["eedf_weighted_L1"],
                    "mean_energy_relative_difference": differences["mean_energy"],
                    "drift_velocity_relative_difference": differences["drift_velocity"],
                    "growth_frequency_relative_difference": differences[
                        "growth_frequency"
                    ],
                    "medium_growth_frequency_s_inv": medium_growth[field],
                    "fine_growth_frequency_s_inv": fine[field][
                        "growth_frequency_s_inv"
                    ],
                    "fine_elapsed_s": fine[field]["elapsed_s"],
                    "formal_target_qualified": int(
                        artifact["decision"]["target_refinement_qualified"]
                    ),
                }
            )
    return rows


def attachment_comparison(
    datasets: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[tuple[int, float], dict[str, Any]]]:
    matrix_path = (
        VALIDATION_ROOT / "attachment_sensitivity" / "propagator_case_matrix.json"
    )
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    p1_attachment = {
        (int(row["mixture_id"]), float(row["E_over_N_Td"])): row
        for row in matrix["cases"]
    }
    two_term_attachment = {
        0: load_table_eedf(
            VALIDATION_ROOT
            / "attachment_sensitivity"
            / "tables"
            / "two_term"
            / "mixture_0000",
            "two_term_attachment",
        ),
        1: load_table_eedf(
            VALIDATION_ROOT
            / "attachment_sensitivity"
            / "tables"
            / "two_term"
            / "mixture_0001",
            "two_term_attachment",
        ),
    }
    no_attachment_growth = {
        "two_term": {
            0: net_growth_rows(VALIDATION_ROOT / "two_term.sqlite", 0),
            2: net_growth_rows(VALIDATION_ROOT / "two_term.sqlite", 2),
        },
        "propagator": {
            0: net_growth_rows(VALIDATION_ROOT / "propagator.sqlite", 0),
            2: net_growth_rows(VALIDATION_ROOT / "propagator.sqlite", 2),
        },
    }
    attachment_growth = {
        0: net_growth_rows(
            VALIDATION_ROOT / "attachment_sensitivity" / "two_term.sqlite", 0
        ),
        1: net_growth_rows(
            VALIDATION_ROOT / "attachment_sensitivity" / "two_term.sqlite", 1
        ),
    }
    rows: list[dict[str, Any]] = []
    for mixture in (MIXTURES[0], MIXTURES[2]):
        key = str(mixture["key"])
        base_id = int(mixture["deterministic_id"])
        attachment_id = int(mixture["attachment_id"])
        for field in FIELDS:
            no_tt = datasets[key]["two_term"].cases[field]
            yes_tt = two_term_attachment[attachment_id].cases[field]
            rows.append(
                {
                    "mixture": key,
                    "solver": "two_term",
                    "E_over_N_Td": field,
                    "status": "passed",
                    "no_attachment_mean_energy_eV": no_tt.reported_mean_energy_eV,
                    "attachment_mean_energy_eV": yes_tt.reported_mean_energy_eV,
                    "no_attachment_net_growth_s_inv": no_attachment_growth[
                        "two_term"
                    ][base_id][field],
                    "attachment_net_growth_s_inv": attachment_growth[attachment_id][
                        field
                    ],
                    **compact_metrics(no_tt, yes_tt),
                }
            )
            no_p1 = datasets[key]["propagator"].cases[field]
            evidence = p1_attachment[(attachment_id, field)]
            if evidence["status"] != "passed":
                rows.append(
                    {
                        "mixture": key,
                        "solver": "propagator",
                        "E_over_N_Td": field,
                        "status": "failed",
                        "error": evidence["error"],
                        "no_attachment_mean_energy_eV": no_p1.reported_mean_energy_eV,
                        "no_attachment_net_growth_s_inv": no_attachment_growth[
                            "propagator"
                        ][base_id][field],
                    }
                )
                continue
            yes_p1 = validated_case(
                solver="propagator_attachment",
                e_over_n_td=field,
                energy_eV=evidence["energy_eV"],
                widths_eV=evidence["energy_widths_eV"],
                density_eV_inv=evidence["eedf"],
                reported_mean_energy_eV=evidence["mean_energy_eV"],
            )
            rows.append(
                {
                    "mixture": key,
                    "solver": "propagator",
                    "E_over_N_Td": field,
                    "status": "passed",
                    "no_attachment_mean_energy_eV": no_p1.reported_mean_energy_eV,
                    "attachment_mean_energy_eV": yes_p1.reported_mean_energy_eV,
                    "no_attachment_net_growth_s_inv": no_attachment_growth[
                        "propagator"
                    ][base_id][field],
                    "attachment_net_growth_s_inv": evidence[
                        "net_ionization_frequency_s_inv"
                    ],
                    **compact_metrics(no_p1, yes_p1),
                }
            )
    return rows, p1_attachment


def probe_mc_attachment() -> list[dict[str, Any]]:
    workflow = load_workflow(
        CONFIG_ROOT / "workflow_propagator_attachment_sensitivity.yaml"
    )
    base = load_config(workflow.base_config_path)
    evidence: list[dict[str, Any]] = []
    for mixture in workflow.mixtures:
        config = _clone_config(
            base,
            mixture=mixture,
            e_over_n_values=(1.0,),
            solver_ids=["monte_carlo"],
            replicate=0,
            seed=2026091499,
        )
        try:
            run(config, write=False)
            evidence.append(
                {"mixture_id": mixture.mixture_id, "status": "unexpectedly_ran"}
            )
        except NotImplementedError as exc:
            evidence.append(
                {
                    "mixture_id": mixture.mixture_id,
                    "status": "unsupported",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )
    return evidence


def source_inventory() -> list[dict[str, Any]]:
    definitions = (
        ("Ar_Biagi_20260914.txt", "Biagi", "Ar", "full_and_common", None),
        ("N2_Biagi_20260914.txt", "Biagi", "N2", "full_and_common", None),
        ("O2_Biagi_20260914.txt", "Biagi", "O2", "full", None),
        (
            "O2_Biagi_binary_with_DA_20260914.txt",
            "Biagi",
            "O2",
            "binary_attachment_sensitivity",
            "O2_Biagi_20260914.txt",
        ),
        (
            "O2_Biagi_no_attachment_20260914.txt",
            "Biagi",
            "O2",
            "three_solver_common",
            "O2_Biagi_20260914.txt",
        ),
        ("Cl2_MuroranIT_20260914.txt", "MuroranIT", "Cl2", "full", None),
        (
            "Cl2_MuroranIT_no_attachment_20260914.txt",
            "MuroranIT",
            "Cl2",
            "three_solver_common",
            "Cl2_MuroranIT_20260914.txt",
        ),
    )
    rows: list[dict[str, Any]] = []
    raw_root = VALIDATION_ROOT / "raw_lxcat"
    keyword = re.compile(
        r"^(ELASTIC|EFFECTIVE|EXCITATION|IONIZATION|ATTACHMENT)$", re.MULTILINE
    )
    for name, database, species, role, parent in definitions:
        path = raw_root / name
        text = path.read_text(encoding="utf-8", errors="replace")
        types: dict[str, int] = {}
        for process_type in keyword.findall(text):
            types[process_type.lower()] = types.get(process_type.lower(), 0) + 1
        rows.append(
            {
                "file": str(path),
                "sha256": sha256(path),
                "database": database,
                "species": species,
                "role": role,
                "parent_file": None if parent is None else str(raw_root / parent),
                "process_count": sum(types.values()),
                "process_types": types,
            }
        )
    return rows


def support_limit(case: EedfCase, retained_mass: float = 1.0 - 1.0e-9) -> float:
    cumulative = np.cumsum(case.density_eV_inv * case.widths_eV)
    index = min(int(np.searchsorted(cumulative, retained_mass)), len(case.widths_eV) - 1)
    return float(case.edges_eV[index + 1])


def plot_attachment_sensitivity(
    datasets: dict[str, dict[str, Any]],
    p1_attachment: dict[tuple[int, float], dict[str, Any]],
) -> list[Path]:
    outputs: list[Path] = []
    for mixture in (MIXTURES[0], MIXTURES[2]):
        key = str(mixture["key"])
        label = str(mixture["label"])
        attachment_id = int(mixture["attachment_id"])
        tt_attachment = load_table_eedf(
            VALIDATION_ROOT
            / "attachment_sensitivity"
            / "tables"
            / "two_term"
            / f"mixture_{attachment_id:04d}",
            "two_term_attachment",
        )
        figure, axes = plt.subplots(2, 2, figsize=(13.5, 8.5))
        for axis, field in zip(axes.flat, FIELDS, strict=True):
            base_tt = datasets[key]["two_term"].cases[field]
            base_p1 = datasets[key]["propagator"].cases[field]
            attach_tt = tt_attachment.cases[field]
            lines: list[tuple[str, EedfCase, str, str]] = [
                ("two-term, no attachment", base_tt, "#D97706", "-"),
                ("two-term, attachment", attach_tt, "#DC2626", "--"),
                ("propagator, no attachment", base_p1, "#4D7C0F", "-"),
            ]
            evidence = p1_attachment[(attachment_id, field)]
            if evidence["status"] == "passed":
                attach_p1 = validated_case(
                    solver="propagator_attachment",
                    e_over_n_td=field,
                    energy_eV=evidence["energy_eV"],
                    widths_eV=evidence["energy_widths_eV"],
                    density_eV_inv=evidence["eedf"],
                    reported_mean_energy_eV=evidence["mean_energy_eV"],
                )
                lines.append(
                    ("propagator, attachment", attach_p1, "#2563EB", "--")
                )
            else:
                axis.text(
                    0.98,
                    0.92,
                    "attachment P1: non-converged",
                    transform=axis.transAxes,
                    ha="right",
                    va="top",
                    color="#991B1B",
                    fontsize=9,
                )
            for line_label, case, color, style in lines:
                axis.step(
                    case.edges_eV[:-1],
                    np.maximum(case.density_eV_inv, 1.0e-300),
                    where="post",
                    color=color,
                    linestyle=style,
                    linewidth=1.7,
                    label=line_label,
                )
            x_limit = max(support_limit(case) for _, case, _, _ in lines)
            axis.set_xlim(0.0, max(2.0, 1.08 * x_limit))
            axis.set_ylim(1.0e-9, 2.0)
            axis.set_yscale("log")
            axis.grid(True, alpha=0.25)
            axis.set_title(f"E/N = {field:g} Td")
            axis.set_xlabel("electron energy, E (eV)")
            axis.set_ylabel("F(E) (eV$^{-1}$)")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        figure.subplots_adjust(top=0.82, hspace=0.34, wspace=0.16)
        figure.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.925),
            ncol=2,
            frameon=False,
        )
        figure.suptitle(
            f"{label}: binary-attachment EEDF sensitivity (MC attachment unsupported)",
            fontsize=15,
            fontweight="bold",
            y=0.985,
        )
        png = OUTPUT / f"{key}_attachment_eedf_sensitivity.png"
        svg = OUTPUT / f"{key}_attachment_eedf_sensitivity.svg"
        figure.savefig(png, dpi=180)
        figure.savefig(svg)
        plt.close(figure)
        outputs.extend((png, svg))
    return outputs


def markdown_report(
    metrics_rows: list[dict[str, Any]],
    mean_rows: list[dict[str, Any]],
    qualification_rows: list[dict[str, Any]],
    refinement_rows: list[dict[str, Any]],
    attachment_rows: list[dict[str, Any]],
    sources: list[dict[str, Any]],
    mc_attachment: list[dict[str, Any]],
) -> str:
    metrics = {
        (row["mixture"], row["E_over_N_Td"], row["pair"]): row
        for row in metrics_rows
    }
    mc_total = len(qualification_rows)
    mc_aggregate_passed = sum(
        int(row["aggregate_quality_passed"]) for row in qualification_rows
    )
    mc_full_passed = sum(
        int(row["full_qualification_passed"]) for row in qualification_rows
    )
    mc_active_passed = sum(
        int(row["active_closure_quality_passed"]) for row in qualification_rows
    )
    active_by_mixture: list[str] = []
    for mixture in MIXTURES:
        passed_fields = [
            float(row["E_over_N_Td"])
            for row in qualification_rows
            if row["mixture"] == mixture["key"]
            and int(row["active_closure_quality_passed"])
        ]
        field_label = (
            "/".join(f"{field:g}" for field in passed_fields) + " Td"
            if passed_fields
            else "なし"
        )
        active_by_mixture.append(f"{mixture['label']}={field_label}")
    high_mc_p1 = [
        float(row["total_variation"])
        for row in metrics_rows
        if float(row["E_over_N_Td"]) == 1000.0
        and row["pair"] == "raw_mc_vs_propagator"
    ]
    high_mc_two_term = [
        float(row["total_variation"])
        for row in metrics_rows
        if float(row["E_over_N_Td"]) == 1000.0
        and row["pair"] == "raw_mc_vs_two_term"
    ]
    p1_closer_count = sum(
        mc_p1 < mc_two_term
        for mc_p1, mc_two_term in zip(
            high_mc_p1,
            high_mc_two_term,
            strict=True,
        )
    )
    lines = [
        "# LXCat混合気体 Swarm EEDF検証（2026-09-14）",
        "",
        "90% Ar + 10%添加気体、300 K、13.3 Pa、DC・B=0、等方散乱閉包、e-e衝突なしで、1/10/100/1000 Tdを比較した。主比較は三手法で同一の二体・非付着衝突集合を用いる。",
        "",
        "## 断面積データ",
        "",
        "Ar・N₂・O₂はLXCat Biagi、Cl₂はLXCat MuroranITの特定complete setを使用した。Ar/N₂は各46過程、O₂はfull 15過程・三手法共通13過程・二体付着感度14過程、Cl₂はfull 21過程・三手法共通20過程である。取得原文はLXCatの再配布制約に従い、ローカルのignored outputsとして保持する。各ファイルのSHA-256と過程内訳は `source_inventory.json` に固定した。",
        "",
        "## 主比較",
        "",
        "| 混合 | E/N [Td] | 平均E 2T/P1/MC [eV] | TV(P1,2T) | TV(MC,2T) | TV(MC,P1) | MC active-closure status |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in mean_rows:
        key = str(row["mixture"])
        field = float(row["E_over_N_Td"])
        mixture_label = next(item["label"] for item in MIXTURES if item["key"] == key)
        p1 = metrics[(key, field, "propagator_vs_two_term")]["total_variation"]
        mc_tt = metrics[(key, field, "raw_mc_vs_two_term")]["total_variation"]
        mc_p1 = metrics[(key, field, "raw_mc_vs_propagator")]["total_variation"]
        qualified = (
            "pass"
            if row["monte_carlo_active_closure_qualified"]
            else "not passed (raw shown)"
        )
        lines.append(
            f"| {mixture_label} | {field:g} | "
            f"{row['two_term_mean_energy_eV']:.4g} / "
            f"{row['propagator_mean_energy_eV']:.4g} / "
            f"{row['monte_carlo_mean_energy_eV']:.4g} | "
            f"{p1:.4f} | {mc_tt:.4f} | {mc_p1:.4f} | {qualified} |"
        )
    maximum = max(metrics_rows, key=lambda row: row["total_variation"])
    p1_max_eedf = max(row["eedf_weighted_L1"] for row in refinement_rows)
    p1_max_mean = max(row["mean_energy_relative_difference"] for row in refinement_rows)
    p1_max_drift = max(row["drift_velocity_relative_difference"] for row in refinement_rows)
    lines.extend(
        [
            "",
            "## 数値品質と判定",
            "",
            f"- 最大の三手法間TV距離は {maximum['mixture']} {maximum['E_over_N_Td']:g} Td の {maximum['pair']} = {maximum['total_variation']:.4f}。",
            "- 二項近似は12/12点で反復・tail判定を完了した。Propagator mediumも12/12点でコア収束した。",
            "- ただしrate-tail警告は残る。二項近似はO₂/N₂の1000 TdとCl₂の全4点、PropagatorはO₂/N₂の1000 TdとCl₂の10/100/1000 Tdで警告。主に高しきい値反応率または1000 Tdのtail確率に関するもので、EEDFのmedium/fine格子収束とは区別している。",
            f"- Propagator fine 12/12点はコア品質を通過し、medium/fineの最大EEDF weighted-L1={p1_max_eedf:.4g}、平均エネルギー差={p1_max_mean:.3%}、ドリフト差={p1_max_drift:.3%}。要求されたEEDF比較としては格子収束している。",
            "- ただし既存の正式target qualificationは3混合ともfalse。1/10 Tdのほぼゼロの正味成長率にも一律1%相対差を課すためであり、このformal failure自体は変更・免除していない。",
            f"- MCの基本aggregate品質は{mc_aggregate_passed}/{mc_total}点、transport等を含むfull qualificationは{mc_full_passed}/{mc_total}点、downstream active-closure複合判定は{mc_active_passed}/{mc_total}点（{', '.join(active_by_mixture)}）。後者はEEDF形状単独の合否ではない。図とsolver-pair指標は判定による置換を行わず、全anchorでraw MC・二項近似・P1を独立表示する。",
            f"- 1000 Tdではraw MC–P1のTV={min(high_mc_p1):.4f}–{max(high_mc_p1):.4f}、raw MC–2Tは{min(high_mc_two_term):.4f}–{max(high_mc_two_term):.4f}。P1の方がMCに近い混合は{p1_closer_count}/{len(high_mc_p1)}件である。",
            "",
            "## 付着過程",
            "",
            "O₂は三体付着を二体断面積として誤用せず、解離性付着だけを保持した14過程で感度計算した。Cl₂は全21過程を使用した。二項近似は8/8点完了、PropagatorはAr/O₂ 4/4点とAr/Cl₂ 3/4点が完了し、Ar/Cl₂ 1 Tdは `no_normalizable_principal_mode` で不合格だった。MCは両混合でattachment loss未実装を明示拒否した。",
            "",
            "## データ上の制約",
            "",
            "- LXCatに原子Clの電子散乱complete setは存在しないため、塩素ケースはAr/Cl₂として明示した。",
            "- Biagi N₂は回転過程を含まず、原典注記の適用範囲は E/N > 1 Td。従って1 Tdは境界参考値である。",
            "- Biagi O₂は回転温度0 K。主三手法比較から付着を除いたため、電気陰性プラズマの生成消滅を表す結果ではなく、共通衝突演算子上のsolver比較である。",
            "- integral cross sectionsを等方散乱として用いており、DCSベースの角度散乱比較ではない。",
            "",
            "## 成果物",
            "",
            "- `solver_pair_metrics.csv`: 全36 solver pair指標",
            "- `mean_energy_comparison.csv`: 全12条件の平均エネルギー",
            "- `monte_carlo_qualification.csv`: MCの全品質証拠",
            "- `propagator_refinement.csv`: medium/fine全12条件",
            "- `attachment_sensitivity.csv`: 付着あり/なし比較",
            "- 各混合の `eedf_physical_comparison.png/.svg`、`eedf_mc_status.json` と本ディレクトリの付着感度図",
            "- `validation_manifest.json`: 入力・出力SHA-256と選択条件",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    datasets = load_primary_datasets()
    metrics_rows, mean_rows, qualification_rows = primary_comparison(datasets)
    refinement_rows = propagator_refinement()
    attachment_rows, p1_attachment = attachment_comparison(datasets)
    mc_attachment = probe_mc_attachment()
    sources = source_inventory()

    generated = {
        "solver_pair_metrics": OUTPUT / "solver_pair_metrics.csv",
        "mean_energy_comparison": OUTPUT / "mean_energy_comparison.csv",
        "monte_carlo_qualification": OUTPUT / "monte_carlo_qualification.csv",
        "propagator_refinement": OUTPUT / "propagator_refinement.csv",
        "attachment_sensitivity": OUTPUT / "attachment_sensitivity.csv",
        "mc_attachment_capability": OUTPUT / "mc_attachment_capability.json",
        "source_inventory": OUTPUT / "source_inventory.json",
        "report": OUTPUT / "README_ja.md",
    }
    write_csv(generated["solver_pair_metrics"], metrics_rows)
    write_csv(generated["mean_energy_comparison"], mean_rows)
    write_csv(generated["monte_carlo_qualification"], qualification_rows)
    write_csv(generated["propagator_refinement"], refinement_rows)
    write_csv(generated["attachment_sensitivity"], attachment_rows)
    generated["mc_attachment_capability"].write_text(
        json.dumps(mc_attachment, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    generated["source_inventory"].write_text(
        json.dumps(sources, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    figures = plot_attachment_sensitivity(datasets, p1_attachment)
    generated["report"].write_text(
        markdown_report(
            metrics_rows,
            mean_rows,
            qualification_rows,
            refinement_rows,
            attachment_rows,
            sources,
            mc_attachment,
        ),
        encoding="utf-8",
    )

    input_paths: list[Path] = []
    input_paths.extend(Path(row["file"]) for row in sources)
    input_paths.extend(sorted(CONFIG_ROOT.glob("*.yaml")))
    input_paths.extend(sorted(TOOLS_ROOT.glob("*.py")))
    input_paths.extend(sorted(TOOLS_ROOT.glob("*.ps1")))
    input_paths.extend(
        [
            VALIDATION_ROOT / "two_term.sqlite",
            VALIDATION_ROOT / "propagator.sqlite",
            VALIDATION_ROOT / "attachment_sensitivity" / "two_term.sqlite",
            VALIDATION_ROOT / "ar_o2" / "monte_carlo.sqlite",
            VALIDATION_ROOT / "ar_n2" / "monte_carlo.sqlite",
            VALIDATION_ROOT / "ar_cl2" / "monte_carlo.sqlite",
            VALIDATION_ROOT
            / "attachment_sensitivity"
            / "propagator_case_matrix.json",
            REPOSITORY_ROOT
            / "docs"
            / "dev"
            / "results"
            / "propagator_p1_deterministic_qualification_20260908.json",
        ]
    )
    for mix in MIXTURES:
        mc_root = VALIDATION_ROOT / str(mix["key"])
        input_paths.extend(
            [
                mc_root / "monte_carlo_analysis.sqlite",
                mc_root / "monte_carlo_analysis.sqlite.projection.json",
            ]
        )
    input_paths.extend(
        VALIDATION_ROOT / "qualifications" / str(mix["target_qualification"])
        for mix in MIXTURES
    )
    input_paths.extend(
        VALIDATION_ROOT / str(mix["key"]) / "eedf_comparison" / "eedf_comparison_manifest.json"
        for mix in MIXTURES
    )
    primary_outputs = [
        VALIDATION_ROOT / str(mix["key"]) / "eedf_comparison" / filename
        for mix in MIXTURES
        for filename in (
            "eedf_physical_comparison.png",
            "eedf_physical_comparison.svg",
            "eedf_difference_metrics.png",
            "eedf_difference_metrics.svg",
            "eedf_comparison_metrics.csv",
            "eedf_mc_status.json",
        )
    ]
    all_outputs = [*generated.values(), *figures, *primary_outputs]
    manifest = {
        "schema": "swarm.validation.lxcat_mixture_eedf.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "conditions": {
                "gas_temperature_K": 300.0,
                "pressure_Pa": 13.3,
                "mixtures": [mix["label"] for mix in MIXTURES],
                "E_over_N_Td": list(FIELDS),
            },
            "solvers": ["two_term", "propagator", "monte_carlo"],
            "common_operator": "attachment excluded for O2 and Cl2",
            "attachment_sensitivity": "two_term and propagator core solves; MC unsupported",
        },
        "lxcat_selection": {
            "Biagi": {"database_id": 6, "Ar_target_id": 37, "N2_target_id": 10, "O2_target_id": 46},
            "MuroranIT": {"database_id": 836, "Cl2_target_id": 95},
            "O2_binary_attachment_exclusion": {
                "excluded_process_id": 66679,
                "excluded": "three-body attachment normalized to gas density",
                "retained_process_id": 66680,
                "retained": "dissociative attachment",
            },
        },
        "quality_summary": {
            "two_term_completed": 12,
            "propagator_medium_completed": 12,
            "propagator_fine_core_quality_passed": sum(
                row["fine_core_quality_passed"] for row in refinement_rows
            ),
            "propagator_task_eedf_refinement_passed": sum(
                row["task_eedf_refinement_passed"] for row in refinement_rows
            ),
            "propagator_formal_target_qualifications_passed": 0,
            "mc_aggregate_quality_passed": sum(
                row["aggregate_quality_passed"] for row in qualification_rows
            ),
            "mc_full_qualification_passed": sum(
                row["full_qualification_passed"] for row in qualification_rows
            ),
            "mc_active_closure_quality_passed": sum(
                row["active_closure_quality_passed"] for row in qualification_rows
            ),
            "attachment_two_term_completed": 8,
            "attachment_propagator_completed": sum(
                row["status"] == "passed"
                for row in attachment_rows
                if row["solver"] == "propagator"
            ),
        },
        "mc_eedf_grid_projection": {
            "reason": "threshold and regular-grid boundary differed only at floating-point resolution",
            "method": "probability-mass-conservative removal of zero-measure cells in analysis copies",
            "source_databases_retained": True,
            "manifests": [
                str(
                    VALIDATION_ROOT
                    / str(mix["key"])
                    / "monte_carlo_analysis.sqlite.projection.json"
                )
                for mix in MIXTURES
            ],
        },
        "numeric_metric_method": "Jensen-Shannon uses log-sum-exp; no EEDF clipping or smoothing",
        "input_files": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in input_paths
        ],
        "output_files": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in all_outputs
        ],
    }
    manifest_path = OUTPUT / "validation_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(manifest_path)


if __name__ == "__main__":
    main()
