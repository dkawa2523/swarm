"""Build the canonical portable-report artifact from reviewed benchmark outputs.

The generated artifact intentionally remains ``partial`` because the clean
COMSOL rerun was blocked by license error -10.  Its visible data-quality status
is ``share_with_caveats`` and all numerical COMSOL comparisons use the archived
14 July 2026 profiles.
"""

from __future__ import annotations

import base64
import csv
import io
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any

from PIL import Image


REPORT_DIR = Path(__file__).resolve().parents[1]
DERIVED_DIR = REPORT_DIR / "data" / "derived"
ARTIFACT_PATH = REPORT_DIR / "artifact.json"
MAX_ARTIFACT_BYTES = 2_900_000

QUANTITY_LABELS_JA = {
    "electron_density": "電子密度",
    "mean_electron_energy": "平均電子エネルギー",
    "electric_potential": "電位",
    "E_over_N": "換算電界 E/N",
    "electron_current_density": "電子伝導電流密度",
    "total_current_density": "電子＋Ar⁺伝導電流密度",
    "excitation_source": "直接励起源 eir2",
    "ionization_source": "直接電離源 eir4",
}


def _read_csv(name: str) -> list[dict[str, str]]:
    with (DERIVED_DIR / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_csv_optional(name: str) -> list[dict[str, str]]:
    path = DERIVED_DIR / name
    if not path.exists():
        return []
    return _read_csv(name)


def _float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _number(row: dict[str, Any], *keys: str, default: float | None = None) -> float:
    """Return the first finite-looking numeric field from a derived row."""

    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return float(value)
    if default is None:
        raise KeyError(f"none of the fields are present: {keys}")
    return default


def _fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}g}"


def _pct(value: float, digits: int = 3) -> str:
    return f"{100.0 * value:.{digits}g}%"


def _typed_rows(name: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw in _read_csv_optional(name):
        row: dict[str, Any] = dict(raw)
        for key, value in raw.items():
            if value == "":
                continue
            if value.lower() in {"true", "false"}:
                row[key] = value.lower() == "true"
                continue
            try:
                numeric = float(value)
                row[key] = numeric if math.isfinite(numeric) else None
            except ValueError:
                pass
        rows.append(row)
    return rows


def _find_row(
    rows: list[dict[str, Any]],
    **conditions: Any,
) -> dict[str, Any] | None:
    for row in rows:
        if all(row.get(key) == value for key, value in conditions.items()):
            return row
    return None


def _read_json_optional(name: str) -> dict[str, Any]:
    path = DERIVED_DIR / name
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _comparison_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw in _read_csv("comparison_metrics.csv"):
        row: dict[str, Any] = dict(raw)
        for key, value in raw.items():
            if key in {"evidence_status", "quantity", "quantity_label", "unit"}:
                continue
            if value != "":
                try:
                    numeric = float(value)
                    row[key] = numeric if math.isfinite(numeric) else None
                except ValueError:
                    pass
        row["common_grid_rows"] = int(raw["common_grid_rows"])
        row["quantity_label_ja"] = QUANTITY_LABELS_JA[raw["quantity"]]
        row["data_quality_status"] = "share_with_caveats"
        rows.append(row)
    return rows


def _runtime_rows() -> list[dict[str, Any]]:
    labels = {
        ("solve_only", "external_swarm"): "run stage・外部Swarm",
        ("solve_only", "builtin_reference"): "run stage・組込み参照",
        ("end_to_end", "external_swarm"): "全工程・外部Swarm",
        ("end_to_end", "builtin_reference"): "全工程・組込み参照",
    }
    rows: list[dict[str, Any]] = []
    for raw in _read_csv("runtime_summary.csv"):
        key = (raw["scope"], raw["path"])
        if key not in labels:
            continue
        rows.append(
            {
                "evidence_status": raw["evidence_status"],
                "data_quality_status": "share_with_caveats",
                "scope": raw["scope"],
                "path": raw["path"],
                "comparison_label_ja": labels[key],
                "repetitions": int(raw["repetitions"]),
                "median_seconds": _float(raw, "median_seconds"),
                "min_seconds": _float(raw, "min_seconds"),
                "max_seconds": _float(raw, "max_seconds"),
            }
        )
    return rows


def _clamp_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw in _read_csv("clamp_metrics.csv"):
        rows.append(
            {
                "evidence_status": raw["evidence_status"],
                "data_quality_status": "share_with_caveats",
                "metric": raw["metric"],
                "predicate": raw["predicate"],
                "threshold": _float(raw, "threshold"),
                "threshold_unit": raw["threshold_unit"],
                "spatial_length_m": _float(raw, "spatial_length_m"),
                "spatial_fraction": _float(raw, "spatial_fraction"),
                "node_fraction": _float(raw, "node_fraction"),
            }
        )
    return rows


def _figure_html(
    filename: str,
    alt: str,
    caption: str,
    *,
    max_height_px: int = 460,
) -> str:
    path = REPORT_DIR / "figures" / filename
    if not path.exists():
        return (
            '<aside style="padding:1rem;border:1px solid #d8dee4;border-radius:.5rem">'
            f"図は解析再実行後に生成される: <code>{filename}</code></aside>"
        )
    # Keep the publication PNG/PDF files untouched.  The portable HTML embeds
    # a bounded WebP preview so the verified payload remains below the
    # delivery tool's 3 MB limit.
    with Image.open(path) as source:
        preview = source.convert("RGB")
        preview.thumbnail((1600, 1200), Image.Resampling.LANCZOS)
        buffer = io.BytesIO()
        preview.save(buffer, format="WEBP", quality=88, method=6)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return (
        '<figure style="margin:0;text-align:center">'
        f'<img src="data:image/webp;base64,{encoded}" alt="{alt}" '
        'style="display:block;max-width:100%;width:auto;height:auto;'
        f'max-height:{max_height_px}px;object-fit:contain;margin:0 auto">'
        f'<figcaption style="margin:.5rem auto 0;max-width:72rem;color:#666">{caption}</figcaption>'
        "</figure>"
    )


def _source(
    source_id: str,
    label: str,
    *,
    path: str | None = None,
    href: str | None = None,
    query: dict[str, Any] | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {"id": source_id, "label": label}
    if path is not None:
        result["path"] = path
    if href is not None:
        result["href"] = href
    if query is not None:
        result["query"] = query
    return result


def build_artifact() -> dict[str, Any]:
    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    comparison = _comparison_rows()
    runtime = _runtime_rows()
    clamp = _clamp_rows()
    regional = _typed_rows("regional_error_attribution.csv")
    current_uniformity = _typed_rows("current_uniformity.csv")
    current_rsd_sensitivity = _typed_rows("current_rsd_sensitivity.csv")
    regime_metrics = _typed_rows("regime_metrics.csv")
    regime_impact = _typed_rows("regime_impact.csv")
    activation_audit = _typed_rows("historical_activation_audit.csv")
    fluid_applicability = _typed_rows("fluid_applicability_diagnostics.csv")
    operating_points = _typed_rows("operating_point_diagnostics.csv")
    townsend_representation = _typed_rows(
        "townsend_source_representation_audit.csv"
    )
    partial_energy = _typed_rows("partial_electron_energy_audit.csv")
    tail_audit = _typed_rows("tail_convergence_audit.csv")
    clamp_impact = _typed_rows("clamp_impact.csv")
    analysis_summary = _read_json_optional("analysis_summary.json")
    clamp_context = analysis_summary.get("clamp_context", {})
    formal_swarm = (
        analysis_summary.get("condition_audit", {}).get("formal_swarm", {})
    )
    by_quantity = {row["quantity"]: row for row in comparison}
    low_field = next(row for row in clamp if row["metric"] == "E_over_N_below_table")
    mean_energy_floor = next(
        row
        for row in clamp
        if row["metric"] == "mean_energy_in_constant_boundary_policy_region"
    )

    def metric(quantity: str, key: str) -> float:
        return float(by_quantity[quantity][key])

    potential_l2 = metric("electric_potential", "relative_L2_spatial_weighted")
    en_l2 = metric("E_over_N", "relative_L2_spatial_weighted")
    density_l2 = metric("electron_density", "relative_L2_spatial_weighted")
    energy_l2 = metric(
        "mean_electron_energy", "relative_L2_spatial_weighted"
    )
    excitation_l2 = metric("excitation_source", "relative_L2_spatial_weighted")
    ionization_l2 = metric("ionization_source", "relative_L2_spatial_weighted")
    density_integral_ratio = metric(
        "electron_density", "integral_ratio_external_over_reference"
    )
    energy_integral_ratio = metric(
        "mean_electron_energy", "integral_ratio_external_over_reference"
    )
    ionization_integral_ratio = metric(
        "ionization_source", "integral_ratio_external_over_reference"
    )

    external_run_stage = _find_row(
        runtime, scope="solve_only", path="external_swarm"
    )
    builtin_run_stage = _find_row(
        runtime, scope="solve_only", path="builtin_reference"
    )
    external_end_to_end = _find_row(
        runtime, scope="end_to_end", path="external_swarm"
    )
    builtin_end_to_end = _find_row(
        runtime, scope="end_to_end", path="builtin_reference"
    )
    run_stage_ratio = (
        _number(builtin_run_stage, "median_seconds")
        / _number(external_run_stage, "median_seconds")
        if external_run_stage and builtin_run_stage
        else float("nan")
    )
    end_to_end_ratio = (
        _number(builtin_end_to_end, "median_seconds")
        / _number(external_end_to_end, "median_seconds")
        if external_end_to_end and builtin_end_to_end
        else float("nan")
    )

    external_high_field = _find_row(
        regime_metrics,
        path="external_swarm",
        regime="E_over_N_above_typical_drift_diffusion_guidance",
    )
    high_field_length_fraction = (
        _number(external_high_field, "spatial_fraction")
        if external_high_field
        else float("nan")
    )
    high_field_threshold = (
        _number(external_high_field, "threshold_Td")
        if external_high_field
        else 500.0
    )
    high_field_ionization = _find_row(
        regime_impact,
        path="external_swarm",
        quantity="ionization_source",
    )
    high_field_ionization_share = (
        _number(high_field_ionization, "absolute_integral_share")
        if high_field_ionization
        else float("nan")
    )

    external_central_current = _find_row(
        current_uniformity,
        path="external_swarm",
        mesh_elements=200.0,
        domain="central_80_percent",
    )
    reference_central_current = _find_row(
        current_uniformity,
        path="builtin_reference",
        mesh_elements=200.0,
        domain="central_80_percent",
    )
    external_full_current = _find_row(
        current_uniformity,
        path="external_swarm",
        mesh_elements=200.0,
        domain="full_domain",
    )

    mean_energy_audit = _find_row(
        activation_audit, quantity="mean_energy_from_E_over_N"
    )
    mobility_audit = _find_row(
        activation_audit, quantity="reduced_mobility_from_mean_energy"
    )
    excitation_audit = _find_row(
        activation_audit,
        quantity="direct_excitation_townsend_from_mean_energy",
    )
    ionization_audit = _find_row(
        activation_audit,
        quantity="direct_ionization_townsend_from_mean_energy",
    )
    external_operating_point = _find_row(operating_points, path="external_swarm")
    reference_operating_point = _find_row(
        operating_points, path="builtin_reference"
    )
    excitation_representation = _find_row(
        townsend_representation, process_type="excitation"
    )
    ionization_representation = _find_row(
        townsend_representation, process_type="ionization"
    )
    external_partial_energy = _find_row(partial_energy, path="external_swarm")
    reference_partial_energy = _find_row(
        partial_energy, path="builtin_reference"
    )
    tail_case_count = len(tail_audit)
    tail_grid_pass_count = sum(
        bool(row.get("formal_tail_grid_gate_pass")) for row in tail_audit
    )
    solver_pass_count = sum(
        bool(row.get("formal_solver_convergence_gate_pass"))
        for row in tail_audit
    )
    formal_case_pass_count = sum(
        bool(row.get("formal_swarm_case_gate_pass")) for row in tail_audit
    )
    max_tail_grid_eV = max(
        (_number(row, "grid_max_eV") for row in tail_audit),
        default=float("nan"),
    )
    tail_solver_failures = [
        row
        for row in tail_audit
        if not bool(row.get("formal_solver_convergence_gate_pass"))
    ]
    external_fluid = _find_row(fluid_applicability, path="external_swarm")
    reference_fluid = _find_row(fluid_applicability, path="builtin_reference")
    floor_current_impact = _find_row(
        clamp_impact,
        region="mean_energy_in_constant_boundary_policy_region",
        quantity="total_current_density",
    )

    def row_pct(row: dict[str, Any] | None, key: str) -> str:
        return _pct(_number(row, key)) if row else "未生成"

    def row_num(row: dict[str, Any] | None, key: str, digits: int = 3) -> str:
        return _fmt(_number(row, key), digits) if row else "未生成"

    def row_int(row: dict[str, Any] | None, key: str) -> str:
        return str(int(_number(row, key))) if row else "未生成"

    headline_metrics = [
        {
            "metric_id": "potential_weighted_l2",
            "value": potential_l2,
            "data_quality_status": "share_with_caveats",
        },
        {
            "metric_id": "en_weighted_l2",
            "value": en_l2,
            "data_quality_status": "share_with_caveats",
        },
        {
            "metric_id": "density_integral_ratio",
            "value": density_integral_ratio,
            "data_quality_status": "share_with_caveats",
        },
        {
            "metric_id": "high_field_length_fraction",
            "value": high_field_length_fraction,
            "data_quality_status": "share_with_caveats",
        },
    ]

    interface_rows = [
        {
            "quantity": "平均電子エネルギー",
            "argument": "E/N",
            "comsol_target": "pes1.enrgXdata/enrgYdata",
            "unit": "eV",
            "written_property": "7 X–Y pairの1組（X/Yの2 property配列）",
            "active_in_LEA": "非活性：平均エネルギーPDEを解く",
            "active_in_LFA": "活性候補",
            "spatial_contract": "LFAでのみ E/N→平均エネルギーを照合",
        },
        {
            "quantity": "換算電子移動度",
            "argument": "平均エネルギー",
            "comsol_target": "pes1.muNXdata/muNYdata",
            "unit": "1/(V m s)",
            "written_property": "書込み・readback対象",
            "active_in_LEA": "活性候補",
            "active_in_LFA": "活性候補",
            "spatial_contract": "平均エネルギー→μNを照合",
        },
        {
            "quantity": "縦方向換算拡散",
            "argument": "平均エネルギー",
            "comsol_target": "pes1.deNXdata/deNYdata",
            "unit": "1/(m s)",
            "written_property": "書込み・readback対象",
            "active_in_LEA": "活性候補",
            "active_in_LFA": "活性候補",
            "spatial_contract": "平均エネルギー→DLNを照合",
        },
        {
            "quantity": "換算エネルギー移動度",
            "argument": "平均エネルギー",
            "comsol_target": "pes1.mueNXdata/mueNYdata",
            "unit": "1/(V m s)",
            "written_property": "書込み・readback対象",
            "active_in_LEA": "活性候補：平均エネルギーPDEで使用",
            "active_in_LFA": "非活性：平均エネルギーPDEなし",
            "spatial_contract": "安定した空間変数未同定、property照合のみ",
        },
        {
            "quantity": "換算エネルギー拡散",
            "argument": "平均エネルギー",
            "comsol_target": "pes1.denNXdata/denNYdata",
            "unit": "1/(m s)",
            "written_property": "書込み・readback対象",
            "active_in_LEA": "活性候補：平均エネルギーPDEで使用",
            "active_in_LFA": "非活性：平均エネルギーPDEなし",
            "spatial_contract": "安定した空間変数未同定、property照合のみ",
        },
        {
            "quantity": "直接励起 eir2 換算Townsend",
            "argument": "平均エネルギー",
            "comsol_target": "eir2.xtownratedata/ytownratedata",
            "unit": "m²",
            "written_property": "書込み・readback対象",
            "active_in_LEA": "活性候補",
            "active_in_LFA": "活性候補",
            "spatial_contract": "平均エネルギー→eir2 α/Nを照合",
        },
        {
            "quantity": "直接電離 eir4 換算Townsend",
            "argument": "平均エネルギー",
            "comsol_target": "eir4.xtownratedata/ytownratedata",
            "unit": "m²",
            "written_property": "書込み・readback対象",
            "active_in_LEA": "活性候補",
            "active_in_LFA": "活性候補",
            "spatial_contract": "平均エネルギー→eir4 α/Nを照合",
        },
    ]
    model_contrast_rows = [
        {
            "row_order": 1,
            "comparison_item": "root電子エネルギー定式化",
            "external_path": "LocalEnergyApproximationE（LEA）",
            "builtin_path": "LocalEnergyApproximationE（LEA）",
            "causal_limit": "両方ともLFAではない",
        },
        {
            "row_order": 2,
            "comparison_item": "電子自由度とclosure",
            "external_path": "平均エネルギーPDE＋外部輸送・Townsend表",
            "builtin_path": "保存済みsol3がEnとplas.F0を含む組込みBoltzmann連成",
            "causal_limit": "同じLEA root propertyでも電子closureは同値でない",
        },
        {
            "row_order": 3,
            "comparison_item": "外部表",
            "external_path": "7 lookup量＝7 X–Y pair＝14 property配列。LEA候補activeは6量",
            "builtin_path": "組込みBoltzmann/EEDFから係数を評価",
            "causal_limit": "written property数をactive closure数と解釈しない",
        },
        {
            "row_order": 4,
            "comparison_item": "反応・電流の表示量",
            "external_path": "eir2直接励起、eir4直接電離、電子＋Ar⁺伝導電流",
            "builtin_path": "同じexport式",
            "causal_limit": "全反応源・変位電流を含む総電流ではない",
        },
        {
            "row_order": 5,
            "comparison_item": "mesh",
            "external_path": "履歴200要素、400は旧summaryのみ",
            "builtin_path": "履歴profileは200要素",
            "causal_limit": "400要素rawがなくgrid convergence未成立",
        },
        {
            "row_order": 6,
            "comparison_item": "要素次数",
            "external_path": "vendor基礎モデルはquadratic",
            "builtin_path": "履歴logにlinear／quadraticの双方",
            "causal_limit": "field別割当の同一性をarchiveで立証できない",
        },
        {
            "row_order": 7,
            "comparison_item": "気体条件",
            "external_path": "旧表300 K・13.3 Pa、profile 293.15 K・13.3322 Pa",
            "builtin_path": "293.15 K・13.3322 Pa",
            "causal_limit": "履歴差に条件差を含む",
        },
        {
            "row_order": 8,
            "comparison_item": "run状態",
            "external_path": "table適用後continuation",
            "builtin_path": "保存済み最終solver sol3",
            "causal_limit": "初期値・保存state効果を分離していない",
        },
    ]

    comparison_sql = (
        "SELECT quantity, CASE quantity "
        "WHEN 'electron_density' THEN '電子密度' "
        "WHEN 'mean_electron_energy' THEN '平均電子エネルギー' "
        "WHEN 'electric_potential' THEN '電位' "
        "WHEN 'E_over_N' THEN '換算電界 E/N' "
        "WHEN 'electron_current_density' THEN '電子伝導電流密度' "
        "WHEN 'total_current_density' THEN '電子＋Ar+伝導電流密度' "
        "WHEN 'excitation_source' THEN '直接励起源 eir2' "
        "WHEN 'ionization_source' THEN '直接電離源 eir4' END AS quantity_label_ja, "
        "relative_L2_spatial_weighted, integral_ratio_external_over_reference, "
        "shape_correlation_spatial_weighted, "
        "relative_L2_trapezoidal_of_squared_samples_sensitivity, "
        "shape_correlation_trapezoidal_nodal_sensitivity, "
        "legacy_L2_unweighted_discrete, "
        "legacy_correlation_unweighted_discrete "
        "FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/comparison_metrics.csv', "
        "header = true)"
    )
    runtime_sql = (
        "SELECT scope, path, repetitions, median_seconds, min_seconds, max_seconds "
        "FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/runtime_summary.csv', "
        "header = true) WHERE scope IN ('solve_only', 'end_to_end')"
    )
    clamp_sql = (
        "SELECT metric, threshold, threshold_unit, spatial_length_m, spatial_fraction, node_fraction "
        "FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/clamp_metrics.csv', "
        "header = true) WHERE metric = 'E_over_N_below_table'"
    )
    interface_sql = (
        "SELECT * FROM (VALUES "
        "('平均電子エネルギー','E/N','pes1.enrgXdata/enrgYdata','eV','written','inactive','active'),"
        "('換算電子移動度','平均エネルギー','pes1.muNXdata/muNYdata','1/(V m s)','written','active','active'),"
        "('縦方向換算拡散','平均エネルギー','pes1.deNXdata/deNYdata','1/(m s)','written','active','active'),"
        "('換算エネルギー移動度','平均エネルギー','pes1.mueNXdata/mueNYdata','1/(V m s)','written','active','inactive'),"
        "('換算エネルギー拡散','平均エネルギー','pes1.denNXdata/denNYdata','1/(m s)','written','active','inactive'),"
        "('直接励起eir2換算Townsend','平均エネルギー','eir2.xtownratedata/ytownratedata','m²','written','active','active'),"
        "('直接電離eir4換算Townsend','平均エネルギー','eir4.xtownratedata/ytownratedata','m²','written','active','active')"
        ") AS closure(quantity,argument,comsol_target,unit,written_property,active_in_LEA,active_in_LFA)"
    )
    model_contrast_sql = (
        "SELECT * FROM (VALUES "
        "(1,'root定式化','LEA','LEA','両方ともLFAではない'),"
        "(2,'電子自由度','平均エネルギーPDE＋外部表','sol3がEn＋plas.F0','closureは同値でない'),"
        "(3,'外部表','7 lookup量＝7 X-Y pair＝14 property配列／LEA active候補6量','組込みBoltzmann','writtenとactiveを分ける'),"
        "(4,'表示量','eir2、eir4、電子＋Ar+伝導電流','同じexport式','全反応源・総電流ではない'),"
        "(5,'mesh','履歴200、400は旧summaryのみ','履歴200','grid convergence未成立'),"
        "(6,'要素次数','vendor基礎モデルはquadratic','履歴logにlinear／quadratic','同一性未立証'),"
        "(7,'気体条件','旧表300 K・13.3 Pa、profile 293.15 K・13.3322 Pa','293.15 K・13.3322 Pa','条件差を含む'),"
        "(8,'run状態','table適用後continuation','保存済みsol3','state効果未分離')"
        ") AS model_contrast(row_order,comparison_item,external_path,builtin_path,causal_limit)"
    )

    sources = [
        _source(
            "comparison_metrics_source",
            "空間重み付き比較指標（2026-07-14履歴profileの再解析）",
            path="reports/comsol_swarm_benchmark_2026/data/derived/comparison_metrics.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": comparison_sql,
                "description": "履歴profileを共通非一様格子へ整列し、区分線形再構成の積を区間ごとに厳密積分して8物理量を比較した派生表。",
                "executed_at": "2026-07-31T07:57:52+09:00",
                "metric_definitions": [
                    "relative_L2_spatial_weighted = 区分線形再構成に対する sqrt(∫(external-reference)^2 dx / ∫reference^2 dx)",
                    "integral_ratio_external_over_reference = ∫external dx / ∫reference dx",
                    "shape_correlation_spatial_weighted = 区分線形再構成の積を厳密積分した空間相関",
                    "trapezoidal_of_squared_samples列は数値積分規約の感度確認",
                ],
                "tables_used": [
                    "reports/comsol_swarm_benchmark_2026/data/archive_2026-07-14/positive_column_results.csv",
                    "reports/comsol_swarm_benchmark_2026/data/archive_2026-07-14/builtin_200V_fresh.csv",
                ],
            },
        ),
        _source(
            "runtime_summary_source",
            "単発run-stage／end-to-end時間集計（2026-07-14）",
            path="reports/comsol_swarm_benchmark_2026/data/derived/runtime_summary.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": runtime_sql,
                "description": "COMSOL run呼出しを囲むrun-stage timerとend-to-end timerを抽出した。内部solverだけの時間ではなく、各経路n=1。",
                "executed_at": "2026-07-31T07:57:52+09:00",
                "filters": ["scope in {solve_only, end_to_end}", "one historical observation per path"],
                "metric_definitions": ["median/min/max are identical when repetitions = 1"],
                "tables_used": [
                    "reports/comsol_swarm_benchmark_2026/data/provisional_runtime.csv"
                ],
            },
        ),
        _source(
            "regional_error_source",
            "区分線形厳密積分による領域別誤差寄与",
            path="reports/comsol_swarm_benchmark_2026/data/derived/regional_error_attribution.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": (
                    "SELECT quantity, region, squared_error_share, "
                    "region_normalized_relative_L2 FROM read_csv_auto("
                    "'reports/comsol_swarm_benchmark_2026/data/derived/"
                    "regional_error_attribution.csv', header=true)"
                ),
                "description": "陰極側10%、中央80%、陽極側10%の幾何学的窓へ二乗誤差を分割。plasma bulk／sheath境界の同定ではない。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "current_rsd_source",
            "伝導電流RSDと中央保持率感度",
            path="reports/comsol_swarm_benchmark_2026/data/derived/current_rsd_sensitivity.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": (
                    "SELECT path, central_retained_fraction, "
                    "piecewise_linear_exact_relative_standard_deviation "
                    "FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/"
                    "data/derived/current_rsd_sensitivity.csv', header=true)"
                ),
                "description": "電子＋Ar+伝導電流の区分線形空間population RSDを中央保持率に対して評価。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "current_uniformity_source",
            "伝導電流の全領域／中央80% RSD",
            path="reports/comsol_swarm_benchmark_2026/data/derived/current_uniformity.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": (
                    "SELECT path, mesh_elements, domain, relative_standard_deviation "
                    "FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/"
                    "data/derived/current_uniformity.csv', header=true)"
                ),
                "description": "raw profileでは区分線形厳密RSD。旧400要素summaryは定義不明の補足値として分離。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "regime_metrics_source",
            "COMSOL drift-diffusion推奨域に対する空間support",
            path="reports/comsol_swarm_benchmark_2026/data/derived/regime_metrics.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": (
                    "SELECT * FROM read_csv_auto('reports/"
                    "comsol_swarm_benchmark_2026/data/derived/regime_metrics.csv',"
                    " header=true)"
                ),
                "description": "E/NがCOMSOLの典型的なdrift-diffusion目安を超える空間長。硬い妥当性境界ではない。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "regime_impact_source",
            "高E/N領域の物理量別絶対積分寄与",
            path="reports/comsol_swarm_benchmark_2026/data/derived/regime_impact.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": (
                    "SELECT path, quantity, threshold_Td, absolute_integral_share "
                    "FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/"
                    "data/derived/regime_impact.csv', header=true)"
                ),
                "description": "field threshold交点と値の零交差で区間分割した絶対積分寄与。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "activation_audit_source",
            "履歴exportと意図したbundle表の数値整合監査",
            path="reports/comsol_swarm_benchmark_2026/data/derived/historical_activation_audit.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": (
                    "SELECT quantity, formulation_activity, relative_error_median, "
                    "relative_error_p95, relative_error_max, "
                    "points_above_1_percent, status FROM read_csv_auto("
                    "'reports/comsol_swarm_benchmark_2026/data/derived/"
                    "historical_activation_audit.csv', header=true)"
                ),
                "description": "意図した2026-07-14 bundleとarchived exportの数値整合のみ。stale-class riskのためactivation証明ではない。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "fluid_applicability_source",
            "weak-ionization・net-flux speed診断",
            path="reports/comsol_swarm_benchmark_2026/data/derived/fluid_applicability_diagnostics.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": "SELECT * FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/fluid_applicability_diagnostics.csv', header=true)",
                "description": "n_e/Nと|J_e|/[e n_e sqrt(2 e mean_energy/m_e)]のoperational diagnostic。hard validity testではない。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "operating_point_source",
            "ballast回路を含むoperating-point診断",
            path="reports/comsol_swarm_benchmark_2026/data/derived/operating_point_diagnostics.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": "SELECT * FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/operating_point_diagnostics.csv', header=true)",
                "description": "source voltage、gap endpoint voltage、ballast drop/current、内部potential maximumをprofileから導出。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "townsend_representation_source",
            "Townsend-flux源とrate-form反実仮想の監査",
            path="reports/comsol_swarm_benchmark_2026/data/derived/townsend_source_representation_audit.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": "SELECT * FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/townsend_source_representation_audit.csv', header=true)",
                "description": "archived Townsend sourceを再構成し、N n_e kとのpostprocess counterfactualを比較。matched rerunではない。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "partial_energy_source",
            "電子field powerとeir2/eir4 threshold lossの部分監査",
            path="reports/comsol_swarm_benchmark_2026/data/derived/partial_electron_energy_audit.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": "SELECT * FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/partial_electron_energy_audit.csv', header=true)",
                "description": "J_e Eとeir2/eir4 threshold-weighted lossのみ。完全なelectron-energy balanceではない。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "tail_audit_source",
            "schema-v2 Swarmのtail・solver convergence監査",
            path="reports/comsol_swarm_benchmark_2026/data/derived/tail_convergence_audit.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": "SELECT * FROM read_csv_auto('reports/comsol_swarm_benchmark_2026/data/derived/tail_convergence_audit.csv', header=true)",
                "description": f"{tail_case_count} E/N点のtail probability、edge/peak、20 keV grid-limit、solver convergenceを別gateで評価。",
                "executed_at": generated_at,
            },
        ),
        _source(
            "clamp_metrics_source",
            "E/N lookup範囲と低電界clamp",
            path="reports/comsol_swarm_benchmark_2026/data/derived/clamp_metrics.csv",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": clamp_sql,
                "description": "2026-07-31 schema-v2表の下限と2026-07-14 COMSOL空間profileを重ねた混合vintage評価。",
                "executed_at": "2026-07-31T07:57:52+09:00",
                "metric_definitions": [
                    f"spatial_fraction = E/Nが{_fmt(_number(low_field, 'threshold'))} "
                    "Td未満となる補間区間長 / 全計算区間長"
                ],
                "tables_used": [
                    "reports/comsol_swarm_benchmark_2026/data/formal_swarm_bundle_solver_qa_20keV/mean_energy_vs_en.csv",
                    "reports/comsol_swarm_benchmark_2026/data/archive_2026-07-14/positive_column_results.csv",
                ],
            },
        ),
        _source(
            "interface_contract_source",
            "COMSOL Java生成、root定式化、7 lookup量／14 property配列の契約",
            path="swarm_workflow/comsol_java.py",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": interface_sql,
                "description": "7 lookup量は7 X–Y pair、すなわち14 COMSOL property配列。LEA候補active 6量／LFA候補active 5量を分離したformulation-aware契約。",
                "executed_at": "2026-07-31T07:57:52+09:00",
                "tables_used": [
                    "swarm_workflow/comsol_mapping.py",
                    "swarm_workflow/comsol_java.py",
                    "swarm_workflow/comsol_verify.py",
                    "Model/maps/positive_column_external.yaml",
                ],
            },
        ),
        _source(
            "model_contrast_source",
            "外部Swarm経路と組込み参照経路のモデル対照",
            path="reports/comsol_swarm_benchmark_2026/manuscript.md",
            query={
                "language": "sql",
                "engine": "DuckDB",
                "sql": model_contrast_sql,
                "description": "入力MPH、closure、mesh、要素次数、気体条件、solve状態の未統制差を表形式にした。",
                "executed_at": generated_at,
                "tables_used": [
                    "reports/comsol_swarm_benchmark_2026/manuscript.md",
                    "reports/comsol_swarm_benchmark_2026/data/archive_2026-07-14/historical_swarm_source_config.yaml",
                ],
            },
        ),
        _source(
            "analysis_summary_source",
            "解析・QA要約",
            path="reports/comsol_swarm_benchmark_2026/data/derived/analysis_summary.json",
            query={
                "language": "json",
                "description": "データvintage、指標定義、clamp、QA、formal gateを記録した解析要約。",
                "executed_at": "2026-07-31T07:57:52+09:00",
            },
        ),
        _source(
            "figure_metadata_source",
            "図契約、出典、hash",
            path="reports/comsol_swarm_benchmark_2026/data/derived/figure_metadata.json",
            query={
                "language": "json",
                "description": "各図のanalytical question、takeaway、fields、palette、source hashを保存。",
                "executed_at": "2026-07-31T07:57:52+09:00",
            },
        ),
        _source(
            "formal_config_source",
            "2026-07-31 solver-QA Swarm候補設定（正式未採用）",
            path="reports/comsol_swarm_benchmark_2026/data/formal_swarm_bundle_solver_qa_20keV/executed_argon_comsol_benchmark_2026.yaml",
            query={
                "language": "yaml",
                "description": (
                    "実行原本とbyte-identical。schema_version 2、two_term、"
                    f"{_fmt(float(formal_swarm.get('gas_temperature_K', float('nan'))))} K、"
                    f"{_fmt(float(formal_swarm.get('pressure_Pa', float('nan'))))} Pa、"
                    f"{_fmt(float(clamp_context.get('table_E_over_N_min_Td', float('nan'))))}–"
                    f"{_fmt(float(clamp_context.get('table_E_over_N_max_Td', float('nan'))))} Tdを定義。"
                ),
                "executed_at": "2026-07-31T07:57:52+09:00",
            },
        ),
        _source(
            "solver_qa_database_source",
            "tail/grid合格・solver残差未完のSwarm database",
            path="reports/comsol_swarm_benchmark_2026/data/formal_swarm_database_solver_qa_20keV.sqlite",
            query={
                "language": "sqlite",
                "description": (
                    f"schema-v2 {tail_case_count}点。tail/grid "
                    f"{tail_grid_pass_count}/{tail_case_count}、solver residual "
                    f"{solver_pass_count}/{tail_case_count}のため正式bundleではない。"
                ),
                "executed_at": generated_at,
            },
        ),
        _source(
            "portable_config_source",
            "pathを報告directory用に変更したschema-v2再現設定",
            path="reports/comsol_swarm_benchmark_2026/repro/argon_comsol_benchmark_2026.yaml",
            query={
                "language": "yaml",
                "description": "実行原本と意味上の条件は同じだが、cross section／output相対pathを変更した再現用複製（SHA-256 f83918d9…29cd）。",
                "executed_at": "2026-07-31T07:57:52+09:00",
            },
        ),
        _source(
            "historical_swarm_condition_source",
            "2026-07-14外部表のSwarm条件証拠",
            path="reports/comsol_swarm_benchmark_2026/data/archive_2026-07-14/historical_swarm_source_config.yaml",
            query={
                "language": "yaml",
                "description": "旧外部表を生成した300 K、13.3 Pa設定。原本とbyte-identical（SHA-256 c1a6d641…9f1）。",
                "executed_at": "2026-07-14T00:00:00+09:00",
            },
        ),
        _source(
            "external_runner_source",
            "外部COMSOL独立反復runner",
            path="reports/comsol_swarm_benchmark_2026/repro/run_external_benchmark.py",
            query={
                "language": "python",
                "description": "run-comsolを独立processで反復し、293.15 K、13.3322 Pa、200/400要素のproperty readback、stage別log、profile、provenance、hash、時間統計をcopy-only archiveへ保存。",
                "executed_at": "2026-07-31T09:00:00+09:00",
            },
        ),
        _source(
            "license_attempt_source",
            "clean再実行のCOMSOL license失敗ログ",
            path="reports/comsol_swarm_benchmark_2026/data/clean_run_attempts/external_run_01_user/apply_stdout.txt",
            query={
                "language": "text",
                "description": "新しいclassは生成されたが、apply開始時にlicense error -10（Product has expired）で停止。",
                "executed_at": "2026-07-31T07:43:10+09:00",
            },
        ),
        _source(
            "compile_archive_source",
            "2026-07-14 compiler誤判定の履歴ログ",
            path="reports/comsol_swarm_benchmark_2026/data/archive_2026-07-14/compile_stdout.txt",
            query={
                "language": "text",
                "description": "return code 0に対して標準出力がCompilation failedを記録した履歴証拠。",
                "executed_at": "2026-07-14T01:34:48+09:00",
            },
        ),
        _source(
            "manuscript_source",
            "和文論文原稿とEnglish abstract",
            path="reports/comsol_swarm_benchmark_2026/manuscript.md",
            query={
                "language": "markdown",
                "description": "問題設定、方法、結果、制約、第三者adapter要件、一次資料を統合した論文原稿。",
                "executed_at": "2026-07-31T07:57:52+09:00",
            },
        ),
        _source(
            "mph_formulation_evidence_source",
            "MPH root formulationとsaved solver自由度の監査証拠",
            path="reports/comsol_swarm_benchmark_2026/SOURCE_INVENTORY.md",
            query={
                "language": "markdown",
                "description": (
                    "positive_column_1d.mph、履歴external-applied MPH、"
                    "Boltzmann MPHのMeanElectronEnergyModel readbackと、"
                    "saved sol3のEn・plas.F0 solve-forを記録。"
                ),
                "executed_at": generated_at,
                "tables_used": [
                    "Model/positive_column_1d.mph",
                    "Model/work/positive_column_external_applied.mph",
                    "Model/positive_column_1d_boltzmann.mph",
                    "Model/generated/boltzmann_study_solver_inspection_stdout.txt",
                ],
            },
        ),
        _source(
            "swarm_transport_source",
            "two_term transport定義とCOMSOL表変換実装",
            path="electron_swarm/solvers/kinetic.py",
            query={
                "language": "python",
                "description": "flux transport metadata、scalar-f0 moment由来のDL=DT、Townsend=k/|drift velocity|の実装根拠。",
                "executed_at": generated_at,
                "tables_used": [
                    "electron_swarm/solvers/kinetic.py",
                    "swarm_workflow/tables.py",
                ],
            },
        ),
        _source(
            "comsol_model_source",
            "COMSOL公式の連成1D DC glow-dischargeモデル",
            href="https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.positive_column_1d/positive_column_1d.html",
        ),
        _source(
            "comsol_interface_source",
            "COMSOL Drift Diffusion interface: LEA/LFA",
            href="https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.02.html",
        ),
        _source(
            "comsol_lookup_source",
            "COMSOL Drift Diffusion Model: lookup-table coefficients",
            href="https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.04.html",
        ),
        _source(
            "comsol_theory_source",
            "COMSOL Drift Diffusion Theory and application regime",
            href="https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.17.html",
        ),
        _source(
            "comsol_guide_source",
            "COMSOL Plasma Module User's Guide 6.4",
            href="https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/PlasmaModuleUsersGuide.pdf",
        ),
        _source(
            "bolsig_paper_source",
            "Hagelaar and Pitchford (2005), BOLSIG+",
            href="https://doi.org/10.1088/0963-0252/14/4/011",
        ),
        _source(
            "lxcat_paper_source",
            "Pancheshnyi et al. (2012), LXCat",
            href="https://doi.org/10.1016/j.chemphys.2011.04.020",
        ),
        _source(
            "dujko_flux_bulk_source",
            "Dujko and White (2008), flux/bulk transport",
            href="https://doi.org/10.1088/1742-6596/133/1/012005",
        ),
        _source(
            "dujko_nonlocal_source",
            "Dujko et al. (2008), nonhydrodynamic transport in cathode fall",
            href="https://doi.org/10.1088/0022-3727/41/24/245205",
        ),
    ]

    manifest_sources = [
        {key: value for key, value in source.items() if key != "query"}
        for source in sources
    ]

    cards = [
        {
            "id": "potential_l2_card",
            "description": "外部Swarmと組込み参照の電位profileに対する空間重み付き相対L2。",
            "dataset": "headline_metrics",
            "sourceId": "comparison_metrics_source",
            "filter": {"metric_id": "potential_weighted_l2"},
            "metrics": [{"label": "電位・重み付きL2", "field": "value", "format": "number"}],
        },
        {
            "id": "en_l2_card",
            "description": "外部Swarmと組込み参照のE/N profileに対する空間重み付き相対L2。",
            "dataset": "headline_metrics",
            "sourceId": "comparison_metrics_source",
            "filter": {"metric_id": "en_weighted_l2"},
            "metrics": [{"label": "E/N・重み付きL2", "field": "value", "format": "number"}],
        },
        {
            "id": "density_ratio_card",
            "description": "外部電子密度の空間積分を組込み参照の積分で除した比。",
            "dataset": "headline_metrics",
            "sourceId": "comparison_metrics_source",
            "filter": {"metric_id": "density_integral_ratio"},
            "metrics": [{"label": "密度積分比 外部/参照", "field": "value", "format": "number"}],
        },
        {
            "id": "high_field_fraction_card",
            "description": "E/NがCOMSOLの典型的drift-diffusion目安を上回る履歴外部profileの空間長割合。",
            "dataset": "headline_metrics",
            "sourceId": "regime_metrics_source",
            "filter": {"metric_id": "high_field_length_fraction"},
            "metrics": [{"label": "高E/N域の長さ割合", "field": "value", "format": "percent"}],
        },
    ]

    charts = [
        {
            "id": "weighted_l2_chart",
            "title": "物理量別の区分線形・空間相対L2",
            "subtitle": "2026年7月14日履歴profile。union grid上の区分線形再構成を区間厳密積分",
            "type": "bar",
            "dataset": "comparison_metrics",
            "sourceId": "comparison_metrics_source",
            "encodings": {
                "x": {"field": "quantity_label_ja", "type": "nominal", "label": "物理量"},
                "y": {
                    "field": "relative_L2_spatial_weighted",
                    "type": "quantitative",
                    "label": "区分線形・空間相対L2",
                },
                "tooltip": [
                    {
                        "field": "shape_correlation_spatial_weighted",
                        "type": "quantitative",
                        "label": "空間重み付き相関",
                    },
                    {
                        "field": "integral_ratio_external_over_reference",
                        "type": "quantitative",
                        "label": "積分比 外部/参照",
                    },
                ],
            },
            "valueFormat": "number",
            "layout": "full",
            "maxRows": 8,
        },
        {
            "id": "integral_ratio_chart",
            "title": "物理量別の空間積分比",
            "subtitle": "外部Swarm / 組込み参照。1が積分量の一致を表す",
            "type": "bar",
            "dataset": "comparison_metrics",
            "sourceId": "comparison_metrics_source",
            "encodings": {
                "x": {"field": "quantity_label_ja", "type": "nominal", "label": "物理量"},
                "y": {
                    "field": "integral_ratio_external_over_reference",
                    "type": "quantitative",
                    "label": "空間積分比 外部/参照",
                },
                "tooltip": [
                    {
                        "field": "relative_L2_spatial_weighted",
                        "type": "quantitative",
                        "label": "空間重み付き相対L2",
                    }
                ],
            },
            "valueFormat": "number",
            "layout": "full",
            "maxRows": 8,
        },
        {
            "id": "runtime_chart",
            "title": "COMSOL run stageとend-to-endの単発時間",
            "subtitle": "run stageは内部solverだけの時間ではない。各経路n=1で分布を表さない",
            "type": "scatter",
            "dataset": "runtime_comparison",
            "sourceId": "runtime_summary_source",
            "encodings": {
                "x": {
                    "field": "comparison_label_ja",
                    "type": "nominal",
                    "label": "scope・経路",
                },
                "y": {
                    "field": "median_seconds",
                    "type": "quantitative",
                    "label": "経過時間",
                    "unit": "s",
                },
                "tooltip": [
                    {"field": "repetitions", "type": "quantitative", "label": "反復数"},
                    {"field": "min_seconds", "type": "quantitative", "label": "最小 s"},
                    {"field": "max_seconds", "type": "quantitative", "label": "最大 s"},
                ],
            },
            "valueFormat": "number",
            "unit": "s",
            "layout": "full",
            "maxRows": 4,
        },
    ]

    tables = [
        {
            "id": "weighted_metrics_table",
            "title": "区分線形・空間比較指標",
            "subtitle": "主列は区分線形再構成の厳密区間積分。台形-of-squaresは感度列",
            "dataset": "comparison_metrics",
            "sourceId": "comparison_metrics_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "relative_L2_spatial_weighted", "direction": "desc"},
            "columns": [
                {"field": "quantity_label_ja", "label": "物理量", "type": "text"},
                {"field": "relative_L2_spatial_weighted", "label": "PL厳密 相対L2", "format": "number"},
                {
                    "field": "relative_L2_trapezoidal_of_squared_samples_sensitivity",
                    "label": "台形-of-squares 感度",
                    "format": "number",
                },
                {
                    "field": "integral_ratio_external_over_reference",
                    "label": "積分比 外部/参照",
                    "format": "number",
                },
                {
                    "field": "shape_correlation_spatial_weighted",
                    "label": "重み付き相関",
                    "format": "number",
                },
            ],
        },
        {
            "id": "closure_contract_table",
            "title": "7 lookup量（7 X–Y pair／14 property配列）とactive set",
            "subtitle": "LEA候補activeは6量、LFA候補activeは5量。property存在≠functional activation",
            "dataset": "interface_contract",
            "sourceId": "interface_contract_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "quantity", "direction": "asc"},
            "columns": [
                {"field": "quantity", "label": "注入量", "type": "text"},
                {"field": "argument", "label": "表引数", "type": "text"},
                {"field": "comsol_target", "label": "COMSOL target", "type": "text"},
                {"field": "unit", "label": "単位", "type": "text"},
                {"field": "written_property", "label": "written/readback", "type": "text"},
                {"field": "active_in_LEA", "label": "LEA", "type": "text"},
                {"field": "active_in_LFA", "label": "LFA", "type": "text"},
                {"field": "spatial_contract", "label": "空間照合契約", "type": "text"},
            ],
        },
        {
            "id": "activation_audit_table",
            "title": "履歴export対 intended bundle の数値監査",
            "subtitle": "archived class provenance未解決のためactivation証明ではない",
            "dataset": "historical_activation_audit",
            "sourceId": "activation_audit_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "quantity_label", "direction": "asc"},
            "columns": [
                {"field": "quantity_label", "label": "照合量", "type": "text"},
                {"field": "formulation_activity", "label": "LEAでの位置づけ", "type": "text"},
                {"field": "relative_error_median", "label": "相対誤差 median", "format": "number"},
                {"field": "relative_error_p95", "label": "相対誤差 p95", "format": "number"},
                {"field": "relative_error_max", "label": "相対誤差 max", "format": "number"},
                {"field": "points_above_1_percent", "label": ">1% 点数", "format": "integer"},
                {"field": "status", "label": "監査結果", "type": "text"},
            ],
        },
        {
            "id": "operating_point_table",
            "title": "ballast回路を含む履歴operating point",
            "subtitle": "同じsource voltageでもgap voltageとcurrentが異なる",
            "dataset": "operating_point_diagnostics",
            "sourceId": "operating_point_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "path_label", "direction": "asc"},
            "columns": [
                {"field": "path_label", "label": "経路", "type": "text"},
                {"field": "source_voltage_V", "label": "source V", "format": "number"},
                {"field": "gap_endpoint_voltage_V", "label": "gap V", "format": "number"},
                {"field": "ballast_drop_V", "label": "ballast drop V", "format": "number"},
                {"field": "ballast_current_A", "label": "ballast current A", "format": "number"},
                {"field": "internal_max_potential_V", "label": "internal max V", "format": "number"},
            ],
        },
        {
            "id": "townsend_representation_table",
            "title": "Townsend源とrate-form反実仮想",
            "subtitle": "postprocess sensitivityのみ。matched COMSOL rerunではない",
            "dataset": "townsend_source_representation",
            "sourceId": "townsend_representation_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "channel_label", "direction": "asc"},
            "columns": [
                {"field": "channel_label", "label": "channel", "type": "text"},
                {"field": "townsend_reconstruction_relative_difference", "label": "Townsend再構成相対差", "format": "number"},
                {"field": "counterfactual_rate_over_observed_townsend", "label": "rate反実仮想/Townsend", "format": "number"},
                {"field": "interpretation", "label": "解釈", "type": "text"},
            ],
        },
        {
            "id": "partial_energy_table",
            "title": "電子energyの部分監査",
            "subtitle": "eir2/eir4 threshold lossのみ。完全収支ではない",
            "dataset": "partial_electron_energy_audit",
            "sourceId": "partial_energy_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "path_label", "direction": "asc"},
            "columns": [
                {"field": "path_label", "label": "経路", "type": "text"},
                {"field": "electron_field_power_signed_W_m2", "label": "∫JₑE dx W/m²", "format": "number"},
                {"field": "direct_eir2_eir4_loss_sum_W_m2", "label": "eir2+eir4 loss W/m²", "format": "number"},
                {"field": "direct_loss_to_field_power_ratio", "label": "部分loss/field power", "format": "number"},
                {"field": "audit_scope", "label": "監査範囲", "type": "text"},
            ],
        },
        {
            "id": "tail_failure_table",
            "title": "Swarm solver convergence未合格点",
            "subtitle": "tail/grid gateは全点合格だが、solver convergence gateは未完",
            "dataset": "tail_solver_failures",
            "sourceId": "tail_audit_source",
            "density": "dense",
            "layout": "full",
            "defaultSort": {"field": "E_over_N_Td", "direction": "asc"},
            "columns": [
                {"field": "E_over_N_Td", "label": "E/N Td", "format": "number"},
                {"field": "iterations", "label": "iterations", "format": "integer"},
                {"field": "residual_L1", "label": "residual L1", "format": "number"},
                {"field": "tail_pass", "label": "tail", "type": "boolean"},
                {"field": "edge_pass", "label": "edge/peak", "type": "boolean"},
                {"field": "formal_solver_convergence_gate_pass", "label": "solver gate", "type": "boolean"},
            ],
        },
        {
            "id": "model_contrast_table",
            "title": "外部Swarm経路と組込み参照経路のモデル対照",
            "subtitle": "観測差をclosure単独へ帰属させないための比較境界",
            "dataset": "model_contrast",
            "sourceId": "model_contrast_source",
            "density": "spacious",
            "layout": "full",
            "defaultSort": {"field": "row_order", "direction": "asc"},
            "columns": [
                {"field": "row_order", "label": "No.", "format": "integer"},
                {"field": "comparison_item", "label": "比較項目", "type": "text"},
                {"field": "external_path", "label": "外部Swarm経路", "type": "text"},
                {"field": "builtin_path", "label": "組込み参照経路", "type": "text"},
                {"field": "causal_limit", "label": "因果解釈の制約", "type": "text"},
            ],
        },
        {
            "id": "legacy_metrics_table",
            "title": "補足：旧離散指標",
            "subtitle": "追跡可能性のためのみ保存。本文結論には用いない",
            "dataset": "comparison_metrics",
            "sourceId": "comparison_metrics_source",
            "density": "dense",
            "layout": "full",
            "defaultSort": {"field": "legacy_L2_unweighted_discrete", "direction": "desc"},
            "columns": [
                {"field": "quantity_label_ja", "label": "物理量", "type": "text"},
                {
                    "field": "legacy_L2_unweighted_discrete",
                    "label": "旧・非重み付き離散L2",
                    "format": "number",
                },
                {
                    "field": "legacy_correlation_unweighted_discrete",
                    "label": "旧・非重み付き相関",
                    "format": "number",
                },
            ],
        },
    ]

    title = (
        "External Swarm Tables in a COMSOL 1D DC Glow Discharge: "
        "A Provisional Formulation-Aware Selected-Channel Hybrid Benchmark"
    )
    blocks = [
        {
            "id": "report_title",
            "type": "markdown",
            "body": (
                f"# {title}\n\n"
                "**和文技術報告・English abstract付き · data quality: `share_with_caveats` · "
                "数値証拠: 2026年7月14日履歴benchmark**"
            ),
        },
        {
            "id": "technical_summary",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## 技術要約\n\n"
                "**最重要の監査結果は、履歴外部MPHと組込み参照MPHのroot電子エネルギー定式化が"
                "`LocalEnergyApproximationE`（LEA）であり、LFAではないことである。** "
                "外部経路は7 lookup量を7 X–Y pair（14 COMSOL property配列）として書込むが、"
                "LEAで候補activeなのは平均エネルギー表を除く6量である。"
                "built-in保存済み`sol3`は平均エネルギー自由度`En`に加えてEEDF自由度`plas.F0`を含む。"
                "さらにCOMSOL側のstepwise/Penning、superelastic、wall、heavy species、Poisson、ballast回路を"
                "保持するため、実態は全electron closure置換ではなくselected-channel hybrid closureである。"
                "したがって、両経路は同じroot propertyでも方程式レベルで同値ではない。\n\n"
                "clean buildを行う修復後runnerは新しいclassを生成したものの、COMSOL apply開始時にlicense error -10 "
                "（`Product has expired`）で停止した。このため、本文のCOMSOL数値は2026年7月14日の履歴profileを"
                "区分線形再構成の厳密空間積分で再解析した暫定技術benchmarkであり、精度優越、実験一致、closure単独の因果差を"
                "主張しない。履歴closure表は300 K・13.3 Pa、COMSOL profileは293.15 K・13.3322 Paであり、"
                "再実行候補設定では一致させたものの実COMSOL runでの確認は未完了である。\n\n"
                f"履歴profileの区分線形・空間相対L2は、電位{_fmt(potential_l2)}、E/N {_fmt(en_l2)}、"
                f"電子密度{_fmt(density_l2)}、平均電子エネルギー{_fmt(energy_l2)}、"
                f"eir2直接励起源{_fmt(excitation_l2)}、eir4直接電離源{_fmt(ionization_l2)}である。"
                f"電子密度積分比は{_fmt(density_integral_ratio)}（外部/参照）。"
                f"単発のCOMSOL run-stage時間比候補は{_fmt(run_stage_ratio)}、"
                f"end-to-end比候補は{_fmt(end_to_end_ratio)}だが、各経路n=1である。"
                f"最新Swarm {tail_case_count}点はtail/grid gateが{tail_grid_pass_count}/{tail_case_count}で合格した一方、"
                f"solver convergenceは{solver_pass_count}/{tail_case_count}、全case gateは"
                f"{formal_case_pass_count}/{tail_case_count}に留まり、正式bundleとは扱わない。"
            ),
        },
        {
            "id": "english_abstract",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## English abstract\n\n"
                "This technical report evaluates an auditable external-swarm table path for the full one-dimensional "
                "argon DC glow-discharge gap in COMSOL. Property inspection shows that the archived external and built-in "
                "models both retain `LocalEnergyApproximationE` (LEA), not a local-field formulation. The interface "
                "represents seven lookup quantities as seven X–Y pairs and fourteen COMSOL "
                "property arrays. Only six quantities are candidate-active in LEA because the E/N-to-mean-energy table is inactive; "
                "under LFA the candidate-active set would instead contain five arrays because the energy-mobility and "
                "energy-diffusivity tables have no mean-energy equation to close. The built-in saved `sol3` additionally "
                "contains `En` and `plas.F0` degrees of freedom. A clean rerun after repairing compiler-error "
                "detection and stale-class rejection was blocked before model application by COMSOL license error -10 "
                "(`Product has expired`). Numerical comparisons therefore use archived 14 July 2026 profiles and are "
                "classified as `share_with_caveats`; the historical Swarm table (300 K, 13.3 Pa) was not generated at the "
                "same gas conditions as the COMSOL profiles (293.15 K, 13.3322 Pa). Exact integration of the piecewise-"
                f"linear reconstructions gave relative L2 differences of {_fmt(potential_l2)} for potential and "
                f"{_fmt(en_l2)} for E/N, but {_fmt(density_l2)} for electron density, {_fmt(energy_l2)} for mean energy, "
                f"{_fmt(excitation_l2)} for the direct eir2 excitation source, and {_fmt(ionization_l2)} for the direct "
                f"eir4 ionization source. The electron-density integral ratio was {_fmt(density_integral_ratio)} "
                f"(external/reference). Single observations suggested COMSOL run-stage and end-to-end time ratios of "
                f"{_fmt(run_stage_ratio)} and {_fmt(end_to_end_ratio)}, respectively; neither is a repeated performance "
                f"estimate. All {tail_grid_pass_count}/{tail_case_count} new Swarm points passed tail/grid checks, but only "
                f"{solver_pass_count}/{tail_case_count} passed the solver-convergence gate. The result supports the utility of an "
                "explicit, provenance-preserving selected-channel hybrid interface, but not equation-level equivalence, "
                "experimental accuracy, or solver superiority."
            ),
        },
        {
            "id": "headline_metrics_block",
            "type": "metric-strip",
            "cardIds": [
                "potential_l2_card",
                "en_l2_card",
                "density_ratio_card",
                "high_field_fraction_card",
            ],
        },
        {
            "id": "profile_finding",
            "type": "markdown",
            "sourceId": "comparison_metrics_source",
            "body": (
                "## 電位・E/Nの形は近いが、密度と平均エネルギーの量的一致は弱い\n\n"
                "状態量profileは、形状相関と振幅を分けて読む必要がある。電位とE/Nは参照形状をよく追う一方、"
                f"電子密度の積分比は{_fmt(density_integral_ratio)}、平均エネルギーは"
                f"{_fmt(energy_integral_ratio)}（外部/参照）である。組込み値はground truthではなく"
                "実装参照であり、差はclosure、反応契約、履歴入力provenance、vendorモデル差を含む。"
            ),
        },
        {
            "id": "profile_figure",
            "type": "html",
            "body": _figure_html(
                "fig03_profile_state.png",
                "電子密度、平均電子エネルギー、電位、E/Nの外部Swarmと組込み参照profile",
                "状態量profile。200 V、13.3322 Pa、2026年7月14日履歴証拠。",
            ),
        },
        {
            "id": "source_current_finding",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## 直接反応源と伝導電流の差は、幾何学的領域を明示して読む\n\n"
                f"eir2直接励起源・eir4直接電離源の区分線形・空間相対L2は"
                f"{_fmt(excitation_l2)}・{_fmt(ionization_l2)}で、電位やE/Nより大きい。"
                "ここで表示する電流は電子伝導電流と電子＋Ar⁺伝導電流であり、変位電流等を含む「全電流」ではない。"
                f"200要素の電子＋Ar⁺伝導電流RSDは中央80%で"
                f"{row_pct(external_central_current, 'relative_standard_deviation')}、全領域で"
                f"{row_pct(external_full_current, 'relative_standard_deviation')}。中央80%は両端の物理長10%を"
                "除いた幾何学的窓であり、bulk/sheath境界を物理的に同定したものではない。"
                "200/400要素の履歴summaryは点数等重みの旧RSDで、400要素raw profileがない。"
                "図では両定義を別panelへ分離し、数値比較しないため、完全なgrid-convergence証拠ではない。"
            ),
        },
        {
            "id": "source_current_figure",
            "type": "html",
            "sourceId": "current_uniformity_source",
            "body": _figure_html(
                "fig04_sources_currents_mesh.png",
                "eir2直接励起源、eir4直接電離源、電子伝導電流、電子＋Ar+伝導電流と200・400要素感度",
                "direct reaction channels・conductive current・mesh感度。中央80%は幾何学的窓。raw profileと旧summaryを相互比較しない。",
            ),
        },
        {
            "id": "spatial_robustness_finding",
            "type": "markdown",
            "sourceId": "regional_error_source",
            "body": (
                "## 差の場所と中央窓の選び方を感度解析した\n\n"
                "二乗誤差を陰極側10%、中央80%、陽極側10%へ分解し、同時に電子＋Ar⁺伝導電流RSDを"
                "中央保持率の関数として評価した。これは境界領域が指標を支配する物理量と、中央領域の振幅差が"
                "支配する物理量を区別するための頑健性確認である。領域名は幾何学的な位置だけを表し、"
                "sheath edgeを独立に同定したものではない。"
            ),
        },
        {
            "id": "spatial_robustness_figure",
            "type": "html",
            "sourceId": "current_rsd_source",
            "body": _figure_html(
                "fig07_spatial_robustness.png",
                "物理量別二乗誤差の陰極側10%、中央80%、陽極側10%寄与と伝導電流RSDの中央保持率感度",
                "領域別誤差帰属とRSD窓感度。区分線形再構成を区間ごとに厳密積分。",
                max_height_px=500,
            ),
        },
        {
            "id": "weighted_metric_explanation",
            "type": "markdown",
            "sourceId": "comparison_metrics_source",
            "body": (
                "## 区分線形・空間L2は局所差、積分比は総量差を表す\n\n"
                "主指標は √[∫(外部−参照)²dx / ∫参照²dx] であり、非一様格子上の局所的な差を評価する。"
                "union grid上の外部・参照profileを区分線形とみなし、その積を各区間で厳密積分する。"
                "従来の台形-of-squared-samplesはcurrent spikeへの感度が大きいため感度列へ格下げした。"
                "値が小さいほど参照に近いが、相関や積分量の一致を保証しない。"
            ),
        },
        {"id": "weighted_l2_block", "type": "chart", "chartId": "weighted_l2_chart"},
        {
            "id": "integral_ratio_explanation",
            "type": "markdown",
            "sourceId": "comparison_metrics_source",
            "body": (
                "積分比は∫外部dx / ∫参照dxで、1が総量の一致を示す。"
                f"電子密度{_fmt(density_integral_ratio)}、平均エネルギー{_fmt(energy_integral_ratio)}、"
                f"eir4直接電離源{_fmt(ionization_integral_ratio)}という結果から、"
                "形状が近い物理量でも振幅差が残ることが分かる。"
            ),
        },
        {"id": "integral_ratio_block", "type": "chart", "chartId": "integral_ratio_chart"},
        {
            "id": "metric_summary_figure",
            "type": "html",
            "sourceId": "comparison_metrics_source",
            "body": _figure_html(
                "fig05_metric_summary.png",
                "8物理量の区分線形・空間相対L2と積分比",
                "主指標を分離表示。L2は局所差、積分比は総量差を表す。",
            ),
        },
        {
            "id": "weighted_metrics_table_explanation",
            "type": "markdown",
            "sourceId": "comparison_metrics_source",
            "body": (
                "次表は図の正確な値と重み付き相関を併記する。符号を持つ電流の積分比は方向規約を保った比であり、"
                "絶対値比ではない。"
            ),
        },
        {"id": "weighted_metrics_table_block", "type": "table", "tableId": "weighted_metrics_table"},
        {
            "id": "runtime_finding",
            "type": "markdown",
            "sourceId": "runtime_summary_source",
            "body": (
                "## runtime差は有望な仮説だが分布ではない\n\n"
                f"履歴値ではCOMSOL run-stageが外部{row_num(external_run_stage, 'median_seconds')} s・"
                f"組込み参照{row_num(builtin_run_stage, 'median_seconds')} s、end-to-endが外部"
                f"{row_num(external_end_to_end, 'median_seconds')} s・組込み参照"
                f"{row_num(builtin_end_to_end, 'median_seconds')} sである。"
                "run-stage timerはCOMSOLの`run`呼出しを囲むもので、内部solverだけの時間ではない。"
                "外部end-to-endはapply・verify・run・exportを含む。ただし全てn=1で中央値・最小・最大が同じため、"
                f"{_fmt(run_stage_ratio)}倍と{_fmt(end_to_end_ratio)}倍を再現可能な性能値として採用しない。"
            ),
        },
        {"id": "runtime_block", "type": "chart", "chartId": "runtime_chart"},
        {
            "id": "runtime_figure",
            "type": "html",
            "sourceId": "runtime_summary_source",
            "body": _figure_html(
                "fig06_runtime_breakdown.png",
                "COMSOL run-stageおよびend-to-end時間と外部経路stage breakdown",
                "各経路n=1。点は観測値であり、分布・中央値推定ではない。",
            ),
        },
        {
            "id": "scope_definitions",
            "type": "markdown",
            "sourceId": "comsol_model_source",
            "body": (
                "## 対象、比較基準、責任境界\n\n"
                "対象はCOMSOL公式Application Libraryの連成1D Ar DC glow-dischargeモデルである。"
                "計算区間はgrounded cathodeからpowered anodeまでのfull gapであり、positive-column部分だけを"
                "切り出した計算ではない。未接続の2D GEC CCPは除外する。"
                "履歴外部経路はCOMSOLのLEAを保持し、Swarm表が電子・エネルギー輸送係数とeir2/eir4 Townsendを供給する。"
                "COMSOLは電子密度・平均エネルギー・Poissonを含む空間PDE、壁・重粒子過程、回路、幾何、mesh、"
                "非線形runを保持する。E/N→平均エネルギー表はLEAではstored/inactiveである。"
                "ゆえに外部経路は、selected transportとdirect reaction channelだけを外部化するhybrid closureである。"
                "built-in Boltzmannは参照実装であり、実験的ground truthではない。"
                "vendorモデル間の要素次数、反応在庫、保存solver状態も異なり得るため、観測差をclosureだけの因果差と"
                "解釈しない。次表は、どの差が未統制かを明示する。"
            ),
        },
        {
            "id": "model_contrast_block",
            "type": "table",
            "tableId": "model_contrast_table",
        },
        {
            "id": "workflow_boundary_explanation",
            "type": "markdown",
            "sourceId": "interface_contract_source",
            "body": (
                "責任境界では7 lookup量を7 X–Y pair、すなわち14 COMSOL property配列として書込むが、"
                "active setはroot定式化で変わる。"
                "LEA候補activeは6、LFA候補activeは5で、7つすべてが同時にactiveになる標準定式化ではない。"
                "これはプラズマモデル全体、wall chemistry、回路、meshを移植するものではない。検証はbundle内容、"
                "root formulation、COMSOL table property、定式化別の空間profileの順に異なる層で行う。"
            ),
        },
        {
            "id": "workflow_boundary_figure",
            "type": "html",
            "body": _figure_html(
                "fig01_workflow_boundary.png",
                "schema v2 SwarmからCOMSOL空間解までのデータフローと責任境界",
                "schema v2 YAML → SQLite → tables → bundle → COMSOL適用・照合 → 空間解。",
            ),
        },
        {
            "id": "lookup_domain_explanation",
            "type": "markdown",
            "sourceId": "analysis_summary_source",
            "body": (
                "## lookup範囲を覆うことと、drift-diffusion適用域は別問題である\n\n"
                f"solver-QA候補のschema-v2表はE/N={_fmt(float(clamp_context.get('table_E_over_N_min_Td', float('nan'))))}"
                f"–{_fmt(float(clamp_context.get('table_E_over_N_max_Td', float('nan'))))} Td、"
                f"履歴COMSOL profileの使用域は{_fmt(float(clamp_context.get('profile_E_over_N_min_Td', float('nan'))))}"
                f"–{_fmt(float(clamp_context.get('profile_E_over_N_max_Td', float('nan'))))} Tdである。"
                f"表下限未満は全長の{_pct(_number(low_field, 'spatial_fraction'))}。"
                f"一方、履歴外部profileではCOMSOLの典型的目安{_fmt(high_field_threshold)} Tdを超える領域が"
                f"空間長の{_pct(high_field_length_fraction)}で、eir4直接電離源の絶対積分の"
                f"{_pct(high_field_ionization_share)}を担う。これは硬い妥当性境界ではないが、"
                "table coverageだけでは局所・drift-diffusion近似の物理妥当性を保証できないことを示す。"
                "図は最新表と履歴profileを重ねた混合vintageである。"
            ),
        },
        {
            "id": "lookup_domain_figure",
            "type": "html",
            "sourceId": "clamp_metrics_source",
            "body": _figure_html(
                "fig02_lookup_domain.png",
                "E/Nと平均エネルギー、輸送係数、励起・電離Townsendのlookup範囲",
                "2026年7月31日schema-v2表と2026年7月14日COMSOL使用範囲。goldはclamp領域。",
            ),
        },
        {
            "id": "activation_regime_finding",
            "type": "markdown",
            "sourceId": "activation_audit_source",
            "body": (
                "## archived exportは輸送・Townsend表と数値整合するが、activation証明ではない\n\n"
                f"E/N→平均エネルギー表との相対誤差medianは"
                f"{row_num(mean_energy_audit, 'relative_error_median')}で、1%超の点は"
                f"{row_int(mean_energy_audit, 'points_above_1_percent')}。これはLEAで平均エネルギー表が"
                "inactiveであることと整合する。対して換算移動度、eir2、eir4は意図した表と数値的に整合した"
                f"（最大相対誤差: {row_num(mobility_audit, 'relative_error_max')}、"
                f"{row_num(excitation_audit, 'relative_error_max')}、"
                f"{row_num(ionization_audit, 'relative_error_max')}）。ただし履歴compileにはstale-class riskがあり、"
                "これらはactivationを証明しない。root property readback、class hash、property readback、"
                "空間照合が同じclean runで揃って初めて正式合格となる。"
            ),
        },
        {
            "id": "activation_regime_figure",
            "type": "html",
            "sourceId": "regime_impact_source",
            "body": _figure_html(
                "fig08_activation_regime.png",
                "履歴export対intended tableの相対誤差と500 Td超領域の物理量別絶対積分寄与",
                "左: activationではなく数値整合監査。右: COMSOLの典型的drift-diffusion目安を超える領域の寄与。",
                max_height_px=500,
            ),
        },
        {
            "id": "activation_audit_table_block",
            "type": "table",
            "tableId": "activation_audit_table",
        },
        {
            "id": "lea_floor_and_tail_audit",
            "type": "markdown",
            "sourceId": "tail_audit_source",
            "body": (
                "## LEAの低平均エネルギーfloorとSwarm収束gate\n\n"
                "履歴外部MPHはLEAなので、低エネルギー監査の主対象はinactiveな"
                "E/N→平均エネルギー表ではなく、PDEが返す平均エネルギーを引数にするtransport／Townsend lookupである。"
                f"最新bundleのconstant boundary policyは平均エネルギー"
                f"{_fmt(_number(mean_energy_floor, 'threshold'))} eV以下で働き、履歴外部profileとの混合vintage overlayでは"
                f"全長の{_pct(_number(mean_energy_floor, 'spatial_fraction'))}、電子＋Ar⁺伝導電流の絶対積分の"
                f"{row_pct(floor_current_impact, 'absolute_integral_share')}に相当する。小さいがゼロではなく、"
                "正式runでは同一vintageで再評価する必要がある。\n\n"
                f"schema-v2 Swarmの{tail_case_count}点はtail probability、edge/peak、20 keV grid-limitをまとめた"
                f"tail/grid gateに{tail_grid_pass_count}/{tail_case_count}で合格し、最大energy gridは"
                f"{_fmt(max_tail_grid_eV)} eVまで自動拡張された。しかしsolver convergenceは"
                f"{solver_pass_count}/{tail_case_count}、tailとsolverを合わせたcase gateも"
                f"{formal_case_pass_count}/{tail_case_count}である。tailが十分小さいことは、反復solverが収束したことを"
                "代替しないため、このSwarm bundleも正式結果へ昇格させない。"
            ),
        },
        {
            "id": "tail_failure_table_block",
            "type": "table",
            "tableId": "tail_failure_table",
        },
        {
            "id": "fluid_applicability_finding",
            "type": "markdown",
            "sourceId": "fluid_applicability_source",
            "body": (
                "## weak ionizationは満たすが、境界局所のscale separationは別監査である\n\n"
                f"最大nₑ/Nは外部{row_num(external_fluid, 'max_electron_to_neutral_ratio')}、参照"
                f"{row_num(reference_fluid, 'max_electron_to_neutral_ratio')}で、両経路ともweakly ionizedである。"
                "一方、|Jₑ|/[e nₑ sqrt(2e mean-energy/mₑ)]はdriftだけでなくdiffusionとboundary fluxを含む"
                "operational ratioであり、COMSOL validity cutoffではない。"
                f"{row_num(external_fluid, 'net_flux_speed_ratio_operational_guide')}を超える空間長割合は外部"
                f"{row_pct(external_fluid, 'spatial_fraction_net_flux_speed_ratio_above_guide')}、参照"
                f"{row_pct(reference_fluid, 'spatial_fraction_net_flux_speed_ratio_above_guide')}で、主にboundary-localである。"
                "したがってpressureと弱電離だけからfluid/local closureの妥当性を確定しない。"
            ),
        },
        {
            "id": "fluid_applicability_figure",
            "type": "html",
            "sourceId": "fluid_applicability_source",
            "body": _figure_html(
                "fig09_fluid_applicability.png",
                "電子密度対中性密度比とnet electron-flux speed対energy speedの空間profile",
                "図9. weak-ionization proxyとboundary-local scale-separation診断。0.1は監査用guideで、COMSOL validity cutoffではない。",
                max_height_px=420,
            ),
        },
        {
            "id": "operating_point_finding",
            "type": "markdown",
            "sourceId": "operating_point_source",
            "body": (
                "## 外部表はballast回路を介して放電operating pointも変える\n\n"
                f"source voltageは両経路で{row_num(external_operating_point, 'source_voltage_V')} Vだが、"
                f"gap endpoint voltageは外部{row_num(external_operating_point, 'gap_endpoint_voltage_V')} V、"
                f"参照{row_num(reference_operating_point, 'gap_endpoint_voltage_V')} V、ballast currentは"
                f"外部{_fmt(1e3 * _number(external_operating_point, 'ballast_current_A'))} mA、参照"
                f"{_fmt(1e3 * _number(reference_operating_point, 'ballast_current_A'))} mAである。"
                "したがってprofile差は固定gap-voltageでclosureだけを交換した応答ではなく、外部表→plasma impedance→"
                "ballast drop/current→field・sourceというcircuit feedbackを含む。内部potential maximumはfield-reversal"
                " proxyに留まり、signed electric fieldをexportしていないため反転位置の確証とはしない。"
            ),
        },
        {
            "id": "operating_point_table_block",
            "type": "table",
            "tableId": "operating_point_table",
        },
        {
            "id": "source_form_and_partial_energy",
            "type": "markdown",
            "sourceId": "townsend_representation_source",
            "body": (
                "## Townsend-vs-rate formとelectron-energyは部分監査として分ける\n\n"
                "履歴eir2/eir4 sourceは(α/N)N|Jₑ|/eのTownsend-flux formであり、同じprofileからの再構成相対差は"
                f"それぞれ{row_num(excitation_representation, 'townsend_reconstruction_relative_difference')}、"
                f"{row_num(ionization_representation, 'townsend_reconstruction_relative_difference')}である。"
                "一方、N nₑ k(mean-energy)というrate-form反実仮想は観測Townsend源の"
                f"{row_num(excitation_representation, 'counterfactual_rate_over_observed_townsend')}倍と"
                f"{row_num(ionization_representation, 'counterfactual_rate_over_observed_townsend')}倍になった。"
                "これはsource representationがoperating pointへ影響し得ることを示すsensitivityで、matched COMSOL rerun"
                "でもrate formの優越証拠でもない。\n\n"
                f"またeir2/eir4 threshold loss合計はsigned electron field powerの外部"
                f"{row_pct(external_partial_energy, 'direct_loss_to_field_power_ratio')}、参照"
                f"{row_pct(reference_partial_energy, 'direct_loss_to_field_power_ratio')}である。"
                "energy flux、elastic/superelastic、stepwise/Penning、secondary emission、boundary項、時間微分を"
                "含まないため、これは完全なelectron-energy balanceでも保存則残差でもない。"
            ),
        },
        {
            "id": "townsend_representation_table_block",
            "type": "table",
            "tableId": "townsend_representation_table",
        },
        {
            "id": "partial_energy_table_block",
            "type": "table",
            "tableId": "partial_energy_table",
        },
        {
            "id": "metric_definitions",
            "type": "markdown",
            "sourceId": "analysis_summary_source",
            "body": (
                "## 指標定義と解析単位\n\n"
                f"外部と参照のprofileは共通区間{_fmt(metric('electric_potential', 'x_min_m'))}"
                f"–{_fmt(metric('electric_potential', 'x_max_m'))} mのunion gridへ補間し、各区間の"
                "区分線形再構成の積を厳密積分した。相関も同じ連続再構成の平均・分散・共分散から求める。"
                "電流RSDは区分線形空間population標準偏差を空間平均の絶対値で除す。"
                "中央80%は両端10%の物理長を除きcut pointを補間挿入する幾何学的窓で、bulk/sheathの物理境界ではない。"
                f"共通grid点数は{int(metric('electric_potential', 'common_grid_rows'))}である。"
                "台形-of-squared-samplesと非重み付き離散指標は感度・履歴追跡用にのみ残す。"
            ),
        },
        {
            "id": "reproducible_method",
            "type": "markdown",
            "sourceId": "formal_config_source",
            "body": (
                "## 再現可能な外部入力経路\n\n"
                "再実行候補inputは`schema_version: 2`、`run.solvers: [{id: two_term}]`、293.15 K、13.3322 Paである。"
                f"ただし現Swarm case gateは{formal_case_pass_count}/{tail_case_count}で、正式input認定前に"
                "未収束点のsolver convergenceを解消しなければならない。"
                "実用経路は schema-v2 YAML → SQLite → canonical tables → COMSOL bundle → MPH適用 → "
                "activation照合 → profile比較。各runはclean class directoryを使い、compilerがreturn code 0でも"
                "失敗文字列を検出し、古いclass・不変MPHを拒否する。Java/class、COMSOL build、電圧列、input MPH、"
                "bundle、mapping、設定hashをprovenanceへ保存する。"
            ),
        },
        {
            "id": "reproduction_commands",
            "type": "markdown",
            "sourceId": "external_runner_source",
            "body": (
                "## 再現コマンド\n\n"
                "repository rootから次の順に実行する。外部runnerは`run-comsol`全体を独立processで3回実行し、"
                "apply・verify・run/solve・export・wall end-to-endを分離して保存する。COMSOL licenseが無効なら"
                "本報告と同じくapply前に明示失敗する。output rootは既存pathを再利用しない。\n\n"
                "```powershell\n"
                "$COMSOL = \"C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics_copy1\\bin\\win64\\comsolbatch.exe\"\n"
                "py -3 -m swarm_workflow.cli sweep examples\\workflow_argon_comsol_benchmark_2026.yaml\n"
                "py -3 -m swarm_workflow.cli build-tables outputs\\comsol_swarm_benchmark_2026\\swarm.sqlite `\n"
                "  --output outputs\\comsol_swarm_benchmark_2026\\tables --source two_term\n"
                "py -3 -m swarm_workflow.cli export-comsol outputs\\comsol_swarm_benchmark_2026\\tables\\mixture_0000 `\n"
                "  --output outputs\\comsol_swarm_benchmark_2026\\bundle\n"
                "py -3 reports\\comsol_swarm_benchmark_2026\\repro\\prepare_low_energy_bundle.py `\n"
                "  outputs\\comsol_swarm_benchmark_2026\\bundle `\n"
                "  outputs\\comsol_swarm_benchmark_2026\\rerun_bundle_low_energy\n"
                "py -3 reports\\comsol_swarm_benchmark_2026\\repro\\run_external_benchmark.py `\n"
                "  --mapping Model\\maps\\positive_column_external.yaml `\n"
                "  --bundle outputs\\comsol_swarm_benchmark_2026\\rerun_bundle_low_energy `\n"
                "  --comsol $COMSOL `\n"
                "  --output-root outputs\\comsol_swarm_benchmark_2026\\external_200elem_runs_clean `\n"
                "  --repetitions 3 --mesh-elements 200 `\n"
                "  --pressure-Pa 13.3322 --gas-temperature-K 293.15\n"
                "py -3 reports\\comsol_swarm_benchmark_2026\\repro\\run_external_benchmark.py `\n"
                "  --mapping Model\\maps\\positive_column_external.yaml `\n"
                "  --bundle outputs\\comsol_swarm_benchmark_2026\\bundle_low_energy_portable `\n"
                "  --comsol $COMSOL `\n"
                "  --output-root outputs\\comsol_swarm_benchmark_2026\\external_400elem_runs_clean `\n"
                "  --repetitions 1 --mesh-elements 400 `\n"
                "  --pressure-Pa 13.3322 --gas-temperature-K 293.15\n"
                "py -3 reports\\comsol_swarm_benchmark_2026\\repro\\run_builtin_benchmark.py `\n"
                "  --input-mph Model\\positive_column_1d_boltzmann.mph --comsol $COMSOL `\n"
                "  --output-root outputs\\comsol_swarm_benchmark_2026\\builtin_200elem_runs `\n"
                "  --repetitions 3 --mesh-elements 200 --voltage-V 200 `\n"
                "  --pressure-Pa 13.3322 --gas-temperature-K 293.15\n"
                "py -3 reports\\comsol_swarm_benchmark_2026\\repro\\benchmark_analysis.py\n"
                "py -3 reports\\comsol_swarm_benchmark_2026\\repro\\build_notebook.py\n"
                "```\n\n"
                "第三者は`data/derived/analysis_summary.json`のformal gatesとsource hashを確認し、"
                "`benchmark_analysis_executed.ipynb`を先頭から再実行してCSV・図を再生成できる。"
            ),
        },
        {
            "id": "closure_contract_explanation",
            "type": "markdown",
            "sourceId": "interface_contract_source",
            "body": (
                "## written propertyとactive closureを分け、root定式化を最初に照合する\n\n"
                "平均エネルギー、換算移動度、縦方向換算拡散、換算エネルギー移動度、換算エネルギー拡散、"
                "eir2直接励起Townsend、eir4直接電離Townsendの7 lookup量について、各X/Yの"
                "14 COMSOL property配列を全行readbackする。"
                "しかし履歴MPHの`MeanElectronEnergyModel=LocalEnergyApproximationE`では、E/N→平均エネルギー表は"
                "stored/inactiveで、候補activeは残る6量である。LFAへ切り替えるなら平均エネルギー表はactiveになる一方、"
                "平均エネルギーPDEが消えるためenergy mobility/diffusivityはinactiveとなり、候補activeは5量である。"
                "LEAの空間照合契約はμN、DLN、eir2 α/N、eir4 α/Nの4量。energy transport 2量は安定した"
                "空間変数が未同定なのでproperty readbackに留める。修復後実runはlicenseで停止したため、実MPHでの"
                "formulation readback・property readback・空間照合は未完了である。"
                "また、eir2/eir4以外のreaction inventoryと回路・wall・heavy-species closureはCOMSOL側に残るため、"
                "「外部Swarmが全放電closureを置換した」とは記述しない。"
            ),
        },
        {
            "id": "mph_formulation_evidence",
            "type": "markdown",
            "sourceId": "mph_formulation_evidence_source",
            "body": (
                "### MPH readbackの物理的意味\n\n"
                "`Model/positive_column_1d.mph`と履歴"
                "`Model/work/positive_column_external_applied.mph`のroot physics propertyはともに"
                "`ElectronProperties/MeanElectronEnergyModel=LocalEnergyApproximationE`である。"
                "`Model/positive_column_1d_boltzmann.mph`もrootはLEAだが、保存済み`sol3/std3`のsolve-forには"
                "`comp1.En`、`comp1.Ne`、`comp1.plas.F0`ほかが含まれる。したがって「external LFA対built-in LEA」"
                "という比較ではなく、「外部表で輸送・Townsendを与えるLEA対、EEDF自由度を含む組込みBoltzmann参照」"
                "として解釈する。"
            ),
        },
        {"id": "closure_contract_block", "type": "table", "tableId": "closure_contract_table"},
        {
            "id": "provenance_and_license",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## 履歴compiler defectは修復したが、licenseが正式計算を遮断した\n\n"
                "2026年7月14日の記録はcompiler return code 0に対して標準出力に`Failed to compile java file`と"
                "`Compilation failed`を含み、実行classが生成Javaと一致した証拠がない。修復後runnerはrunごとの"
                "classを生成しSHA-256を保存したが、2026年7月31日のapplyはCOMSOL license error -10 "
                "（`Product has expired`）でmodel apply前に停止した。よって古いMPHを新結果として採用せず、"
                "正式COMSOL数値と3反復時間は欠測のままである。"
            ),
        },
        {
            "id": "external_input_utility",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## 外部Swarm入力の有用性は交換可能性と監査可能性にある\n\n"
                "外部表はCOMSOL空間PDEの責任範囲を保ったまま、定式化に応じたelectron closure係数を交換し、"
                "単位、有効範囲、process、quality、"
                "入力hashをbundle manifestで追跡できる。同じ表を複数の空間条件へ再利用できる可能性があり、"
                "局所closureが妥当な範囲ではrun時間短縮候補になる。一方、非局所性、時間履歴、e-e衝突、磁場、"
                "state-resolved populationなどを表が表現しない場合、`feature_policy`に従ってfail・skip・"
                "明示approximationとし、要求物理を黙って無視してはならない。"
                "本実装はselected-channel hybridなので、再利用時も回路と残存COMSOL chemistryを含む"
                "operating-point変化を追跡し、外部表単独の効果と呼ばない。"
            ),
        },
        {
            "id": "swarm_transport_conventions",
            "type": "markdown",
            "sourceId": "swarm_transport_source",
            "body": (
                "### 本Swarm表の規約上の制約\n\n"
                "現`two_term` transport metadataはflux定義である。換算Townsendはrate coefficientを"
                "電子drift speedの絶対値で除して構成するため、COMSOLのflux-Townsend源形式と整合する。"
                "一方、現在の縦方向拡散と横方向拡散は同じscalar-f0 momentから評価され、DL=DTである。"
                "したがって「縦方向換算拡散」というCOMSOL property名へ接続しても、独立なanisotropic diffusion"
                "計算を行ったとは主張しない。第三者adapterではflux/bulkとlongitudinal/transverseの定義を必須metadata"
                "として保持する必要がある。"
            ),
        },
        {
            "id": "third_party_adapter",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## BOLSIG+・MCIG等を直結するにはsourceとsolver modeを分離したadapterが必要\n\n"
                "現製品では第三者Swarm出力はbenchmark referenceとして読めてもCOMSOL bundleへ直結しない。将来adapterは"
                "source名・version・入力hashを`two_term`、`multi_term`、`monte_carlo`というcanonical solver idと"
                "分離し、flux/bulk、縦/横拡散、単位付き輸送係数、過程別rate/Townsend、EEDF規約、混合比・温度・圧力、"
                "有効範囲、外挿policy、不確かさ、feature treatmentを保持する必要がある。必須量の欠損や規約不明は"
                "明示的に失敗させる。公開schemaとadapter実装自体は本稿の範囲外である。"
            ),
        },
        {
            "id": "primary_sources",
            "type": "markdown",
            "body": (
                "## 物理記述の一次資料\n\n"
                "- [COMSOL 6.4 official coupled DC Glow Discharge, 1D model]"
                "(https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.positive_column_1d/positive_column_1d.html)\n"
                "- [COMSOL Plasma Module User's Guide 6.4]"
                "(https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/PlasmaModuleUsersGuide.pdf)\n"
                "- [COMSOL Drift Diffusion interface: LEA/LFA]"
                "(https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.02.html)\n"
                "- [COMSOL Drift Diffusion Model: lookup tables]"
                "(https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.04.html)\n"
                "- [COMSOL Drift Diffusion Theory]"
                "(https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.17.html)\n"
                "- [Hagelaar and Pitchford (2005), BOLSIG+ method]"
                "(https://doi.org/10.1088/0963-0252/14/4/011)\n"
                "- [Pancheshnyi et al. (2012), LXCat]"
                "(https://doi.org/10.1016/j.chemphys.2011.04.020)\n"
                "- [Dujko and White (2008), flux/bulk transport]"
                "(https://doi.org/10.1088/1742-6596/133/1/012005)\n"
                "- [Dujko et al. (2008), cathode-fall nonhydrodynamic transport]"
                "(https://doi.org/10.1088/0022-3727/41/24/245205)\n\n"
                "積分断面積だけから完全な微分散乱は決まらない。本Swarmの`two_term`表は局所的なtwo-term"
                "Boltzmann計算から作るが、COMSOL履歴空間モデルのroot定式化はLEAである。この2つの「local」を"
                "混同しない。"
                "DCSに基づくexact multi-term solverとは記述しない。"
            ),
        },
        {
            "id": "limitations",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## 制約、不確かさ、頑健性\n\n"
                "本結果は1D Ar、DC、無磁場、履歴profile一組に限定される。clean COMSOL solve、独立3反復、"
                "400要素raw profile、formulation-aware activation照合がなく、履歴class provenanceにも欠陥がある。"
                "7 lookup量（14 property配列）の存在は7 active closureを意味せず、履歴LEAでは平均エネルギー表がinactiveである。"
                f"最新Swarmはtail/grid {tail_grid_pass_count}/{tail_case_count}に対しsolver convergence "
                f"{solver_pass_count}/{tail_case_count}で、全点正式合格ではない。最新Swarm表と"
                "履歴COMSOL profileを重ねたclamp評価はvintageが異なる。履歴Swarm表の300 K・13.3 Paと"
                "COMSOLの293.15 K・13.3322 Paも一致しない。built-in pathとの反応形式・保存solver状態・"
                "要素次数も完全同一ではなく、selected-channel hybrid closureとballast回路feedbackを含む。"
                "Townsend-vs-rateとelectron energyはpostprocess／部分監査に留まる。実験data、cross-section uncertainty、wall chemistry・metastable感度を"
                "評価していないため、自然に対するvalidationでもsolver優越の証明でもない。"
            ),
        },
        {
            "id": "legacy_metrics_explanation",
            "type": "markdown",
            "sourceId": "comparison_metrics_source",
            "body": (
                "## 補足：旧非重み付き離散指標\n\n"
                "旧artifactとの追跡可能性のため、点数を等重みとした離散L2と相関を次表に残す。非一様meshでは"
                "空間測度を正しく表さないため、本文の主張・順位・結論には用いない。"
            ),
        },
        {"id": "legacy_metrics_block", "type": "table", "tableId": "legacy_metrics_table"},
        {
            "id": "recommended_next_steps",
            "type": "markdown",
            "sourceId": "manuscript_source",
            "body": (
                "## 次に必要な検証\n\n"
                "1. 有効なCOMSOL license下で外部系・組込み参照を独立processで各3回実行し、中央値と範囲を保存する。\n"
                "2. 200要素を主結果、200/400要素を感度確認として両方のraw profileを保存する。\n"
                f"3. Swarm {tail_case_count}点すべてでtail/gridとsolver convergenceを合格させ、未収束点を含むbundleを正式入力にしない。\n"
                "4. Java/class、build、電圧列、MPH、bundle、config hashが揃ったrunだけを正式採用する。\n"
                "5. root formulationをreadbackし、7 X–Y pair／14 property配列と定式化別active setを照合する。LEAでは4空間量を照合し、energy transport 2量は変数契約を追加するか制約として残す。\n"
                "6. 同一vintageでmean-energy floor、高E/N、fluid diagnostic、circuit operating point、Townsend/rate form、完全energy balanceを再評価する。"
            ),
        },
        {
            "id": "further_questions",
            "type": "markdown",
            "body": (
                "## 今後の問い\n\n"
                "- equation-levelの反応契約とvendorモデル差を揃えたとき、profile差のどこまでが残るか。\n"
                "- n≥3の独立JVM・license checkoutで、速度比はどの範囲に収まるか。\n"
                "- 低電界および低平均エネルギーfloorは同一runの積分量へどの程度影響するか。\n"
                "- 第三者adapterの最小必須contractは、どのsolver familyまで共通化できるか。"
            ),
        },
    ]

    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "description": (
                "COMSOLのfull 1D Ar DC glow-dischargeに対するselected-channel hybrid Swarm表入力を、"
                "LEA/LFA active set、circuit feedback、fluid premise、tail/solver convergenceまで監査し、"
                "履歴比較と再現性制約を整理した和文技術報告。English abstractを含む。"
            ),
            "generatedAt": generated_at,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": manifest_sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "partial",
            "datasets": {
                "headline_metrics": headline_metrics,
                "comparison_metrics": comparison,
                "runtime_comparison": runtime,
                "clamp_metrics": clamp,
                "regional_error_attribution": regional,
                "current_uniformity": current_uniformity,
                "current_rsd_sensitivity": current_rsd_sensitivity,
                "regime_metrics": regime_metrics,
                "regime_impact": regime_impact,
                "fluid_applicability_diagnostics": fluid_applicability,
                "operating_point_diagnostics": operating_points,
                "townsend_source_representation": townsend_representation,
                "partial_electron_energy_audit": partial_energy,
                "tail_convergence_audit": tail_audit,
                "tail_solver_failures": tail_solver_failures,
                "clamp_impact": clamp_impact,
                "historical_activation_audit": activation_audit,
                "interface_contract": interface_rows,
                "model_contrast": model_contrast_rows,
                "report_status": [
                    {
                        "data_quality_status": "share_with_caveats",
                        "formal_result_available": False,
                        "formal_swarm_input_available": False,
                        "tail_grid_gate_passed": tail_grid_pass_count,
                        "solver_convergence_gate_passed": solver_pass_count,
                        "swarm_cases_total": tail_case_count,
                        "license_error": -10,
                    }
                ],
            },
            "accessIssues": [
                {
                    "id": "formal_rerun_license_blocked",
                    "scope": "正式COMSOL再実行",
                    "sourceId": "license_attempt_source",
                    "message": "COMSOL license error -10（Product has expired）によりmodel apply前で停止した。",
                },
                {
                    "id": "historical_compile_provenance_defect",
                    "scope": "2026-07-14履歴結果",
                    "sourceId": "compile_archive_source",
                    "message": "compiler return code 0でも標準出力にcompile失敗があり、stale classを排除できない。",
                },
                {
                    "id": "runtime_repetitions_missing",
                    "dataset": "runtime_comparison",
                    "sourceId": "runtime_summary_source",
                    "message": "外部・組込み参照とも各scope n=1で、中央値と範囲を推定できない。",
                },
                {
                    "id": "formulation_aware_activation_incomplete",
                    "dataset": "interface_contract",
                    "sourceId": "interface_contract_source",
                    "message": "7 lookup量＝7 X-Y pair＝14 property配列だが、property存在≠functional activation。履歴LEA候補activeは6量、LFAなら5量。",
                },
                {
                    "id": "formal_swarm_solver_convergence_incomplete",
                    "dataset": "tail_convergence_audit",
                    "sourceId": "tail_audit_source",
                    "message": f"tail/grid gateは{tail_grid_pass_count}/{tail_case_count}だが、solver convergenceと全case gateは{solver_pass_count}/{tail_case_count}。formal inputとして未採用。",
                },
                {
                    "id": "historical_activation_provenance_unresolved",
                    "dataset": "historical_activation_audit",
                    "sourceId": "activation_audit_source",
                    "message": "履歴exportとintended tablesの数値整合はactivation証明ではない。compile logのstale-class riskが残る。",
                },
                {
                    "id": "raw_400_mesh_profile_missing",
                    "scope": "mesh感度",
                    "sourceId": "analysis_summary_source",
                    "message": "400要素はsummary値のみでraw profileがなく、独立再計算できない。",
                },
            ],
        },
        "sources": sources,
    }


def main() -> None:
    artifact = build_artifact()
    encoded = json.dumps(artifact, ensure_ascii=False, indent=2).encode("utf-8")
    if len(encoded) > MAX_ARTIFACT_BYTES:
        raise RuntimeError(
            f"artifact payload is {len(encoded):,} bytes, above {MAX_ARTIFACT_BYTES:,}"
        )
    ARTIFACT_PATH.write_bytes(encoded + b"\n")
    print(
        json.dumps(
            {
                "artifact": str(ARTIFACT_PATH),
                "bytes": len(encoded) + 1,
                "blocks": len(artifact["manifest"]["blocks"]),
                "charts": len(artifact["manifest"]["charts"]),
                "tables": len(artifact["manifest"]["tables"]),
                "embedded_figures": 9,
                "data_quality_status": "share_with_caveats",
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
