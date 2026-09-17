"""Read-only, like-for-like exports from saved GEC-ICP COMSOL solutions.

This module deliberately does not run a study or save a model.  It prepares a
small Java program that reopens an existing MPH with ``ModelUtil.loadCopy`` and
exports the same physical fields at an exact common time and at that model's
terminal saved time.  Keeping this separate from :mod:`java` prevents result
comparison concerns from growing the apply/solve generator.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
from typing import Any

from swarm_workflow._io import write_json
from swarm_workflow.comsol.java import java_path, java_string
from swarm_workflow.comsol.runtime import (
    ComsolExecutionSummary,
    ComsolRunContext,
    execute_generated_comsol_java,
)

from .contracts import GecIcpMapping, GecIcpModelSpec


COMPARISON_EXPORT_CLASS = "SwarmGecIcpComparisonExport"

FIELD_EXPRESSIONS = (
    "plas.ne",
    "e_const*plas.ebar",
    "plas.n_wArs",
    "plas.n_wAr_1p",
    "mf.Qrh",
)
FIELD_UNITS = ("1/m^3", "eV", "1/m^3", "1/m^3", "W/m^3")
FIELD_DESCRIPTIONS = (
    "electron density",
    "mean electron energy",
    "argon metastable number density",
    "argon ion number density",
    "RF resistive power deposition density",
)

VOLUME_EXPRESSIONS = (
    "1",
    "plas.ne",
    "plas.ne*e_const*plas.ebar",
    "plas.n_wArs",
    "plas.n_wAr_1p",
    "mf.Qrh",
)
VOLUME_UNITS = ("m^3", "1", "eV", "1", "1", "W")
VOLUME_DESCRIPTIONS = (
    "axisymmetric plasma volume",
    "electron inventory",
    "electron mean-energy inventory",
    "argon metastable inventory",
    "argon ion inventory",
    "plasma RF resistive deposition",
)

EXPECTED_EXPORT_NAMES = (
    "solution_times.csv",
    "fields_common_time.csv",
    "volume_common_time.csv",
    "coil_power_common_time.csv",
    "fields_terminal.csv",
    "volume_terminal.csv",
    "coil_power_terminal.csv",
    "volume_time_series.csv",
    "coil_power_time_series.csv",
)


class GecIcpComparisonExportError(RuntimeError):
    """Raised when a saved-solution comparison export cannot be prepared."""


@dataclass(frozen=True, slots=True)
class GecIcpComparisonExportPlan:
    """Immutable preparation record for one saved GEC-ICP solution."""

    case_id: str
    input_mph: Path
    input_mph_sha256: str
    output_directory: Path
    common_time_s: float
    java_path: Path
    manifest_path: Path
    expected_outputs: tuple[Path, ...]


def prepare_gec_icp_comparison_export(
    mapping: GecIcpMapping,
    *,
    case_id: str,
    input_mph: str | Path,
    output_directory: str | Path,
    common_time_s: float = 1.0e-3,
) -> GecIcpComparisonExportPlan:
    """Prepare, but do not execute, a read-only saved-solution export."""

    source = Path(input_mph).resolve()
    output = Path(output_directory).resolve()
    _validate_inputs(case_id, source, output, common_time_s)
    output.mkdir(parents=True, exist_ok=True)
    java_file = output / f"{COMPARISON_EXPORT_CLASS}.java"
    manifest = output / "comparison_export_manifest.json"
    expected = tuple(output / name for name in EXPECTED_EXPORT_NAMES)
    source_sha256 = _sha256_file(source)
    java_file.write_text(
        generate_comparison_export_java(
            mapping.model,
            input_mph=source,
            output_directory=output,
            common_time_s=common_time_s,
        ),
        encoding="utf-8",
    )
    write_json(
        manifest,
        {
            "schema_version": 1,
            "stage": "prepare-gec-icp-comparison-export",
            "status": "prepared",
            "case_id": case_id,
            "read_only": True,
            "input_mph": {
                "path": str(source),
                "sha256": source_sha256,
                "size_bytes": source.stat().st_size,
            },
            "model_mapping": {
                "path": str(mapping.path.resolve()),
                "sha256": _sha256_file(mapping.path.resolve()),
                "solution": mapping.model.solution,
                "dataset": mapping.model.dataset,
            },
            "time_selection": {
                "common": {
                    "mode": "exact_transient_interpolation",
                    "time_s": common_time_s,
                },
                "terminal": {"mode": "last_saved_solution"},
            },
            "field_expressions": [
                {
                    "expression": expression,
                    "unit": unit,
                    "description": description,
                }
                for expression, unit, description in zip(
                    FIELD_EXPRESSIONS,
                    FIELD_UNITS,
                    FIELD_DESCRIPTIONS,
                    strict=True,
                )
            ],
            "volume_expressions": [
                {
                    "expression": expression,
                    "unit": unit,
                    "description": description,
                }
                for expression, unit, description in zip(
                    VOLUME_EXPRESSIONS,
                    VOLUME_UNITS,
                    VOLUME_DESCRIPTIONS,
                    strict=True,
                )
            ],
            "generated_java": str(java_file),
            "expected_outputs": [str(path) for path in expected],
            "mutation_contract": (
                "loadCopy_only; no study.run, clearSolutionData, or model.save"
            ),
        },
    )
    return GecIcpComparisonExportPlan(
        case_id=case_id,
        input_mph=source,
        input_mph_sha256=source_sha256,
        output_directory=output,
        common_time_s=common_time_s,
        java_path=java_file,
        manifest_path=manifest,
        expected_outputs=expected,
    )


def execute_gec_icp_comparison_export(
    mapping: GecIcpMapping,
    plan: GecIcpComparisonExportPlan,
    *,
    comsol_executable: str | Path | None = None,
) -> ComsolExecutionSummary:
    """Execute one prepared read-only export and verify fresh evidence.

    The source MPH and model mapping are rebound immediately before execution,
    every expected CSV is cleared from the case-local output directory, and the
    source MPH hash is checked again afterwards.  This makes a successful
    manifest evidence of the saved solution that was actually read rather than
    evidence left by an older export.
    """

    _validate_execution_binding(mapping, plan)
    for output in plan.expected_outputs:
        resolved = output.resolve()
        if resolved.parent != plan.output_directory.resolve():
            raise GecIcpComparisonExportError(
                f"comparison output escaped its case directory: {resolved}"
            )
        resolved.unlink(missing_ok=True)

    context = ComsolRunContext(
        # COMSOL's file-system policy permits writes below the batch working
        # directory.  Use the narrow common ancestor of the immutable model,
        # mapping, and case-local export directory rather than the map folder.
        root=_comparison_execution_root(mapping, plan),
        mapping_path=mapping.path.resolve(),
        input_mph=plan.input_mph,
        # The comparison program never writes an MPH.  The runtime field is
        # intentionally bound to the same read-only source and freshness checks
        # remain disabled.
        output_mph=plan.input_mph,
        log_path=plan.output_directory / "logs",
        bundle_path=None,
    )
    try:
        summary = execute_generated_comsol_java(
            context,
            plan.java_path,
            operation=f"gec_icp_comparison_{plan.case_id}",
            comsol_executable=comsol_executable,
            study="saved_solution_read_only",
            require_fresh_output=False,
        )
        _materialize_solution_times(plan, summary.stdout_paths[-1])
        _validate_completed_export(plan)
    except Exception as exc:
        _update_execution_manifest(
            plan,
            status="failed",
            execution=None,
            error=f"{type(exc).__name__}: {exc}",
        )
        raise

    _update_execution_manifest(
        plan,
        status="completed",
        execution=summary,
        error=None,
    )
    return summary


def _materialize_solution_times(
    plan: GecIcpComparisonExportPlan,
    stdout_path: Path,
) -> None:
    prefix = "SWARM_GEC_ICP_SOLUTION_TIME\t"
    indexed: dict[int, float] = {}
    for line in stdout_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(prefix):
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise GecIcpComparisonExportError(
                "malformed GEC ICP solution-time sentinel"
            )
        try:
            index = int(fields[1])
            time_s = float(fields[2])
        except ValueError as exc:
            raise GecIcpComparisonExportError(
                "non-numeric GEC ICP solution-time sentinel"
            ) from exc
        if index < 1 or not math.isfinite(time_s) or index in indexed:
            raise GecIcpComparisonExportError(
                "invalid or duplicate GEC ICP solution-time sentinel"
            )
        indexed[index] = time_s
    expected_indices = list(range(1, len(indexed) + 1))
    if not indexed or sorted(indexed) != expected_indices:
        raise GecIcpComparisonExportError(
            "GEC ICP solution-time sentinels are incomplete"
        )
    times = [indexed[index] for index in expected_indices]
    if any(right <= left for left, right in zip(times, times[1:], strict=False)):
        raise GecIcpComparisonExportError(
            "GEC ICP saved solution times are not strictly increasing"
        )
    if not times[0] <= plan.common_time_s <= times[-1]:
        raise GecIcpComparisonExportError(
            "common comparison time is outside saved solution support"
        )
    payload = ["solution_index,time_s"]
    payload.extend(
        f"{index},{time_s:.17g}"
        for index, time_s in zip(expected_indices, times, strict=True)
    )
    (plan.output_directory / "solution_times.csv").write_text(
        "\n".join(payload) + "\n",
        encoding="utf-8",
    )


def _comparison_execution_root(
    mapping: GecIcpMapping,
    plan: GecIcpComparisonExportPlan,
) -> Path:
    try:
        return Path(
            os.path.commonpath(
                (
                    mapping.path.resolve(),
                    plan.input_mph.resolve(),
                    plan.output_directory.resolve(),
                )
            )
        )
    except ValueError as exc:
        raise GecIcpComparisonExportError(
            "comparison mapping, MPH, and output must share an execution root"
        ) from exc


def _validate_execution_binding(
    mapping: GecIcpMapping,
    plan: GecIcpComparisonExportPlan,
) -> None:
    if not plan.manifest_path.is_file() or not plan.java_path.is_file():
        raise GecIcpComparisonExportError(
            "comparison export must be prepared before execution"
        )
    manifest = json.loads(plan.manifest_path.read_text(encoding="utf-8"))
    expected_mapping_hash = manifest.get("model_mapping", {}).get("sha256")
    current_mapping_hash = _sha256_file(mapping.path.resolve())
    if expected_mapping_hash != current_mapping_hash:
        raise GecIcpComparisonExportError(
            "GEC ICP model mapping changed after comparison preparation"
        )
    current_source_hash = _sha256_file(plan.input_mph)
    if current_source_hash != plan.input_mph_sha256:
        raise GecIcpComparisonExportError(
            "comparison input MPH changed after preparation"
        )


def _validate_completed_export(plan: GecIcpComparisonExportPlan) -> None:
    current_source_hash = _sha256_file(plan.input_mph)
    if current_source_hash != plan.input_mph_sha256:
        raise GecIcpComparisonExportError(
            "comparison input MPH changed during read-only export"
        )
    missing = [
        str(path)
        for path in plan.expected_outputs
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise GecIcpComparisonExportError(
            "COMSOL comparison export did not produce fresh nonempty outputs: "
            + ", ".join(missing)
        )


def _update_execution_manifest(
    plan: GecIcpComparisonExportPlan,
    *,
    status: str,
    execution: ComsolExecutionSummary | None,
    error: str | None,
) -> None:
    payload: dict[str, Any] = json.loads(plan.manifest_path.read_text(encoding="utf-8"))
    payload["stage"] = "execute-gec-icp-comparison-export"
    payload["status"] = status
    payload["input_mph"]["post_execution_sha256"] = _sha256_file(plan.input_mph)
    payload["outputs"] = [
        {
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256_file(path),
        }
        for path in plan.expected_outputs
        if path.is_file()
    ]
    payload["error"] = error
    payload["execution"] = (
        {
            "operation": execution.operation,
            "log_dir": str(execution.log_dir),
            "command_json": str(execution.command_json),
            "result_json": str(execution.result_json),
            "provenance_json": str(execution.provenance_json),
            "return_codes": list(execution.return_codes),
            "total_time_s": execution.total_time_s,
            "comsol_version": execution.comsol_version,
            "comsol_build": execution.comsol_build,
        }
        if execution is not None
        else None
    )
    write_json(plan.manifest_path, payload)


def generate_comparison_export_java(
    model: GecIcpModelSpec,
    *,
    input_mph: Path,
    output_directory: Path,
    common_time_s: float,
) -> str:
    """Render the read-only Java exporter for one saved solution."""

    output = output_directory.resolve()
    model_name = "swarmGecIcpComparisonExport"
    lines = [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "import java.util.Locale;",
        "",
        f"public class {COMPARISON_EXPORT_CLASS} {{",
        "  public static void main(String[] args) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string(model_name)
        + ", "
        + java_string(java_path(input_mph.resolve()))
        + ");",
        f"    double[] savedTimes = model.sol({java_string(model.solution)}).getPVals();",
        "    if (savedTimes.length == 0) {",
        '      throw new IllegalStateException("saved solution has no time values");',
        "    }",
        "    double minimumTime = savedTimes[0];",
        "    double maximumTime = savedTimes[0];",
        "    for (double time : savedTimes) {",
        "      if (!Double.isFinite(time)) {",
        '        throw new IllegalStateException("saved solution time is nonfinite");',
        "      }",
        "      minimumTime = Math.min(minimumTime, time);",
        "      maximumTime = Math.max(maximumTime, time);",
        "    }",
        f"    double commonTime = {common_time_s:.17e};",
        "    double timeTolerance = 1.0e-12*Math.max(1.0, Math.abs(maximumTime));",
        "    if (commonTime < minimumTime-timeTolerance ||",
        "        commonTime > maximumTime+timeTolerance) {",
        '      throw new IllegalStateException("common comparison time is outside saved support");',
        "    }",
        "    emitSolutionTimes(savedTimes);",
        f"    int[] plasmaDomains = model.component({java_string(model.component)})",
        f"        .physics({java_string(model.plasma_physics)}).selection().entities(2);",
        "    if (plasmaDomains.length == 0) {",
        '      throw new IllegalStateException("plasma domain selection is empty");',
        "    }",
    ]
    lines.extend(
        _field_export_lines(
            model.dataset,
            "swCmpFieldsCommon",
            output / "fields_common_time.csv",
            selection_mode="interp",
            common_time_s=common_time_s,
        )
    )
    lines.extend(
        _numerical_export_lines(
            model.dataset,
            "swCmpVolumeCommon",
            "IntSurface",
            VOLUME_EXPRESSIONS,
            VOLUME_UNITS,
            VOLUME_DESCRIPTIONS,
            output / "volume_common_time.csv",
            selection_mode="interp",
            common_time_s=common_time_s,
            selection="plasmaDomains",
            volume=True,
        )
    )
    lines.extend(
        _numerical_export_lines(
            model.dataset,
            "swCmpCoilCommon",
            "EvalGlobal",
            ("mf.PCoil_1",),
            ("W",),
            ("configured coil-power readback",),
            output / "coil_power_common_time.csv",
            selection_mode="interp",
            common_time_s=common_time_s,
        )
    )
    lines.extend(
        _field_export_lines(
            model.dataset,
            "swCmpFieldsTerminal",
            output / "fields_terminal.csv",
            selection_mode="last",
        )
    )
    lines.extend(
        _numerical_export_lines(
            model.dataset,
            "swCmpVolumeTerminal",
            "IntSurface",
            VOLUME_EXPRESSIONS,
            VOLUME_UNITS,
            VOLUME_DESCRIPTIONS,
            output / "volume_terminal.csv",
            selection_mode="last",
            selection="plasmaDomains",
            volume=True,
        )
    )
    lines.extend(
        _numerical_export_lines(
            model.dataset,
            "swCmpCoilTerminal",
            "EvalGlobal",
            ("mf.PCoil_1",),
            ("W",),
            ("configured coil-power readback",),
            output / "coil_power_terminal.csv",
            selection_mode="last",
        )
    )
    lines.extend(
        _numerical_export_lines(
            model.dataset,
            "swCmpVolumeSeries",
            "IntSurface",
            VOLUME_EXPRESSIONS,
            VOLUME_UNITS,
            VOLUME_DESCRIPTIONS,
            output / "volume_time_series.csv",
            selection_mode="all",
            selection="plasmaDomains",
            volume=True,
        )
    )
    lines.extend(
        _numerical_export_lines(
            model.dataset,
            "swCmpCoilSeries",
            "EvalGlobal",
            ("mf.PCoil_1",),
            ("W",),
            ("configured coil-power readback",),
            output / "coil_power_time_series.csv",
            selection_mode="all",
        )
    )
    lines.extend(
        [
            "    // Read-only export: no study execution, solution clearing, or model save.",
            f"    ModelUtil.remove({java_string(model_name)});",
            "  }",
            "",
            "  private static void emitSolutionTimes(double[] times) {",
            "      for (int index = 0; index < times.length; index++) {",
            '        System.out.printf(Locale.ROOT, "SWARM_GEC_ICP_SOLUTION_TIME\\t%d\\t%.17g%n",',
            "            index+1, times[index]);",
            "      }",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def _field_export_lines(
    dataset: str,
    tag: str,
    output: Path,
    *,
    selection_mode: str,
    common_time_s: float | None = None,
) -> list[str]:
    qtag = java_string(tag)
    expressions = ", ".join(java_string(value) for value in FIELD_EXPRESSIONS)
    units = ", ".join(java_string(value) for value in FIELD_UNITS)
    descriptions = ", ".join(java_string(value) for value in FIELD_DESCRIPTIONS)
    lines = [
        f"    try {{ model.result().export().remove({qtag}); }} catch (Exception ignored) {{}}",
        f'    model.result().export().create({qtag}, "Data");',
        f'    model.result().export({qtag}).set("data", {java_string(dataset)});',
        f'    model.result().export({qtag}).set("expr", new String[]{{{expressions}}});',
        f'    model.result().export({qtag}).set("unit", new String[]{{{units}}});',
        f'    model.result().export({qtag}).set("descr", new String[]{{{descriptions}}});',
        f'    model.result().export({qtag}).set("innerinput", {java_string(selection_mode)});',
    ]
    lines.extend(_time_value_lines("export", qtag, selection_mode, common_time_s))
    lines.extend(
        [
            f'    model.result().export({qtag}).set("header", "on");',
            f'    model.result().export({qtag}).set("fullprec", "on");',
            f'    model.result().export({qtag}).set("ifexists", "overwrite");',
            f'    model.result().export({qtag}).set("includecoords", true);',
            f'    model.result().export({qtag}).set("includenan", false);',
            f'    model.result().export({qtag}).set("sort", "on");',
            f'    model.result().export({qtag}).set("filename", '
            + java_string(java_path(output))
            + ");",
            f"    model.result().export({qtag}).run();",
        ]
    )
    return lines


def _numerical_export_lines(
    dataset: str,
    tag: str,
    kind: str,
    expressions: tuple[str, ...],
    units: tuple[str, ...],
    descriptions: tuple[str, ...],
    output: Path,
    *,
    selection_mode: str,
    common_time_s: float | None = None,
    selection: str | None = None,
    volume: bool = False,
) -> list[str]:
    if not (len(expressions) == len(units) == len(descriptions)):
        raise GecIcpComparisonExportError(
            "comparison expressions, units, and descriptions must align"
        )
    qtag = java_string(tag)
    table_tag = java_string(f"{tag}Table")
    export_tag = java_string(f"{tag}Export")
    expression_values = ", ".join(java_string(value) for value in expressions)
    unit_values = ", ".join(java_string(value) for value in units)
    description_values = ", ".join(java_string(value) for value in descriptions)
    lines = [
        f"    try {{ model.result().numerical().remove({qtag}); }} catch (Exception ignored) {{}}",
        f"    try {{ model.result().table().remove({table_tag}); }} catch (Exception ignored) {{}}",
        f"    try {{ model.result().export().remove({export_tag}); }} catch (Exception ignored) {{}}",
        f'    model.result().table().create({table_tag}, "Table");',
        f"    model.result().numerical().create({qtag}, {java_string(kind)});",
        f'    model.result().numerical({qtag}).set("data", {java_string(dataset)});',
        f'    model.result().numerical({qtag}).set("expr", new String[]{{{expression_values}}});',
        f'    model.result().numerical({qtag}).set("unit", new String[]{{{unit_values}}});',
        f'    model.result().numerical({qtag}).set("descr", new String[]{{{description_values}}});',
        f'    model.result().numerical({qtag}).set("innerinput", {java_string(selection_mode)});',
    ]
    lines.extend(_time_value_lines("numerical", qtag, selection_mode, common_time_s))
    if selection is not None:
        lines.append(
            f"    model.result().numerical({qtag}).selection().set({selection});"
        )
    if volume:
        lines.append(f'    model.result().numerical({qtag}).set("intvolume", true);')
    lines.extend(
        [
            f'    model.result().numerical({qtag}).set("table", {table_tag});',
            f"    model.result().numerical({qtag}).setResult();",
            f'    model.result().export().create({export_tag}, "Table");',
            f'    model.result().export({export_tag}).set("table", {table_tag});',
            f'    model.result().export({export_tag}).set("header", "on");',
            f'    model.result().export({export_tag}).set("prec", "full");',
            f'    model.result().export({export_tag}).set("ifexists", "overwrite");',
            f'    model.result().export({export_tag}).set("filename", '
            + java_string(java_path(output))
            + ");",
            f"    model.result().export({export_tag}).run();",
        ]
    )
    return lines


def _time_value_lines(
    feature_kind: str,
    tag_literal: str,
    selection_mode: str,
    common_time_s: float | None,
) -> list[str]:
    if selection_mode != "interp":
        return []
    if common_time_s is None:
        raise GecIcpComparisonExportError("interpolated export requires common_time_s")
    owner = (
        f"model.result().export({tag_literal})"
        if feature_kind == "export"
        else f"model.result().numerical({tag_literal})"
    )
    return [f'    {owner}.set("t", new double[]{{{common_time_s:.17e}}});']


def _validate_inputs(
    case_id: str,
    source: Path,
    output: Path,
    common_time_s: float,
) -> None:
    if re.fullmatch(r"[a-z][a-z0-9_]*", case_id) is None:
        raise GecIcpComparisonExportError("case_id must match [a-z][a-z0-9_]*")
    if not source.is_file() or source.suffix.lower() != ".mph":
        raise GecIcpComparisonExportError(
            f"comparison input must be an existing MPH: {source}"
        )
    if output == source or output in source.parents:
        raise GecIcpComparisonExportError(
            "comparison output directory must not be the source MPH"
        )
    if not math.isfinite(common_time_s) or common_time_s < 0.0:
        raise GecIcpComparisonExportError(
            "common comparison time must be finite and nonnegative"
        )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = (
    "COMPARISON_EXPORT_CLASS",
    "EXPECTED_EXPORT_NAMES",
    "FIELD_DESCRIPTIONS",
    "FIELD_EXPRESSIONS",
    "FIELD_UNITS",
    "GecIcpComparisonExportError",
    "GecIcpComparisonExportPlan",
    "VOLUME_DESCRIPTIONS",
    "VOLUME_EXPRESSIONS",
    "VOLUME_UNITS",
    "execute_gec_icp_comparison_export",
    "generate_comparison_export_java",
    "prepare_gec_icp_comparison_export",
)
