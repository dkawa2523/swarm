"""Validated result readers and comparison metrics for GEC CCP plots."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

from swarm_workflow.comsol.models.gec_ccp.plots import contracts as _contracts
from swarm_workflow.comsol.models.gec_ccp.validation.bundle_context import (
    GEC_OPTIONAL_PROVENANCE_HASH_KEYS,
    GEC_PROVENANCE_HASH_KEYS,
)

def _validated_plot_provenance(
    bundle: Path,
    bundle_manifest: dict[str, Any],
    results: Path,
) -> dict[str, Any]:
    plan_path = results / "gec_ccp_plan.json"
    status_path = results / "gec_ccp_run_status.json"
    if not plan_path.exists() or not status_path.exists():
        raise _contracts.GecCcpPlotError(
            "COMSOL result plotting requires gec_ccp_plan.json and "
            "gec_ccp_run_status.json from the same run"
        )
    try:
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise _contracts.GecCcpPlotError("COMSOL run provenance JSON is unreadable") from exc
    if plan.get("status") != "ready":
        raise _contracts.GecCcpPlotError("GEC CCP plan is not ready")
    plan_sha256 = hashlib.sha256(plan_path.read_bytes()).hexdigest()
    if status.get("plan_sha256") != plan_sha256:
        raise _contracts.GecCcpPlotError(
            "GEC CCP status does not belong to the supplied plan"
        )
    model = plan.get("model")
    input_mph = model.get("input_mph") if isinstance(model, dict) else None
    if (
        not isinstance(input_mph, dict)
        or status.get("input_mph_sha256") != input_mph.get("sha256")
    ):
        raise _contracts.GecCcpPlotError("GEC CCP input MPH provenance is inconsistent")
    generated = plan.get("generated_java")
    expected_java_hashes = (
        {
            name: metadata.get("sha256")
            for name, metadata in generated.items()
            if isinstance(metadata, dict)
        }
        if isinstance(generated, dict)
        else {}
    )
    if (
        not expected_java_hashes
        or status.get("generated_java_sha256") != expected_java_hashes
        or status.get("bundle_artifacts_verified") is not True
    ):
        raise _contracts.GecCcpPlotError("GEC CCP generated-input provenance is inconsistent")
    if status.get("solve_status") != "completed":
        raise _contracts.GecCcpPlotError("COMSOL result provenance is not a completed solve")
    quality_status = status.get("quality_status")
    if quality_status != "passed" or status.get("quality_accepted") is not True:
        raise _contracts.GecCcpPlotError(
            "COMSOL result failed physics-quality acceptance and cannot be "
            "plotted as an accepted result"
        )
    if status.get("status") != "completed":
        raise _contracts.GecCcpPlotError("COMSOL result provenance is not accepted")
    result_artifacts = _validated_result_artifacts(plan, status, results)
    runtime = status.get("comsol_runtime")
    if not isinstance(runtime, dict) or not all(
        isinstance(runtime.get(name), str) and runtime[name]
        for name in ("version", "build", "executable")
    ):
        raise _contracts.GecCcpPlotError("COMSOL runtime identity is missing")
    output_mph = status.get("output_mph")
    try:
        output_mph_path = Path(str(output_mph["path"])).resolve()
        output_mph_sha256 = hashlib.sha256(
            output_mph_path.read_bytes()
        ).hexdigest()
    except (KeyError, TypeError, OSError) as exc:
        raise _contracts.GecCcpPlotError("COMSOL output MPH provenance is invalid") from exc
    if output_mph.get("sha256") != output_mph_sha256:
        raise _contracts.GecCcpPlotError("COMSOL output MPH has changed since the solve")
    plan_bundle = plan.get("bundle")
    closure = plan.get("closure")
    if not isinstance(plan_bundle, dict) or not isinstance(closure, dict):
        raise _contracts.GecCcpPlotError("GEC CCP plan lacks bundle or closure metadata")
    try:
        planned_bundle = Path(str(plan_bundle["path"])).resolve()
    except (KeyError, TypeError, OSError) as exc:
        raise _contracts.GecCcpPlotError("GEC CCP plan bundle path is invalid") from exc
    if planned_bundle != bundle:
        raise _contracts.GecCcpPlotError(
            f"plot bundle differs from solved bundle: {bundle} != {planned_bundle}"
        )
    manifest_path = bundle / "manifest.json"
    if (
        not manifest_path.is_file()
        or plan_bundle.get("manifest_sha256")
        != hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    ):
        raise _contracts.GecCcpPlotError(
            "plot bundle manifest differs from the solved bundle manifest"
        )
    if plan_bundle.get("source") != bundle_manifest.get("source"):
        raise _contracts.GecCcpPlotError("plot bundle solver source differs from solved bundle")
    bundle_hashes = bundle_manifest.get("hashes")
    required_hashes = set(GEC_PROVENANCE_HASH_KEYS)
    allowed_hashes = required_hashes | set(GEC_OPTIONAL_PROVENANCE_HASH_KEYS)
    if (
        not isinstance(bundle_hashes, dict)
        or not required_hashes.issubset(bundle_hashes)
        or not set(bundle_hashes).issubset(allowed_hashes)
        or plan_bundle.get("hashes")
        != {name: bundle_hashes[name] for name in GEC_PROVENANCE_HASH_KEYS}
    ):
        raise _contracts.GecCcpPlotError("plot bundle hashes differ from solved bundle")
    quality_thresholds = bundle_manifest.get("quality_thresholds")
    if (
        not isinstance(quality_thresholds, dict)
        or plan_bundle.get("quality_thresholds") != quality_thresholds
    ):
        raise _contracts.GecCcpPlotError(
            "plot bundle quality policy differs from solved bundle"
        )
    electron_transport = closure.get("electron_transport")
    reaction_model = closure.get("reaction_model")
    if electron_transport not in {
        "comsol",
        "swarm_mobility_einstein",
        _contracts.GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
        _contracts.CANONICAL_RESTRICTED_TRANSPORT_CLOSURE,
    } or reaction_model not in {
        "comsol_eedf",
        "external_rates",
        *_contracts.FUNCTION_EEDF_REACTION_MODELS,
    }:
        raise _contracts.GecCcpPlotError("GEC CCP plan has invalid closure axes")
    source_names = {
        "two_term": "two-term Swarm",
        "monte_carlo": "Monte Carlo Swarm",
        "multi_term": "multi-term Swarm",
    }
    source_label = source_names.get(
        str(plan_bundle.get("source")), str(plan_bundle.get("source"))
    )
    transport_labels = {
        "comsol": "COMSOL transport",
        "swarm_mobility_einstein": "Swarm mobility + COMSOL Einstein diffusion",
        _contracts.GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE: (
            "Swarm mobility/energy transport + COMSOL Einstein particle "
            "diffusion"
        ),
        _contracts.CANONICAL_RESTRICTED_TRANSPORT_CLOSURE: (
            "restricted COMSOL Specify All particle/energy transport"
        ),
    }
    reaction_labels = {
        "comsol_eedf": "COMSOL EEDF reactions",
        "external_rates": "imported Swarm rates",
        "function_eedf": "Function-EEDF + COMSOL cross-section integration",
        "function_eedf_preintegrated_inelastic": (
            "elastic Function-EEDF cross-section + excitation/ionization "
            "preintegrated Swarm rates"
        ),
    }
    external_label = (
        f"{source_label}: {transport_labels[electron_transport]}; "
        f"{reaction_labels[reaction_model]}"
    )
    return {
        "plan": str(plan_path),
        "status": str(status_path),
        "quality_status": quality_status,
        "result_role": status.get("result_role"),
        "planned_result_role": plan.get("result_role"),
        "physical_target_accepted": status.get("physical_target_accepted"),
        "plan_sha256": plan_sha256,
        "bundle": plan_bundle,
        "closure": {
            "electron_transport": electron_transport,
            "reaction_model": reaction_model,
        },
        "include_builtin_reference": bool(
            isinstance(plan.get("solve"), dict)
            and plan["solve"].get("include_builtin_reference") is True
        ),
        "results": result_artifacts,
        "external_label": external_label,
    }


def _validated_result_artifacts(
    plan: dict[str, Any],
    status: dict[str, Any],
    results: Path,
) -> list[dict[str, Any]]:
    expected_raw = plan.get("expected_results")
    artifacts_raw = status.get("results")
    if not isinstance(expected_raw, list) or not expected_raw:
        raise _contracts.GecCcpPlotError("GEC CCP plan lacks its expected result set")
    if not isinstance(artifacts_raw, list) or not artifacts_raw:
        raise _contracts.GecCcpPlotError(
            "COMSOL status lacks canonical result-file provenance"
        )

    root = results.resolve()
    expected: set[str] = set()
    for raw_path in expected_raw:
        if not isinstance(raw_path, str) or not raw_path:
            raise _contracts.GecCcpPlotError("GEC CCP expected result set is invalid")
        candidate = Path(raw_path)
        resolved = (
            candidate.resolve()
            if candidate.is_absolute()
            else (root / candidate).resolve()
        )
        try:
            relative = resolved.relative_to(root).as_posix()
        except ValueError as exc:
            raise _contracts.GecCcpPlotError(
                "GEC CCP expected result is outside the result directory"
            ) from exc
        if relative in expected:
            raise _contracts.GecCcpPlotError("GEC CCP expected result set contains duplicates")
        expected.add(relative)

    artifacts: dict[str, dict[str, Any]] = {}
    for metadata in artifacts_raw:
        if not isinstance(metadata, dict) or set(metadata) != {
            "path",
            "size_bytes",
            "sha256",
        }:
            raise _contracts.GecCcpPlotError(
                "COMSOL result-file provenance entry is invalid"
            )
        relative = metadata["path"]
        size_bytes = metadata["size_bytes"]
        digest = metadata["sha256"]
        if (
            not isinstance(relative, str)
            or not relative
            or not isinstance(size_bytes, int)
            or isinstance(size_bytes, bool)
            or size_bytes < 0
            or not isinstance(digest, str)
            or re.fullmatch(r"[0-9a-f]{64}", digest) is None
        ):
            raise _contracts.GecCcpPlotError(
                "COMSOL result-file provenance entry is invalid"
            )
        resolved = (root / Path(relative)).resolve()
        try:
            canonical_relative = resolved.relative_to(root).as_posix()
        except ValueError as exc:
            raise _contracts.GecCcpPlotError(
                "COMSOL result-file provenance escapes the result directory"
            ) from exc
        if relative != canonical_relative or relative in artifacts:
            raise _contracts.GecCcpPlotError(
                "COMSOL result-file provenance path is not canonical"
            )
        artifacts[relative] = metadata

    if set(artifacts) != expected:
        raise _contracts.GecCcpPlotError(
            "COMSOL result-file provenance set differs from the plan"
        )
    for relative in sorted(expected):
        path = root / Path(relative)
        try:
            data = path.read_bytes()
        except OSError as exc:
            raise _contracts.GecCcpPlotError(
                f"COMSOL result file is missing or unreadable: {relative}"
            ) from exc
        metadata = artifacts[relative]
        if (
            len(data) != metadata["size_bytes"]
            or hashlib.sha256(data).hexdigest() != metadata["sha256"]
        ):
            raise _contracts.GecCcpPlotError(
                f"COMSOL result file has changed since the solve: {relative}"
            )
    return [artifacts[path] for path in sorted(artifacts)]


def _comparison_input_mph_sha256(provenance: dict[str, Any]) -> str:
    try:
        plan = json.loads(Path(provenance["plan"]).read_text(encoding="utf-8"))
        value = plan["model"]["input_mph"]["sha256"]
    except (KeyError, TypeError, OSError, json.JSONDecodeError) as exc:
        raise _contracts.GecCcpPlotError("comparison input MPH provenance is invalid") from exc
    if not isinstance(value, str) or len(value) != 64:
        raise _contracts.GecCcpPlotError("comparison input MPH hash is invalid")
    return value


def _save(
    fig: Any,
    stem: Path,
    plt: Any,
    *,
    include_svg: bool = True,
) -> list[Path]:
    paths = [stem.with_suffix(".png")]
    fig.savefig(paths[0], dpi=180, bbox_inches="tight")
    if include_svg:
        paths.append(stem.with_suffix(".svg"))
        fig.savefig(paths[-1], bbox_inches="tight")
    plt.close(fig)
    return paths


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise _contracts.GecCcpPlotError(f"missing plot input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _numeric_rows(path: Path) -> list[dict[str, float]]:
    rows = _rows(path)
    parsed: list[dict[str, float]] = []
    for row in rows:
        parsed_row: dict[str, float] = {}
        for key, value in row.items():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number):
                parsed_row[key] = number
        parsed.append(parsed_row)
    return parsed


def _comsol_numeric_matrix(path: Path) -> Any:
    try:
        import numpy as np
    except ImportError as exc:
        raise _contracts.GecCcpPlotError("numpy is required") from exc
    lines = [
        line
        for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("%")
    ]
    if not lines:
        raise _contracts.GecCcpPlotError(f"COMSOL CSV is empty: {path}")
    delimiter = "," if "," in lines[0] else None
    data = np.genfromtxt(lines, delimiter=delimiter, invalid_raise=False)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    data = data[~np.isnan(data).all(axis=1)]
    return data


def _comsol_spatial_data(path: Path) -> tuple[Any, Any]:
    """Return the varying spatial coordinate and result columns.

    COMSOL exports a 2D cut line as ``R,Z,expressions...`` even when one
    coordinate is constant.  Select and sort by the coordinate with the
    largest span instead of assuming that the first numeric column is the
    line coordinate.
    """
    try:
        import numpy as np
    except ImportError as exc:
        raise _contracts.GecCcpPlotError("numpy is required") from exc
    data = _comsol_numeric_matrix(path)
    dimension = _comsol_export_dimension(path)
    if dimension < 1 or data.shape[1] <= dimension:
        raise _contracts.GecCcpPlotError(
            f"COMSOL spatial CSV has invalid dimension metadata: {path}"
        )
    coordinate_columns = data[:, :dimension]
    spans = np.ptp(coordinate_columns, axis=0)
    coordinate_index = int(np.argmax(spans))
    coordinate = coordinate_columns[:, coordinate_index]
    order = np.argsort(coordinate, kind="stable")
    return coordinate[order], data[order, dimension:]


def _comsol_export_dimension(path: Path) -> int:
    for line in path.read_text(
        encoding="utf-8-sig", errors="replace"
    ).splitlines():
        match = re.match(r"\s*%\s*Dimension\s*,\s*(\d+)\s*$", line)
        if match:
            return int(match.group(1))
    return 1


def _spatial_comparison_metrics(
    baseline_path: Path,
    external_path: Path,
) -> dict[str, Any]:
    import numpy as np

    baseline_coordinate, baseline = _comsol_spatial_data(baseline_path)
    external_coordinate, external = _comsol_spatial_data(external_path)
    field_names = (
        "electron_density_m3",
        "electron_temperature_eV",
        "potential_V",
        "ionization_source_m3_s",
        "absorbed_power_W_m3",
    )
    count = min(baseline.shape[1], external.shape[1], len(field_names))
    low = max(float(baseline_coordinate.min()), float(external_coordinate.min()))
    high = min(float(baseline_coordinate.max()), float(external_coordinate.max()))
    mask = (baseline_coordinate >= low) & (baseline_coordinate <= high)
    coordinate = baseline_coordinate[mask]
    metrics: dict[str, Any] = {
        "coordinate_overlap_m": [low, high],
        "samples": int(coordinate.size),
        "fields": {},
    }
    for index, name in enumerate(field_names[:count]):
        reference = baseline[mask, index]
        candidate = np.interp(
            coordinate,
            external_coordinate,
            external[:, index],
        )
        reference_norm = float(np.linalg.norm(reference))
        reference_abs_peak = float(np.max(np.abs(reference)))
        candidate_abs_peak = float(np.max(np.abs(candidate)))
        metrics["fields"][name] = {
            "baseline_min": float(np.min(reference)),
            "baseline_max": float(np.max(reference)),
            "external_min": float(np.min(candidate)),
            "external_max": float(np.max(candidate)),
            "external_to_baseline_abs_peak_ratio": (
                candidate_abs_peak / reference_abs_peak
                if reference_abs_peak > 0.0
                else None
            ),
            "relative_l2_difference": (
                float(np.linalg.norm(candidate - reference)) / reference_norm
                if reference_norm > 0.0
                else None
            ),
        }
    return metrics


def _waveform_comparison_metrics(
    baseline_path: Path,
    external_path: Path,
) -> dict[str, Any]:
    import numpy as np

    waveforms: dict[str, dict[str, float]] = {}
    for name, path in (("baseline", baseline_path), ("external", external_path)):
        phase, voltage, current = _comsol_waveform_data(path)
        order = np.argsort(phase, kind="stable")
        phase, voltage, current = phase[order], voltage[order], current[order]
        if (
            phase.size < 2
            or not np.all(np.isfinite((phase, voltage, current)))
            or np.any(np.diff(phase) <= 0.0)
        ):
            raise _contracts.GecCcpPlotError(
                "COMSOL waveform metrics require finite values at distinct RF phases"
            )
        # Integrate on each export's own grid. Equal sample weights double-count
        # periodic endpoints and bias nonuniform grids; resampling can lose peaks.
        span = float(phase[-1] - phase[0])
        waveforms[name] = {
            "voltage_peak_to_peak_V": float(np.ptp(voltage)),
            "voltage_mean_V": float(np.trapezoid(voltage, phase) / span),
            "current_rms_A": float(
                np.sqrt(np.trapezoid(np.square(current), phase) / span)
            ),
        }
    baseline_voltage_pp = waveforms["baseline"]["voltage_peak_to_peak_V"]
    external_voltage_pp = waveforms["external"]["voltage_peak_to_peak_V"]
    baseline_current_rms = waveforms["baseline"]["current_rms_A"]
    external_current_rms = waveforms["external"]["current_rms_A"]
    return {
        "averaging_method": "trapezoidal_over_each_exported_rf_period",
        "voltage_peak_to_peak_V": {
            "baseline": baseline_voltage_pp,
            "external": external_voltage_pp,
            "external_to_baseline_ratio": (
                external_voltage_pp / baseline_voltage_pp
                if baseline_voltage_pp > 0.0
                else None
            ),
        },
        "voltage_mean_V": {
            "baseline": waveforms["baseline"]["voltage_mean_V"],
            "external": waveforms["external"]["voltage_mean_V"],
        },
        "current_rms_A": {
            "baseline": baseline_current_rms,
            "external": external_current_rms,
            "external_to_baseline_ratio": (
                external_current_rms / baseline_current_rms
                if baseline_current_rms > 0.0
                else None
            ),
        },
    }


def _comsol_waveform_data(path: Path) -> tuple[Any, Any, Any]:
    """Return normalized RF phase, terminal voltage, and terminal current.

    A converted COMSOL time dataset is exported in row-wise form.  A native
    periodic dataset is exported in wide form, with one V/I column pair per
    phase and the global terminal values repeated at every spatial node.
    Supporting both forms avoids a redundant periodic-to-time conversion.
    """
    import numpy as np

    data = _comsol_numeric_matrix(path)
    if data.shape[1] < 3:
        raise _contracts.GecCcpPlotError(
            "COMSOL waveform CSV needs coordinate, voltage, current"
        )
    names = _comsol_column_names(path, data.shape[1])
    voltage_indices = [
        index for index, name in enumerate(names)
        if name.startswith("ptp.mct1.V")
    ]
    current_indices = [
        index for index, name in enumerate(names)
        if name.startswith("ptp.mct1.I")
    ]
    if len(voltage_indices) == len(current_indices) == 1:
        return (
            _normalized_coordinate(data[:, 0]),
            data[:, voltage_indices[0]],
            data[:, current_indices[0]],
        )
    if len(voltage_indices) < 2 or len(voltage_indices) != len(current_indices):
        raise _contracts.GecCcpPlotError(
            "COMSOL waveform CSV has unmatched phase-resolved V/I columns"
        )

    def phase_columns(indices: list[int], expression: str) -> dict[float, int]:
        columns: dict[float, int] = {}
        for index in indices:
            match = re.search(r"\s@\st=([^\s]+)\s*$", names[index])
            if match is None:
                raise _contracts.GecCcpPlotError(
                    f"COMSOL wide waveform column lacks phase time: {names[index]}"
                )
            try:
                time_s = float(match.group(1))
            except ValueError as exc:
                raise _contracts.GecCcpPlotError(
                    f"COMSOL wide waveform has invalid phase time: {names[index]}"
                ) from exc
            if time_s in columns:
                raise _contracts.GecCcpPlotError(
                    f"COMSOL wide waveform repeats {expression} at t={time_s}"
                )
            columns[time_s] = index
        return columns

    voltage_columns = phase_columns(voltage_indices, "voltage")
    current_columns = phase_columns(current_indices, "current")
    if voltage_columns.keys() != current_columns.keys():
        raise _contracts.GecCcpPlotError(
            "COMSOL wide waveform V/I columns use different phase times"
        )
    times = np.asarray(sorted(voltage_columns))
    voltage = np.asarray([data[0, voltage_columns[time]] for time in times])
    current = np.asarray([data[0, current_columns[time]] for time in times])
    for time, index in voltage_columns.items():
        if not np.allclose(data[:, index], data[0, index], rtol=1.0e-10, atol=1.0e-12):
            raise _contracts.GecCcpPlotError(
                f"COMSOL terminal voltage is not global at t={time}"
            )
    for time, index in current_columns.items():
        if not np.allclose(data[:, index], data[0, index], rtol=1.0e-10, atol=1.0e-12):
            raise _contracts.GecCcpPlotError(
                f"COMSOL terminal current is not global at t={time}"
            )
    return _normalized_coordinate(times), voltage, current


def _comsol_radial_midplane_data(
    path: Path,
    expected_center_height_cm: float,
) -> tuple[Any, Any, float]:
    """Return a sorted horizontal COMSOL cut at the expected gap midpoint."""
    import numpy as np

    data = _comsol_numeric_matrix(path)
    dimension = _comsol_export_dimension(path)
    if dimension != 2 or data.shape[1] <= dimension:
        raise _contracts.GecCcpPlotError(
            f"COMSOL radial-midplane CSV is not a 2D cut export: {path}"
        )
    coordinates = np.asarray(data[:, :2], dtype=float)
    fields = np.asarray(data[:, dimension:], dtype=float)
    valid_coordinates = np.isfinite(coordinates).all(axis=1)
    coordinates = coordinates[valid_coordinates]
    fields = fields[valid_coordinates, :]
    if coordinates.shape[0] < 2:
        raise _contracts.GecCcpPlotError(f"COMSOL radial cut has too few points: {path}")
    radius_cm = coordinates[:, 0] * 100.0
    axial_cm = coordinates[:, 1] * 100.0
    axial_span = float(np.ptp(axial_cm))
    if axial_span > 1.0e-7:
        raise _contracts.GecCcpPlotError(
            f"COMSOL radial cut is not at constant height: {path}"
        )
    actual_height_cm = float(np.mean(axial_cm))
    if not math.isclose(
        actual_height_cm,
        expected_center_height_cm,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        raise _contracts.GecCcpPlotError(
            "COMSOL radial cut is not at the GEC gap midpoint: "
            f"expected {expected_center_height_cm:.6g} cm, "
            f"found {actual_height_cm:.6g} cm in {path}"
        )
    if float(np.ptp(radius_cm)) <= 0.0:
        raise _contracts.GecCcpPlotError(f"COMSOL radial cut has no radial span: {path}")
    order = np.argsort(radius_cm, kind="stable")
    return radius_cm[order], fields[order, :], actual_height_cm


def _comsol_domain_data(path: Path) -> tuple[Any, Any]:
    import numpy as np

    data = _comsol_numeric_matrix(path)
    dimension = _comsol_export_dimension(path)
    if dimension != 2 or data.shape[1] <= dimension:
        raise _contracts.GecCcpPlotError(f"COMSOL domain CSV is not a 2D field export: {path}")
    coordinates = np.asarray(data[:, :2], dtype=float)
    fields = np.asarray(data[:, dimension:], dtype=float)
    valid = np.isfinite(coordinates).all(axis=1)
    return coordinates[valid], fields[valid, :]


def _comsol_field_matrix(path: Path, expression: str) -> Any:
    data = _comsol_numeric_matrix(path)
    names = _comsol_column_names(path, data.shape[1])
    indices = [
        index
        for index, name in enumerate(names)
        if name == expression or name.startswith(f"{expression} (")
    ]
    if not indices:
        raise _contracts.GecCcpPlotError(
            f"COMSOL CSV does not contain {expression}: {path}"
        )
    return data[:, indices]


def _comsol_sorted_phase_field(path: Path, expression: str) -> tuple[Any, Any]:
    """Return a phase-resolved field sorted along its varying cut coordinate."""
    import numpy as np

    data = _comsol_numeric_matrix(path)
    dimension = _comsol_export_dimension(path)
    if dimension < 1 or data.shape[1] <= dimension:
        raise _contracts.GecCcpPlotError(
            f"COMSOL phase CSV has invalid dimension metadata: {path}"
        )
    spans = np.ptp(data[:, :dimension], axis=0)
    coordinate = data[:, int(np.argmax(spans))]
    order = np.argsort(coordinate, kind="stable")
    field = _comsol_field_matrix(path, expression)
    return coordinate[order], field[order, :]


def _comsol_column_names(path: Path, expected_columns: int) -> list[str]:
    candidates: list[list[str]] = []
    for line in path.read_text(
        encoding="utf-8-sig", errors="replace"
    ).splitlines():
        if not line.lstrip().startswith("%"):
            continue
        payload = line.lstrip()[1:].strip()
        fields = next(csv.reader([payload]))
        if len(fields) == expected_columns:
            candidates.append([field.strip() for field in fields])
    if not candidates:
        raise _contracts.GecCcpPlotError(f"COMSOL CSV column header is absent: {path}")
    return candidates[-1]


def _nearest(values: list[float], target: float) -> float:
    return min(values, key=lambda value: abs(math.log(value / target)))


def _positive(value: object) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number > 0.0


def _normalized_coordinate(values: Any) -> Any:
    low = values.min()
    span = values.max() - low
    return (values - low) / span if span > 0.0 else values * 0.0
