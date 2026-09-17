"""Immutable solver-qualification evidence used by downstream workflows."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import sqlite3
from typing import Any, Mapping

from .propagator_source import (
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
    PROPAGATOR_SOURCE_FINGERPRINT_SCHEMA,
    propagator_qualification_source_fingerprint,
)


PROPAGATOR_CORE_QUALIFICATION_SCHEMA = (
    "swarm.propagator_p1_deterministic_qualification.v1"
)
PROPAGATOR_CORE_QUALIFICATION_FILE = "propagator_p1_core_qualification.json"
PROPAGATOR_TARGET_QUALIFICATION_SCHEMA = "swarm.propagator_target_qualification.v1"
PROPAGATOR_TARGET_QUALIFICATION_FILE = "propagator_target_qualification.json"
PROPAGATOR_TARGET_QUALIFIER_FINGERPRINT_SCHEMA = (
    "swarm.propagator_target_qualifier_source_tree_sha256.v1"
)
PROPAGATOR_TARGET_QUALIFIER_FILES = (
    "swarm_workflow/quality/solver.py",
    "tools/qualify_propagator_target.py",
)


class SolverQualificationError(ValueError):
    """Raised when solver evidence is stale, incomplete, or altered."""


@dataclass(frozen=True, slots=True)
class SolverQualification:
    """A validated immutable qualification artifact."""

    source_path: Path
    sha256: str
    schema: str
    decision: str
    implementation_fingerprint_schema: str
    implementation_fingerprint_sha256: str

    def identity_entry(self) -> dict[str, Any]:
        return {
            "sha256": self.sha256,
            "schema": self.schema,
            "decision": self.decision,
            "implementation_fingerprint": {
                "schema": self.implementation_fingerprint_schema,
                "sha256": self.implementation_fingerprint_sha256,
            },
        }

    def manifest_entry(self, *, file: str) -> dict[str, Any]:
        return {"file": file, **self.identity_entry()}


@dataclass(frozen=True, slots=True)
class PropagatorTargetQualification:
    """Validated refinement evidence for a bounded Propagator target range."""

    source_path: Path
    sha256: str
    schema: str
    decision: str
    core_qualification_sha256: str
    medium_evidence_sha256: str
    target: str
    fields_Td: tuple[float, ...]
    operating_range_bracket_Td: tuple[float, float]
    table_support_cap_Td: float
    mixture_id: int
    medium_grid: tuple[int, int]
    fine_grid: tuple[int, int]

    def manifest_entry(self, *, file: str) -> dict[str, Any]:
        return {
            "file": file,
            "sha256": self.sha256,
            "schema": self.schema,
            "decision": self.decision,
            "core_qualification_sha256": self.core_qualification_sha256,
            "medium_evidence_sha256": self.medium_evidence_sha256,
            "fields_Td": list(self.fields_Td),
            "medium_grid": list(self.medium_grid),
            "fine_grid": list(self.fine_grid),
            "target": self.target,
            "operating_range_bracket_Td": list(self.operating_range_bracket_Td),
            "table_support_cap_Td": self.table_support_cap_Td,
            "mixture_id": self.mixture_id,
        }


@dataclass(frozen=True, slots=True)
class PropagatorTargetRequirement:
    """Downstream-owned scope required from generic target evidence."""

    target: str
    fields_Td: tuple[float, ...]
    operating_range_bracket_Td: tuple[float, float]
    table_support_cap_Td: float
    mixture_id: int
    medium_grid: tuple[int, int]
    fine_grid: tuple[int, int]
    max_scalar_refinement_limit: float
    max_eedf_weighted_L1_limit: float


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def propagator_target_qualifier_fingerprint(
    repository_root: str | Path,
) -> dict[str, Any]:
    """Hash the generic target-refinement algorithm and its validator."""

    root = Path(repository_root).resolve()
    files = list(PROPAGATOR_TARGET_QUALIFIER_FILES)
    return {
        "schema": PROPAGATOR_TARGET_QUALIFIER_FINGERPRINT_SCHEMA,
        "sha256": _source_tree_digest(root, files),
        "file_count": len(files),
        "files": files,
    }


def validate_propagator_core_qualification(
    path: str | Path,
    *,
    repository_root: str | Path | None = None,
) -> SolverQualification:
    """Validate a passing P1 artifact against the current source and inputs."""

    artifact_path = Path(path).resolve()
    try:
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SolverQualificationError(
            f"cannot read Propagator qualification artifact: {artifact_path}"
        ) from exc
    if not isinstance(payload, dict):
        raise SolverQualificationError(
            "Propagator qualification artifact must be a JSON object"
        )
    if payload.get("schema") != PROPAGATOR_CORE_QUALIFICATION_SCHEMA:
        raise SolverQualificationError("unsupported Propagator qualification schema")
    decision = payload.get("decision")
    if not isinstance(decision, dict):
        raise SolverQualificationError(
            "Propagator qualification artifact lacks its decision"
        )
    if decision.get("p1_deterministic_core_qualified") is not True:
        raise SolverQualificationError(
            "Propagator deterministic P1 core is not qualified"
        )
    blockers = decision.get("blocking_gates")
    if blockers != []:
        raise SolverQualificationError(
            "Propagator qualification contains blocking gates"
        )
    scope = payload.get("scope")
    if not isinstance(scope, dict) or (
        scope.get("phase") != "P0_P1"
        or scope.get("quick") is not False
        or scope.get("included_solver") != "propagator"
        or scope.get("field_model") != "homogeneous_dc_B0"
    ):
        raise SolverQualificationError(
            "Propagator qualification is not a full deterministic P0/P1 run"
        )
    component_decisions = (
        "case_quality_passed",
        "medium_fine_refinement_passed",
        "mixed_gas_angular_refinement_passed",
        "energy_ceiling_independence_passed",
        "inelastic_weak_transfer_passed",
    )
    if any(decision.get(name) is not True for name in component_decisions):
        raise SolverQualificationError(
            "Propagator qualification component decision failed"
        )

    environment = payload.get("environment")
    fingerprint = (
        environment.get("implementation_fingerprint")
        if isinstance(environment, dict)
        else None
    )
    if not isinstance(fingerprint, dict):
        raise SolverQualificationError(
            "Propagator qualification lacks an implementation fingerprint"
        )
    expected_digest = fingerprint.get("sha256")
    fingerprint_schema = fingerprint.get("schema")
    if (
        fingerprint_schema != PROPAGATOR_SOURCE_FINGERPRINT_SCHEMA
        or not isinstance(expected_digest, str)
    ):
        raise SolverQualificationError(
            "Propagator implementation fingerprint is malformed"
        )
    root = (
        Path(repository_root).resolve()
        if repository_root is not None
        else Path(__file__).resolve().parents[2]
    )
    actual_fingerprint = propagator_qualification_source_fingerprint(root)
    if fingerprint != actual_fingerprint:
        raise SolverQualificationError(
            "Propagator implementation changed after qualification: "
            f"expected {expected_digest}, found {actual_fingerprint['sha256']}"
        )

    input_hashes = environment.get("input_sha256")
    if not isinstance(input_hashes, dict) or not input_hashes:
        raise SolverQualificationError("Propagator qualification lacks input hashes")
    _validate_input_hashes(root, input_hashes)
    return SolverQualification(
        source_path=artifact_path,
        sha256=file_sha256(artifact_path),
        schema=PROPAGATOR_CORE_QUALIFICATION_SCHEMA,
        decision="p1_deterministic_core_qualified",
        implementation_fingerprint_schema=fingerprint_schema,
        implementation_fingerprint_sha256=expected_digest,
    )


def validate_propagator_target_qualification(
    path: str | Path,
    *,
    core_qualification: SolverQualification,
    database_path: str | Path | None = None,
    provenance_hashes: Mapping[str, Any] | None = None,
    requirement: PropagatorTargetRequirement | None = None,
    allow_self_described_generic: bool = False,
) -> PropagatorTargetQualification:
    """Validate independent refinement evidence for one target range.

    ``allow_self_described_generic`` is limited to copying already validated
    evidence into a generic table bundle.  A consuming model adapter must
    still provide its own explicit ``requirement``.
    """

    artifact_path = Path(path).resolve()
    try:
        payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SolverQualificationError(
            f"cannot read Propagator target qualification artifact: {artifact_path}"
        ) from exc
    if not isinstance(payload, dict):
        raise SolverQualificationError(
            "Propagator target qualification artifact must be a JSON object"
        )
    schema = payload.get("schema")
    if schema != PROPAGATOR_TARGET_QUALIFICATION_SCHEMA:
        raise SolverQualificationError(
            "unsupported Propagator target qualification schema"
        )
    environment = payload.get("environment")
    target_fingerprint = (
        environment.get("target_qualifier_fingerprint")
        if isinstance(environment, dict)
        else None
    )
    current_target_fingerprint = propagator_target_qualifier_fingerprint(
        Path(__file__).resolve().parents[2]
    )
    if target_fingerprint != current_target_fingerprint:
        raise SolverQualificationError(
            "Propagator target qualification algorithm changed after qualification"
        )

    scope = payload.get("scope")
    if not isinstance(scope, dict):
        raise SolverQualificationError(
            "Propagator target qualification lacks its scope"
        )
    target = scope.get("target")
    if not isinstance(target, str) or not target:
        raise SolverQualificationError(
            "Propagator target qualification has malformed target"
        )
    if scope.get("included_solver") != "propagator":
        raise SolverQualificationError(
            "Propagator target qualification has an unexpected solver"
        )
    if scope.get("field_model") != "homogeneous_dc_B0":
        raise SolverQualificationError(
            "Propagator target qualification has an unexpected field model"
        )
    fields = _strictly_increasing_positive_tuple(
        scope.get("fields_Td"), "target fields"
    )
    bracket_values = _strictly_increasing_positive_tuple(
        scope.get("operating_range_bracket_Td"),
        "operating range bracket",
        allow_equal=True,
    )
    if len(bracket_values) != 2 or any(value not in fields for value in bracket_values):
        raise SolverQualificationError(
            "Propagator target qualification operating bracket is inconsistent"
        )
    bracket = (bracket_values[0], bracket_values[1])
    table_support_cap = _positive_float(
        scope.get("table_support_cap_Td"), "table support cap"
    )
    if table_support_cap not in fields or table_support_cap < bracket[1]:
        raise SolverQualificationError(
            "Propagator target qualification support cap is inconsistent"
        )
    medium_grid = _integer_pair(scope.get("medium_grid"), "medium grid")
    fine_grid = _integer_pair(scope.get("fine_grid"), "fine grid")
    if (
        min(*medium_grid, *fine_grid) < 1
        or medium_grid[1] % 2
        or fine_grid[1] % 2
        or fine_grid[0] < medium_grid[0]
        or fine_grid[1] < medium_grid[1]
        or fine_grid == medium_grid
    ):
        raise SolverQualificationError(
            "Propagator target qualification uses an unexpected refinement pair"
        )
    mixture_id = scope.get("mixture_id")
    if (
        isinstance(mixture_id, bool)
        or not isinstance(mixture_id, int)
        or mixture_id < 0
    ):
        raise SolverQualificationError(
            "Propagator target qualification has malformed mixture id"
        )
    scalar_limit = _positive_float(
        scope.get("scalar_refinement_limit"), "scalar refinement limit"
    )
    eedf_limit = _positive_float(
        scope.get("eedf_weighted_L1_limit"), "EEDF refinement limit"
    )
    if requirement is None and not allow_self_described_generic:
        raise SolverQualificationError(
            "Propagator target qualification requires an explicit target requirement"
        )
    if requirement is not None:
        actual_scope = (
            target,
            fields,
            bracket,
            table_support_cap,
            mixture_id,
            medium_grid,
            fine_grid,
        )
        required_scope = (
            requirement.target,
            requirement.fields_Td,
            requirement.operating_range_bracket_Td,
            requirement.table_support_cap_Td,
            requirement.mixture_id,
            requirement.medium_grid,
            requirement.fine_grid,
        )
        if actual_scope != required_scope:
            raise SolverQualificationError(
                "Propagator target qualification does not match the required scope"
            )
        if (
            scalar_limit > requirement.max_scalar_refinement_limit
            or eedf_limit > requirement.max_eedf_weighted_L1_limit
        ):
            raise SolverQualificationError(
                "Propagator target qualification tolerances exceed the required scope"
            )

    decision = payload.get("decision")
    if not isinstance(decision, dict):
        raise SolverQualificationError(
            "Propagator target qualification lacks its decision"
        )
    if decision.get("target_refinement_qualified") is not True:
        raise SolverQualificationError("Propagator target refinement is not qualified")
    if (
        decision.get("fine_case_quality_passed") is not True
        or decision.get("medium_fine_refinement_passed") is not True
    ):
        raise SolverQualificationError(
            "Propagator target qualification component decision failed"
        )
    if decision.get("blocking_gates") != []:
        raise SolverQualificationError(
            "Propagator target qualification contains blocking gates"
        )

    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        raise SolverQualificationError(
            "Propagator target qualification lacks its inputs"
        )
    core_entry = inputs.get("core_qualification")
    expected_core_entry = core_qualification.identity_entry()
    if core_entry != expected_core_entry:
        raise SolverQualificationError(
            "Propagator target qualification is not bound to the supplied core evidence"
        )
    if (
        inputs.get(PROPAGATOR_SOLVER_SOURCE_METADATA_KEY)
        != core_qualification.implementation_fingerprint_sha256
    ):
        raise SolverQualificationError(
            "Propagator target medium-grid source differs from the qualified core"
        )
    medium_digest = inputs.get("medium_evidence_sha256")
    if not isinstance(medium_digest, str) or len(medium_digest) != 64:
        raise SolverQualificationError(
            "Propagator target qualification lacks its medium-grid evidence digest"
        )
    if provenance_hashes is not None:
        if not isinstance(provenance_hashes, Mapping):
            raise SolverQualificationError(
                "downstream Propagator provenance hashes are malformed"
            )
        for name in (
            "base_config_sha256",
            "workflow_config_sha256",
            "cross_sections_sha256",
            PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
        ):
            if inputs.get(name) != provenance_hashes.get(name):
                raise SolverQualificationError(
                    "Propagator target qualification provenance differs from "
                    f"the downstream manifest: {name}"
                )

    fine_runs = payload.get("fine_runs")
    refinement = payload.get("medium_fine_refinement")
    _validate_target_fine_runs(fine_runs, fields, fine_grid)
    _validate_target_refinement(
        refinement,
        fields,
        scalar_limit=scalar_limit,
        eedf_limit=eedf_limit,
    )
    if database_path is not None:
        actual_digest, metadata = propagator_target_medium_evidence(
            database_path,
            fields,
            expected_grid=medium_grid,
            mixture_id=mixture_id,
        )
        if actual_digest != medium_digest:
            raise SolverQualificationError(
                "Propagator target medium-grid evidence changed after qualification"
            )
        hash_pairs = {
            "base_config_sha256": "base_config_sha256",
            "workflow_config_sha256": "workflow_config_sha256",
            "cross_sections_sha256": "cross_sections_sha256",
            PROPAGATOR_SOLVER_SOURCE_METADATA_KEY: (
                PROPAGATOR_SOLVER_SOURCE_METADATA_KEY
            ),
        }
        for input_name, metadata_name in hash_pairs.items():
            recorded = inputs.get(input_name)
            current = metadata.get(metadata_name)
            if not isinstance(recorded, str) or recorded != current:
                raise SolverQualificationError(
                    "Propagator target qualification provenance differs from "
                    f"the workflow database: {input_name}"
                )

    return PropagatorTargetQualification(
        source_path=artifact_path,
        sha256=file_sha256(artifact_path),
        schema=schema,
        decision="target_refinement_qualified",
        core_qualification_sha256=core_qualification.sha256,
        medium_evidence_sha256=medium_digest,
        target=target,
        fields_Td=fields,
        operating_range_bracket_Td=bracket,
        table_support_cap_Td=table_support_cap,
        mixture_id=mixture_id,
        medium_grid=medium_grid,
        fine_grid=fine_grid,
    )


def propagator_target_medium_evidence(
    database_path: str | Path,
    fields_Td: tuple[float, ...],
    *,
    expected_grid: tuple[int, int] | None = None,
    mixture_id: int = 0,
) -> tuple[str, dict[str, str]]:
    """Digest only the qualified medium-grid values used by the target gate."""

    rows: list[dict[str, Any]] = []
    try:
        connection = sqlite3.connect(Path(database_path))
        connection.row_factory = sqlite3.Row
        metadata = {
            str(row[0]): str(row[1])
            for row in connection.execute(
                "SELECT key, value FROM metadata ORDER BY key"
            ).fetchall()
        }
        for field in fields_Td:
            case = connection.execute(
                "SELECT mean_energy_eV, drift_velocity_m_s, diagnostics_json "
                "FROM cases WHERE mixture_id=? AND solver='propagator' "
                "AND e_over_n_Td=? AND replicate=0",
                (mixture_id, field),
            ).fetchone()
            quality = connection.execute(
                "SELECT passed FROM aggregate_quality WHERE mixture_id=? "
                "AND solver='propagator' AND e_over_n_Td=?",
                (mixture_id, field),
            ).fetchone()
            bins = connection.execute(
                "SELECT energy_eV, energy_width_eV, eedf FROM eedf_bins "
                "WHERE mixture_id=? AND solver='propagator' "
                "AND e_over_n_Td=? AND replicate=0 ORDER BY bin_index",
                (mixture_id, field),
            ).fetchall()
            if case is None or quality is None or int(quality[0]) != 1 or not bins:
                raise SolverQualificationError(
                    f"workflow database lacks qualified Propagator evidence at {field:g} Td"
                )
            diagnostics = json.loads(str(case["diagnostics_json"]))
            propagator = diagnostics.get("propagator")
            if not isinstance(propagator, dict):
                raise SolverQualificationError(
                    f"workflow database lacks Propagator diagnostics at {field:g} Td"
                )
            diagnostic_names = (
                "energy_cells",
                "polar_cells",
                "energy_max_eV",
                "growth_frequency_s_inv",
                "iterations",
                "operator_residual_L1",
                "number_balance_residual",
                "negative_population_mass",
                "tail_probability",
                "tail_probability_target",
                "outer_acceleration_flux_fraction",
                "estimated_peak_memory_bytes",
            )
            try:
                selected_diagnostics = {
                    name: propagator[name] for name in diagnostic_names
                }
            except KeyError as exc:
                raise SolverQualificationError(
                    f"workflow database has incomplete Propagator diagnostics at {field:g} Td"
                ) from exc
            if expected_grid is not None and (
                int(selected_diagnostics["energy_cells"]) < expected_grid[0]
                or int(selected_diagnostics["polar_cells"]) != expected_grid[1]
            ):
                raise SolverQualificationError(
                    "workflow database Propagator grid differs from the "
                    f"qualified target grid at {field:g} Td"
                )
            rows.append(
                {
                    "E_over_N_Td": field,
                    "mean_energy_eV": float(case["mean_energy_eV"]),
                    "drift_velocity_m_s": float(case["drift_velocity_m_s"]),
                    "diagnostics": selected_diagnostics,
                    "eedf": [
                        [float(item[0]), float(item[1]), float(item[2])]
                        for item in bins
                    ],
                }
            )
    except (sqlite3.Error, OSError, json.JSONDecodeError) as exc:
        raise SolverQualificationError(
            f"cannot read Propagator target database: {database_path}"
        ) from exc
    finally:
        if "connection" in locals():
            connection.close()
    encoded = json.dumps(
        rows,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return sha256(encoded).hexdigest(), metadata


def validate_copied_qualification(
    directory: str | Path,
    entry: Mapping[str, Any],
) -> Path:
    """Validate a qualification copy named by a table or bundle manifest."""

    file = entry.get("file")
    expected = entry.get("sha256")
    if (
        not isinstance(file, str)
        or not file
        or Path(file).name != file
        or not isinstance(expected, str)
    ):
        raise SolverQualificationError(
            "solver qualification manifest entry is malformed"
        )
    path = Path(directory) / file
    if not path.is_file() or file_sha256(path) != expected:
        raise SolverQualificationError(
            f"solver qualification copy is missing or altered: {file}"
        )
    return path


def _strictly_increasing_positive_tuple(
    value: Any,
    label: str,
    *,
    allow_equal: bool = False,
) -> tuple[float, ...]:
    if not isinstance(value, list) or not value:
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        )
    try:
        result = tuple(float(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        ) from exc
    if any(not math.isfinite(item) or item <= 0 for item in result):
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        )
    if any(
        right < left if allow_equal else right <= left
        for left, right in zip(result, result[1:])
    ):
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        )
    return result


def _integer_pair(value: Any, label: str) -> tuple[int, int]:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
    ):
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        )
    return int(value[0]), int(value[1])


def _positive_float(value: Any, label: str) -> float:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        ) from exc
    if not 0.0 < numeric < float("inf"):
        raise SolverQualificationError(
            f"Propagator target qualification has malformed {label}"
        )
    return numeric


def _validate_target_fine_runs(
    value: Any,
    fields: tuple[float, ...],
    fine_grid: tuple[int, int],
) -> None:
    if not isinstance(value, list) or len(value) != len(fields):
        raise SolverQualificationError(
            "Propagator target qualification has incomplete fine runs"
        )
    try:
        indexed = {
            float(row["E_over_N_Td"]): row for row in value if isinstance(row, dict)
        }
    except (KeyError, TypeError, ValueError) as exc:
        raise SolverQualificationError(
            "Propagator target qualification has malformed fine-run fields"
        ) from exc
    if tuple(sorted(indexed)) != tuple(sorted(fields)):
        raise SolverQualificationError(
            "Propagator target qualification fine-run fields are inconsistent"
        )
    for field in fields:
        row = indexed[field]
        if (
            row.get("status") != "ok"
            or row.get("quality_gates_passed") is not True
            or int(row.get("energy_cells", -1)) < fine_grid[0]
            or int(row.get("polar_cells", -1)) != fine_grid[1]
        ):
            raise SolverQualificationError(
                f"Propagator target fine run failed at {field:g} Td"
            )


def _validate_target_refinement(
    value: Any,
    fields: tuple[float, ...],
    *,
    scalar_limit: float,
    eedf_limit: float,
) -> None:
    if not isinstance(value, list) or len(value) != len(fields):
        raise SolverQualificationError(
            "Propagator target qualification has incomplete refinement comparisons"
        )
    try:
        indexed = {
            float(row["E_over_N_Td"]): row for row in value if isinstance(row, dict)
        }
    except (KeyError, TypeError, ValueError) as exc:
        raise SolverQualificationError(
            "Propagator target qualification has malformed refinement fields"
        ) from exc
    if tuple(sorted(indexed)) != tuple(sorted(fields)):
        raise SolverQualificationError(
            "Propagator target refinement fields are inconsistent"
        )
    for field in fields:
        row = indexed[field]
        differences = row.get("relative_differences")
        gated = row.get("gated_scalars")
        try:
            scalar_max = max(float(differences[name]) for name in gated)
            eedf = float(row["eedf_weighted_L1"])
        except (KeyError, TypeError, ValueError) as exc:
            raise SolverQualificationError(
                f"Propagator target refinement is malformed at {field:g} Td"
            ) from exc
        if (
            row.get("status") != "available"
            or row.get("passed") is not True
            or scalar_max > scalar_limit
            or eedf > eedf_limit
        ):
            raise SolverQualificationError(
                f"Propagator target refinement failed at {field:g} Td"
            )


def _source_tree_digest(root: Path, files: list[str]) -> str:
    digest = sha256()
    for relative_text in files:
        relative = Path(relative_text)
        if relative.is_absolute() or ".." in relative.parts:
            raise SolverQualificationError(
                "Propagator fingerprint contains an unsafe source path"
            )
        source = (root / relative).resolve()
        if not source.is_relative_to(root) or not source.is_file():
            raise SolverQualificationError(
                f"Propagator fingerprint source is unavailable: {relative_text}"
            )
        canonical = relative.as_posix().encode("utf-8")
        digest.update(canonical + b"\0" + source.read_bytes() + b"\0")
    return digest.hexdigest()


def _validate_input_hashes(root: Path, hashes: Mapping[str, Any]) -> None:
    for relative_text, expected in hashes.items():
        if not isinstance(relative_text, str) or not isinstance(expected, str):
            raise SolverQualificationError(
                "Propagator qualification input hashes are malformed"
            )
        relative = Path(relative_text)
        if relative.is_absolute() or ".." in relative.parts:
            raise SolverQualificationError(
                "Propagator qualification contains an unsafe input path"
            )
        source = (root / relative).resolve()
        if (
            not source.is_relative_to(root)
            or not source.is_file()
            or file_sha256(source) != expected
        ):
            raise SolverQualificationError(
                f"Propagator qualification input changed: {relative_text}"
            )
