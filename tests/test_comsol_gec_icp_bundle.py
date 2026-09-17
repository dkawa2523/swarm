from __future__ import annotations

import ast
import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest

import electron_swarm.solvers.monte_carlo.evidence as mc_evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)

from swarm_workflow.comsol.input.export.contracts import (
    COMSOL_FUNCTION_EEDF_TABLE,
)
from swarm_workflow.comsol.input.function_eedf import (
    COMSOL_EEDF_COLUMNS,
    ComsolEedfImportContract,
    piecewise_linear_weighted_moments,
    read_comsol_function_eedf_grid,
)
from swarm_workflow.selection import MC_QUALIFICATION_TABLE
from swarm_workflow.comsol.models.gec_icp import (
    GecIcpBundleError,
    validate_gec_icp_bundle,
)
from swarm_workflow.quality.solver import (
    PROPAGATOR_CORE_QUALIFICATION_FILE,
    validate_propagator_core_qualification,
)
from swarm_workflow.quality.table import quality_table_schema
from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
    SOURCE_CHOICES,
)


ROOT = Path(__file__).resolve().parents[1]
CORE_QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_p1_deterministic_qualification_20260908.json"
)
ICP_TARGET_QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_gec_icp_mean_energy_guarded_qualification_20260913.json"
)
PRESSURE_PA = 2.66644
TEMPERATURE_K = 300.0
BOLTZMANN_J_K = 1.380649e-23
CANONICAL_ARGON_CROSS_SECTIONS = (
    ROOT / "examples" / "cross_sections" / "argon_application_library.csv"
)


@pytest.mark.parametrize("source", SOURCE_CHOICES)
def test_gec_icp_accepts_each_qualified_function_eedf_source(
    tmp_path: Path, source: str
) -> None:
    bundle = _write_bundle(tmp_path, source)

    evidence = validate_gec_icp_bundle(
        bundle,
        source,
        expected_transport_definition="flux",
        expected_mc_qualification_profile=(
            MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            if source == "monte_carlo"
            else None
        ),
    )

    assert evidence["source"] == source
    assert evidence["transport_definition"] == "flux"
    assert evidence["pressure_Pa"] == pytest.approx(PRESSURE_PA)
    assert evidence["temperature_K"] == pytest.approx(TEMPERATURE_K)
    assert evidence["mean_energy_support_eV"] == pytest.approx((1.0, 2.0))
    assert Path(evidence["table_paths"]["function_eedf"]).is_file()
    assert Path(evidence["table_paths"]["transport"]).name == (
        "transport_vs_mean_energy.csv"
    )
    assert Path(evidence["table_paths"]["quality"]).name == "quality.csv"
    assert evidence["manifest_sha256"] == _digest(bundle / "manifest.json")
    assert evidence["artifact_verification"]["verified_files"] == (
        4 if source == "monte_carlo" else 3
    )
    if source == "propagator":
        assert evidence["provenance"]["solver_qualification"]["decision"] == (
            "p1_deterministic_core_qualified"
        )
    if source == "monte_carlo":
        assert (
            evidence["provenance"]["mc_qualification"]["qualification_profile"]
            == MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
        )


def test_gec_icp_rejects_a_changed_manifest_listed_file(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, "two_term")
    with (bundle / "transport_vs_mean_energy.csv").open(
        "a", encoding="utf-8"
    ) as stream:
        stream.write("\n")

    with pytest.raises(GecIcpBundleError, match="SHA-256 mismatch"):
        validate_gec_icp_bundle(bundle, "two_term")


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("field", "homogeneous DC"),
        ("mixture", "must be Ar"),
        ("pressure", "pressure does not match"),
        ("temperature", "temperature does not match"),
    ),
)
def test_gec_icp_rejects_wrong_physical_context(
    tmp_path: Path, mutation: str, message: str
) -> None:
    bundle = _write_bundle(tmp_path, "two_term")
    manifest = _manifest(bundle)
    if mutation == "field":
        manifest["physical_context"]["field"]["type"] = "rf_periodic"
    elif mutation == "mixture":
        manifest["mixture"]["species"][0]["species"] = "Ne"
    elif mutation == "pressure":
        manifest["physical_context"]["gas_number_density_m3"] *= 2.0
    else:
        manifest["physical_context"]["gas_temperature_K"] = 400.0
    _write_manifest(bundle, manifest)

    with pytest.raises(GecIcpBundleError, match=message):
        validate_gec_icp_bundle(bundle, "two_term")


def test_gec_icp_rejects_source_mismatch(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, "two_term")

    with pytest.raises(GecIcpBundleError, match="source mismatch"):
        validate_gec_icp_bundle(bundle, "propagator")


def test_gec_icp_rejects_transport_definition_mismatch(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, "two_term")

    with pytest.raises(GecIcpBundleError, match="transport_definition mismatch"):
        validate_gec_icp_bundle(
            bundle,
            "two_term",
            expected_transport_definition="different_flux_contract",
        )


def test_gec_icp_rejects_noncanonical_cross_section_provenance(
    tmp_path: Path,
) -> None:
    bundle = _write_bundle(tmp_path, "two_term")
    manifest = _manifest(bundle)
    manifest["hashes"]["cross_sections_sha256"] = "f" * 64
    _write_manifest(bundle, manifest)

    with pytest.raises(
        GecIcpBundleError,
        match="does not match the repository canonical",
    ):
        validate_gec_icp_bundle(bundle, "two_term")


def test_gec_icp_rejects_noncanonical_active_function_eedf(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, "two_term")
    manifest = _manifest(bundle)
    manifest["tables"][COMSOL_FUNCTION_EEDF_TABLE]["canonical_comsol_input"] = False
    _write_manifest(bundle, manifest)

    with pytest.raises(GecIcpBundleError, match="canonical active Function-EEDF"):
        validate_gec_icp_bundle(bundle, "two_term")


def test_gec_icp_rejects_unaccepted_quality_row(tmp_path: Path) -> None:
    bundle = _write_bundle(tmp_path, "two_term")
    quality_columns = quality_table_schema("two_term").columns
    _write_csv(
        bundle / "quality.csv",
        quality_columns,
        _quality_rows("two_term", accepted=(False, True)),
    )
    _refresh_table_digest(bundle, "quality.csv")

    with pytest.raises(GecIcpBundleError, match="unaccepted quality points"):
        validate_gec_icp_bundle(bundle, "two_term")


def test_gec_icp_requires_passing_propagator_core_qualification(
    tmp_path: Path,
) -> None:
    bundle = _write_bundle(tmp_path, "propagator")
    manifest = _manifest(bundle)
    manifest["solver_qualification"] = None
    _write_manifest(bundle, manifest)

    with pytest.raises(GecIcpBundleError, match="lacks core solver qualification"):
        validate_gec_icp_bundle(bundle, "propagator")


def test_gec_icp_accepts_required_generic_icp_target_qualification(
    tmp_path: Path,
) -> None:
    bundle = _write_bundle(tmp_path, "propagator")

    evidence = validate_gec_icp_bundle(
        bundle, "propagator", expected_transport_definition="flux"
    )

    target = evidence["provenance"]["target_qualification"]
    assert target["schema"] == "swarm.propagator_target_qualification.v1"
    assert target["target"] == "argon_gec_icp_restricted_local_mean_energy"
    assert target["fields_Td"] == [1500.0, 2500.0, 2744.2936035]


def test_gec_icp_requires_generic_icp_target_qualification(
    tmp_path: Path,
) -> None:
    bundle = _write_bundle(tmp_path, "propagator")
    manifest = _manifest(bundle)
    manifest["target_qualification"] = None
    _write_manifest(bundle, manifest)

    with pytest.raises(
        GecIcpBundleError,
        match="lacks required generic ICP target qualification",
    ):
        validate_gec_icp_bundle(bundle, "propagator")


def test_gec_icp_rejects_noncanonical_target_qualification_schema(
    tmp_path: Path,
) -> None:
    bundle = _write_bundle(tmp_path, "propagator")
    _add_target_qualification(bundle, mutation="schema")

    with pytest.raises(GecIcpBundleError, match="generic propagator target"):
        validate_gec_icp_bundle(bundle, "propagator")


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ("provenance", "provenance differs"),
        ("failed", "not qualified"),
        ("target", "required scope"),
    ),
)
def test_gec_icp_rejects_invalid_generic_target_qualification(
    tmp_path: Path, mutation: str, message: str
) -> None:
    bundle = _write_bundle(tmp_path, "propagator")
    _add_target_qualification(bundle, mutation=mutation)

    with pytest.raises(GecIcpBundleError, match=message):
        validate_gec_icp_bundle(bundle, "propagator")


@pytest.mark.parametrize("failure", ("profile", "row"))
def test_gec_icp_requires_passing_restricted_lmea_mc_qualification(
    tmp_path: Path, failure: str
) -> None:
    bundle = _write_bundle(tmp_path, "monte_carlo")
    if failure == "profile":
        manifest = _manifest(bundle)
        manifest["source_policy"]["qualification_profile"] = "full_transport"
        _write_manifest(bundle, manifest)
    else:
        columns = (
            "E_over_N_Td",
            "qualification_profile",
            "active_closure_quality_passed",
        )
        _write_csv(
            bundle / MC_QUALIFICATION_TABLE,
            columns,
            (
                (1.0, MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA, 1),
                (2.0, MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA, 0),
            ),
        )
        _refresh_table_digest(bundle, MC_QUALIFICATION_TABLE)

    with pytest.raises(GecIcpBundleError, match="restricted-LMEA"):
        validate_gec_icp_bundle(bundle, "monte_carlo")


def test_gec_icp_bundle_owner_does_not_import_gec_ccp() -> None:
    path = ROOT / "swarm_workflow" / "comsol" / "models" / "gec_icp" / "bundle.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }

    assert not any("gec_ccp" in name for name in imported)


def _write_bundle(root: Path, source: str) -> Path:
    bundle = root / source
    bundle.mkdir()
    active_name, active_metadata = _write_function_eedf(bundle, source)
    transport_columns = (
        "mean_energy_eV",
        "E_over_N_Td",
        "reduced_mobility_m2_V_s_m3",
    )
    _write_csv(
        bundle / "transport_vs_mean_energy.csv",
        transport_columns,
        ((1.0, 1.0, 1.0e24), (2.0, 2.0, 9.0e23)),
    )
    quality_columns = quality_table_schema(source).columns
    quality_rows = _quality_rows(source, accepted=(True, True))
    _write_csv(bundle / "quality.csv", quality_columns, quality_rows)
    tables = {
        active_name: active_metadata,
        "transport_vs_mean_energy.csv": {
            "argument": "mean_energy_eV",
            "artifact_role": "canonical_comsol_coefficient_input",
            "columns": list(transport_columns),
            "sha256": _digest(bundle / "transport_vs_mean_energy.csv"),
        },
        "quality.csv": {
            "argument": "E_over_N_Td",
            "artifact_role": "swarm_quality_audit",
            "columns": list(quality_columns),
            "sha256": _digest(bundle / "quality.csv"),
        },
    }
    hashes: dict[str, str] = {
        "workflow_config_sha256": "c" * 64,
        "base_config_sha256": "a" * 64,
        "cross_sections_sha256": _canonical_cross_sections_hash(),
    }
    mc_sampling_plan = None
    mc_qualification = None
    source_policy: dict[str, object] = {
        "source": source,
        "field_type": "dc",
        "rf_frequency_Hz": None,
        "postprocess": "none",
        "transport_definition": "flux",
    }
    if source == "monte_carlo":
        sampling_entries = [
            {
                "e_over_n_Td": field,
                "particles": 128,
                "warmup_collisions": 16,
                "max_collisions": 32,
                "tail_max_collisions": 32,
                "replicas": 4,
                "transport_correlation_lag_barriers": 4,
                "transport_estimator": "single_field",
            }
            for field in (1.0, 2.0)
        ]
        sampling_json = json.dumps(
            sampling_entries, sort_keys=True, separators=(",", ":")
        )
        sampling_digest = sha256(sampling_json.encode("utf-8")).hexdigest()
        hashes.update(
            {
                "mc_sampling_plan_json": sampling_json,
                "mc_sampling_plan_sha256": sampling_digest,
                "mc_eedf_estimator_schema_version": (
                    mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
                ),
                "mc_solver_source_sha256": monte_carlo_source_sha256(),
            }
        )
        mc_sampling_plan = {
            "status": "validated_against_mc_cases",
            "sha256": sampling_digest,
            "entries": sampling_entries,
        }
        source_policy.update(
            {
                "qualification_profile": (
                    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
                ),
                "qualification_outputs": [
                    "function_eedf",
                    "reduced_mobility",
                    "direct_elastic_energy_loss",
                ],
            }
        )
        qualification_columns = (
            "E_over_N_Td",
            "qualification_profile",
            "active_closure_quality_passed",
        )
        _write_csv(
            bundle / MC_QUALIFICATION_TABLE,
            qualification_columns,
            (
                (1.0, MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA, 1),
                (2.0, MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA, 1),
            ),
        )
        tables[MC_QUALIFICATION_TABLE] = {
            "argument": "E_over_N_Td",
            "artifact_role": "all_planned_mc_anchor_qualification",
            "columns": list(qualification_columns),
            "sha256": _digest(bundle / MC_QUALIFICATION_TABLE),
        }
        mc_qualification = {
            "file": MC_QUALIFICATION_TABLE,
            "all_planned_anchors": 2,
            "qualified_anchors": 2,
            "coefficient_tables_available": True,
            "table_build_failure": None,
        }

    solver_qualification = None
    if source == "propagator":
        qualification = validate_propagator_core_qualification(CORE_QUALIFICATION)
        hashes["propagator_solver_source_sha256"] = (
            qualification.implementation_fingerprint_sha256
        )
        target = bundle / PROPAGATOR_CORE_QUALIFICATION_FILE
        target.write_bytes(CORE_QUALIFICATION.read_bytes())
        solver_qualification = qualification.manifest_entry(file=target.name)

    manifest = {
        "format_version": 1,
        "stage": "export-comsol",
        "status": "ok",
        "source": source,
        "physical_context": {
            "schema": "swarm_physical_context.v1",
            "gas_temperature_K": TEMPERATURE_K,
            "gas_number_density_m3": (PRESSURE_PA / (BOLTZMANN_J_K * TEMPERATURE_K)),
            "field": {"type": "dc", "magnetic_field": {"enabled": False}},
        },
        "mixture": {
            "mixture_id": 0,
            "species": [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}],
        },
        "hashes": hashes,
        "mc_sampling_plan": mc_sampling_plan,
        "mc_qualification": mc_qualification,
        "solver_qualification": solver_qualification,
        "target_qualification": None,
        "source_policy": source_policy,
        "valid_ranges": {"mean_energy_eV": [1.0, 2.0]},
        "monotonicity": {"mean_energy_strictly_monotonic": True},
        "quality_summary": {"passed": True, "failed_points": 0, "total_points": 2},
        "missing_required_coefficients": [],
        "eedf_artifacts": {
            "comsol_input": {
                "file": active_name,
                "sha256": active_metadata["sha256"],
                "role": active_metadata["artifact_role"],
            }
        },
        "tables": tables,
    }
    _write_manifest(bundle, manifest)
    if source == "propagator":
        _add_target_qualification(bundle)
    return bundle


def _write_function_eedf(bundle: Path, source: str) -> tuple[str, dict[str, object]]:
    energies = np.linspace(0.0, 40.0, 321)
    means = np.asarray([1.0, 2.0])
    values = np.asarray([_normalized_exponential(energies, mean) for mean in means])
    active_name = COMSOL_FUNCTION_EEDF_TABLE
    active_path = bundle / active_name
    _write_csv(
        active_path,
        COMSOL_EEDF_COLUMNS,
        (
            (energy, mean, values[mean_index, energy_index])
            for mean_index, mean in enumerate(means)
            for energy_index, energy in enumerate(energies)
        ),
    )
    grid = read_comsol_function_eedf_grid(active_path)
    contract = ComsolEedfImportContract()
    return active_name, {
        "argument": "electron_energy_eV,mean_energy_eV",
        "argument_order": list(COMSOL_EEDF_COLUMNS[:2]),
        "artifact_role": "canonical_comsol_function_eedf_input",
        "canonical_comsol_input": True,
        "columns": list(COMSOL_EEDF_COLUMNS),
        "format": contract.format,
        "representation": contract.representation,
        "grid_shape": [len(means), len(energies)],
        "mean_energy_grid_points": len(means),
        "energy_grid_points": len(energies),
        "grid_axis_order": list(COMSOL_EEDF_COLUMNS[:2]),
        "row_order": contract.row_order,
        "mean_energy_range_eV": [float(means[0]), float(means[-1])],
        "projected_normalization_error_max": grid.normalization_error_max,
        "projected_mean_energy_relative_error_max": (
            grid.mean_energy_relative_error_max
        ),
        "projected_nonnegative_minimum": grid.minimum_value,
        "comsol_import": contract.import_settings(),
        "sha256": _digest(active_path),
    }


def _add_target_qualification(bundle: Path, *, mutation: str | None = None) -> None:
    payload = json.loads(ICP_TARGET_QUALIFICATION.read_text(encoding="utf-8"))
    manifest = _manifest(bundle)
    payload["inputs"]["core_qualification"] = {
        key: value
        for key, value in manifest["solver_qualification"].items()
        if key != "file"
    }
    for name in (
        "base_config_sha256",
        "workflow_config_sha256",
        "cross_sections_sha256",
        "propagator_solver_source_sha256",
    ):
        payload["inputs"][name] = manifest["hashes"][name]
    if mutation == "schema":
        payload["schema"] = "swarm.propagator_gec_target_qualification.v1"
    elif mutation == "provenance":
        payload["inputs"]["workflow_config_sha256"] = "d" * 64
    elif mutation == "failed":
        payload["decision"]["target_refinement_qualified"] = False
        payload["decision"]["blocking_gates"] = ["fixture_failure"]
    elif mutation == "target":
        payload["scope"]["target"] = "different_icp_model"
    target_path = bundle / "propagator_target_qualification.json"
    target_path.write_text(json.dumps(payload), encoding="utf-8")
    manifest = _manifest(bundle)
    manifest["target_qualification"] = _target_manifest_entry(target_path, payload)
    _write_manifest(bundle, manifest)


def _target_manifest_entry(path: Path, payload: dict[str, object]) -> dict[str, object]:
    scope = payload["scope"]
    inputs = payload["inputs"]
    return {
        "file": path.name,
        "sha256": _digest(path),
        "schema": payload["schema"],
        "decision": "target_refinement_qualified",
        "core_qualification_sha256": inputs["core_qualification"]["sha256"],
        "medium_evidence_sha256": inputs["medium_evidence_sha256"],
        "fields_Td": scope["fields_Td"],
        "medium_grid": scope["medium_grid"],
        "fine_grid": scope["fine_grid"],
        "target": scope["target"],
        "operating_range_bracket_Td": scope["operating_range_bracket_Td"],
        "table_support_cap_Td": scope["table_support_cap_Td"],
        "mixture_id": scope["mixture_id"],
    }


def _normalized_exponential(energies: np.ndarray, requested_mean: float) -> np.ndarray:
    lower = 1.0e-3
    upper = 100.0
    for _ in range(100):
        decay = 0.5 * (lower + upper)
        values = np.exp(-decay * energies)
        normalization, first_moment = piecewise_linear_weighted_moments(
            energies, values
        )
        if first_moment / normalization > requested_mean:
            lower = decay
        else:
            upper = decay
    result = np.exp(-0.5 * (lower + upper) * energies)
    normalization, _ = piecewise_linear_weighted_moments(energies, result)
    return result / normalization


def _quality_rows(
    source: str, *, accepted: tuple[bool, bool]
) -> tuple[tuple[object, ...], ...]:
    columns = quality_table_schema(source).columns
    rows = []
    for field, passed in zip((1.0, 2.0), accepted, strict=True):
        values: dict[str, object] = {name: "" for name in columns}
        values.update(
            {
                "E_over_N_Td": field,
                "E_over_N_V_m2": field * 1.0e-21,
                "passed": int(passed if source != "monte_carlo" else False),
                "aggregate_quality_passed": int(passed),
                "failure_reasons_json": "[]",
                "aggregate_failure_reasons_json": "[]",
                "eedf_normalization_error": 0.0,
                "valid_replicates": 4 if source == "monte_carlo" else 1,
                "uncertainty_available": int(source == "monte_carlo"),
                "required_rate_min_process_peak_fraction": 0.0,
                "quality_policy_reevaluated": 0,
                "solver_diagnostics_available": 1,
                "solver_diagnostics_passed": 1,
                "quality_source": "fixture",
            }
        )
        if source == "monte_carlo":
            values.update(
                {
                    "qualification_profile": (
                        MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
                    ),
                    "active_closure_quality_passed": int(passed),
                    "active_closure_failure_reasons_json": "[]",
                }
            )
        rows.append(tuple(values[name] for name in columns))
    return tuple(rows)


def _write_csv(
    path: Path,
    columns: tuple[str, ...],
    rows: object,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def _manifest(bundle: Path) -> dict[str, object]:
    return json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))


def _write_manifest(bundle: Path, manifest: dict[str, object]) -> None:
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _refresh_table_digest(bundle: Path, name: str) -> None:
    manifest = _manifest(bundle)
    manifest["tables"][name]["sha256"] = _digest(bundle / name)
    _write_manifest(bundle, manifest)


def _canonical_cross_sections_hash() -> str:
    file_digest = _digest(CANONICAL_ARGON_CROSS_SECTIONS)
    digest = sha256()
    digest.update(str(CANONICAL_ARGON_CROSS_SECTIONS.resolve()).encode("utf-8"))
    digest.update(b"\0")
    digest.update(file_digest.encode("ascii"))
    digest.update(b"\0")
    return digest.hexdigest()


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()
