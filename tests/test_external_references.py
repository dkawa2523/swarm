from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.references import load_bolsig_reference, load_mcig_reference
from electron_swarm.references.common import (
    ExternalReferenceConfig,
    ReferenceCaseResult,
    eepf_to_eedf,
    parse_external_reference_configs,
    reference_comparison_metrics,
    reference_to_swarm_case,
)
from electron_swarm.references.runner import run_external_reference_command
from tests.product_helpers import base_product_config, write_config
from tools.benchmark_ar_eedf_consistency import run_benchmark
from tools.benchmark_ar_external_references import (
    _classify_bolsig_mismatch,
    _classify_mcig_mismatch,
    _mc_confidence_status,
    load_benchmark_config,
    run_benchmark as run_external_benchmark,
)
from tools import benchmark_ar_bolsig_mcig_triage as triage

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures" / "references"


def test_bolsig_canonical_csv_parse_normalizes_eedf() -> None:
    [case] = load_bolsig_reference(
        ExternalReferenceConfig(
            id="bolsig_plus",
            path=FIXTURES / "bolsig_reference.csv",
            format="electron_swarm_reference_csv",
        )
    )
    assert case.reference_id == "bolsig_plus"
    assert case.e_over_n_Td == 50.0
    widths = np.array([1.0, 1.0, 1.0])
    assert np.sum(case.eedf_eV_inv * widths) == pytest.approx(1.0)
    assert case.rates["ionization"] == pytest.approx(1.0e-15)


def test_bolsig_text_eepf_to_eedf_conversion() -> None:
    [case] = load_bolsig_reference(
        ExternalReferenceConfig(
            id="bolsig_plus",
            path=FIXTURES / "bolsig_eepf.txt",
            format="bolsig_text",
            eedf_convention="eepf",
        )
    )
    expected = eepf_to_eedf(
        case.energy_eV,
        np.array([0.8485281374, 0.2449489743, 0.0632455532]),
    )
    expected = expected / np.sum(expected)
    assert case.eedf_eV_inv == pytest.approx(expected)
    assert case.metadata["eedf_convention_input"] == "eepf_eV_m32"


def test_mcig_canonical_csv_marks_uncertainty_unavailable() -> None:
    [case] = load_mcig_reference(
        ExternalReferenceConfig(
            id="mcig",
            path=FIXTURES / "mcig_reference.csv",
            format="electron_swarm_reference_csv",
            uncertainty="unavailable",
        )
    )
    assert case.reference_id == "mcig"
    assert case.metadata["uncertainty_unavailable"] is True
    assert "scalar_ci95" in case.metadata


def test_reference_missing_file_and_unsupported_format_fail_clearly(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="not found"):
        load_bolsig_reference(
            ExternalReferenceConfig(
                id="bolsig_plus",
                path=tmp_path / "missing.csv",
                format="electron_swarm_reference_csv",
            )
        )
    data = base_product_config(tmp_path)
    data["references"] = {
        "external": [
            {
                "id": "bolsig_plus",
                "path": "missing.csv",
                "format": "unsupported",
            }
        ]
    }
    with pytest.raises(ValueError, match="references.external\\[0\\].format"):
        parse_external_reference_configs(data, base=tmp_path)


def test_external_reference_command_generates_output(tmp_path: Path) -> None:
    script = tmp_path / "write_reference.py"
    script.write_text(
        "\n".join(
            [
                "from pathlib import Path",
                "import sys",
                "Path(sys.argv[1]).write_text(",
                "    'case_id,E_over_N_Td,energy_eV,eedf_eV_inv\\n'",
                "    'cmd_0000,50,0.5,0.6\\n'",
                "    'cmd_0000,50,1.5,0.4\\n',",
                "    encoding='utf-8',",
                ")",
            ]
        ),
        encoding="utf-8",
    )
    output = tmp_path / "generated_bolsig.csv"
    command = f'"{sys.executable}" "{script}" "{{output}}"'
    path = run_external_reference_command(
        reference_id="bolsig_plus",
        command=command,
        output=output,
    )
    assert path == output.resolve()
    assert output.exists()
    assert output.with_suffix(output.suffix + ".log").exists()


def test_reference_comparison_metrics_against_swarm_case() -> None:
    [reference] = load_bolsig_reference(
        ExternalReferenceConfig(
            id="bolsig_plus",
            path=FIXTURES / "bolsig_reference.csv",
            format="electron_swarm_reference_csv",
        )
    )
    candidate = reference_to_swarm_case(reference)
    candidate.solver = "two_term"
    metrics = reference_comparison_metrics(reference, candidate)
    assert metrics["eedf_relative_l1"] == pytest.approx(0.0)
    assert metrics["mobility_relative_difference"] == pytest.approx(0.0)
    assert metrics["diffusion_T_relative_difference"] == pytest.approx(0.0)


def test_ar_benchmark_reports_missing_optional_external_reference(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["run"]["case_prefix"] = "ar_eedf"
    data["references"] = {
        "external": [
            {
                "id": "bolsig_plus",
                "path": "missing_bolsig.csv",
                "format": "electron_swarm_reference_csv",
            }
        ]
    }
    config_path = write_config(tmp_path, data, "ar_reference_missing.yaml")
    paths = run_benchmark(config_path)
    ref_failures = pd.read_csv(tmp_path / "ar_reference_failure_analysis.csv")
    assert "missing_external_reference" in set(ref_failures["category"])
    assert any(path.name == "ar_reference_report.md" for path in paths)


def test_ar_bolsig_plus_equivalence_config_parses() -> None:
    cfg, references = load_benchmark_config(
        ROOT / "configs" / "benchmarks" / "ar_bolsig_plus_equivalence.yaml"
    )
    assert cfg.schema_version == 2
    assert [item.id for item in cfg.run.solvers] == ["two_term", "multi_term"]
    assert references[0].id == "bolsig_plus"


def test_ar_mcig_reference_config_parses() -> None:
    cfg, references = load_benchmark_config(
        ROOT / "configs" / "benchmarks" / "ar_mcig_reference.yaml"
    )
    assert cfg.schema_version == 2
    assert [item.id for item in cfg.run.solvers] == ["two_term", "multi_term", "monte_carlo"]
    assert references[0].id == "mcig"


def test_bolsig_equivalence_missing_reference_skips(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["run"]["case_prefix"] = "ar_bolsig_equiv"
    data["references"] = {
        "external": [
            {
                "id": "bolsig_plus",
                "path": "missing_bolsig.csv",
                "format": "electron_swarm_reference_csv",
            }
        ]
    }
    paths = run_external_benchmark(write_config(tmp_path, data, "bolsig_missing.yaml"))
    failures = pd.read_csv(tmp_path / "prod_failure_analysis.csv")
    report = (tmp_path / "prod_report.md").read_text(encoding="utf-8")
    assert "missing_external_reference" in set(failures["category"])
    assert "SKIP_EXTERNAL_REFERENCE" in report
    assert any(path.name == "prod_report.md" for path in paths)


def test_bolsig_equivalence_synthetic_reference_passes(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["run"]["case_prefix"] = "ar_bolsig_equiv"
    seed_data = base_product_config(tmp_path, ["two_term"])
    seed_data["run"]["case_prefix"] = "ar_bolsig_equiv"
    cfg = load_config(write_config(tmp_path, seed_data, "seed.yaml"))
    [two_term] = run(cfg, write=False).cases
    reference_path = tmp_path / "synthetic_bolsig.csv"
    pd.DataFrame(
        {
            "case_id": two_term.case_id,
            "E_over_N_Td": two_term.e_over_n_Td,
            "energy_eV": two_term.energy_eV,
            "eedf_eV_inv": two_term.eedf,
            "mean_energy_eV": two_term.mean_energy_eV,
            "drift_velocity_m_s": two_term.drift_velocity_m_s,
            "mobility_m2_V_s": two_term.mobility_m2_V_s,
            "diffusion_L_m2_s": two_term.diffusion_L_m2_s,
            "diffusion_T_m2_s": two_term.diffusion_T_m2_s,
            "net_ionization_frequency_s": two_term.net_ionization_frequency_s,
        }
    ).to_csv(reference_path, index=False)
    data["references"] = {
        "external": [
            {
                "id": "bolsig_plus",
                "path": reference_path.as_posix(),
                "format": "electron_swarm_reference_csv",
            }
        ]
    }
    run_external_benchmark(write_config(tmp_path, data, "bolsig_pass.yaml"))
    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    two_row = summary[summary["candidate_solver"] == "two_term"].iloc[0]
    assert two_row["status"] == "PASS"
    assert two_row["eedf_relative_l1"] < 1.0e-12


def test_bolsig_convention_mismatch_classification() -> None:
    rows = _classify_bolsig_mismatch(
        {
            "mean_energy_relative_difference": 0.2,
            "eedf_relative_l1": 0.2,
            "major_rate_relative_difference": 0.0,
            "log_tail_error": 0.0,
            "normalization_error_reference": 0.0,
        },
        candidate_solver="two_term",
        eover=50.0,
    )
    categories = {row["category"] for row in rows}
    assert "bolsig_convention_mismatch" in categories
    assert "energy_grid_mismatch" in categories


def test_mcig_reference_missing_skips(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["run"]["case_prefix"] = "ar_mcig_ref"
    data["references"] = {
        "external": [
            {
                "id": "mcig",
                "path": "missing_mcig.csv",
                "format": "electron_swarm_reference_csv",
                "angular_model": "unknown",
            }
        ]
    }
    run_external_benchmark(
        write_config(tmp_path, data, "mcig_missing.yaml"),
        reference="mcig",
    )
    failures = pd.read_csv(tmp_path / "prod_failure_analysis.csv")
    report = (tmp_path / "prod_report.md").read_text(encoding="utf-8")
    assert "missing_external_reference" in set(failures["category"])
    assert "SKIP_EXTERNAL_REFERENCE" in report


def test_mcig_angular_mismatch_classification() -> None:
    [reference] = load_mcig_reference(
        ExternalReferenceConfig(
            id="mcig",
            path=FIXTURES / "mcig_reference.csv",
            format="electron_swarm_reference_csv",
            angular_model="mcig_default",
        )
    )
    rows = _classify_mcig_mismatch(
        {"eedf_relative_l1": 0.2, "major_rate_relative_difference": 0.0},
        candidate_solver="multi_term",
        candidate_method="pn_closure_direct",
        reference_case=reference,
        eover=50.0,
        angular_model_status="mismatch",
        angular_evidence="reference=mcig_default; candidate=isotropic",
        confidence_status="unknown",
    )
    assert {row["category"] for row in rows} == {"angular_model_mismatch"}


def test_mcig_confidence_interval_status() -> None:
    [reference] = load_mcig_reference(
        ExternalReferenceConfig(
            id="mcig",
            path=FIXTURES / "mcig_reference.csv",
            format="electron_swarm_reference_csv",
            uncertainty="reported",
            angular_model="isotropic",
        )
    )
    candidate = reference_to_swarm_case(reference)
    candidate.solver = "two_term"
    status, coverage = _mc_confidence_status(reference, candidate, {})
    assert status == "pass_within_mc_uncertainty"
    assert coverage == pytest.approx(1.0)


def test_mcig_synthetic_reference_comparison_degraded_for_unknown_angular(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["run"]["case_prefix"] = "ar_mcig_ref"
    seed_cfg = load_config(write_config(tmp_path, data, "seed_mcig.yaml"))
    [two_term] = run(seed_cfg, write=False).cases
    reference_path = tmp_path / "synthetic_mcig.csv"
    pd.DataFrame(
        {
            "case_id": two_term.case_id,
            "E_over_N_Td": two_term.e_over_n_Td,
            "energy_eV": two_term.energy_eV,
            "eedf_eV_inv": two_term.eedf,
            "mean_energy_eV": two_term.mean_energy_eV,
            "drift_velocity_m_s": two_term.drift_velocity_m_s,
            "mobility_m2_V_s": two_term.mobility_m2_V_s,
            "diffusion_L_m2_s": two_term.diffusion_L_m2_s,
            "diffusion_T_m2_s": two_term.diffusion_T_m2_s,
        }
    ).to_csv(reference_path, index=False)
    data["references"] = {
        "external": [
            {
                "id": "mcig",
                "path": reference_path.as_posix(),
                "format": "electron_swarm_reference_csv",
                "angular_model": "unknown",
            }
        ]
    }
    run_external_benchmark(write_config(tmp_path, data, "mcig_degraded.yaml"), reference="mcig")
    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    assert summary.loc[0, "status"] == "DEGRADED"
    assert summary.loc[0, "angular_model_status"] == "unknown"


def _triage_swarm_case(
    solver: str,
    eedf: list[float],
    *,
    lmax: int | None = None,
    angular_model: str = "isotropic",
    rate: float = 1.0,
    mean_energy: float = 1.2,
    drift: float = 1000.0,
) -> SwarmCaseResult:
    energy = np.array([0.5, 1.5, 2.5], dtype=float)
    values = np.array(eedf, dtype=float)
    metadata = {"angular_model": angular_model}
    if solver == "multi_term":
        metadata.update({"solver_method": "pn_closure_direct", "lmax": lmax or 1})
    return SwarmCaseResult(
        solver=solver,
        case_id="triage_0000",
        e_over_n_Td=50.0,
        mean_energy_eV=mean_energy,
        drift_velocity_m_s=drift,
        mobility_m2_V_s=2.0,
        reduced_mobility_m2_V_s_m3=0.0,
        diffusion_L_m2_s=0.1,
        diffusion_T_m2_s=0.2,
        reduced_diffusion_L_m2_s_m3=0.0,
        reduced_diffusion_T_m2_s_m3=0.0,
        net_ionization_frequency_s=rate,
        effective_townsend_m2=0.0,
        energy_eV=energy,
        eedf=values,
        eepf=values / np.sqrt(np.maximum(energy, 1.0e-30)),
        energy_widths_eV=np.ones_like(energy),
        rates=[
            RateResult(
                solver=solver,
                case_id="triage_0000",
                e_over_n_Td=50.0,
                species="Ar",
                process="ionization",
                process_type="IONIZATION",
                threshold_eV=15.0,
                rate_coefficient_m3_s=rate,
                mixture_weighted_rate_m3_s=rate,
            )
        ],
        metadata=metadata,
        schema_version="2",
    )


def _triage_reference(
    reference_id: str,
    eedf: list[float],
    *,
    angular_model: str = "isotropic",
    rate: float = 1.0,
    mean_energy: float = 1.2,
    drift: float = 1000.0,
    ci95: float | None = None,
) -> ReferenceCaseResult:
    metadata = {"angular_model": angular_model}
    if ci95 is not None:
        metadata["scalar_ci95"] = {
            "mean_energy_eV": ci95,
            "drift_velocity_m_s": ci95 * 1000.0,
            "mobility_m2_V_s": ci95,
            "diffusion_L_m2_s": ci95,
            "diffusion_T_m2_s": ci95,
        }
    return ReferenceCaseResult(
        reference_id=reference_id,
        case_id="triage_0000",
        e_over_n_Td=50.0,
        energy_eV=np.array([0.5, 1.5, 2.5], dtype=float),
        eedf_eV_inv=np.array(eedf, dtype=float),
        rates={"ionization": rate},
        scalars={
            "mean_energy_eV": mean_energy,
            "drift_velocity_m_s": drift,
            "mobility_m2_V_s": 2.0,
            "diffusion_L_m2_s": 0.1,
            "diffusion_T_m2_s": 0.2,
            "net_ionization_frequency_s": rate,
        },
        metadata=metadata,
    )


def _run_synthetic_triage(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    references: list[ReferenceCaseResult],
    solver_cases: list[SwarmCaseResult],
) -> pd.DataFrame:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["run"]["case_prefix"] = "triage"
    config_path = write_config(tmp_path, data, "triage.yaml")
    monkeypatch.setattr(triage, "load_selected_references", lambda *args, **kwargs: (references, []))
    monkeypatch.setattr(triage, "run_solver_variants", lambda cfg: (solver_cases, []))
    triage.run_triage(config_path)
    return pd.read_csv(tmp_path / "ar_triage_failure_analysis.csv")


def test_triage_classifies_multi_term_lmax1_code_regression(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failures = _run_synthetic_triage(
        tmp_path,
        monkeypatch,
        references=[],
        solver_cases=[
            _triage_swarm_case("two_term", [0.6, 0.3, 0.1]),
            _triage_swarm_case("multi_term", [0.2, 0.3, 0.5], lmax=1),
        ],
    )
    assert "code_regression_multi_term_lmax1" in set(failures["category"])


def test_triage_classifies_two_term_bolsig_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failures = _run_synthetic_triage(
        tmp_path,
        monkeypatch,
        references=[_triage_reference("bolsig_plus", [0.1, 0.1, 0.8])],
        solver_cases=[
            _triage_swarm_case("two_term", [0.6, 0.3, 0.1]),
            _triage_swarm_case("multi_term", [0.6, 0.3, 0.1], lmax=1),
        ],
    )
    assert "code_or_bolsig_input_mismatch" in set(failures["category"])


def test_triage_classifies_mc_uncertainty_limited(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failures = _run_synthetic_triage(
        tmp_path,
        monkeypatch,
        references=[_triage_reference("mcig", [0.1, 0.1, 0.8], ci95=1.0)],
        solver_cases=[
            _triage_swarm_case("two_term", [0.6, 0.3, 0.1]),
            _triage_swarm_case("multi_term", [0.6, 0.3, 0.1], lmax=1),
            _triage_swarm_case("monte_carlo", [0.6, 0.3, 0.1]),
        ],
    )
    assert "mc_uncertainty_limited" in set(failures["category"])


def test_triage_classifies_angular_model_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failures = _run_synthetic_triage(
        tmp_path,
        monkeypatch,
        references=[_triage_reference("mcig", [0.6, 0.3, 0.1], angular_model="mcig_default")],
        solver_cases=[
            _triage_swarm_case("two_term", [0.6, 0.3, 0.1], angular_model="isotropic"),
            _triage_swarm_case("multi_term", [0.6, 0.3, 0.1], lmax=1, angular_model="isotropic"),
        ],
    )
    assert "angular_model_mismatch" in set(failures["category"])


def test_triage_classifies_rate_only_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failures = _run_synthetic_triage(
        tmp_path,
        monkeypatch,
        references=[_triage_reference("mcig", [0.6, 0.3, 0.1], rate=100.0)],
        solver_cases=[
            _triage_swarm_case("two_term", [0.6, 0.3, 0.1], rate=1.0),
            _triage_swarm_case("multi_term", [0.6, 0.3, 0.1], lmax=1, rate=1.0),
        ],
    )
    assert "rate_convolution_or_cross_section_projection_issue" in set(failures["category"])


def test_triage_classifies_tail_only_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    failures = _run_synthetic_triage(
        tmp_path,
        monkeypatch,
        references=[_triage_reference("mcig", [0.9, 0.099999, 1.0e-6])],
        solver_cases=[
            _triage_swarm_case("two_term", [0.9, 0.05, 0.05]),
            _triage_swarm_case("multi_term", [0.9, 0.05, 0.05], lmax=1),
        ],
    )
    assert "tail_boundary_or_sampling_issue" in set(failures["category"])
