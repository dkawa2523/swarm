from __future__ import annotations

import tomllib
from pathlib import Path

import pandas as pd
import pytest

import electron_swarm
from electron_swarm import load_config, run
from electron_swarm.core.result_metadata import PRODUCT_CASE_METADATA_KEYS


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted((ROOT / "examples").glob("*.yaml"))
PRODUCT_EXAMPLES = [
    path for path in EXAMPLES if not path.name.startswith("workflow_")
]
LIGHT_EXAMPLES = [
    ROOT / "examples" / "two_term.yaml",
    ROOT / "examples" / "multi_term_direct_lmax1.yaml",
    ROOT / "examples" / "multi_term_direct_lmax4.yaml",
    ROOT / "examples" / "pn_dcs_moment_table.yaml",
    ROOT / "examples" / "ar_o2_base.yaml",
    ROOT / "examples" / "argon_comsol_base.yaml",
]
MC_EXAMPLES = [
    ROOT / "examples" / "compare_three_solvers.yaml",
    ROOT / "examples" / "magnetic_mc.yaml",
]


def test_examples_are_schema_v2_configs() -> None:
    required_examples = {
        "two_term.yaml",
        "multi_term_direct_lmax1.yaml",
        "multi_term_direct_lmax4.yaml",
        "compare_three_solvers.yaml",
        "magnetic_mc.yaml",
        "pn_dcs_moment_table.yaml",
        "ar_o2_base.yaml",
        "argon_comsol_base.yaml",
    }
    assert required_examples <= {path.name for path in PRODUCT_EXAMPLES}
    for path in PRODUCT_EXAMPLES:
        cfg = load_config(path)
        assert cfg.schema_version == 2
        assert cfg.run.solvers


@pytest.mark.parametrize("path", LIGHT_EXAMPLES, ids=lambda p: p.name)
def test_light_examples_run_without_writing(path: Path) -> None:
    result = run(load_config(path), write=False)
    assert result.cases
    assert all(case.schema_version == "2" for case in result.cases)


def test_pn_dcs_example_runs_with_moment_table() -> None:
    result = run(load_config(ROOT / "examples" / "pn_dcs_moment_table.yaml"), write=False)
    [case] = result.cases
    assert case.metadata["solver_method"] == "pn_dcs"
    assert case.metadata["angular_model"] == "moment_table"
    assert case.metadata["angular_moment_source"] == "moment_table"
    assert case.metadata["ordinary_integral_xs_closure"] is False


def test_readme_quickstart_target_runs() -> None:
    result = run(load_config(ROOT / "examples" / "two_term.yaml"), write=False)
    assert len(result.cases) == 1
    assert result.cases[0].solver == "two_term"


def test_direct_lmax1_example_runs_as_gated_product_path() -> None:
    result = run(load_config(ROOT / "examples" / "multi_term_direct_lmax1.yaml"), write=False)
    [case] = result.cases
    assert case.solver == "multi_term"
    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["direct_pn_operator"] is True
    assert case.metadata["ordinary_integral_xs_closure"] is True
    assert case.metadata["exact_dcs_based"] is False


def test_direct_lmax4_example_runs_in_limited_product_scope() -> None:
    result = run(load_config(ROOT / "examples" / "multi_term_direct_lmax4.yaml"), write=False)
    [case] = result.cases
    assert case.solver == "multi_term"
    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["lmax"] == 4
    assert case.metadata["transport_definition"] == "f0_gradient_reconstruction"
    assert set(case.metadata) <= PRODUCT_CASE_METADATA_KEYS
    assert case.metadata["exact_dcs_based"] is False


@pytest.mark.mc
@pytest.mark.parametrize("path", MC_EXAMPLES, ids=lambda p: p.name)
def test_mc_examples_run_without_writing(path: Path) -> None:
    result = run(load_config(path), write=False)
    assert result.cases


def test_release_metadata_is_product_shaped() -> None:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = pyproject["project"]
    assert project["name"] == "electron-swarm"
    assert project["version"] == electron_swarm.__version__
    assert project["readme"] == "README.md"
    assert project["scripts"]["electron-swarm"] == "electron_swarm.runner:main"
    assert project["scripts"]["swarm-workflow"] == "swarm_workflow.cli:main"

    dependencies = set(project["dependencies"])
    assert {"json5", "molmass", "matplotlib"}.isdisjoint(dependencies)
    assert "matplotlib" in set(project["optional-dependencies"]["plot"])
    assert "pytest" in set(project["optional-dependencies"]["dev"])
    assert pyproject["tool"]["setuptools"]["packages"]["find"]["include"] == [
        "electron_swarm*",
        "swarm_workflow*",
    ]


def test_legacy_top_level_runtime_projects_are_removed() -> None:
    for name in ["swarm_mc", "swarm_comsol_exporter", "run_sweep.py"]:
        assert not (ROOT / name).exists()


def test_summary_csv_uses_stable_metadata_allowlist(tmp_path: Path) -> None:
    cfg = load_config(ROOT / "examples" / "two_term.yaml")
    cfg.output.directory = tmp_path
    run(cfg, write=True)
    summary = pd.read_csv(tmp_path / f"{cfg.output.base_name}_summary.csv")
    row = summary.iloc[0]
    assert row["solver_method"] == "native_sg"
    assert "meta_grid_n_cells" not in summary.columns
    assert "meta_transport_definition" in summary.columns
    assert "meta_tail_rate_fraction_max" not in summary.columns


def test_public_docs_do_not_advertise_obsolete_api() -> None:
    forbidden = [
        "run.mode",
        "boltzmann_two_term",
        "multiterm_boltzmann",
        "allow_experimental_operator",
        "output.compatibility",
        "summary_boltzmann",
        "summary_multiterm",
        "Ar_N2_sweep",
        "swarm_mc",
    ]
    docs = [ROOT / "README.md"]
    docs.extend((ROOT / "docs").rglob("*.md"))
    docs.extend(EXAMPLES)
    docs.extend((ROOT / "configs").rglob("*.yaml"))
    for path in docs:
        text = path.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{token!r} appears in {path}"
