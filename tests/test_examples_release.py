from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

import electron_swarm
from electron_swarm import load_config, run


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted((ROOT / "examples").glob("*.yaml"))
LIGHT_EXAMPLES = [
    ROOT / "examples" / "two_term.yaml",
    ROOT / "examples" / "multi_term_surrogate.yaml",
    ROOT / "examples" / "pn_dcs_moment_table.yaml",
]
MC_EXAMPLES = [
    ROOT / "examples" / "compare_three_solvers.yaml",
    ROOT / "examples" / "magnetic_mc.yaml",
]


def test_examples_are_schema_v2_configs() -> None:
    assert {path.name for path in EXAMPLES} == {
        "two_term.yaml",
        "multi_term_surrogate.yaml",
        "compare_three_solvers.yaml",
        "magnetic_mc.yaml",
        "pn_dcs_moment_table.yaml",
    }
    for path in EXAMPLES:
        cfg = load_config(path)
        assert cfg.schema_version == 2
        assert cfg.run.solvers


@pytest.mark.parametrize("path", LIGHT_EXAMPLES, ids=lambda p: p.name)
def test_light_examples_run_without_writing(path: Path) -> None:
    result = run(load_config(path), write=False)
    assert result.cases
    assert all(case.schema_version == "2" for case in result.cases)


def test_readme_quickstart_target_runs() -> None:
    result = run(load_config(ROOT / "examples" / "two_term.yaml"), write=False)
    assert len(result.cases) == 1
    assert result.cases[0].solver == "two_term"


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

    dependencies = set(project["dependencies"])
    assert {"json5", "molmass", "matplotlib"}.isdisjoint(dependencies)
    assert "matplotlib" in set(project["optional-dependencies"]["plot"])
    assert "pytest" in set(project["optional-dependencies"]["dev"])
    assert pyproject["tool"]["setuptools"]["packages"]["find"]["include"] == [
        "electron_swarm*"
    ]


def test_legacy_top_level_runtime_projects_are_removed() -> None:
    for name in ["swarm_mc", "swarm_comsol_exporter", "run_sweep.py"]:
        assert not (ROOT / name).exists()


def test_case_summary_dict_uses_stable_metadata_allowlist() -> None:
    result = run(load_config(ROOT / "examples" / "two_term.yaml"), write=False)
    summary = result.cases[0].summary_dict()
    assert summary["meta_solver_method"] == "native_sg"
    assert "meta_grid_n_cells" not in summary
    assert "meta_transport_definition" in summary
    assert "meta_tail_probability" in summary


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
