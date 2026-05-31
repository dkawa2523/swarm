from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from electron_swarm import load_config
from benchmarks.run_product_benchmarks import main as benchmark_main


@pytest.mark.regression
def test_product_benchmark_runner_writes_summary(tmp_path: Path) -> None:
    output = tmp_path / "benchmark_summary.csv"
    rc = benchmark_main(
        [
            "--config",
            "configs/benchmarks/ar_simple.yaml",
            "--output",
            str(output),
            "--fail-fast",
        ]
    )

    assert rc == 0
    summary = pd.read_csv(output)
    assert list(summary["status"]) == ["ok"]
    assert summary.loc[0, "cases"] > 0
    assert "two_term" in summary.loc[0, "solvers"]
    assert "multi_term" in summary.loc[0, "solvers"]


def test_ar_eedf_consistency_configs_parse() -> None:
    for path in [
        Path("configs/benchmarks/ar_bolsig_eedf_consistency.yaml"),
        Path("configs/benchmarks/ar_bolsig_eedf_consistency_mc_manual.yaml"),
        Path("configs/benchmarks/ar_bolsig_plus_equivalence.yaml"),
        Path("configs/benchmarks/ar_mcig_reference.yaml"),
        Path("configs/benchmarks/ar_bolsig_mcig_triage.yaml"),
    ]:
        cfg = load_config(path)
        assert cfg.schema_version == 2
        assert [item.id for item in cfg.run.solvers][:2] == ["two_term", "multi_term"]
