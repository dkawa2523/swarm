from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

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
