import json
from pathlib import Path

import pandas as pd
import yaml

from swarm_comsol_export.converter import export_to_comsol


def test_export_minimal(tmp_path: Path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()

    # transport.csv
    pd.DataFrame({
        "epsbar_eV": [1.0, 2.0, 3.0],
        "mu_m2_V_s": [0.1, 0.2, 0.3],
        "De_m2_s": [1.0e-2, 2.0e-2, 3.0e-2],
    }).to_csv(run_dir / "transport.csv", index=False)

    # rates.csv
    pd.DataFrame({
        "epsbar_eV": [1.0, 2.0, 3.0],
        "k_ion": [1e-15, 2e-15, 3e-15],
        "k_attach": [5e-16, 6e-16, 7e-16],
    }).to_csv(run_dir / "rates.csv", index=False)

    # eedf.csv stacked
    rows = []
    for epsbar in [1.0, 2.0, 3.0]:
        for eps in [0.0, 1.0, 2.0]:
            rows.append({"eps_eV": eps, "epsbar_eV": epsbar, "f_eV_m32": 1.0})
    pd.DataFrame(rows).to_csv(run_dir / "eedf.csv", index=False)

    cfg = {
        "input": {
            "transport_csv": "transport.csv",
            "rates_csv": "rates.csv",
            "eedf": {
                "mode": "stacked_csv",
                "stacked_csv": "eedf.csv",
                "columns": {
                    "eps_eV": ["eps_eV"],
                    "epsbar_eV": ["epsbar_eV"],
                    "f": ["f_eV_m32"],
                },
            },
            "columns": {
                "epsbar_eV": ["epsbar_eV"],
                "mu": ["mu_m2_V_s"],
                "De": ["De_m2_s"],
                "rate_prefix": "k_",
            },
        },
        "processing": {"sort_by_epsbar": True},
        "output": {"delimiter": ",", "float_format": "%.6e"},
    }

    out_dir = tmp_path / "comsol"
    paths = export_to_comsol(run_dir, out_dir, cfg["input"], cfg["processing"], cfg["output"])

    assert paths.transport_csv.exists()
    assert paths.rates_csv.exists()
    assert paths.eedf_csv.exists()
    assert paths.report_json.exists()

    rep = json.loads(paths.report_json.read_text())
    assert rep["transport"]["rows"] == 3
