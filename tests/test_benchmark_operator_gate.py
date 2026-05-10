import builtins
import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "benchmark_operator_gate.py"


def _load_gate_module():
    spec = importlib.util.spec_from_file_location("benchmark_operator_gate", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_benchmark_gate_bolos_skip_is_optional(monkeypatch):
    gate = _load_gate_module()
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "bolos":
            raise ImportError("test forces missing BOLOS")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    scenario = gate.LoadedScenario("dummy", None, None, 50.0)

    rows = gate.run_bolos_gate([scenario], require_bolos=False)
    assert rows[0]["status"] == "skipped_optional"

    rows = gate.run_bolos_gate([scenario], require_bolos=True)
    assert rows[0]["status"] == "FAIL"


def test_benchmark_gate_lmax1_matches_native_two_term():
    gate = _load_gate_module()
    scenario = gate._load_scenarios(quick=True)[:1]

    rows = gate.run_lmax1_gate(scenario)

    assert rows
    assert {row["status"] for row in rows} == {"ok"}
    assert max(float(row["relerr"]) for row in rows) <= gate.LMAX1_TOLERANCE


def test_benchmark_gate_json_uses_fixed_columns(tmp_path: Path):
    gate = _load_gate_module()
    path = tmp_path / "gate.json"
    rows = [
        gate._row(
            "scenario",
            "candidate",
            "reference",
            "metric",
            1.0,
            1.0,
            0.1,
            "ok",
            "note",
        )
    ]

    gate.write_json(rows, path)

    data = json.loads(path.read_text(encoding="utf-8"))
    assert list(data[0]) == list(gate.REPORT_COLUMNS)
    assert data[0]["status"] == "ok"
