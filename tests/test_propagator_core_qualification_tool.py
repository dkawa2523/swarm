from __future__ import annotations

from pathlib import Path
import json

import pytest

from tools import qualify_propagator_p1_deterministic as qualifier


def test_source_guard_rejects_changed_input(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "input.yaml"
    source.write_text("before\n", encoding="utf-8")
    monkeypatch.setattr(qualifier, "ROOT", tmp_path)
    monkeypatch.setattr(
        qualifier,
        "propagator_qualification_source_fingerprint",
        lambda _root: {"schema": "test", "sha256": "a" * 64},
    )
    inputs = (source,)
    starting_state = qualifier._qualification_source_state(inputs)

    source.write_text("after\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="changed during execution"):
        qualifier._require_unchanged_qualification_sources(starting_state, inputs)


def test_source_guard_rejects_changed_implementation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "input.yaml"
    source.write_text("unchanged\n", encoding="utf-8")
    monkeypatch.setattr(qualifier, "ROOT", tmp_path)
    fingerprints = iter(
        (
            {"schema": "test", "sha256": "a" * 64},
            {"schema": "test", "sha256": "b" * 64},
        )
    )
    monkeypatch.setattr(
        qualifier,
        "propagator_qualification_source_fingerprint",
        lambda _root: next(fingerprints),
    )
    inputs = (source,)
    starting_state = qualifier._qualification_source_state(inputs)

    with pytest.raises(RuntimeError, match="changed during execution"):
        qualifier._require_unchanged_qualification_sources(starting_state, inputs)


def test_quick_run_uses_separate_output_and_cannot_overwrite_formal_artifact(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    quick_output = tmp_path / "quick.json"
    formal_output = tmp_path / "formal.json"
    monkeypatch.setattr(qualifier, "QUICK_QUALIFICATION_OUTPUT", quick_output)
    monkeypatch.setattr(qualifier, "FULL_QUALIFICATION_OUTPUT", formal_output)
    monkeypatch.setattr(
        qualifier,
        "run_qualification",
        lambda **_kwargs: {
            "decision": {"p1_deterministic_core_qualified": True}
        },
    )

    assert qualifier.main(["--quick"]) == 0
    assert json.loads(quick_output.read_text(encoding="utf-8"))["decision"][
        "p1_deterministic_core_qualified"
    ]

    with pytest.raises(SystemExit):
        qualifier.main(["--quick", "--output", str(formal_output)])
