"""Canonical workflow quality-policy types, codecs, and provenance checks."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from hashlib import sha256
import json
import math


QUALITY_THRESHOLDS_METADATA_KEY = "quality_thresholds_json"


@dataclass(frozen=True, slots=True)
class RequiredRateRse:
    """Reaction-rate gates that apply independently of major-rate ranking."""

    process_type: dict[str, float | None] = field(
        default_factory=lambda: {"excitation": None, "ionization": None}
    )
    process: dict[str, float | None] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class QualityThresholds:
    mobility_rse: float = 0.02
    diffusion_rse: float = 0.05
    major_rate_rse: float = 0.05
    major_rate_fraction: float = 0.01
    required_rate_min_process_peak_fraction: float = 0.0
    eedf_normalization_error: float = 1.0e-6
    required_rate_rse: RequiredRateRse = field(default_factory=RequiredRateRse)


def quality_thresholds_payload(thresholds: QualityThresholds) -> dict[str, object]:
    """Return the complete validated quality policy as a JSON-ready mapping."""

    return asdict(thresholds)


def quality_thresholds_json(thresholds: QualityThresholds) -> str:
    """Serialize the quality policy in its canonical database representation."""

    return json.dumps(
        quality_thresholds_payload(thresholds),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )

def quality_thresholds_sha256(thresholds: QualityThresholds) -> str:
    return sha256(quality_thresholds_json(thresholds).encode("utf-8")).hexdigest()


def quality_policy_provenance(
    source: QualityThresholds,
    evaluation: QualityThresholds,
    *,
    source_json: str | None = None,
) -> dict[str, object]:
    if source_json is None:
        source_json = quality_thresholds_json(source)
    else:
        try:
            source_payload = json.loads(source_json)
            parsed_source = parse_quality_thresholds(source_payload)
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("source quality policy JSON is invalid") from exc
        if parsed_source != source or source_json != quality_thresholds_json(source):
            raise ValueError("source quality policy JSON is noncanonical")
    evaluation_json = quality_thresholds_json(evaluation)
    return {
        "source": json.loads(source_json),
        "source_sha256": sha256(source_json.encode("utf-8")).hexdigest(),
        "evaluation": quality_thresholds_payload(evaluation),
        "evaluation_sha256": quality_thresholds_sha256(evaluation),
        "quality_only_reevaluation": (
            quality_thresholds_json(source) != evaluation_json
        ),
    }


def parse_quality_thresholds(raw: object) -> QualityThresholds:
    if raw is None:
        return QualityThresholds()
    if not isinstance(raw, dict):
        raise ValueError("workflow quality must be a mapping")
    allowed = {
        "mobility_rse",
        "diffusion_rse",
        "major_rate_rse",
        "major_rate_fraction",
        "required_rate_min_process_peak_fraction",
        "eedf_normalization_error",
        "required_rate_rse",
    }
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"Unsupported workflow quality fields: {sorted(unknown)}")
    values: dict[str, object] = {}
    for key, value in raw.items():
        if key == "required_rate_rse":
            values[key] = _parse_required_rate_rse(value)
            continue
        number = float(value)
        if not math.isfinite(number) or number < 0.0:
            raise ValueError(f"workflow quality.{key} must be finite and nonnegative")
        if key in {
            "major_rate_fraction",
            "required_rate_min_process_peak_fraction",
        } and number > 1.0:
            raise ValueError(f"workflow quality.{key} must be in [0, 1]")
        values[key] = number
    return QualityThresholds(**values)


def _parse_required_rate_rse(raw: object) -> RequiredRateRse:
    if raw is None:
        return RequiredRateRse()
    if not isinstance(raw, dict):
        raise ValueError("workflow quality.required_rate_rse must be a mapping")
    unknown = set(raw) - {"process", "process_type"}
    if unknown:
        raise ValueError(
            "Unsupported workflow quality.required_rate_rse fields: "
            f"{sorted(unknown)}"
        )

    def parse_selector(name: str) -> dict[str, float | None]:
        selector = raw.get(name, {})
        if not isinstance(selector, dict):
            raise ValueError(
                f"workflow quality.required_rate_rse.{name} must be a mapping"
            )
        parsed: dict[str, float | None] = {}
        normalized_names: set[str] = set()
        for raw_name, raw_limit in selector.items():
            gate_name = str(raw_name).strip()
            normalized = gate_name.casefold()
            if not gate_name:
                raise ValueError(
                    f"workflow quality.required_rate_rse.{name} names must be nonempty"
                )
            if normalized in normalized_names:
                raise ValueError(
                    f"workflow quality.required_rate_rse.{name} contains duplicate "
                    f"case-insensitive name {gate_name!r}"
                )
            normalized_names.add(normalized)
            if raw_limit is None:
                parsed[gate_name] = None
                continue
            limit = float(raw_limit)
            if not math.isfinite(limit) or limit < 0.0:
                raise ValueError(
                    f"workflow quality.required_rate_rse.{name}.{gate_name} "
                    "must be null or finite and nonnegative"
                )
            parsed[gate_name] = limit
        return parsed

    defaults = RequiredRateRse()
    return RequiredRateRse(
        process_type=(
            parse_selector("process_type")
            if "process_type" in raw
            else dict(defaults.process_type)
        ),
        process=(
            parse_selector("process")
            if "process" in raw
            else dict(defaults.process)
        ),
    )
