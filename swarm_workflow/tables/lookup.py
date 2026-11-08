"""Read and interpolate qualified DC coefficients at fixed gas conditions."""

from __future__ import annotations

from bisect import bisect_left
from collections.abc import Mapping
from copy import deepcopy
import csv
from dataclasses import dataclass, field
from hashlib import sha256
from io import StringIO
import json
import math
from pathlib import Path
from typing import Any

from electron_swarm.core.physical_context import PHYSICAL_CONTEXT_SCHEMA

from . import contracts


_FIELD = "E_over_N_Td"
_RATE = "rate_coefficient_m3_s"
_ELASTIC_LOSS = "elastic_energy_loss_rate_coefficient_eV_m3_s"
_TRANSPORT = frozenset(contracts.CASE_COLUMNS[2:])
_RESTRICTED_MC_QUANTITIES = frozenset(
    ("mean_energy_eV", "reduced_mobility_m2_V_s_m3", _ELASTIC_LOSS)
)


class CoefficientLookupError(ValueError):
    """A coefficient or its fixed physical context is unavailable."""


@dataclass(frozen=True, slots=True)
class Reaction:
    species: str
    process: str
    process_type: str
    threshold_eV: float | None


@dataclass(frozen=True, slots=True)
class CoefficientSeries:
    quantity: str
    unit: str
    source: str
    reaction: Reaction | None
    fields_Td: tuple[float, ...]
    values: tuple[float | None, ...]

    def at(self, e_over_n_Td: float) -> float:
        """Return an anchor or a linear interpolation between qualified neighbors."""
        target = _finite(e_over_n_Td, _FIELD)
        index = bisect_left(self.fields_Td, target)
        if index < len(self.fields_Td) and target == self.fields_Td[index]:
            value = self.values[index]
            if value is not None:
                return value
        elif 0 < index < len(self.fields_Td):
            lower, upper = self.values[index - 1 : index + 1]
            if lower is not None and upper is not None:
                left, right = self.fields_Td[index - 1 : index + 1]
                fraction = (target - left) / (right - left)
                return (1 - fraction) * lower + fraction * upper
        raise CoefficientLookupError(
            f"{self.quantity} has no qualified value or bracket at {target:g} Td"
        )


@dataclass(frozen=True, slots=True)
class CoefficientTables:
    directory: Path
    source: str
    gas_temperature_K: float
    gas_number_density_m3: float
    mixture: tuple[tuple[str, float], ...]
    manifest_sha256: str
    _manifest: dict[str, Any] = field(repr=False)

    @property
    def source_policy(self) -> dict[str, Any]:
        """Describe the producer's closure, transport definition, and evidence."""
        return deepcopy(self._manifest["source_policy"])

    @property
    def reactions(self) -> tuple[Reaction, ...]:
        rows = self._read_table("rates_vs_en.csv", (_FIELD, _RATE))
        return tuple(dict.fromkeys(_reaction(row) for row in rows))

    def series(
        self, quantity: str, *, reaction: Reaction | None = None
    ) -> CoefficientSeries:
        """Select a quantity by its canonical column and an exact reaction key."""
        if (quantity == _RATE) != (reaction is not None):
            raise CoefficientLookupError(
                "reaction selection is required only for rate_coefficient_m3_s"
            )
        policy = self._manifest["source_policy"]
        profile = policy.get("qualification_profile")
        if (
            self.source == contracts.MC_SOLVER
            and profile == contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            and quantity not in _RESTRICTED_MC_QUANTITIES
        ):
            raise CoefficientLookupError(
                f"{quantity} is outside the restricted MC qualification"
            )

        if quantity == _RATE:
            filename = "rates_vs_en.csv"
        elif quantity == _ELASTIC_LOSS:
            filename = contracts.ELASTIC_ENERGY_LOSS_TABLE
        elif quantity in _TRANSPORT:
            filename = "transport_vs_en.csv"
        else:
            raise CoefficientLookupError(f"unknown coefficient quantity {quantity!r}")
        rows = self._read_table(filename, (_FIELD, quantity))
        if reaction is not None:
            rows = [row for row in rows if _reaction(row) == reaction]
        fields, qualified = self._support(quantity)
        values: dict[float, float | None] = {}
        for row in rows:
            point = _finite(row[_FIELD], _FIELD)
            if point not in qualified or point in values:
                raise CoefficientLookupError(
                    f"unexpected or duplicate {quantity} point at {point:g} Td"
                )
            raw = row[quantity]
            value = _finite(raw, quantity) if raw != "" else None
            if quantity == _RATE and value is not None:
                if value < 0:
                    raise CoefficientLookupError(
                        "reaction rate coefficient must be nonnegative"
                    )
                if self.source == contracts.MC_SOLVER and value == 0:
                    # A zero-count MC estimate is censored evidence.
                    value = None
            values[point] = value
        if all(value is None for value in values.values()):
            raise CoefficientLookupError(
                f"no qualified values for {quantity} and reaction {reaction}"
            )
        return CoefficientSeries(
            quantity,
            contracts.UNITS[quantity],
            self.source,
            reaction,
            fields,
            tuple(values.get(point) for point in fields),
        )

    def _support(self, quantity: str) -> tuple[tuple[float, ...], set[float]]:
        if self.source == contracts.MC_SOLVER:
            component = (
                "eedf_rate" if quantity in {_RATE, _ELASTIC_LOSS} else "transport"
            )
            policy = self._manifest["source_policy"]
            eligible = policy.get(f"{component}_eligible_E_over_N_Td")
            excluded = policy.get(f"{component}_excluded_E_over_N_Td")
            if (
                not isinstance(eligible, list)
                or not eligible
                or not isinstance(excluded, dict)
            ):
                raise CoefficientLookupError(f"missing MC {component} anchor evidence")
            qualified = {_finite(point, _FIELD) for point in eligible}
            rejected = {_finite(point, _FIELD) for point in excluded}
            if len(qualified) != len(eligible) or qualified & rejected:
                raise CoefficientLookupError(
                    f"inconsistent MC {component} anchor evidence"
                )
            return tuple(sorted(qualified | rejected)), qualified

        rows = self._read_table("quality.csv", (_FIELD,), required=("passed",))
        fields = {_finite(row[_FIELD], _FIELD) for row in rows}
        if (
            len(fields) != len(rows)
            or not rows
            or any(row["passed"] != "1" for row in rows)
        ):
            raise CoefficientLookupError(
                "deterministic coefficient anchors must pass numerical quality"
            )
        return tuple(sorted(fields)), fields

    def _read_table(
        self, name: str, quantities: tuple[str, ...], *, required: tuple[str, ...] = ()
    ) -> list[dict[str, str]]:
        record = self._manifest["tables"].get(name)
        if not isinstance(record, dict) or record.get("argument") != _FIELD:
            raise CoefficientLookupError(f"{name} is unavailable on the Td axis")
        units = record.get("units")
        if not isinstance(units, dict) or any(
            units.get(key) != contracts.UNITS[key] for key in quantities
        ):
            raise CoefficientLookupError(f"{name} has missing or incompatible units")
        try:
            content = (self.directory / name).read_bytes()
            if sha256(content).hexdigest() != record.get("sha256"):
                raise CoefficientLookupError(
                    f"{name} SHA-256 differs from the manifest"
                )
            reader = csv.DictReader(StringIO(content.decode("utf-8-sig"), newline=""))
            columns = reader.fieldnames or []
            if (
                columns != record.get("columns")
                or len(set(columns)) != len(columns)
                or not set((*quantities, *required)).issubset(columns)
            ):
                raise CoefficientLookupError(
                    f"{name} columns differ from the declared coefficient table"
                )
            rows = list(reader)
        except (OSError, UnicodeError, csv.Error) as exc:
            raise CoefficientLookupError(f"cannot read {name}") from exc
        if any(None in row or None in row.values() for row in rows):
            raise CoefficientLookupError(f"{name} has a row with the wrong width")
        return rows


def load_coefficients(
    directory: str | Path,
    *,
    gas_temperature_K: float,
    gas_number_density_m3: float,
    mole_fractions: Mapping[str, float],
) -> CoefficientTables:
    """Bind a build-tables mixture to the caller's actual fixed gas state."""
    root = Path(directory).resolve()
    try:
        content = (root / "manifest.json").read_bytes()
        manifest = json.loads(content)
    except (OSError, ValueError) as exc:
        raise CoefficientLookupError(
            f"cannot read the mixture manifest in {root}"
        ) from exc
    if (
        not isinstance(manifest, dict)
        or manifest.get("stage") != "build-tables"
        or manifest.get("format_version") != contracts.FORMAT_VERSION
        or manifest.get("source") not in contracts.SOURCE_CHOICES
        or not isinstance(manifest.get("tables"), dict)
    ):
        raise CoefficientLookupError("choose a current build-tables mixture directory")
    context = manifest.get("physical_context")
    policy = manifest.get("source_policy")
    if (
        not isinstance(context, dict)
        or context.get("schema") != PHYSICAL_CONTEXT_SCHEMA
        or not isinstance(policy, dict)
        or policy.get("source") != manifest["source"]
        or policy.get("field_type") != "dc"
        or policy.get("postprocess") != "none"
    ):
        raise CoefficientLookupError(
            "coefficient lookup requires unmodified DC tables with physical context"
        )
    source = manifest["source"]
    if (
        source == contracts.MC_SOLVER
        and policy.get("qualification_profile")
        not in contracts.MC_QUALIFICATION_PROFILES
    ):
        raise CoefficientLookupError("unknown MC coefficient qualification profile")

    temperature = _finite(context.get("gas_temperature_K"), "table gas temperature")
    density = _finite(context.get("gas_number_density_m3"), "table gas number density")
    for label, requested, expected in (
        ("gas temperature in K", gas_temperature_K, temperature),
        ("gas number density in m^-3", gas_number_density_m3, density),
    ):
        if expected <= 0 or not math.isclose(
            _finite(requested, label), expected, rel_tol=1e-12
        ):
            raise CoefficientLookupError(
                f"{label} differs from table value {expected:.17g}"
            )
    try:
        species = manifest["mixture"]["species"]
        fractions = {item["species"]: item["fraction"] for item in species}
    except (KeyError, TypeError) as exc:
        raise CoefficientLookupError("manifest lacks a fixed gas mixture") from exc
    if len(fractions) != len(species):
        raise CoefficientLookupError("table gas mixture contains duplicate species")
    mixture = _active_fractions(fractions)
    requested = _active_fractions(mole_fractions)
    if requested.keys() != mixture.keys() or any(
        not math.isclose(value, mixture[name], rel_tol=1e-12)
        for name, value in requested.items()
    ):
        raise CoefficientLookupError("gas mole fractions differ from the table mixture")
    return CoefficientTables(
        root,
        source,
        temperature,
        density,
        tuple(sorted(mixture.items())),
        sha256(content).hexdigest(),
        manifest,
    )


def _reaction(row: dict[str, str]) -> Reaction:
    try:
        species, process, process_type, threshold = (
            row[key] for key in ("species", "process", "process_type", "threshold_eV")
        )
    except KeyError as exc:
        raise CoefficientLookupError(
            "rate table lacks the complete reaction identity"
        ) from exc
    if not species or not process or not process_type:
        raise CoefficientLookupError("rate table contains an empty reaction identity")
    return Reaction(
        species,
        process,
        process_type,
        _finite(threshold, "threshold_eV") if threshold else None,
    )


def _finite(value: object, label: str) -> float:
    try:
        if isinstance(value, bool):
            raise ValueError
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise CoefficientLookupError(f"{label} must be finite") from exc
    if not math.isfinite(number):
        raise CoefficientLookupError(f"{label} must be finite")
    return number


def _active_fractions(fractions: Mapping[str, float]) -> dict[str, float]:
    if not isinstance(fractions, Mapping):
        raise CoefficientLookupError(
            "gas mole fractions must map species names to fractions"
        )
    values = {key: _finite(value, "mole fraction") for key, value in fractions.items()}
    if any(
        not isinstance(key, str) or not key or not 0 <= value <= 1
        for key, value in values.items()
    ) or not math.isclose(math.fsum(values.values()), 1.0, abs_tol=1e-12, rel_tol=0):
        raise CoefficientLookupError(
            "gas mole fractions must be nonnegative and sum to one"
        )
    return {key: value for key, value in values.items() if value > 0}
