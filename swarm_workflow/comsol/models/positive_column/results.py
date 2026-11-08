"""Read native field CSVs and hash-bound positive-column result summaries."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from io import StringIO
import json
import math
from pathlib import Path
import re
from tempfile import NamedTemporaryFile
from typing import Any, TYPE_CHECKING
from xml.etree import ElementTree
from zipfile import BadZipFile, ZipFile

from .contracts import PositiveColumnWorkflowError
from .result_units import canonical_unit

if TYPE_CHECKING:
    from .config import ComsolModelMapping, ResultOutput

RESULT_SCHEMA = "swarm.positive_column_result.v1"
FIELD_CSV_MARKER = "SWARM_POSITIVE_COLUMN_FIELDS_CSV"
IDENTITY_COLUMNS = (
    "dataset",
    "solution_index",
    "point_index",
    "time_s",
    "parameter_name",
    "parameter_unit",
    "parameter_value",
)
_TEXT = {"dataset", "parameter_name", "parameter_unit"}
_NONNEGATIVE = {"electron_density", "excitation_source", "ionization_source"}


@dataclass(frozen=True, slots=True)
class PositiveColumnResults:
    rows: tuple[dict[str, Any], ...]
    units: dict[str, str]
    metadata: dict[str, Any]


def read_saved_argon_context(
    model: Path, mapping: ComsolModelMapping, selection: ResultOutput
) -> dict[str, Any] | None:
    """Recognize a saved positive-column model with constant pure-Ar gas data."""
    try:
        with ZipFile(model) as archive:
            root = ElementTree.fromstring(archive.read("dmodel.xml"))
    except (BadZipFile, KeyError, ElementTree.ParseError):
        return None
    physics = [
        node for node in root.iter("Physics")
        if node.get("op") == "ColdPlasma" and node.get("tag") == mapping.model.physics
    ]
    datasets = [
        node for node in root.iter("DatasetFeature")
        if node.get("op") == "Solution" and node.get("tag") == selection.dataset
    ]
    if len(physics) != 1 or len(datasets) != 1:
        return None
    component = datasets[0].find("propertyValue[@name='p:comp']")
    if component is None or component.get("value") != mapping.model.component:
        return None

    features = list(physics[0].iter("PhysicsFeature"))
    gas = [node for node in features if node.get("op") == "PlasmaEsModel"]
    contacts = [
        node for node in features
        if node.get("tag") == mapping.run.voltage_feature and node.get("op") == "MetalContact"
    ]
    if len(gas) != 1 or gas[0].get("tag") != mapping.closure.feature or len(contacts) != 1:
        return None

    def parameter(node: ElementTree.Element, name: str) -> str:
        entries = [entry for entry in node.findall("param") if entry.get("param") == name]
        encoded = entries[0].get("value", "") if len(entries) == 1 else ""
        value = re.fullmatch(r"1\|1,'([^']+)'", encoded)
        return value[1] if value is not None else ""

    species = [node for node in features if node.get("op") == "Species"]
    names = {parameter(node, "specName"): node for node in species}
    if len(species) != 4 or set(names) != {"e", "Ar", "Ars", "Ar_1p"}:
        return None
    if any(
        "DISABLED" in node.findtext("entityFlags", default="")
        for node in (physics[0], datasets[0], gas[0], contacts[0], *species)
    ):
        return None
    for name in ("Ar", "Ars", "Ar_1p"):
        if (
            parameter(names[name], "PresetSpeciesData") != "Ar"
            or parameter(names[name], "FromMassConstraint") != ("1" if name == "Ar" else "0")
        ):
            return None

    fields = {item.name: item.expression.strip() for item in selection.fields}
    pressure = parameter(gas[0], "pA")
    voltage = parameter(contacts[0], mapping.run.voltage_property)
    if (
        parameter(gas[0], "T_src") != "userdef"
        or not pressure or fields.get("gas_pressure") != pressure
        or not voltage or fields.get("applied_voltage") != voltage
    ):
        return None
    temperature_expression = parameter(gas[0], "T")
    mass_expression = parameter(names["Ar"], "M")
    if not temperature_expression.endswith("[K]") or not mass_expression.endswith("[kg/mol]"):
        return None
    try:
        temperature = float(temperature_expression.removesuffix("[K]"))
        mass = float(mass_expression.removesuffix("[kg/mol]"))
    except ValueError:
        return None
    if not all(math.isfinite(value) and value > 0 for value in (temperature, mass)):
        return None
    return {
        "gas_temperature_K": temperature,
        "gas_mixture": [{"species": "Ar", "fraction": 1.0, "molar_mass_kg_mol": mass}],
    }


def copy_comsol_csv(stdout: Path, output: Path, *, marker: str) -> bool:
    """Publish a complete CSV block from one COMSOL execution."""
    if not stdout.is_file():
        return False
    begin = f"{marker}_BEGIN".encode("ascii")
    end = f"{marker}_END".encode("ascii")
    temporary = None
    try:
        with stdout.open("rb") as source:
            if not any(line.rstrip(b"\r\n") == begin for line in source):
                return False
            output.parent.mkdir(parents=True, exist_ok=True)
            with NamedTemporaryFile(dir=output.parent, prefix=".csv-", delete=False) as stream:
                temporary = Path(stream.name)
                for line in source:
                    if line.rstrip(b"\r\n") == end:
                        break
                    stream.write(line)
                else:
                    return False
        temporary.replace(output)
        return True
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def read_mapped_results(
    path: str | Path, selection: ResultOutput
) -> PositiveColumnResults:
    return read_result_csv(
        path,
        units={field.name: field.unit for field in selection.fields},
        dataset=selection.dataset,
        solutions=selection.solutions,
    )


def read_result_csv(
    path: str | Path,
    *,
    units: dict[str, str],
    dataset: str,
    solutions: str,
    content: bytes | None = None,
) -> PositiveColumnResults:
    """Normalize declared numeric fields present; consumers require their own subset."""
    path = Path(path)
    content = path.read_bytes() if content is None else content
    reader = csv.DictReader(StringIO(content.decode("utf-8-sig"), newline=""))
    names = reader.fieldnames or []
    if not names or len(set(names)) != len(names):
        raise PositiveColumnWorkflowError(
            "result CSV has missing or duplicate column names"
        )
    if not set(IDENTITY_COLUMNS) <= set(names):
        raise PositiveColumnWorkflowError(
            "result CSV is missing saved-solution identity columns"
        )
    fields = [name for name in names if name not in IDENTITY_COLUMNS]
    if not set(fields) <= set(units):
        raise PositiveColumnWorkflowError(
            "result CSV has fields without declared units"
        )
    try:
        resolved = {name: canonical_unit(name, units[name]) for name in fields}
    except ValueError as error:
        raise PositiveColumnWorkflowError(str(error)) from error
    rows = []
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise PositiveColumnWorkflowError(
                "result CSV row width differs from its header"
            )
        parsed = {}
        for name, value in row.items():
            if name in _TEXT:
                parsed[name] = value or None
                continue
            if name in {"solution_index", "point_index"}:
                try:
                    number = int(value)
                except ValueError as error:
                    raise PositiveColumnWorkflowError(
                        "result CSV has invalid solution or point indices"
                    ) from error
                if number < (1 if name == "solution_index" else 0):
                    raise PositiveColumnWorkflowError(
                        "result CSV has invalid solution or point indices"
                    )
                parsed[name] = number
                continue
            if not value and name in {"time_s", "parameter_value"}:
                parsed[name] = None
                continue
            try:
                number = float(value) * resolved.get(name, (None, 1.0))[1]
            except ValueError as error:
                raise PositiveColumnWorkflowError(
                    f"result CSV field {name!r} is not numeric"
                ) from error
            if not math.isfinite(number):
                raise PositiveColumnWorkflowError(
                    f"result CSV field {name!r} is non-finite"
                )
            if name in _NONNEGATIVE and number < 0:
                raise PositiveColumnWorkflowError(
                    f"result CSV field {name!r} must be nonnegative"
                )
            parsed[name] = number
        rows.append(parsed)
    if not rows:
        raise PositiveColumnWorkflowError("positive-column result CSV is empty")
    _validate_selection(rows, dataset, solutions)
    return PositiveColumnResults(
        tuple(rows),
        {name: value[0] for name, value in resolved.items()},
        {
            "csv_sha256": sha256(content).hexdigest(),
            "source_units": {name: units[name] for name in fields},
        },
    )


def _validate_selection(
    rows: list[dict[str, Any]], dataset: str, solutions: str
) -> None:
    groups, contexts = {}, {}
    for row in rows:
        if row["dataset"] != dataset:
            raise PositiveColumnWorkflowError(
                "result CSV dataset differs from the selected dataset"
            )
        solution, point = row["solution_index"], row["point_index"]
        points = groups.setdefault(solution, set())
        if point in points:
            raise PositiveColumnWorkflowError(
                "result CSV has duplicate solution/point rows"
            )
        points.add(point)
        context = tuple(
            row[name]
            for name in (
                "time_s",
                "parameter_name",
                "parameter_unit",
                "parameter_value",
            )
        )
        if contexts.setdefault(solution, context) != context:
            raise PositiveColumnWorkflowError(
                "result CSV mixes parameter/time values within one solution; time_s must be constant"
            )
        if row["parameter_name"] is not None and row["parameter_value"] is None:
            raise PositiveColumnWorkflowError(
                "result CSV has a missing solution parameter value"
            )
        if row["parameter_name"] == "t" and row["time_s"] is None:
            raise PositiveColumnWorkflowError(
                "transient result CSV is missing physical time_s"
            )
        if row["parameter_name"] is None and (
            row["parameter_unit"] is not None or row["parameter_value"] is not None
        ):
            raise PositiveColumnWorkflowError(
                "result CSV has parameter data without a parameter name"
            )
    if solutions == "last" and len(groups) != 1:
        raise PositiveColumnWorkflowError(
            "last-solution export must contain one saved solution"
        )
    if solutions not in {"all", "last"}:
        raise PositiveColumnWorkflowError("result selection must be all or last")
    first = next(iter(groups.values()))
    if first != set(range(len(first))) or any(
        points != first for points in groups.values()
    ):
        raise PositiveColumnWorkflowError(
            "result CSV has incomplete spatial points across saved solutions"
        )


def read_positive_column_results(summary_path: str | Path) -> PositiveColumnResults:
    """Restore canonical field values and the original dataset/solution identities.

    A unitless CSV is not a result artifact. The summary and CSV can be moved
    together; the CSV filename is relative and its actual bytes are hash checked.
    """
    path = Path(summary_path).resolve()
    if path.suffix.lower() != ".json":
        raise PositiveColumnWorkflowError(
            "a bare CSV needs its unit and solution metadata. "
            "Pass positive_column_run_summary.json or positive_column_export.json"
        )
    encoded = path.read_bytes()
    metadata = json.loads(encoded)
    if not isinstance(metadata, dict):
        raise PositiveColumnWorkflowError("result summary must be a JSON object")
    if (
        metadata.get("artifact_schema") != RESULT_SCHEMA
        or metadata.get("status") != "completed"
    ):
        raise PositiveColumnWorkflowError(
            "expected a completed positive-column result artifact; re-export obsolete results with the current adapter"
        )
    filename = metadata.get("result_csv")
    if not isinstance(filename, str) or Path(filename).name != filename:
        raise PositiveColumnWorkflowError(
            "result CSV must be a filename inside its summary directory"
        )
    content = (path.parent / filename).read_bytes()
    if sha256(content).hexdigest() != metadata.get("result_csv_sha256"):
        raise PositiveColumnWorkflowError("result CSV is incomplete or altered")
    selection = metadata.get("result_selection")
    if (
        not isinstance(selection, dict)
        or not isinstance(selection.get("field_units"), dict)
        or not isinstance(selection.get("dataset"), str)
        or selection.get("solutions") not in ("all", "last")
    ):
        raise PositiveColumnWorkflowError(
            "result summary is missing units or saved-solution selection"
        )
    result = read_result_csv(
        path.parent / filename,
        content=content,
        units=selection["field_units"],
        dataset=selection["dataset"],
        solutions=selection["solutions"],
    )
    if set(result.units) != set(selection["field_units"]):
        raise PositiveColumnWorkflowError("result CSV omitted configured output fields")
    return PositiveColumnResults(
        result.rows,
        result.units,
        {
            **metadata,
            "summary_path": str(path),
            "summary_sha256": sha256(encoded).hexdigest(),
            "source_units": result.metadata["source_units"],
        },
    )
