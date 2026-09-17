"""Licensed post-apply readback for the GEC-ICP closure contract."""

from __future__ import annotations

from hashlib import sha256
import math
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json
from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.java import java_path, java_string

from .java import TRANSPORT_TABLE, read_mobility_table
from .run_contracts import (
    GecIcpClosureReadbackPlan,
    GecIcpRunMapping,
    GecIcpWorkflowError,
)


READBACK_CLASS = "SwarmGecIcpClosureReadback"
_CONTRACT_MARKER = "SWARM_ICP_CLOSURE_CONTRACT"
_TABLE_MARKER = "SWARM_ICP_MOBILITY_TABLE"
_VALUE_MARKER = "SWARM_ICP_MOBILITY_VALUE"
_SCHEMA = "swarm.gec_icp_closure_readback.v1"
_ANCHOR_TOLERANCE = 1.0e-10
_MEAN_ENERGY_IDENTITY_TOLERANCE = 1.0e-10
_MOBILITY_TAG = "sw_icp_log_muN"
_VARIABLE_TAG = "swIcpClosureVars"
_VARIABLE_NAME = "sw_icp_logeps"


def prepare_gec_icp_closure_readback(
    mapping: GecIcpRunMapping,
    *,
    output_directory: Path,
    write_java: bool,
) -> GecIcpClosureReadbackPlan:
    """Create the immutable audit inputs and optional generated Java source."""

    transport_path = mapping.bundle.path / TRANSPORT_TABLE
    energies, _ = read_mobility_table(transport_path)
    output_directory.mkdir(parents=True, exist_ok=True)
    plan = GecIcpClosureReadbackPlan(
        mapping=mapping,
        java_path=output_directory / f"{READBACK_CLASS}.java",
        report_path=output_directory / "closure_readback.json",
        transport_path=transport_path,
        transport_sha256=_file_sha256(transport_path),
        anchor_count=len(energies),
    )
    if write_java:
        plan.java_path.write_text(
            generate_closure_readback_java(plan),
            encoding="utf-8",
        )
    return plan


def generate_closure_readback_java(plan: GecIcpClosureReadbackPlan) -> str:
    """Read the saved output MPH and emit a provenance-bound sentinel stream."""

    if _file_sha256(plan.transport_path) != plan.transport_sha256:
        raise GecIcpWorkflowError(
            "GEC ICP mobility table changed after closure-readback planning"
        )
    energies, _ = read_mobility_table(plan.transport_path)
    if len(energies) != plan.anchor_count:
        raise GecIcpWorkflowError(
            "GEC ICP mobility anchor count changed after readback planning"
        )
    mapping = plan.mapping
    model = mapping.model_mapping.model
    reactions = mapping.model_mapping.reactions
    component = java_string(model.component)
    physics = java_string(model.plasma_physics)
    plasma_feature = java_string(model.plasma_feature)
    dataset = java_string(model.dataset)
    model_path = java_string(java_path(mapping.run.output_mph))
    transport_path = java_string(java_path(plan.transport_path))
    log_energies = ", ".join(f"{math.log(value):.17g}" for value in energies)
    support = MeanEnergyArgument(energies[0], energies[-1])
    log_argument = support.expression("En-Ne")
    log_ebar_argument = support.expression("log(plas.ebar/1[V])")
    relative_expression = "abs(exp(En-Ne)*1[V]-plas.ebar)/max(abs(plas.ebar),1e-30[V])"
    support_expression = f"abs(({log_argument})-({log_ebar_argument}))"
    physics_target = f"model.component({component}).physics({physics})"
    feature_target = f"{physics_target}.feature({plasma_feature})"
    lines = [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "",
        f"public class {READBACK_CLASS} {{",
        "  private static void contract(String key, String value) {",
        f'    System.out.println("{_CONTRACT_MARKER}\\t" + key + "\\t"',
        '        + (value == null ? "" : value));',
        "  }",
        "",
        "  private static String joined(String[] values) {",
        '    return values == null ? "" : String.join("|", values);',
        "  }",
        "",
        "  private static double maxFinite(double[] values) {",
        "    if (values == null || values.length == 0) {",
        '      throw new IllegalStateException("empty mean-energy identity result");',
        "    }",
        "    double maximum = 0.0;",
        "    for (double value : values) {",
        "      if (!Double.isFinite(value)) {",
        '        throw new IllegalStateException("nonfinite mean-energy identity result");',
        "      }",
        "      maximum = Math.max(maximum, Math.abs(value));",
        "    }",
        "    return maximum;",
        "  }",
        "",
        "  public static void main(String[] ignored) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string("swarmGecIcpClosureReadback")
        + ", "
        + model_path
        + ");",
        f'    contract("model_path", {model_path});',
        f'    contract("transport_path", {transport_path});',
        '    contract("transport_sha256", ' + java_string(plan.transport_sha256) + ");",
        f'    contract("source", {java_string(mapping.bundle.expected_source)});',
        f'    contract("anchor_count", Integer.toString({plan.anchor_count}));',
        '    contract("ElectronProperties.ReducedProps", Boolean.toString('
        f'{physics_target}.prop("ElectronProperties").getBoolean("ReducedProps")));',
        '    contract("ElectronProperties.TensorElectronProps", Boolean.toString('
        f'{physics_target}.prop("ElectronProperties").getBoolean("TensorElectronProps")));',
        '    contract("ElectronProperties.IncludeThermalDiffusion", Boolean.toString('
        f'{physics_target}.prop("ElectronProperties").getBoolean("IncludeThermalDiffusion")));',
        '    contract("ElectronProperties.MeanElectronEnergyModel", '
        f'{physics_target}.prop("ElectronProperties")'
        '.getString("MeanElectronEnergyModel"));',
        '    contract("pes1.SpecifyElectronDensityAndEnergy", '
        f'{feature_target}.getString("SpecifyElectronDensityAndEnergy"));',
        f'    contract("pes1.muN", joined({feature_target}.getStringArray("muN")));',
        f'    contract("{_VARIABLE_TAG}.{_VARIABLE_NAME}", '
        f'model.component({component}).variable("{_VARIABLE_TAG}")'
        f'.get("{_VARIABLE_NAME}"));',
        f'    FunctionFeature mobility = model.func("{_MOBILITY_TAG}");',
        f'    contract("{_MOBILITY_TAG}.source", mobility.getString("source"));',
        f'    contract("{_MOBILITY_TAG}.argunit", joined(mobility.getStringArray("argunit")));',
        f'    contract("{_MOBILITY_TAG}.fununit", joined(mobility.getStringArray("fununit")));',
        f'    contract("{_MOBILITY_TAG}.interp", mobility.getString("interp"));',
        f'    contract("{_MOBILITY_TAG}.extrap", mobility.getString("extrap"));',
        '    String[][] mobilityTable = mobility.getStringMatrix("table");',
        f'    contract("{_MOBILITY_TAG}.table_rows", '
        "Integer.toString(mobilityTable == null ? 0 : mobilityTable.length));",
        "    if (mobilityTable != null) {",
        "      for (int index = 0; index < mobilityTable.length; index++) {",
        "        String[] row = mobilityTable[index];",
        "        if (row == null || row.length != 2) {",
        '          throw new IllegalStateException("malformed mobility table row");',
        "        }",
        f'        System.out.println("{_TABLE_MARKER}\\t" + index + "\\t"',
        '            + row[0] + "\\t" + row[1]);',
        "      }",
        "    }",
    ]
    for role, tag in (
        ("elastic", reactions.elastic),
        ("excitation", reactions.excitation),
        ("superelastic", reactions.superelastic),
        ("ionization", reactions.ionization),
        ("stepwise_ionization", reactions.stepwise_ionization),
    ):
        target = f"{physics_target}.feature({java_string(tag)})"
        prefix = f"reaction.{role}.{tag}"
        lines.extend(
            [
                f'    contract("{prefix}.SpecifyReactionUsing", '
                f'{target}.getString("SpecifyReactionUsing"));',
                f'    contract("{prefix}.eedf", {target}.getString("eedf"));',
                f'    contract("{prefix}.UseTownsendTwoTermsBoltzmann", '
                f'Boolean.toString({target}.getBoolean("UseTownsendTwoTermsBoltzmann")));',
            ]
        )
    lines.extend(
        [
            f"    double[] logEnergies = new double[]{{{log_energies}}};",
            "    for (int index = 0; index < logEnergies.length; index++) {",
            "      double logEnergy = logEnergies[index];",
            '      model.param().set("sw_icp_readback_value", '
            f'"{_MOBILITY_TAG}(" + Double.toString(logEnergy) + ")");',
            '      double value = model.param().evaluate("sw_icp_readback_value");',
            "      if (!Double.isFinite(value)) {",
            '        throw new IllegalStateException("nonfinite mobility anchor value");',
            "      }",
            f'      System.out.println("{_VALUE_MARKER}\\t" + index + "\\t"',
            '          + Double.toString(logEnergy) + "\\t" + Double.toString(value));',
            "    }",
            f"    int[] plasmaDomains = {physics_target}.selection().entities(2);",
            "    if (plasmaDomains.length == 0) {",
            '      throw new IllegalStateException("plasma domain selection is empty");',
            "    }",
            '    String identityTag = "swIcpClosureMeanIdentity";',
            "    try { model.result().numerical().remove(identityTag); }",
            "    catch (Exception missing) {}",
            '    model.result().numerical().create(identityTag, "MaxSurface");',
            f'    model.result().numerical(identityTag).set("data", {dataset});',
            '    model.result().numerical(identityTag).set("expr", new String[]{'
            + java_string(relative_expression)
            + ", "
            + java_string(support_expression)
            + "});",
            '    model.result().numerical(identityTag).set("unit", '
            'new String[]{"1", "1"});',
            "    model.result().numerical(identityTag).selection().set(plasmaDomains);",
            "    double[][] identity = model.result().numerical(identityTag).getReal();",
            "    if (identity == null || identity.length != 2) {",
            '      throw new IllegalStateException("incomplete mean-energy identity result");',
            "    }",
            '    contract("saved_solution.mean_energy_relative_error_max", '
            "Double.toString(maxFinite(identity[0])));",
            '    contract("saved_solution.support_log_difference_max", '
            "Double.toString(maxFinite(identity[1])));",
            '    ModelUtil.remove("swarmGecIcpClosureReadback");',
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def audit_gec_icp_closure_readback(
    plan: GecIcpClosureReadbackPlan,
    stdout_path: str | Path,
) -> dict[str, Any]:
    """Parse one sentinel stream, evaluate every gate, and save canonical JSON."""

    failures: list[str] = []
    contract: dict[str, str] = {}
    table_rows: dict[int, tuple[float, float]] = {}
    values: dict[int, tuple[float, float]] = {}
    try:
        contract, table_rows, values = _parse_readback_log(Path(stdout_path))
        energies, mobilities = read_mobility_table(plan.transport_path)
    except (OSError, UnicodeError, ValueError, GecIcpWorkflowError) as exc:
        failures.append(f"invalid_readback_log:{exc}")
        energies, mobilities = [], []

    current_sha256 = (
        _file_sha256(plan.transport_path) if plan.transport_path.is_file() else None
    )
    expected_contract = _expected_contract(plan, energies)
    expected_keys = set(expected_contract) | {
        "saved_solution.mean_energy_relative_error_max",
        "saved_solution.support_log_difference_max",
    }
    missing = sorted(expected_keys - set(contract))
    extra = sorted(set(contract) - expected_keys)
    if missing:
        failures.append("missing_contract_keys:" + ",".join(missing))
    if extra:
        failures.append("unexpected_contract_keys:" + ",".join(extra))
    mismatches = sorted(
        key
        for key, expected in expected_contract.items()
        if key in contract and contract[key] != expected
    )
    if mismatches:
        failures.append("contract_value_mismatch:" + ",".join(mismatches))
    if current_sha256 != plan.transport_sha256:
        failures.append("transport_sha256_changed")
    if len(energies) != plan.anchor_count:
        failures.append("transport_anchor_count_changed")

    expected_indices = set(range(len(energies)))
    if set(table_rows) != expected_indices:
        failures.append("mobility_table_indices_incomplete")
    if set(values) != expected_indices:
        failures.append("mobility_value_indices_incomplete")

    anchor_rows: list[dict[str, Any]] = []
    maximum_table_error = 0.0
    maximum_value_error = 0.0
    for index, (energy, mobility) in enumerate(zip(energies, mobilities, strict=True)):
        expected_log_energy = math.log(energy)
        expected_log_mobility = math.log(mobility)
        table = table_rows.get(index)
        evaluated = values.get(index)
        table_errors = (
            [
                abs(table[0] - expected_log_energy),
                abs(table[1] - expected_log_mobility),
            ]
            if table is not None
            else [math.inf, math.inf]
        )
        value_errors = (
            [
                abs(evaluated[0] - expected_log_energy),
                abs(evaluated[1] - expected_log_mobility),
            ]
            if evaluated is not None
            else [math.inf, math.inf]
        )
        maximum_table_error = max(maximum_table_error, *table_errors)
        maximum_value_error = max(maximum_value_error, *value_errors)
        anchor_rows.append(
            {
                "index": index,
                "mean_energy_eV": energy,
                "expected_log_mean_energy": expected_log_energy,
                "expected_log_reduced_mobility": expected_log_mobility,
                "table_log_mean_energy": table[0] if table else None,
                "table_log_reduced_mobility": table[1] if table else None,
                "evaluated_log_mean_energy": evaluated[0] if evaluated else None,
                "evaluated_log_reduced_mobility": evaluated[1] if evaluated else None,
            }
        )
    if (
        not math.isfinite(maximum_table_error)
        or maximum_table_error > _ANCHOR_TOLERANCE
    ):
        failures.append("mobility_table_anchor_mismatch")
    if (
        not math.isfinite(maximum_value_error)
        or maximum_value_error > _ANCHOR_TOLERANCE
    ):
        failures.append("mobility_evaluated_anchor_mismatch")

    relative_error = _optional_nonnegative_float(
        contract.get("saved_solution.mean_energy_relative_error_max")
    )
    support_log_error = _optional_nonnegative_float(
        contract.get("saved_solution.support_log_difference_max")
    )
    if relative_error is None or relative_error > _MEAN_ENERGY_IDENTITY_TOLERANCE:
        failures.append("saved_solution_mean_energy_identity_mismatch")
    if support_log_error is None or support_log_error > _MEAN_ENERGY_IDENTITY_TOLERANCE:
        failures.append("saved_solution_support_argument_mismatch")

    unique_failures = list(dict.fromkeys(failures))
    gates = {
        "saved_contract": not any(
            reason.startswith(
                (
                    "invalid_readback_log",
                    "missing_contract_keys",
                    "unexpected_contract_keys",
                    "contract_value_mismatch",
                )
            )
            for reason in unique_failures
        ),
        "transport_provenance": not any(
            reason in {"transport_sha256_changed", "transport_anchor_count_changed"}
            for reason in unique_failures
        ),
        "mobility_table_anchors": "mobility_table_anchor_mismatch"
        not in unique_failures
        and "mobility_table_indices_incomplete" not in unique_failures,
        "mobility_evaluated_anchors": "mobility_evaluated_anchor_mismatch"
        not in unique_failures
        and "mobility_value_indices_incomplete" not in unique_failures,
        "saved_solution_mean_energy_identity": not any(
            reason
            in {
                "saved_solution_mean_energy_identity_mismatch",
                "saved_solution_support_argument_mismatch",
            }
            for reason in unique_failures
        ),
    }
    passed = not unique_failures and all(gates.values())
    report = {
        "schema": _SCHEMA,
        "stage": "audit-gec-icp-closure-readback",
        "status": "passed" if passed else "failed",
        "passed": passed,
        "source": plan.mapping.bundle.expected_source,
        "output_mph": str(plan.mapping.run.output_mph),
        "transport": {
            "path": str(plan.transport_path),
            "planned_sha256": plan.transport_sha256,
            "current_sha256": current_sha256,
            "planned_anchor_count": plan.anchor_count,
            "current_anchor_count": len(energies),
        },
        "contract": {
            "expected": expected_contract,
            "readback": contract,
        },
        "mobility_anchor_audit": {
            "absolute_tolerance": _ANCHOR_TOLERANCE,
            "maximum_table_log_error": maximum_table_error,
            "maximum_evaluated_log_error": maximum_value_error,
            "rows": anchor_rows,
        },
        "saved_solution_mean_energy_identity": {
            "relative_error_max": relative_error,
            "support_log_difference_max": support_log_error,
            "limit": _MEAN_ENERGY_IDENTITY_TOLERANCE,
        },
        "gates": gates,
        "failure_reasons": unique_failures,
    }
    write_json(plan.report_path, report)
    return report


def _expected_contract(
    plan: GecIcpClosureReadbackPlan,
    energies: list[float],
) -> dict[str, str]:
    mapping = plan.mapping
    reactions = mapping.model_mapping.reactions
    support_expression = (
        MeanEnergyArgument(energies[0], energies[-1]).expression("En-Ne")
        if len(energies) >= 2
        else ""
    )
    mobility = "exp(sw_icp_log_muN(sw_icp_logeps))*1[1/(V*m*s)]"
    zero = "0[1/(V*m*s)]"
    # COMSOL serializes the isotropic 3-D tensor as the six symmetric
    # components xx, xy, yy, xz, yz, zz after reopening the MPH.
    tensor = [mobility, zero, mobility, zero, zero, mobility]
    expected = {
        "model_path": java_path(mapping.run.output_mph),
        "transport_path": java_path(plan.transport_path),
        "transport_sha256": plan.transport_sha256,
        "source": mapping.bundle.expected_source,
        "anchor_count": str(plan.anchor_count),
        "ElectronProperties.ReducedProps": "true",
        "ElectronProperties.TensorElectronProps": "false",
        "ElectronProperties.IncludeThermalDiffusion": "false",
        "ElectronProperties.MeanElectronEnergyModel": "LocalEnergyApproximationE",
        "pes1.SpecifyElectronDensityAndEnergy": "SpecifyMueOnly",
        "pes1.muN": "|".join(tensor),
        f"{_VARIABLE_TAG}.{_VARIABLE_NAME}": support_expression,
        f"{_MOBILITY_TAG}.source": "table",
        f"{_MOBILITY_TAG}.argunit": "1",
        f"{_MOBILITY_TAG}.fununit": "1",
        f"{_MOBILITY_TAG}.interp": "piecewisecubic",
        f"{_MOBILITY_TAG}.extrap": "const",
        f"{_MOBILITY_TAG}.table_rows": str(plan.anchor_count),
    }
    for role, tag in (
        ("elastic", reactions.elastic),
        ("excitation", reactions.excitation),
        ("superelastic", reactions.superelastic),
        ("ionization", reactions.ionization),
        ("stepwise_ionization", reactions.stepwise_ionization),
    ):
        prefix = f"reaction.{role}.{tag}"
        expected[f"{prefix}.SpecifyReactionUsing"] = "UseCrossSectionData"
        expected[f"{prefix}.eedf"] = "FromPhysicsInterfaceProperty"
        expected[f"{prefix}.UseTownsendTwoTermsBoltzmann"] = "false"
    return expected


def _parse_readback_log(
    path: Path,
) -> tuple[
    dict[str, str],
    dict[int, tuple[float, float]],
    dict[int, tuple[float, float]],
]:
    if not path.is_file():
        raise GecIcpWorkflowError(f"missing GEC ICP closure readback log: {path}")
    contract: dict[str, str] = {}
    tables: dict[int, tuple[float, float]] = {}
    values: dict[int, tuple[float, float]] = {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        if line.startswith(_CONTRACT_MARKER + "\t"):
            parts = line.split("\t", 2)
            if len(parts) != 3 or not parts[1] or parts[1] in contract:
                raise GecIcpWorkflowError(
                    "duplicate or malformed closure contract sentinel"
                )
            contract[parts[1]] = parts[2]
        elif line.startswith(_TABLE_MARKER + "\t"):
            index, x_value, y_value = _indexed_numeric_sentinel(line, _TABLE_MARKER)
            if index in tables:
                raise GecIcpWorkflowError("duplicate mobility-table sentinel index")
            tables[index] = (x_value, y_value)
        elif line.startswith(_VALUE_MARKER + "\t"):
            index, x_value, y_value = _indexed_numeric_sentinel(line, _VALUE_MARKER)
            if index in values:
                raise GecIcpWorkflowError("duplicate mobility-value sentinel index")
            values[index] = (x_value, y_value)
    if not contract:
        raise GecIcpWorkflowError("closure readback log contains no contract sentinels")
    return contract, tables, values


def _indexed_numeric_sentinel(line: str, marker: str) -> tuple[int, float, float]:
    parts = line.split("\t")
    if len(parts) != 4:
        raise GecIcpWorkflowError(f"malformed {marker} sentinel")
    try:
        index = int(parts[1])
        x_value = float(parts[2])
        y_value = float(parts[3])
    except ValueError as exc:
        raise GecIcpWorkflowError(f"nonnumeric {marker} sentinel") from exc
    if index < 0 or not math.isfinite(x_value) or not math.isfinite(y_value):
        raise GecIcpWorkflowError(f"invalid {marker} sentinel")
    return index, x_value, y_value


def _optional_nonnegative_float(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed) or parsed < 0.0:
        return None
    return parsed


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


__all__ = [
    "READBACK_CLASS",
    "audit_gec_icp_closure_readback",
    "generate_closure_readback_java",
    "prepare_gec_icp_closure_readback",
]
