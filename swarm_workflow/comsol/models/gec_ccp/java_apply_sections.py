"""Cohesive section renderers for the GEC-CCP model-application source."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any

from swarm_workflow.comsol.java import java_path, java_string
import swarm_workflow.comsol.models.gec_ccp.closure_arguments as closure_arguments
import swarm_workflow.comsol.models.gec_ccp.java_support as java_support
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _active_rate_table,
    _active_transport_functions,
    _reaction_uses_preintegrated_rate,
    _thermal_diffusion_enabled,
    _transport_property_tensors,
    _uses_direct_mc_rates,
    _uses_external_elastic_energy_loss,
    _uses_external_rates,
    _uses_function_eedf,
    _uses_hybrid_einstein_transport,
    _uses_transport,
    _validate_restricted_closure_policy,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    AVOGADRO_PER_MOL,
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.data import _lookup, _read_csv
from swarm_workflow.comsol.models.gec_ccp.validation.bundle_guards import (
    GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    GEC_ELASTIC_ENERGY_LOSS_TABLE,
)
from swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency import (
    _active_closure_mean_energy_support,
    _build_dense_function_eedf_rate_closure,
)


@dataclass(frozen=True)
class ApplyJavaContext:
    """Validated inputs shared by every application-source section."""

    mapping: GecCcpMapping
    input_mph: Path
    transport_rows: list[dict[str, Any]]
    rate_rows: list[dict[str, Any]]
    uses_transport: bool
    uses_rates: bool
    uses_function_eedf: bool
    uses_external_elastic_loss: bool
    target: str


def build_apply_java_context(
    mapping: GecCcpMapping,
    *,
    input_mph: Path | None,
    preintegrated_rate_rows: list[dict[str, Any]] | None,
) -> ApplyJavaContext:
    """Load and validate all table inputs before rendering Java."""

    uses_transport = _uses_transport(mapping.closure)
    uses_rates = _uses_external_rates(mapping.closure)
    uses_function_eedf = _uses_function_eedf(mapping.closure)
    uses_external_elastic_loss = _uses_external_elastic_energy_loss(mapping.closure)
    transport = mapping.bundle.path / "transport_vs_mean_energy.csv"
    rows = _read_csv(transport) if uses_transport else []
    _validate_restricted_closure_policy(
        mapping,
        source=mapping.bundle.expected_source,
        transport_rows=rows if uses_transport else None,
    )
    rate_rows: list[dict[str, Any]] = []
    if uses_rates:
        if mapping.closure.reaction_model == FUNCTION_EEDF_PREINTEGRATED_INELASTIC:
            if preintegrated_rate_rows is None:
                preintegrated_rate_rows, _ = _build_dense_function_eedf_rate_closure(
                    mapping,
                    input_mph=input_mph or mapping.model.input_mph,
                )
            rate_rows = preintegrated_rate_rows
        else:
            rate_rows = _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
    target = (
        f"model.component({java_string(mapping.model.component)})"
        f".physics({java_string(mapping.model.physics)})"
        f".feature({java_string(mapping.model.plasma_feature)})"
    )
    return ApplyJavaContext(
        mapping=mapping,
        input_mph=input_mph or mapping.model.input_mph,
        transport_rows=rows,
        rate_rows=rate_rows,
        uses_transport=uses_transport,
        uses_rates=uses_rates,
        uses_function_eedf=uses_function_eedf,
        uses_external_elastic_loss=uses_external_elastic_loss,
        target=target,
    )


def render_stabilization(context: ApplyJavaContext) -> list[str]:
    mapping = context.mapping
    component = java_string(mapping.model.component)
    physics = java_string(mapping.model.physics)
    return [
        "    // Stabilization is an explicit model input, not a hidden",
        "    // convergence aid. The schema default is OFF/OFF.",
        f"    model.component({component}).physics({physics})"
        '.prop("Stabilization").set("SourceStabilization", '
        + java_support.java_boolean(mapping.run.source_stabilization)
        + ");",
        f"    model.component({component}).physics({physics})"
        '.prop("Stabilization").set("ReactionSourceStabilization", '
        + java_support.java_boolean(mapping.run.reaction_source_stabilization)
        + ");",
    ]


def render_transport_mode(context: ApplyJavaContext) -> list[str]:
    if not context.uses_transport:
        return []
    mapping = context.mapping
    component = java_string(mapping.model.component)
    physics = java_string(mapping.model.physics)
    lines = [
        "    // Swarm exports reduced transport coefficients. Use",
        "    // COMSOL's native reduced-property contract so COMSOL",
        "    // owns the neutral-density conversion and Jacobian.",
        f"    model.component({component}).physics({physics})"
        '.prop("ElectronProperties").set("ReducedProps", true);',
        "    // Manual Swarm tensors are independent inputs. Do not",
        "    // ask COMSOL to synthesize magnetized transport. The",
        "    // thermal-diffusion option is an explicit closure choice.",
        f"    model.component({component}).physics({physics})"
        '.prop("ElectronProperties").set("TensorElectronProps", false);',
        f"    model.component({component}).physics({physics})"
        '.prop("ElectronProperties").set("IncludeThermalDiffusion", '
        + java_support.java_boolean(_thermal_diffusion_enabled(mapping.closure))
        + ");",
    ]
    if mapping.closure.electron_transport == "swarm_mobility_einstein":
        lines.extend(
            [
                "    // Import Swarm electron mobility while retaining",
                "    // COMSOL's mean-energy-dependent Einstein diffusion",
                "    // and energy-flux closure.",
                f"    {context.target}.set("
                '"SpecifyElectronDensityAndEnergy", "SpecifyMueOnly");',
            ]
        )
    elif _uses_hybrid_einstein_transport(mapping.closure):
        lines.extend(
            [
                "    // Qualified standard-local-energy hybrid: import",
                "    // mobility and solver-native energy transport, but",
                "    // evaluate particle diffusion",
                "    // from DeN=muN*Te inside COMSOL. External De is not",
                "    // created or bound by this closure.",
                f"    {context.target}.set("
                '"SpecifyElectronDensityAndEnergy", "SpecifyAll");',
            ]
        )
    else:
        lines.extend(
            [
                "    // Preserve COMSOL's solved mean energy while importing",
                "    // the complete local Swarm particle/energy transport.",
                f"    {context.target}.set("
                '"SpecifyElectronDensityAndEnergy", "SpecifyAll");',
            ]
        )
    return lines


def render_closure_arguments(context: ApplyJavaContext) -> list[str]:
    if not (
        context.uses_transport
        or context.uses_rates
        or context.uses_external_elastic_loss
    ):
        return []
    mapping = context.mapping
    closure_support = _active_closure_mean_energy_support(
        mapping,
        preintegrated_rate_rows=(context.rate_rows if context.uses_rates else None),
    )
    if not closure_support["passed"]:
        raise GecCcpWorkflowError(
            "active transport/rate tables have no usable common mean-"
            f"energy support: {closure_support.get('reason')}"
        )
    component = java_string(mapping.model.component)
    lines = [
        "    // Evaluate smooth closure arguments in one variable node.",
        "    // COMSOL defines ne=exp(Ne_per) and en=exp(En_per).",
        "    // Their log ratio is evaluated directly so nonphysical",
        "    // Newton trials cannot overflow exp(En)/exp(Ne).",
        f"    try {{ model.component({component})"
        '.variable("swClosureVars"); } catch (Exception ex) {',
        f'      model.component({component}).variable().create("swClosureVars");',
        "    }",
    ]
    if mapping.bundle.expected_source == "monte_carlo":
        lines.extend(_render_mc_closure_arguments(context, closure_support))
    else:
        support_minimum, support_maximum = closure_support[
            "common_intersection_mean_energy_eV"
        ]
        log_energy_expression = closure_arguments.smooth_log_energy_argument(
            "En_per-Ne_per",
            minimum_eV=support_minimum,
            maximum_eV=support_maximum,
        )
        lines.append(
            f"    model.component({component})"
            '.variable("swClosureVars").set("sw_logeps", '
            + java_string(log_energy_expression)
            + ");"
        )
    return lines


def _render_mc_closure_arguments(
    context: ApplyJavaContext,
    closure_support: dict[str, Any],
) -> list[str]:
    mapping = context.mapping
    component = java_string(mapping.model.component)
    lines: list[str] = []
    transport_support = closure_support.get("transport_intersection_mean_energy_eV")
    if context.uses_transport:
        if not isinstance(transport_support, list):
            raise GecCcpWorkflowError("Monte Carlo transport support is unavailable")
        transport_minimum, transport_maximum = closure_arguments.closure_argument_range(
            mapping, transport_support
        )
        transport_argument = closure_arguments.smooth_log_energy_argument(
            "En_per-Ne_per",
            minimum_eV=transport_minimum,
            maximum_eV=transport_maximum,
        )
        lines.extend(
            [
                "    // Transport and reaction tables retain independent",
                "    // supports; no reaction zero anchor truncates transport.",
                f"    model.component({component})"
                '.variable("swClosureVars").set('
                '"sw_logeps_transport", ' + java_string(transport_argument) + ");",
            ]
        )
    rate_supports = closure_support.get("rate_supports_mean_energy_eV", {})
    for reaction in mapping.reactions:
        if not (
            context.uses_rates
            and _reaction_uses_preintegrated_rate(mapping.closure, reaction)
        ):
            continue
        rate_support = rate_supports.get(reaction.process_type)
        if not isinstance(rate_support, list):
            raise GecCcpWorkflowError(
                f"Monte Carlo rate support is unavailable for {reaction.process_type}"
            )
        rate_minimum, rate_maximum = closure_arguments.closure_argument_range(
            mapping, rate_support
        )
        rate_argument = closure_arguments.smooth_log_energy_argument(
            "En_per-Ne_per",
            minimum_eV=rate_minimum,
            maximum_eV=rate_maximum,
        )
        lines.append(
            f"    model.component({component})"
            '.variable("swClosureVars").set('
            + java_string(f"sw_logeps_rate_{reaction.process_type}")
            + ", "
            + java_string(rate_argument)
            + ");"
        )
    if context.uses_external_elastic_loss:
        elastic_support = closure_support.get(
            "elastic_energy_loss_support_mean_energy_eV"
        )
        if not isinstance(elastic_support, list):
            raise GecCcpWorkflowError(
                "Monte Carlo elastic energy-loss support is unavailable"
            )
        elastic_minimum, elastic_maximum = closure_arguments.closure_argument_range(
            mapping, elastic_support
        )
        elastic_argument = closure_arguments.smooth_log_energy_argument(
            "En_per-Ne_per",
            minimum_eV=elastic_minimum,
            maximum_eV=elastic_maximum,
        )
        lines.append(
            f"    model.component({component})"
            '.variable("swClosureVars").set("sw_logeps_el", '
            + java_string(elastic_argument)
            + ");"
        )
    return lines


def render_transport_functions(context: ApplyJavaContext) -> list[str]:
    mapping = context.mapping
    lines: list[str] = []
    for tag, _, column, _ in _active_transport_functions(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        x_values, y_values = _lookup(context.transport_rows, "mean_energy_eV", column)
        if any(value <= 0.0 for value in y_values):
            raise GecCcpWorkflowError(
                f"differentiable log closure requires positive values in {column}"
            )
        function_tag = f"sw_log_{tag.removeprefix('sw_')}"
        lines.extend(
            java_support.inline_interpolation_lines(
                function_tag,
                [math.log(value) for value in x_values],
                [math.log(value) for value in y_values],
            )
        )
    return lines


def render_elastic_interpolation(context: ApplyJavaContext) -> list[str]:
    if not context.uses_external_elastic_loss:
        return []
    mapping = context.mapping
    elastic_rows = _read_csv(mapping.bundle.path / GEC_ELASTIC_ENERGY_LOSS_TABLE)
    x_values, y_values = _lookup(
        elastic_rows,
        "mean_energy_eV",
        GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    )
    if any(value <= 0.0 for value in y_values):
        raise GecCcpWorkflowError(
            "external elastic energy-loss log closure requires positive coefficients"
        )
    return java_support.inline_interpolation_lines(
        "sw_logKel",
        [math.log(value) for value in x_values],
        [math.log(value) for value in y_values],
    )


def render_transport_bindings(context: ApplyJavaContext) -> list[str]:
    # Bind requested reduced properties directly. COMSOL owns the
    # neutral-density conversion through ElectronProperties.ReducedProps.
    # Mobility-only deterministic closures bind only muN. Full two-term
    # transport is isotropic after its L=T audit; MC retains L/T anisotropy.
    lines: list[str] = []
    mapping = context.mapping
    for property_name, _, tensor_values, _ in _transport_property_tensors(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        tensor = ", ".join(java_string(value) for value in tensor_values)
        lines.append(
            f"    {context.target}.set({java_string(property_name)}, "
            f"new String[]{{{tensor}}});"
        )
    return lines


def render_reaction_bindings(context: ApplyJavaContext) -> list[str]:
    mapping = context.mapping
    lines: list[str] = []
    for reaction in (
        reaction
        for reaction in mapping.reactions
        if context.uses_rates
        and _reaction_uses_preintegrated_rate(mapping.closure, reaction)
    ):
        selected, rate_metadata = _active_rate_table(
            mapping,
            context.rate_rows,
            process_type=reaction.process_type,
        )
        x_values, y_values = _lookup(
            selected, "mean_energy_eV", "rate_coefficient_m3_s"
        )
        reaction_target = (
            f"model.component({java_string(mapping.model.component)})"
            f".physics({java_string(mapping.model.physics)})"
            f".feature({java_string(reaction.feature)})"
        )
        if _uses_direct_mc_rates(mapping):
            function_tag = f"sw_knorm_{reaction.process_type}_e"
            rate_scale = float(rate_metadata["normalization_rate_m3_s"])
            lines.extend(
                java_support.inline_interpolation_lines(
                    function_tag,
                    [math.log(value) for value in x_values],
                    [value / rate_scale for value in y_values],
                )
            )
            particle_rate_value = (
                f"{rate_scale:.17e}*{function_tag}("
                f"sw_logeps_rate_{reaction.process_type})"
            )
        else:
            function_tag = f"sw_logk_{reaction.process_type}_e"
            lines.extend(
                java_support.inline_interpolation_lines(
                    function_tag,
                    [math.log(value) for value in x_values],
                    [math.log(value) for value in y_values],
                )
            )
            rate_argument = (
                f"sw_logeps_rate_{reaction.process_type}"
                if mapping.bundle.expected_source == "monte_carlo"
                else "sw_logeps"
            )
            particle_rate_value = f"exp({function_tag}({rate_argument}))"
        particle_rate = f"{particle_rate_value}*1[m^3/s]"
        lines.extend(
            [
                f"    // {reaction.name}: preserve the Electron Impact Reaction",
                "    // stoichiometry and energy loss, while replacing its rate.",
                "    // This active table is the sole owner of the reaction",
                "    // rate; its Swarm/EEDF provenance is recorded in the plan.",
                f'    {reaction_target}.set("SpecifyReactionUsing", "RateConstant");',
                f'    {reaction_target}.set("RateConstantForm", "UseRate");',
                f'    {reaction_target}.set("UseTownsendTwoTermsBoltzmann", false);',
                "    // Swarm k is per target particle [m^3/s]; COMSOL kf is",
                "    // molar [m^3/(mol*s)], hence the explicit Avogadro factor.",
                f'    {reaction_target}.set("kf", '
                + java_string(f"{AVOGADRO_PER_MOL}*{particle_rate}")
                + ");",
            ]
        )
    return lines


def render_external_elastic_owner(context: ApplyJavaContext) -> list[str]:
    if not context.uses_external_elastic_loss:
        return []
    mapping = context.mapping
    component = java_string(mapping.model.component)
    physics_tag = java_string(mapping.model.physics)
    physics_target = f"model.component({component}).physics({physics_tag})"
    elastic_argument = (
        "sw_logeps_el"
        if mapping.bundle.expected_source == "monte_carlo"
        else "sw_logeps"
    )
    return [
        "    // The solver-native elastic energy moment is the sole",
        "    // elastic owner in the electron energy equation. Disable",
        "    // only eir1; inelastic particle and energy sources remain",
        "    // owned by their original Electron Impact Reaction nodes.",
        f'    {physics_target}.feature("eir1").active(false);',
        f'    try {{ {physics_target}.feature().remove("swElLoss"); }}',
        "    catch (Exception missing) {}",
        f'    {physics_target}.create("swElLoss", "GeneralPowerDeposition", 2);',
        f'    {physics_target}.feature("swElLoss").selection().set(new int[]{{1}});',
        f'    {physics_target}.feature("swElLoss").set("Qgen", '
        + java_string(
            f"-ptp.ne*ptp.n_wAr*exp(sw_logKel({elastic_argument}))*1[eV*m^3/s]"
        )
        + ");",
    ]


def render_function_eedf(context: ApplyJavaContext) -> list[str]:
    if not context.uses_function_eedf:
        return []
    mapping = context.mapping
    spec = mapping.closure.function_eedf
    if spec is None:  # Parser validation should make this unreachable.
        raise GecCcpWorkflowError("missing Function-EEDF specification")
    component = java_string(mapping.model.component)
    physics = (
        f"model.component({component}).physics({java_string(mapping.model.physics)})"
    )
    function_tag = java_string(spec.function_tag)
    function_file = java_string(java_path(mapping.bundle.path / spec.table))
    function_name_table = "new String[][]{new String[]{" + function_tag + ', "1"}}'
    common_function_lines = [
        "    // COMSOL accepts an arbitrary plasma EEDF only from a",
        "    // genuine component-local two-argument Interpolation.",
        "    // The imported physical f0(E,<E>) is linear (C0).",
        f"    try {{ model.component({component}).func().remove("
        f"{function_tag}); }} catch (Exception missing) {{}}",
        f"    model.component({component}).func().create("
        f'{function_tag}, "Interpolation");',
        f'    model.component({component}).func({function_tag}).set("source", "file");',
    ]
    lines = _render_spreadsheet_function_eedf(
        component=component,
        physics=physics,
        function_tag=function_tag,
        function_file=function_file,
        function_name_table=function_name_table,
        common_function_lines=common_function_lines,
    )
    for reaction in (
        reaction
        for reaction in mapping.reactions
        if not _reaction_uses_preintegrated_rate(mapping.closure, reaction)
    ):
        reaction_target = f"{physics}.feature({java_string(reaction.feature)})"
        lines.extend(
            [
                f'    {reaction_target}.set("SpecifyReactionUsing", '
                '"UseCrossSectionData");',
                f'    {reaction_target}.set("eedf", "FromPhysicsInterfaceProperty");',
                f'    {reaction_target}.set("UseTownsendTwoTermsBoltzmann", false);',
            ]
        )
    return lines


def _render_spreadsheet_function_eedf(
    *,
    component: str,
    physics: str,
    function_tag: str,
    function_file: str,
    function_name_table: str,
    common_function_lines: list[str],
) -> list[str]:
    return [
        *common_function_lines,
        f"    model.component({component}).func({function_tag})"
        f'.set("filename", {function_file});',
        f'    model.component({component}).func({function_tag}).set("nargs", 2);',
        f"    model.component({component}).func({function_tag})"
        '.set("struct", "spreadsheet");',
        f"    model.component({component}).func({function_tag})"
        '.set("scaledata", "auto");',
        f"    model.component({component}).func({function_tag})"
        f'.set("funcnametable", {function_name_table});',
        f"    model.component({component}).func({function_tag})"
        '.set("argunit", "eV,eV");',
        f'    model.component({component}).func({function_tag}).set("fununit", "1");',
        f"    model.component({component}).func({function_tag})"
        '.set("interp", "linear");',
        f"    model.component({component}).func({function_tag})"
        '.set("extrap", "const");',
        f"    model.component({component}).func({function_tag}).importData();",
        f'    {physics}.prop("EEDFSettings").set("eedf", {function_tag});',
        "    // Cross-section reactions use this active external EEDF.",
        "    // Preintegrated reactions retain their positive smooth kf.",
    ]


def render_apply_footer(context: ApplyJavaContext) -> list[str]:
    return [
        "    model.save("
        + java_string(java_path(context.mapping.model.output_mph))
        + ");",
        '    ModelUtil.remove("swarmGecApply");',
        "  }",
        "",
        "  private static void replaceInterpolation(Model model, String tag) {",
        "    try { model.func().remove(tag); } catch (Exception ignored) {}",
        '    model.func().create(tag, "Interpolation");',
        "  }",
        "",
        "}",
        "",
    ]
