"""Section renderers for GEC-CCP result-export Java sources."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from swarm_workflow.comsol.java import java_string
import swarm_workflow.comsol.models.gec_ccp.java_support as java_support
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _transport_audit_components,
    _uses_external_elastic_energy_loss,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import GecCcpMapping


@dataclass(frozen=True, slots=True)
class ExportJavaContext:
    mapping: GecCcpMapping
    output_directory: Path
    period_dataset: str
    waveform_dataset: str
    include_external_closure_audit: bool
    external_elastic_loss: bool
    elastic_phase_expression: str
    elastic_period_average_expression: str
    actual_transport_expressions: tuple[str, ...]
    transport_audit_units: tuple[str, ...]
    full_domain_transport_audit: bool
    active_rate_expressions: tuple[str, ...]
    closure_phase_expressions: tuple[str, ...]
    closure_phase_units: tuple[str, ...]


def build_export_java_context(
    mapping: GecCcpMapping,
    *,
    output_directory: Path,
    period_dataset: str,
    waveform_dataset: str,
    include_external_closure_audit: bool,
) -> ExportJavaContext:
    external_elastic_loss = bool(
        include_external_closure_audit
        and _uses_external_elastic_energy_loss(mapping.closure)
    )
    elastic_argument = (
        "sw_logeps_el"
        if mapping.bundle.expected_source == "monte_carlo"
        else "sw_logeps"
    )
    # The converted physical-time dataset must reconstruct the constitutive
    # expression instead of evaluating periodic Ne_per/En_per feature state.
    elastic_phase_expression = (
        f"-ptp.ne*ptp.n_wAr*exp(sw_logKel({elastic_argument}))*1[eV*m^3/s]"
    )
    elastic_period_average_expression = (
        "ptp.xdintop_ptp3((" + elastic_phase_expression + ")/ptp.xdim)"
    )
    transport_audit = (
        _transport_audit_components(
            mapping.closure, source=mapping.bundle.expected_source
        )
        if include_external_closure_audit
        else ()
    )
    actual_transport_expr = tuple(item[2] for item in transport_audit)
    transport_audit_units = tuple(item[4] for item in transport_audit)
    full_domain_transport_audit = bool(
        include_external_closure_audit
        and (
            mapping.results.role == "physical_target"
            or mapping.run.support_policy == "qualified_guard_sensitivity"
        )
    )
    active_rate_expr = (
        ("ptp.kf_2", "ptp.kf_3")
        if external_elastic_loss
        else ("ptp.kf_1", "ptp.kf_2", "ptp.kf_3")
    )
    closure_phase_expr = (
        "ptp.ebar",
        "ptp.Nn",
        *(
            actual_transport_expr
            if include_external_closure_audit
            else ("ptp.muerr", "ptp.Derr", "ptp.muenrr", "ptp.Denrr")
        ),
        *active_rate_expr,
    )
    closure_phase_units = (
        "V",
        "1/m^3",
        *(
            transport_audit_units
            if include_external_closure_audit
            else ("m^2/(V*s)", "m^2/s", "m^2/(V*s)", "m^2/s")
        ),
        *("m^3/(s*mol)" for _ in active_rate_expr),
    )
    return ExportJavaContext(
        mapping=mapping,
        output_directory=output_directory,
        period_dataset=period_dataset,
        waveform_dataset=waveform_dataset,
        include_external_closure_audit=include_external_closure_audit,
        external_elastic_loss=external_elastic_loss,
        elastic_phase_expression=elastic_phase_expression,
        elastic_period_average_expression=elastic_period_average_expression,
        actual_transport_expressions=actual_transport_expr,
        transport_audit_units=transport_audit_units,
        full_domain_transport_audit=full_domain_transport_audit,
        active_rate_expressions=active_rate_expr,
        closure_phase_expressions=closure_phase_expr,
        closure_phase_units=closure_phase_units,
    )


def render_elastic_domain_assertion(mapping: GecCcpMapping) -> list[str]:
    if not _uses_external_elastic_energy_loss(mapping.closure):
        return []
    physics = (
        f"model.component({java_string(mapping.model.component)})"
        f".physics({java_string(mapping.model.physics)})"
    )
    return [
        "    // Fail closed unless elastic loss acts on plasma domain 1 only.",
        f'    int[] swElLossDomains = {physics}.feature("swElLoss")'
        ".selection().entities(2);",
        "    if (swElLossDomains.length != 1 || swElLossDomains[0] != 1) {",
        '      throw new IllegalStateException("swElLoss must select domain 1 only");',
        "    }",
    ]


def render_field_exports(context: ExportJavaContext) -> list[str]:
    mapping = context.mapping
    output = context.output_directory
    period_expr = ("ptp.neav", "ptp.Teav", "ptp.Vav", "ptp.Re_av")
    period_units = ("1/m^3", "V", "V", "1/(m^3*s)")
    phase_expr = ("ptp.ne", "ptp.Te", "V")
    phase_units = ("1/m^3", "V", "V")
    lines = [
        f"    model.result().dataset({java_string(mapping.run.axis_dataset)})"
        f'.set("data", {java_string(context.period_dataset)});',
        f"    model.result().dataset({java_string(mapping.run.radial_dataset)})"
        f'.set("data", {java_string(context.period_dataset)});',
    ]
    for tag, dataset, filename in (
        ("swDomainAvg", context.period_dataset, "domain_period_average.csv"),
        ("swAxisAvg", mapping.run.axis_dataset, "axis_period_average.csv"),
        ("swRadialAvg", mapping.run.radial_dataset, "radial_period_average.csv"),
    ):
        lines.extend(
            java_support.data_export_lines(
                tag, dataset, period_expr, period_units, output / filename
            )
        )
    lines.append(
        f"    model.result().dataset({java_string(mapping.run.axis_dataset)})"
        f'.set("data", {java_string(mapping.run.phase_dataset)});'
    )
    lines.extend(
        java_support.data_export_lines(
            "swAxisPhase",
            mapping.run.axis_dataset,
            phase_expr,
            phase_units,
            output / "axis_phase_resolved.csv",
        )
    )
    lines.extend(_render_domain_phase_export(context))
    lines.extend(
        java_support.data_export_lines(
            "swClosurePhase",
            mapping.run.axis_dataset,
            context.closure_phase_expressions,
            context.closure_phase_units,
            output / "closure_phase.csv",
        )
    )
    lines.append(
        f"    model.result().dataset({java_string(mapping.run.radial_dataset)})"
        f'.set("data", {java_string(mapping.run.phase_dataset)});'
    )
    lines.extend(
        java_support.data_export_lines(
            "swClosurePhaseRadial",
            mapping.run.radial_dataset,
            context.closure_phase_expressions,
            context.closure_phase_units,
            output / "closure_phase_radial.csv",
        )
    )
    lines.extend(
        java_support.data_export_lines(
            "swWaveform",
            context.waveform_dataset,
            ("ptp.mct1.V", "ptp.mct1.I"),
            ("V", "A"),
            output / "electrode_waveform.csv",
        )
    )
    return lines


def _render_domain_phase_export(context: ExportJavaContext) -> list[str]:
    mapping = context.mapping
    return java_support.data_export_lines(
        "swDomainPhaseClosure",
        mapping.run.phase_dataset,
        (
            "ptp.ne",
            "ptp.ebar",
            "ptp.Nn",
            "ptp.Er",
            "ptp.Ez",
            *(
                context.actual_transport_expressions
                if context.full_domain_transport_audit
                else ()
            ),
            "ptp.wAr_1p",
            *(
                ("ptp.n_wAr", context.elastic_phase_expression)
                if context.external_elastic_loss
                else ()
            ),
            *context.active_rate_expressions,
        ),
        (
            "1/m^3",
            "V",
            "1/m^3",
            "V/m",
            "V/m",
            *(
                context.transport_audit_units
                if context.full_domain_transport_audit
                else ()
            ),
            "1",
            *(("1/m^3", "W/m^3") if context.external_elastic_loss else ()),
            *("m^3/(s*mol)" for _ in context.active_rate_expressions),
        ),
        context.output_directory / "domain_phase_closure.csv",
    )


def render_conservation_exports(context: ExportJavaContext) -> list[str]:
    lines = java_support.numerical_table_export_lines(
        "swConservationVolume",
        "IntSurface",
        context.period_dataset,
        (
            "ptp.Re_av",
            "N_A_const*ptp.R_wAr_1p_av",
            "ptp.Pcap_av",
            "ptp.xdintop_ptp3(ptp.Sen/ptp.xdim)",
            (
                context.elastic_period_average_expression
                if context.external_elastic_loss
                else "0[W/m^3]"
            ),
            "1",
        ),
        ("1/s", "1/s", "W", "W", "W", "m^3"),
        (
            "electron volume source",
            "argon ion volume source",
            "electron absorbed power",
            "electron collisional signed power",
            "external elastic signed power",
            "axisymmetric plasma volume",
        ),
        (1,),
        "intvolume",
        context.output_directory / "conservation_volume.csv",
    )
    lines.extend(_render_wall_conservation_export(context))
    lines.extend(
        java_support.numerical_table_export_lines(
            "swConservationTerminalPower",
            "EvalGlobal",
            context.period_dataset,
            (
                "ptp.mct1.PowerT",
                "P0",
                "ptp.mct1.PowerT-P0",
                "ptp.mct1.Ipa",
            ),
            ("W", "W", "W", "A"),
            (
                "terminal power",
                "prescribed power",
                "terminal power constraint residual",
                "period-averaged terminal current",
            ),
            None,
            None,
            context.output_directory / "conservation_terminal_power.csv",
        )
    )
    return lines


def _render_wall_conservation_export(context: ExportJavaContext) -> list[str]:
    return java_support.numerical_table_export_lines(
        "swConservationWall",
        "IntLine",
        context.period_dataset,
        (
            "ptp.xdintop_ptp3((ptp.ne_bnd+ptp.tne_bnd)/ptp.xdim)",
            "N_A_const*ptp.Rsurf_1_wAr_1p",
            "e_const*ptp.xdintop_ptp3((ptp.en_bnd+ptp.ten_bnd)/ptp.xdim)",
            (
                "e_const*ptp.xdintop_ptp3((5/6)*(1-ptp.re)*ptp.en*"
                "ptp.nue_th/((1+ptp.re)*ptp.xdim))"
            ),
            "e_const*ptp.xdintop_ptp3(ptp.seflux/ptp.xdim)",
            "e_const*ptp.xdintop_ptp3(ptp.ten_bnd/ptp.xdim)",
            "nr*down(ptp.gflux_ne_avr)+nz*down(ptp.gflux_ne_avz)",
            (
                "e_const*(nr*ptp.xdintop_ptp3("
                "down(ptp.gflux_enr)/ptp.xdim)+"
                "nz*ptp.xdintop_ptp3("
                "down(ptp.gflux_enz)/ptp.xdim))"
            ),
        ),
        ("1/s", "1/s", "W", "W", "W", "W", "1/s", "W"),
        (
            "signed electron particle weak boundary term",
            "argon ion wall weak term",
            "signed electron energy weak boundary term",
            "random-motion electron energy wall loss",
            "secondary-electron energy wall input",
            "thermionic-electron energy wall input",
            "strong-trace electron particle flux diagnostic",
            "strong-trace electron energy flux diagnostic",
        ),
        tuple(range(2, 13)),
        "intsurface",
        context.output_directory / "conservation_wall.csv",
    )


def render_export_footer(class_name: str) -> list[str]:
    return [
        "    ModelUtil.remove(" + java_string(class_name) + ");",
        "  }",
        "}",
        "",
    ]
