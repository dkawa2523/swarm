"""Physics-feature section parser."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, cast

from electron_swarm.core.config import (
    AngularScatteringConfig,
    ElectronElectronConfig,
    EnergyGridPolicyConfig,
    FieldConfig,
    FiniteKConfig,
    IonizationConfig,
    MagneticFieldConfig,
    MomentTableConfig,
    PhysicsConfig,
    TimeDependentFieldConfig,
)
from electron_swarm.core.config_sections.common import (
    as_mapping_section,
    as_path,
    bool_field,
    float_field,
    integer_field,
    reject_unknown_fields,
    strict_float,
    string_value,
    validate_literal,
)


def _parse_finite_k(raw: dict[str, Any]) -> FiniteKConfig:
    enabled = bool_field(raw, "enabled", False, "physics.finite_k.enabled")
    if raw.get("k_m_inv") is None:
        return FiniteKConfig(enabled=enabled, k_m_inv=None)
    return FiniteKConfig(
        enabled=enabled,
        k_m_inv=strict_float(raw["k_m_inv"], "physics.finite_k.k_m_inv"),
    )


def _parse_angular_scattering(
    raw: dict[str, Any], base: Path
) -> AngularScatteringConfig:
    angular_name = string_value(
        raw.get("model", "isotropic"),
        "physics.angular_scattering.model",
    )
    if angular_name == "dcs_table" or "dcs_table" in raw:
        raise ValueError(
            "physics.angular_scattering.dcs_table input is not implemented; "
            "use model=moment_table"
        )
    reject_unknown_fields(
        raw,
        {"model", "higher_moment_closure", "moment_table"},
        "physics.angular_scattering",
    )

    closure_defaults = {
        "isotropic": "zero",
        "momentum_power": "power",
        "maxent_p1": "maxent",
        "moment_table": "table",
    }
    if angular_name not in closure_defaults:
        raise ValueError(
            "physics.angular_scattering.model must be isotropic, momentum_power, "
            "maxent_p1, or moment_table"
        )
    moment_table_raw = as_mapping_section(
        raw, "moment_table", "physics.angular_scattering.moment_table"
    )
    higher_moment_closure = cast(
        Literal["zero", "power", "maxent", "table"],
        validate_literal(
            string_value(
                raw.get("higher_moment_closure", closure_defaults[angular_name]),
                "physics.angular_scattering.higher_moment_closure",
            ),
            {"zero", "power", "maxent", "table"},
            "physics.angular_scattering.higher_moment_closure",
        ),
    )
    moment_table: MomentTableConfig | None = None
    if angular_name == "moment_table" or moment_table_raw:
        reject_unknown_fields(
            moment_table_raw,
            {"path", "format", "provenance", "extrapolation"},
            "physics.angular_scattering.moment_table",
        )
        if "path" not in moment_table_raw:
            raise ValueError("physics.angular_scattering.moment_table.path is required")
        moment_table = MomentTableConfig(
            path=as_path(moment_table_raw["path"], base)
            or Path(moment_table_raw["path"]),
            format=cast(
                Literal["normalized_legendre_moments"],
                validate_literal(
                    string_value(
                        moment_table_raw.get(
                            "format", "normalized_legendre_moments"
                        ),
                        "physics.angular_scattering.moment_table.format",
                    ),
                    {"normalized_legendre_moments"},
                    "physics.angular_scattering.moment_table.format",
                ),
            ),
            provenance=cast(
                Literal["dcs_derived", "model_derived", "unknown"],
                validate_literal(
                    string_value(
                        moment_table_raw.get("provenance", "unknown"),
                        "physics.angular_scattering.moment_table.provenance",
                    ),
                    {"dcs_derived", "model_derived", "unknown"},
                    "physics.angular_scattering.moment_table.provenance",
                ),
            ),
            extrapolation=cast(
                Literal["error"],
                validate_literal(
                    string_value(
                        moment_table_raw.get("extrapolation", "error"),
                        "physics.angular_scattering.moment_table.extrapolation",
                    ),
                    {"error"},
                    "physics.angular_scattering.moment_table.extrapolation",
                ),
            ),
        )

    return AngularScatteringConfig(
        model=cast(
            Literal["isotropic", "momentum_power", "maxent_p1", "moment_table"],
            angular_name,
        ),
        higher_moment_closure=higher_moment_closure,
        moment_table=moment_table,
    )


def _parse_electron_electron(raw: dict[str, Any]) -> ElectronElectronConfig:
    reject_unknown_fields(
        raw,
        {
            "enabled",
            "model",
            "strength_model",
            "relaxation_fraction",
            "conserve_mean_energy",
            "fallback_temperature_eV",
        },
        "physics.electron_electron",
    )
    enabled = bool_field(raw, "enabled", False, "physics.electron_electron.enabled")
    model = validate_literal(
        string_value(
            raw.get("model", "none"),
            "physics.electron_electron.model",
        ),
        {"none", "relaxation_postprocess", "fp_energy"},
        "physics.electron_electron.model",
    )

    relaxation_fraction = strict_float(
        raw.get("relaxation_fraction", 0.05),
        "physics.electron_electron.relaxation_fraction",
    )
    fallback_temperature_eV = strict_float(
        raw.get("fallback_temperature_eV", 2.0),
        "physics.electron_electron.fallback_temperature_eV",
    )
    return ElectronElectronConfig(
        enabled=enabled,
        model=cast(Literal["none", "relaxation_postprocess", "fp_energy"], model),
        strength_model=cast(
            Literal["simple_relaxation", "density_based"],
            validate_literal(
                string_value(
                    raw.get("strength_model", "simple_relaxation"),
                    "physics.electron_electron.strength_model",
                ),
                {"simple_relaxation", "density_based"},
                "physics.electron_electron.strength_model",
            ),
        ),
        relaxation_fraction=relaxation_fraction,
        conserve_mean_energy=bool_field(
            raw,
            "conserve_mean_energy",
            True,
            "physics.electron_electron.conserve_mean_energy",
        ),
        fallback_temperature_eV=fallback_temperature_eV,
    )


def _parse_ionization(raw: dict[str, Any]) -> IonizationConfig:
    reject_unknown_fields(
        raw,
        {"energy_sharing", "secondary_electron_energy_eV"},
        "physics.ionization",
    )
    energy_sharing = cast(
        Literal["equal", "primary_secondary", "loss_only"],
        validate_literal(
            string_value(
                raw.get("energy_sharing", "equal"),
                "physics.ionization.energy_sharing",
            ),
            {"equal", "primary_secondary", "loss_only"},
            "physics.ionization.energy_sharing",
        ),
    )
    secondary_electron_energy_eV = strict_float(
        raw.get("secondary_electron_energy_eV", 0.0),
        "physics.ionization.secondary_electron_energy_eV",
    )
    return IonizationConfig(
        energy_sharing=energy_sharing,
        secondary_electron_energy_eV=secondary_electron_energy_eV,
    )


def _parse_energy_grid_policy(raw: dict[str, Any]) -> EnergyGridPolicyConfig:
    reject_unknown_fields(
        raw,
        {
            "adaptive",
            "threshold_refinement",
            "tail_probability_target",
            "tail_metrics",
            "tail_threshold_eV",
            "tail_rate_warning_fraction",
            "max_eV_limit",
        },
        "physics.energy_grid_policy",
    )
    tail_threshold = (
        strict_float(
            raw["tail_threshold_eV"],
            "physics.energy_grid_policy.tail_threshold_eV",
        )
        if raw.get("tail_threshold_eV") is not None
        else None
    )
    tail_rate_warning_fraction = strict_float(
        raw.get("tail_rate_warning_fraction", 0.05),
        "physics.energy_grid_policy.tail_rate_warning_fraction",
    )
    tail_metrics = raw.get("tail_metrics", True)
    if not isinstance(tail_metrics, bool):
        raise ValueError("physics.energy_grid_policy.tail_metrics must be a boolean")
    return EnergyGridPolicyConfig(
        adaptive=bool_field(
            raw, "adaptive", True, "physics.energy_grid_policy.adaptive"
        ),
        threshold_refinement=bool_field(
            raw,
            "threshold_refinement",
            True,
            "physics.energy_grid_policy.threshold_refinement",
        ),
        tail_probability_target=strict_float(
            raw.get("tail_probability_target", 1.0e-8),
            "physics.energy_grid_policy.tail_probability_target",
        ),
        tail_metrics=tail_metrics,
        tail_threshold_eV=tail_threshold,
        tail_rate_warning_fraction=tail_rate_warning_fraction,
        max_eV_limit=strict_float(
            raw.get("max_eV_limit", 2000.0),
            "physics.energy_grid_policy.max_eV_limit",
        ),
    )


def _parse_field(raw: dict[str, Any]) -> FieldConfig:
    reject_unknown_fields(
        raw,
        {"type", "magnetic_field", "time_dependent"},
        "physics.field",
    )
    magnetic_raw = as_mapping_section(
        raw, "magnetic_field", "physics.field.magnetic_field"
    )
    reject_unknown_fields(
        magnetic_raw,
        {"enabled", "B_T", "angle_EB_deg"},
        "physics.field.magnetic_field",
    )
    magnetic_B_T = float_field(
        magnetic_raw, "B_T", 0.0, "physics.field.magnetic_field.B_T"
    )
    magnetic_angle = float_field(
        magnetic_raw,
        "angle_EB_deg",
        0.0,
        "physics.field.magnetic_field.angle_EB_deg",
    )

    field_type = cast(
        Literal["dc", "rf", "time_dependent"],
        validate_literal(
            string_value(raw.get("type", "dc"), "physics.field.type"),
            {"dc", "rf", "time_dependent"},
            "physics.field.type",
        ),
    )
    time_raw = as_mapping_section(
        raw,
        "time_dependent",
        "physics.field.time_dependent",
    )
    reject_unknown_fields(
        time_raw,
        {
            "waveform",
            "frequency_Hz",
            "amplitude_definition",
            "momentum_response",
            "phase_steps",
            "max_periods",
            "periodic_tolerance",
        },
        "physics.field.time_dependent",
    )
    frequency = (
        None
        if time_raw.get("frequency_Hz") is None
        else float_field(
            time_raw,
            "frequency_Hz",
            0.0,
            "physics.field.time_dependent.frequency_Hz",
        )
    )
    phase_steps = integer_field(
        time_raw,
        "phase_steps",
        48,
        "physics.field.time_dependent.phase_steps",
    )
    max_periods = integer_field(
        time_raw,
        "max_periods",
        2000,
        "physics.field.time_dependent.max_periods",
    )
    periodic_tolerance = strict_float(
        time_raw.get("periodic_tolerance", 1.0e-7),
        "physics.field.time_dependent.periodic_tolerance",
    )
    return FieldConfig(
        type=field_type,
        magnetic_field=MagneticFieldConfig(
            enabled=bool_field(
                magnetic_raw,
                "enabled",
                False,
                "physics.field.magnetic_field.enabled",
            ),
            B_T=magnetic_B_T,
            angle_EB_deg=magnetic_angle,
        ),
        time_dependent=TimeDependentFieldConfig(
            waveform=cast(
                Literal["sinusoidal"],
                validate_literal(
                    string_value(
                        time_raw.get("waveform", "sinusoidal"),
                        "physics.field.time_dependent.waveform",
                    ),
                    {"sinusoidal"},
                    "physics.field.time_dependent.waveform",
                ),
            ),
            frequency_Hz=frequency,
            amplitude_definition=cast(
                Literal["rms"],
                validate_literal(
                    string_value(
                        time_raw.get("amplitude_definition", "rms"),
                        "physics.field.time_dependent.amplitude_definition",
                    ),
                    {"rms"},
                    "physics.field.time_dependent.amplitude_definition",
                ),
            ),
            momentum_response=cast(
                Literal["instantaneous"],
                validate_literal(
                    string_value(
                        time_raw.get("momentum_response", "instantaneous"),
                        "physics.field.time_dependent.momentum_response",
                    ),
                    {"instantaneous"},
                    "physics.field.time_dependent.momentum_response",
                ),
            ),
            phase_steps=phase_steps,
            max_periods=max_periods,
            periodic_tolerance=periodic_tolerance,
        ),
    )


def parse_physics(raw: dict[str, Any], base: Path) -> PhysicsConfig:
    phys_raw = as_mapping_section(raw, "physics", "physics")
    reject_unknown_fields(
        phys_raw,
        {
            "field",
            "angular_scattering",
            "electron_electron",
            "ionization",
            "energy_grid_policy",
            "finite_k",
        },
        "physics",
    )
    return PhysicsConfig(
        field=_parse_field(as_mapping_section(phys_raw, "field", "physics.field")),
        angular_scattering=_parse_angular_scattering(
            as_mapping_section(
                phys_raw,
                "angular_scattering",
                "physics.angular_scattering",
            ),
            base,
        ),
        electron_electron=_parse_electron_electron(
            as_mapping_section(
                phys_raw,
                "electron_electron",
                "physics.electron_electron",
            )
        ),
        ionization=_parse_ionization(
            as_mapping_section(phys_raw, "ionization", "physics.ionization")
        ),
        energy_grid_policy=_parse_energy_grid_policy(
            as_mapping_section(
                phys_raw,
                "energy_grid_policy",
                "physics.energy_grid_policy",
            )
        ),
        finite_k=_parse_finite_k(
            as_mapping_section(phys_raw, "finite_k", "physics.finite_k")
        ),
    )
