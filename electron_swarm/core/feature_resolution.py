"""Value-dependent physics feature resolution for solver descriptors."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from electron_swarm.core.capabilities import (
    ModelFidelity,
    SolverCapabilities,
    SupportLevel,
)

if TYPE_CHECKING:  # pragma: no cover
    from electron_swarm.core.cross_sections import CrossSectionSet


@dataclass(frozen=True, slots=True)
class FeatureTreatment:
    feature: str
    requested: bool
    support: SupportLevel
    treatment: str
    fidelity: ModelFidelity | None = None
    assumption: str | None = None
    reason: str | None = None


def _requested_magnetic_field(config: object) -> bool:
    magnetic_field = config.physics.field.magnetic_field
    return bool(magnetic_field.enabled and magnetic_field.B_T != 0.0)


def _requested_rf_field(config: object) -> bool:
    return str(config.physics.field.type) in {"rf", "time_dependent"}


def _requested_electron_electron(config: object) -> bool:
    return bool(config.physics.electron_electron.enabled)


def _requested_ionization_source(config: object) -> bool:
    return str(config.physics.ionization.energy_sharing) != "equal"


def _requested_finite_k(config: object) -> bool:
    feature = getattr(config.physics, "finite_k", None)
    return bool(feature and feature.enabled)


def _tail_resolution(
    solver: str,
    config: object,
    level: SupportLevel,
) -> tuple[bool, str]:
    if solver == "monte_carlo":
        requested = config.solvers.monte_carlo.tail_max_collisions is not None
        return requested, "configured" if requested else "disabled"
    requested = bool(config.physics.energy_grid_policy.adaptive)
    return requested, level.value if requested else "none"


def _electron_electron_treatment(config: object) -> str:
    ee = config.physics.electron_electron
    if not ee.enabled:
        return "none"
    if getattr(ee, "strength_model", "simple_relaxation") == "density_based":
        raise NotImplementedError(
            "electron_electron strength_model='density_based' is not implemented; "
            "use strength_model='simple_relaxation'"
        )
    if ee.model == "relaxation_postprocess":
        return "relaxation_postprocess"
    if ee.model == "fp_energy":
        return "fp_energy"
    return "unsupported"


def _rf_resolution(
    solver: str,
    config: object,
    level: SupportLevel,
) -> tuple[SupportLevel, str]:
    field_type = str(config.physics.field.type)
    if field_type == "dc":
        return level, "none"
    if solver == "two_term" and field_type == "time_dependent":
        return (
            SupportLevel.APPROXIMATE,
            "time_periodic_f0:instantaneous_f1:sinusoidal_rms",
        )
    if field_type == "rf":
        return SupportLevel.UNSUPPORTED, "high_frequency_stationary_rf_not_implemented"
    return SupportLevel.UNSUPPORTED, "time_dependent_field_not_implemented"


def _ionization_resolution(
    solver: str,
    config: object,
    caps: SolverCapabilities,
) -> tuple[SupportLevel, str]:
    model = str(config.physics.ionization.energy_sharing)
    if solver == "multi_term":
        if model != "equal":
            return SupportLevel.UNSUPPORTED, "unsupported"
        method = str(config.solvers.multi_term.method)
        treatment = (
            "rate_convolution"
            if method in {"pn_closure_direct", "pn_dcs"}
            else "unsupported"
        )
        return caps.ionization_source, treatment
    if solver == "monte_carlo":
        return SupportLevel.APPROXIMATE, model
    if solver == "propagator" and model == "loss_only":
        return SupportLevel.APPROXIMATE, "loss_only_single_daughter"
    return caps.ionization_source, model


def _angular_resolution(
    solver: str,
    config: object,
    caps: SolverCapabilities,
    cross_sections: "CrossSectionSet | None",
) -> tuple[
    SupportLevel,
    str,
    ModelFidelity,
    str,
    str | None,
]:
    angular = config.physics.angular_scattering
    fidelity = _angular_model_fidelity(angular)
    assumption = f"{angular.model}:{angular.higher_moment_closure}"
    if solver in {"monte_carlo", "propagator"}:
        if angular.model not in {"isotropic", "maxent_p1"}:
            reason = (
                f"monte_carlo has no product MC sampler for angular model "
                f"{angular.model!r}"
                if solver == "monte_carlo"
                else f"propagator has no particle transition kernel for "
                f"{angular.model!r}"
            )
            return (
                SupportLevel.UNSUPPORTED,
                (
                    f"same_as_physics:{angular.model}:unsupported_sampler"
                    if solver == "monte_carlo"
                    else f"same_as_physics:{angular.model}:unsupported_kernel"
                ),
                fidelity,
                assumption,
                reason,
            )
        if cross_sections is not None:
            from electron_swarm.core.scattering import (
                resolve_particle_scattering,
            )

            try:
                resolve_particle_scattering(
                    cross_sections,
                    angular_model=angular.model,
                    active_species=(
                        component.species
                        for component in config.conditions.gas_mixture
                        if component.fraction > 0.0
                    ),
                )
            except (NotImplementedError, ValueError) as exc:
                return (
                    SupportLevel.UNSUPPORTED,
                    (
                        "maxent_p1:missing_explicit_total_or_momentum_xs"
                        if angular.model == "maxent_p1"
                        else "isotropic:missing_physical_scattering_xs"
                    ),
                    fidelity,
                    assumption,
                    str(exc),
                )
        suffix = (
            "cell_kernel"
            if solver == "propagator"
            else "sampler_supported"
        )
        return (
            caps.angular_scattering,
            f"same_as_physics:{angular.model}:{suffix}",
            fidelity,
            assumption,
            None,
        )
    if solver == "multi_term" and config.solvers.multi_term.method == "pn_dcs":
        if angular.model != "moment_table":
            return (
                SupportLevel.UNSUPPORTED,
                f"same_as_physics:{angular.model}:unsupported_pn_dcs_source",
                fidelity,
                assumption,
                "multi_term method 'pn_dcs' requires angular model "
                "'moment_table'",
            )
        return (
            caps.angular_scattering,
            f"{angular.model}:pn_dcs_moment_table",
            fidelity,
            assumption,
            None,
        )
    if solver == "multi_term" and angular.model == "moment_table":
        return (
            SupportLevel.UNSUPPORTED,
            "moment_table:unsupported_pn_closure_direct_source",
            fidelity,
            assumption,
            "multi_term method 'pn_closure_direct' uses ordinary integral "
            "cross-section closure; use method 'pn_dcs' for moment_table",
        )
    if solver == "two_term" and angular.model == "moment_table":
        return (
            SupportLevel.UNSUPPORTED,
            "moment_table:unsupported_two_term_source",
            fidelity,
            assumption,
            "two_term uses the momentum-transfer integral and cannot consume "
            "a higher-moment table; use multi_term method 'pn_dcs'",
        )
    operator = (
        "two_term_l1_operator"
        if solver == "two_term"
        else "pn_closure_direct"
    )
    return (
        caps.angular_scattering,
        f"same_as_physics:{angular.model}:{operator}",
        fidelity,
        assumption,
        None,
    )


def _angular_model_fidelity(angular: object) -> ModelFidelity:
    if str(angular.model) != "moment_table":
        return ModelFidelity.INTEGRAL_XS_CLOSURE
    table = getattr(angular, "moment_table", None)
    provenance = str(getattr(table, "provenance", "unknown"))
    return {
        "dcs_derived": ModelFidelity.DCS_DERIVED_MOMENTS,
        "model_derived": ModelFidelity.MODEL_DERIVED_MOMENTS,
    }.get(provenance, ModelFidelity.UNKNOWN_MOMENT_PROVENANCE)


def resolve_solver_features(
    solver: str,
    config: object,
    cross_sections: "CrossSectionSet | None",
    caps: SolverCapabilities,
) -> dict[str, FeatureTreatment]:
    """Resolve one solver's requested physics using values and XS inventory."""

    ee_treatment = _electron_electron_treatment(config)
    ee_requested = _requested_electron_electron(config)
    ee_level = caps.electron_electron
    if ee_requested and ee_treatment == "unsupported":
        ee_level = SupportLevel.UNSUPPORTED

    (
        angular_level,
        angular_treatment,
        angular_fidelity,
        angular_assumption,
        angular_reason,
    ) = _angular_resolution(
        solver,
        config,
        caps,
        cross_sections,
    )
    ionization_level, ionization_treatment = _ionization_resolution(
        solver,
        config,
        caps,
    )
    magnetic_requested = _requested_magnetic_field(config)
    rf_requested = _requested_rf_field(config)
    rf_level, rf_treatment = _rf_resolution(solver, config, caps.rf_field)
    tail_requested, tail_treatment = _tail_resolution(
        solver,
        config,
        caps.tail_refinement,
    )
    return {
        "angular_scattering": FeatureTreatment(
            feature="angular_scattering",
            requested=True,
            support=angular_level,
            treatment=angular_treatment,
            fidelity=angular_fidelity,
            assumption=angular_assumption,
            reason=angular_reason,
        ),
        "ionization_source": FeatureTreatment(
            feature="ionization_source",
            requested=_requested_ionization_source(config),
            support=ionization_level,
            treatment=ionization_treatment,
        ),
        "electron_electron": FeatureTreatment(
            feature="electron_electron",
            requested=ee_requested,
            support=ee_level,
            treatment=ee_treatment if ee_treatment != "none" else "none",
        ),
        "magnetic_field": FeatureTreatment(
            feature="magnetic_field",
            requested=magnetic_requested,
            support=caps.magnetic_field,
            treatment=(
                "boris_lorentz_push"
                if solver == "monte_carlo" and magnetic_requested
                else "none"
            ),
        ),
        "rf_field": FeatureTreatment(
            feature="rf_field",
            requested=rf_requested,
            support=rf_level,
            treatment=rf_treatment if rf_requested else "none",
        ),
        "tail_refinement": FeatureTreatment(
            feature="tail_refinement",
            requested=tail_requested,
            support=caps.tail_refinement,
            treatment=tail_treatment,
        ),
        "finite_k": FeatureTreatment(
            feature="finite_k",
            requested=_requested_finite_k(config),
            support=SupportLevel.UNSUPPORTED,
            treatment="none",
        ),
    }


def two_term_feature_resolver(
    config: object,
    cross_sections: "CrossSectionSet | None",
    caps: SolverCapabilities,
) -> dict[str, FeatureTreatment]:
    return resolve_solver_features("two_term", config, cross_sections, caps)


def multi_term_feature_resolver(
    config: object,
    cross_sections: "CrossSectionSet | None",
    caps: SolverCapabilities,
) -> dict[str, FeatureTreatment]:
    return resolve_solver_features("multi_term", config, cross_sections, caps)


def monte_carlo_feature_resolver(
    config: object,
    cross_sections: "CrossSectionSet | None",
    caps: SolverCapabilities,
) -> dict[str, FeatureTreatment]:
    return resolve_solver_features("monte_carlo", config, cross_sections, caps)


def propagator_feature_resolver(
    config: object,
    cross_sections: "CrossSectionSet | None",
    caps: SolverCapabilities,
) -> dict[str, FeatureTreatment]:
    return resolve_solver_features("propagator", config, cross_sections, caps)
