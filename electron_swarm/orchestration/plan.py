"""Build product solve plans from requested solvers and physics features."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

from electron_swarm.core.capabilities import (
    SolverCapabilities,
    SupportLevel,
    get_solver_capabilities,
)

DEGRADED_SUPPORT = {SupportLevel.APPROXIMATE}
DEGRADED_SUPPORT_VALUES = {level.value for level in DEGRADED_SUPPORT}
POLICY_FEATURE_ORDER = (
    "ionization_source",
    "electron_electron",
    "magnetic_field",
    "rf_field",
    "finite_k",
    "angular_scattering",
    "tail_refinement",
)
FEATURE_ORDER = (
    "angular_scattering",
    "ionization_source",
    "electron_electron",
    "magnetic_field",
    "rf_field",
    "tail_refinement",
    "finite_k",
)


@dataclass(frozen=True, slots=True)
class FeatureTreatment:
    feature: str
    requested: bool
    support: str
    treatment: str
    reason: str | None = None


@dataclass(frozen=True, slots=True)
class SolverPlanItem:
    solver: str
    requested: bool
    runnable: bool
    skipped: bool = False
    degraded: bool = False
    skip_reason: str | None = None
    warnings: tuple[str, ...] = ()
    feature_treatments: dict[str, FeatureTreatment] = field(default_factory=dict)

    def treatment(self, feature: str, default: str = "none") -> str:
        treatment = self.feature_treatments.get(feature)
        if treatment is not None:
            return treatment.treatment
        return default

    def to_metadata(self) -> dict[str, Any]:
        row: dict[str, Any] = {
            "solver": self.solver,
            "requested": self.requested,
            "runnable": self.runnable,
            "skipped": self.skipped,
            "degraded": self.degraded,
            "skip_reason": self.skip_reason,
            "warnings": "; ".join(self.warnings),
        }
        row.update(
            {
                f"effective_{feature}": self.feature_treatments[feature].treatment
                for feature in FEATURE_ORDER
                if feature in self.feature_treatments
            }
        )
        return row


def _requested_magnetic_field(config: object) -> bool:
    return bool(config.physics.field.magnetic_field.enabled)


def _requested_rf_field(config: object) -> bool:
    return str(config.physics.field.type) in {"rf", "time_dependent"}


def _requested_electron_electron(config: object) -> bool:
    return bool(config.physics.electron_electron.enabled)


def _requested_ionization_source(config: object) -> bool:
    ionization = config.physics.ionization
    return str(ionization.energy_sharing) != "equal"


def _default_ionization_treatment(solver: str, config: object) -> str:
    model = str(config.physics.ionization.energy_sharing)
    if solver == "multi_term":
        if model != "equal":
            return "unsupported"
        method = str(config.solvers.multi_term.method)
        return (
            "rate_convolution"
            if method in {"pn_closure_direct", "pn_dcs"}
            else "unsupported"
        )
    if solver == "monte_carlo":
        return model
    return model


def _angular_treatment(solver: str, config: object, level: SupportLevel) -> str:
    angular = config.physics.angular_scattering
    if solver == "monte_carlo":
        if angular.model == "isotropic":
            return "same_as_physics:isotropic:sampler_supported"
        if angular.model == "maxent_p1":
            return "same_as_physics:maxent_p1:requires_total_and_momentum_xs"
        return f"same_as_physics:{angular.model}:unsupported_sampler"
    if solver == "multi_term" and config.solvers.multi_term.method == "pn_dcs":
        return f"{angular.model}:pn_dcs_moment_table"
    return f"{angular.model}:{angular.higher_moment_closure}:{level.value}"


def _mc_same_as_physics_sampler_unsupported(config: object) -> bool:
    return config.physics.angular_scattering.model not in {"isotropic", "maxent_p1"}


def _magnetic_treatment(solver: str, config: object) -> str | None:
    if solver != "monte_carlo":
        return None
    return "boris_lorentz_push"


def _requested_finite_k(config: object) -> bool:
    return bool(getattr(config.physics, "finite_k", None) and config.physics.finite_k.enabled)


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


def _unsupported_message(
    solver: str,
    feature: str,
    treatment: str,
    config: object,
) -> str:
    if solver == "monte_carlo" and feature == "angular_scattering" and (
        treatment.endswith("unsupported_sampler")
    ):
        model = config.physics.angular_scattering.model
        return f"monte_carlo has no product MC sampler for angular model {model!r}"
    return f"{solver} does not support requested {feature.replace('_', ' ')}"


def _skip_treatment(feature: str, treatment: str) -> str:
    if feature == "magnetic_field":
        return "skipped"
    if "unsupported" in treatment:
        return treatment
    return "unsupported"


def _apply_requested_feature_policy(
    *,
    solver: str,
    treatment: FeatureTreatment,
    unsupported_policy: str,
    degraded_policy: str,
    warnings: list[str],
    config: object,
) -> tuple[FeatureTreatment, bool, bool, str | None]:
    if treatment.support == SupportLevel.UNSUPPORTED.value:
        msg = _unsupported_message(
            solver,
            treatment.feature,
            treatment.treatment,
            config,
        )
        if unsupported_policy == "fail":
            raise ValueError(msg)
        if unsupported_policy == "skip_solver":
            return (
                replace(
                    treatment,
                    treatment=_skip_treatment(treatment.feature, treatment.treatment),
                    reason=msg,
                ),
                False,
                False,
                msg,
            )
        raise ValueError(f"{msg}; unsupported policy {unsupported_policy!r} is invalid")
    if treatment.support in DEGRADED_SUPPORT_VALUES:
        msg = (
            f"{solver} handles requested {treatment.feature} "
            f"as {treatment.treatment}"
        )
        if degraded_policy == "fail":
            raise ValueError(msg)
        warnings.append(msg)
        return treatment, True, True, None
    return treatment, True, False, None


def _solver_feature_treatments(
    *,
    solver: str,
    config: object,
    caps: SolverCapabilities,
    ionization_level: SupportLevel,
    ee_treatment: str,
) -> dict[str, FeatureTreatment]:
    ee_requested = _requested_electron_electron(config)
    ee_level = caps.electron_electron
    if ee_requested and ee_treatment == "unsupported":
        ee_level = SupportLevel.UNSUPPORTED

    angular_support = caps.angular_scattering
    angular_treatment = _angular_treatment(solver, config, caps.angular_scattering)
    if solver == "monte_carlo" and _mc_same_as_physics_sampler_unsupported(config):
        angular_support = SupportLevel.UNSUPPORTED

    magnetic_requested = _requested_magnetic_field(config)
    tail_requested = bool(config.physics.energy_grid_policy.adaptive)
    return {
        "angular_scattering": FeatureTreatment(
            feature="angular_scattering",
            requested=True,
            support=angular_support.value,
            treatment=angular_treatment,
        ),
        "ionization_source": FeatureTreatment(
            feature="ionization_source",
            requested=_requested_ionization_source(config),
            support=ionization_level.value,
            treatment=_default_ionization_treatment(solver, config),
        ),
        "electron_electron": FeatureTreatment(
            feature="electron_electron",
            requested=ee_requested,
            support=ee_level.value,
            treatment=ee_treatment if ee_treatment != "none" else "none",
        ),
        "magnetic_field": FeatureTreatment(
            feature="magnetic_field",
            requested=magnetic_requested,
            support=caps.magnetic_field.value,
            treatment=_magnetic_treatment(solver, config)
            if magnetic_requested
            else "none",
        ),
        "rf_field": FeatureTreatment(
            feature="rf_field",
            requested=_requested_rf_field(config),
            support=SupportLevel.UNSUPPORTED.value,
            treatment="none",
        ),
        "tail_refinement": FeatureTreatment(
            feature="tail_refinement",
            requested=tail_requested,
            support=caps.tail_refinement.value,
            treatment=caps.tail_refinement.value if tail_requested else "none",
        ),
        "finite_k": FeatureTreatment(
            feature="finite_k",
            requested=_requested_finite_k(config),
            support=SupportLevel.UNSUPPORTED.value,
            treatment="none",
        ),
    }


def build_solve_plan(config: object) -> list[SolverPlanItem]:
    plan: list[SolverPlanItem] = []
    unsupported_policy = config.feature_policy.unsupported
    degraded_policy = config.feature_policy.degraded
    ee_treatment = _electron_electron_treatment(config)

    for item in config.run.solvers:
        solver = item.id
        if not getattr(item, "enabled", True):
            plan.append(
                SolverPlanItem(
                    solver=solver,
                    requested=True,
                    runnable=False,
                    skipped=True,
                    skip_reason="solver disabled",
                    warnings=("solver disabled",),
                )
            )
            continue
        caps = get_solver_capabilities(solver)
        warnings: list[str] = []
        ionization_level = caps.ionization_source
        if (
            solver == "monte_carlo"
            and str(config.physics.ionization.energy_sharing)
            in {"equal", "primary_secondary", "loss_only"}
        ):
            ionization_level = SupportLevel.APPROXIMATE
        treatments = _solver_feature_treatments(
            solver=solver,
            config=config,
            caps=caps,
            ionization_level=ionization_level,
            ee_treatment=ee_treatment,
        )
        degraded = False
        runnable = True
        skip_reason: str | None = None

        for feature in POLICY_FEATURE_ORDER:
            treatment = treatments[feature]
            if not treatment.requested:
                continue
            treatment, runnable, feature_degraded, skip_reason = (
                _apply_requested_feature_policy(
                    solver=solver,
                    treatment=treatment,
                    unsupported_policy=unsupported_policy,
                    degraded_policy=degraded_policy,
                    warnings=warnings,
                    config=config,
                )
            )
            treatments[feature] = treatment
            degraded = degraded or feature_degraded
            if not runnable:
                break

        plan.append(
            SolverPlanItem(
                solver=solver,
                requested=True,
                runnable=runnable,
                skipped=not runnable,
                degraded=degraded,
                skip_reason=skip_reason,
                warnings=tuple(warnings),
                feature_treatments=treatments,
            )
        )
    return plan


def solver_plan_metadata(plan: list[SolverPlanItem]) -> list[dict[str, Any]]:
    return [item.to_metadata() for item in plan]
