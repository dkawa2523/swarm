"""Build product solve plans from requested solvers and physics features."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from electron_swarm.core.capabilities import SupportLevel, get_solver_capabilities

DEGRADED_SUPPORT = {SupportLevel.APPROXIMATE}
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
    effective_physics: dict[str, str] = field(default_factory=dict)
    feature_treatments: dict[str, FeatureTreatment] = field(default_factory=dict)

    def treatment(self, feature: str, default: str = "none") -> str:
        treatment = self.feature_treatments.get(feature)
        if treatment is not None:
            return treatment.treatment
        return self.effective_physics.get(feature, default)

    def to_metadata(self) -> dict[str, Any]:
        treatments = self.feature_treatments or {
            key: FeatureTreatment(
                feature=key,
                requested=True,
                support="unknown",
                treatment=value,
            )
            for key, value in self.effective_physics.items()
        }
        return {
            "solver": self.solver,
            "requested": self.requested,
            "runnable": self.runnable,
            "skipped": self.skipped,
            "degraded": self.degraded,
            "skip_reason": self.skip_reason,
            "warnings": "; ".join(self.warnings),
            **{
                f"effective_{key}": treatment.treatment
                for key, treatment in treatments.items()
            },
        }


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


def _handle_requested_feature(
    *,
    solver: str,
    feature: str,
    level: SupportLevel,
    unsupported_policy: str,
    degraded_policy: str,
    effective: dict[str, str],
    warnings: list[str],
    effective_value: str | None = None,
) -> tuple[bool, bool, str | None]:
    if level == SupportLevel.UNSUPPORTED:
        msg = f"{solver} does not support requested {feature.replace('_', ' ')}"
        if unsupported_policy == "fail":
            effective[feature] = "unsupported"
            raise ValueError(msg)
        if unsupported_policy == "skip_solver":
            effective[feature] = "skipped" if feature == "magnetic_field" else "unsupported"
            return False, False, msg
        effective[feature] = "unsupported"
        raise ValueError(f"{msg}; unsupported policy {unsupported_policy!r} is invalid")
    treatment = effective_value or level.value
    effective[feature] = treatment
    if level in DEGRADED_SUPPORT:
        msg = f"{solver} handles requested {feature} as {treatment}"
        if degraded_policy == "fail":
            raise ValueError(msg)
        warnings.append(msg)
        return True, True, None
    return True, False, None


def _feature_treatments(
    *,
    effective: dict[str, str],
    requested: dict[str, bool],
    support: dict[str, str],
    skip_reason: str | None,
) -> dict[str, FeatureTreatment]:
    treatments: dict[str, FeatureTreatment] = {}
    for feature in FEATURE_ORDER:
        value = effective.get(feature, "none")
        reason = None
        if requested.get(feature, False) and (
            value in {"unsupported", "skipped"} or "unsupported" in value
        ):
            reason = skip_reason
        treatments[feature] = FeatureTreatment(
            feature=feature,
            requested=bool(requested.get(feature, False)),
            support=support.get(feature, "unknown"),
            treatment=value,
            reason=reason,
        )
    return treatments


def build_solve_plan(config: object) -> list[SolverPlanItem]:
    plan: list[SolverPlanItem] = []
    unsupported_policy = config.feature_policy.unsupported
    degraded_policy = config.feature_policy.degraded
    b_requested = _requested_magnetic_field(config)
    rf_requested = _requested_rf_field(config)
    ee_requested = _requested_electron_electron(config)
    ee_treatment = _electron_electron_treatment(config)
    ionization_requested = _requested_ionization_source(config)
    finite_k_requested = _requested_finite_k(config)

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
                    effective_physics={},
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
        requested_features = {
            "angular_scattering": True,
            "ionization_source": ionization_requested,
            "electron_electron": ee_requested,
            "magnetic_field": b_requested,
            "rf_field": rf_requested,
            "tail_refinement": bool(config.physics.energy_grid_policy.adaptive),
            "finite_k": finite_k_requested,
        }
        support_by_feature = {
            "angular_scattering": caps.angular_scattering.value,
            "ionization_source": ionization_level.value,
            "electron_electron": caps.electron_electron.value,
            "magnetic_field": caps.magnetic_field.value,
            "rf_field": SupportLevel.UNSUPPORTED.value,
            "tail_refinement": caps.tail_refinement.value,
            "finite_k": SupportLevel.UNSUPPORTED.value,
        }
        angular_effective = _angular_treatment(
            solver,
            config,
            caps.angular_scattering,
        )
        effective: dict[str, str] = {
            "angular_scattering": angular_effective,
            "electron_electron": "none",
            "magnetic_field": "none",
            "rf_field": "none",
            "tail_refinement": caps.tail_refinement.value
            if config.physics.energy_grid_policy.adaptive
                else "none",
            "ionization_source": _default_ionization_treatment(solver, config),
            "finite_k": "none",
        }
        degraded = False
        runnable = True
        skip_reason: str | None = None

        if solver == "monte_carlo" and _mc_same_as_physics_sampler_unsupported(config):
            model = config.physics.angular_scattering.model
            msg = (
                "monte_carlo has no product MC sampler for "
                f"angular model {model!r}"
            )
            effective["angular_scattering"] = f"same_as_physics:{model}:unsupported_sampler"
            if unsupported_policy == "fail":
                raise ValueError(msg)
            if unsupported_policy == "skip_solver":
                runnable = False
                skip_reason = msg
            else:
                raise ValueError(msg)

        ee_level = caps.electron_electron
        if ee_requested and ee_treatment == "unsupported":
            ee_level = SupportLevel.UNSUPPORTED
            support_by_feature["electron_electron"] = SupportLevel.UNSUPPORTED.value
        if runnable:
            for requested, feature, level, effective_value in (
                (
                    ionization_requested,
                    "ionization_source",
                    ionization_level,
                    config.physics.ionization.energy_sharing,
                ),
                (
                    ee_requested,
                    "electron_electron",
                    ee_level,
                    ee_treatment if ee_treatment != "none" else None,
                ),
                (
                    b_requested,
                    "magnetic_field",
                    caps.magnetic_field,
                    _magnetic_treatment(solver, config),
                ),
                (rf_requested, "rf_field", SupportLevel.UNSUPPORTED, None),
                (finite_k_requested, "finite_k", SupportLevel.UNSUPPORTED, None),
            ):
                if not requested:
                    continue
                runnable, feature_degraded, skip_reason = _handle_requested_feature(
                    solver=solver,
                    feature=feature,
                    level=level,
                    unsupported_policy=unsupported_policy,
                    degraded_policy=degraded_policy,
                    effective=effective,
                    warnings=warnings,
                    effective_value=effective_value,
                )
                degraded = degraded or feature_degraded
                if not runnable:
                    break

        if runnable:
            for requested, feature, level, effective_value in (
                (
                    True,
                    "angular_scattering",
                    caps.angular_scattering,
                    angular_effective,
                ),
                (
                    bool(config.physics.energy_grid_policy.adaptive),
                    "tail_refinement",
                    caps.tail_refinement,
                    caps.tail_refinement.value,
                ),
            ):
                if not requested:
                    effective[feature] = "none" if feature == "tail_refinement" else effective[feature]
                    continue
                runnable, feature_degraded, skip_reason = _handle_requested_feature(
                    solver=solver,
                    feature=feature,
                    level=level,
                    unsupported_policy=unsupported_policy,
                    degraded_policy=degraded_policy,
                    effective=effective,
                    warnings=warnings,
                    effective_value=effective_value,
                )
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
                effective_physics=effective,
                feature_treatments=_feature_treatments(
                    effective=effective,
                    requested=requested_features,
                    support=support_by_feature,
                    skip_reason=skip_reason,
                ),
            )
        )
    return plan


def solver_plan_metadata(plan: list[SolverPlanItem]) -> list[dict[str, Any]]:
    return [item.to_metadata() for item in plan]
