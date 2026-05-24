"""Build product solve plans from requested solvers and physics features."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from electron_swarm.core.capabilities import SupportLevel, get_solver_capabilities

DEGRADED_SUPPORT = {
    SupportLevel.APPROXIMATE,
    SupportLevel.POSTPROCESS,
    SupportLevel.SURROGATE,
    SupportLevel.DIAGNOSTIC,
}


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

    def to_metadata(self) -> dict[str, Any]:
        caps = get_solver_capabilities(self.solver)
        return {
            "solver": self.solver,
            "requested": self.requested,
            "runnable": self.runnable,
            "skipped": self.skipped,
            "degraded": self.degraded,
            "skip_reason": self.skip_reason,
            "warnings": "; ".join(self.warnings),
            **{f"effective_{key}": value for key, value in self.effective_physics.items()},
            "capability_electron_neutral": caps.electron_neutral.value,
            "capability_angular_scattering": caps.angular_scattering.value,
            "capability_ionization_source": caps.ionization_source.value,
            "capability_electron_electron": caps.electron_electron.value,
            "capability_magnetic_field": caps.magnetic_field.value,
            "capability_tail_refinement": caps.tail_refinement.value,
            "capability_bulk_transport": caps.bulk_transport.value,
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
        return "surrogate_rate_only" if model == "equal" else "unsupported"
    if solver == "monte_carlo":
        return "external_adapter" if model == "equal" else "unsupported"
    return model


def _angular_treatment(solver: str, config: object, level: SupportLevel) -> str:
    angular = config.physics.angular_scattering
    if solver == "monte_carlo":
        mode = config.solvers.monte_carlo.angular_scattering
        if mode == "same_as_physics":
            if angular.model in {"isotropic", "maxent_p1"}:
                return f"same_as_physics:{angular.model}:sampler_supported"
            return f"same_as_physics:{angular.model}:unsupported_sampler"
        return f"external_adapter:{level.value}"
    return f"{angular.model}:{angular.higher_moment_closure}:{level.value}"


def _mc_same_as_physics_sampler_unsupported(config: object) -> bool:
    return (
        config.solvers.monte_carlo.angular_scattering == "same_as_physics"
        and config.physics.angular_scattering.model not in {"isotropic", "maxent_p1"}
    )


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
    allow_unsupported_fallback: bool,
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
        if feature == "magnetic_field" and allow_unsupported_fallback:
            effective[feature] = "ignored_fallback"
            warnings.append(
                f"{solver} ignores requested magnetic_field as explicit fallback"
            )
            return True, True, None
        effective[feature] = "unsupported"
        raise ValueError(
            f"{msg}; unsupported features cannot be silently approximated"
        )
    treatment = effective_value or level.value
    effective[feature] = treatment
    if level in DEGRADED_SUPPORT:
        msg = f"{solver} handles requested {feature} as {treatment}"
        if degraded_policy == "fail":
            raise ValueError(msg)
        warnings.append(msg)
        return True, True, None
    return True, False, None


def build_solve_plan(config: object) -> list[SolverPlanItem]:
    plan: list[SolverPlanItem] = []
    unsupported_policy = config.feature_policy.unsupported
    degraded_policy = config.feature_policy.degraded
    allow_unsupported_fallback = config.feature_policy.allow_unsupported_fallback
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
            "bulk_transport": "none",
            "finite_k": "none",
        }
        degraded = False
        runnable = True
        skip_reason: str | None = None

        if solver == "monte_carlo" and _mc_same_as_physics_sampler_unsupported(config):
            model = config.physics.angular_scattering.model
            msg = (
                "monte_carlo same_as_physics has no product MC sampler for "
                f"angular model {model!r}"
            )
            effective["angular_scattering"] = f"same_as_physics:{model}:unsupported_sampler"
            if unsupported_policy == "fail":
                raise ValueError(msg)
            if unsupported_policy == "skip_solver":
                runnable = False
                skip_reason = msg
            elif config.solvers.monte_carlo.backend == "internal":
                raise ValueError(
                    f"{msg}; internal monte_carlo cannot use an external "
                    "metadata-validation fallback"
                )
            elif not allow_unsupported_fallback:
                raise ValueError(
                    f"{msg}; unsupported features cannot be silently approximated"
                )
            else:
                effective["angular_scattering"] = (
                    f"same_as_physics:{model}:external_metadata_validation_fallback"
                )
                angular_effective = effective["angular_scattering"]
                degraded = True
                warnings.append(
                    "monte_carlo same_as_physics uses external metadata-validation "
                    "fallback; no product sampler is available"
                )

        ee_level = caps.electron_electron
        if ee_requested and ee_treatment == "unsupported":
            ee_level = SupportLevel.UNSUPPORTED
        if runnable:
            for requested, feature, level, effective_value in (
                (
                    ionization_requested,
                    "ionization_source",
                    caps.ionization_source,
                    config.physics.ionization.energy_sharing,
                ),
                (
                    ee_requested,
                    "electron_electron",
                    ee_level,
                    ee_treatment if ee_treatment != "none" else None,
                ),
                (b_requested, "magnetic_field", caps.magnetic_field, None),
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
                    allow_unsupported_fallback=allow_unsupported_fallback,
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
                    allow_unsupported_fallback=allow_unsupported_fallback,
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
            )
        )
    return plan


def solver_plan_metadata(plan: list[SolverPlanItem]) -> list[dict[str, Any]]:
    return [item.to_metadata() for item in plan]
