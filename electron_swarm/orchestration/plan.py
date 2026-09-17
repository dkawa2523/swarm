"""Build product solve plans from requested solvers and physics features."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING, Any

from electron_swarm.core.capabilities import SupportLevel
from electron_swarm.core.feature_resolution import FeatureTreatment
from electron_swarm.core.solver_registry import solver_descriptor

if TYPE_CHECKING:  # pragma: no cover
    from electron_swarm.core.cross_sections import CrossSectionSet


DEGRADED_SUPPORT = {SupportLevel.APPROXIMATE}
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

    def fidelity(self, feature: str, default: str = "") -> str:
        treatment = self.feature_treatments.get(feature)
        if treatment is not None and treatment.fidelity is not None:
            return treatment.fidelity.value
        return default

    def assumption(self, feature: str, default: str = "") -> str:
        treatment = self.feature_treatments.get(feature)
        if treatment is not None and treatment.assumption is not None:
            return treatment.assumption
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
        angular = self.feature_treatments.get("angular_scattering")
        if angular is not None:
            row.update(
                {
                    "angular_scattering_support": angular.support.value,
                    "angular_scattering_fidelity": self.fidelity(
                        "angular_scattering"
                    ),
                    "angular_scattering_assumption": self.assumption(
                        "angular_scattering"
                    ),
                }
            )
        return row


def _unsupported_message(
    solver: str,
    treatment: FeatureTreatment,
) -> str:
    if treatment.reason:
        return treatment.reason
    return (
        f"{solver} does not support requested "
        f"{treatment.feature.replace('_', ' ')}"
    )


def _skip_treatment(feature: str, treatment: str) -> str:
    if feature == "magnetic_field":
        return "skipped"
    if "unsupported" in treatment or "missing_" in treatment:
        return treatment
    return "unsupported"


def _apply_requested_feature_policy(
    *,
    solver: str,
    treatment: FeatureTreatment,
    unsupported_policy: str,
    degraded_policy: str,
    warnings: list[str],
) -> tuple[FeatureTreatment, bool, bool, str | None]:
    if treatment.support == SupportLevel.UNSUPPORTED:
        msg = _unsupported_message(solver, treatment)
        if unsupported_policy == "fail":
            raise ValueError(msg)
        if unsupported_policy == "skip_solver":
            return (
                replace(
                    treatment,
                    treatment=_skip_treatment(
                        treatment.feature,
                        treatment.treatment,
                    ),
                    reason=msg,
                ),
                False,
                False,
                msg,
            )
        raise ValueError(
            f"{msg}; unsupported policy {unsupported_policy!r} is invalid"
        )
    if treatment.support in DEGRADED_SUPPORT:
        msg = (
            f"{solver} handles requested {treatment.feature} "
            f"as {treatment.treatment}"
        )
        if degraded_policy == "fail":
            raise ValueError(msg)
        warnings.append(msg)
        return treatment, True, True, None
    return treatment, True, False, None


def build_solve_plan(
    config: object,
    cross_sections: "CrossSectionSet | None" = None,
) -> list[SolverPlanItem]:
    """Resolve the complete finite solver plan before executing any solver."""

    validate_unique_solver_ids(config)

    plan: list[SolverPlanItem] = []
    unsupported_policy = config.feature_policy.unsupported
    degraded_policy = config.feature_policy.degraded

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

        descriptor = solver_descriptor(solver)
        treatments = descriptor.feature_resolver(
            config,
            cross_sections,
            descriptor.capabilities,
        )
        warnings: list[str] = []
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


def validate_unique_solver_ids(config: object) -> None:
    """Reject duplicate canonical solver requests before orchestration splits."""

    solver_ids = [item.id for item in config.run.solvers]
    duplicate_solver_ids = sorted(
        solver_id
        for solver_id in set(solver_ids)
        if solver_ids.count(solver_id) > 1
    )
    if duplicate_solver_ids:
        raise ValueError(
            "run.solvers contains duplicate solver ids: "
            f"{duplicate_solver_ids}; each canonical solver id may appear only once"
        )


def solver_plan_metadata(plan: list[SolverPlanItem]) -> list[dict[str, Any]]:
    return [item.to_metadata() for item in plan]
