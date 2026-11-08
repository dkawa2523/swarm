"""Publication and comparison choices, independent of numerical model defaults."""

from dataclasses import dataclass, field as dc_field
from pathlib import Path
from typing import Literal
from importlib.util import find_spec

from electron_swarm.chemistry.catalog import PathwayScreen
from .solver_ids import CANONICAL_SOLVER_IDS, SolverId

FigureKind = Literal["eedf", "transport", "reaction_rates", "energy_angle", "rf_phase"]
FIGURE_KINDS = ("eedf", "transport", "reaction_rates", "energy_angle", "rf_phase")
_FIGURE_SOURCES = {
    "eedf": frozenset(CANONICAL_SOLVER_IDS),
    "transport": frozenset(CANONICAL_SOLVER_IDS),
    "reaction_rates": frozenset(CANONICAL_SOLVER_IDS),
    "energy_angle": frozenset({"propagator"}),
    "rf_phase": frozenset({"two_term"}),
}


@dataclass(frozen=True, slots=True)
class RateExportConfig:
    solver: SolverId
    quantities: tuple[str, ...] = ("reaction_rates",)
    pathways: PathwayScreen | None = None

    def __post_init__(self) -> None:
        if self.solver not in CANONICAL_SOLVER_IDS:
            raise ValueError("output.export.solver must be a canonical solver id")
        supported = {"reaction_rates", "reaction_pathways"}
        if (
            not self.quantities
            or any(
                not isinstance(item, str) or item not in supported
                for item in self.quantities
            )
            or len(set(self.quantities)) != len(self.quantities)
        ):
            raise ValueError(
                "output.export.quantities must select unique reaction_rates or reaction_pathways"
            )
        if "reaction_pathways" in self.quantities:
            if not isinstance(self.pathways, PathwayScreen):
                raise ValueError(
                    "reaction_pathways requires output.export.pathways with a catalog"
                )
        elif self.pathways is not None:
            raise ValueError(
                "output.export.pathways requires reaction_pathways in quantities"
            )


@dataclass(slots=True)
class ComparisonConfig:
    enabled: bool = False
    reference_solver: SolverId | None = None
    candidate_solvers: list[SolverId] = dc_field(default_factory=list)
    compare_eedf: bool = True
    required: bool = False


@dataclass(slots=True)
class OutputConfig:
    directory: Path | str = Path("outputs")
    base_name: str = "swarm"
    export: RateExportConfig | None = None
    figures: tuple[FigureKind, ...] = ()

    def __post_init__(self) -> None:
        self.validate()
        self.directory = Path(self.directory)
        self.figures = tuple(self.figures)

    def validate(self) -> None:
        if (
            not isinstance(self.figures, (list, tuple))
            or any(
                not isinstance(kind, str) or kind not in FIGURE_KINDS
                for kind in self.figures
            )
            or len(set(self.figures)) != len(self.figures)
        ):
            raise ValueError(
                f"output.figures must be a unique list from {list(FIGURE_KINDS)}"
            )
        if self.figures and find_spec("matplotlib") is None:
            raise ValueError(
                "figure output requires matplotlib: pip install 'electron-swarm[plot]'"
            )
        if self.export is not None and not isinstance(self.export, RateExportConfig):
            raise ValueError("output.export must be a RateExportConfig")
        if (
            not isinstance(self.directory, (str, Path))
            or not str(self.directory).strip()
            or "\0" in str(self.directory)
        ):
            raise ValueError("output.directory must be a nonempty filesystem path")
        if (
            not isinstance(self.base_name, str)
            or not self.base_name.strip()
            or self.base_name in {".", ".."}
            or any(char in self.base_name for char in '/\\:\0<>"|?*')
            or any(ord(char) < 32 for char in self.base_name)
        ):
            raise ValueError("output.base_name must be a nonempty filename component")

    def export_plan(self) -> dict | None:
        if self.export is None:
            return None
        plan = {
            "solver": self.export.solver,
            "quantities": list(self.export.quantities),
            "independent_variable": "E_over_N_Td",
            "coefficient_unit": "m^3/s",
            "data_role": "coefficient_points",
            "directory": str(self.directory),
            "csv": f"{self.base_name}_reaction_rates.csv",
        }
        if self.export.pathways is not None:
            plan["pathways"] = self.export.pathways.plan()
            plan["pathways_csv"] = f"{self.base_name}_reaction_pathways.csv"
        return plan

    def figure_plan(self, *, runnable, field_type: str) -> dict | None:
        if not self.figures:
            return None
        sources = {}
        for kind in self.figures:
            eligible = set(runnable) & _FIGURE_SOURCES[kind]
            if not eligible:
                raise ValueError(f"requested {kind} figure has no runnable source")
            if kind == "rf_phase" and field_type != "time_dependent":
                raise ValueError("rf_phase figures require a time-dependent field")
            if kind == "reaction_rates" and field_type != "dc":
                raise NotImplementedError(
                    "reaction_rates figures currently require DC E/N; use rf_phase for RF electron response"
                )
            sources[kind] = sorted(eligible)
        return {
            "requested": list(self.figures),
            "formats": ["png", "svg"],
            "sources": sources,
            "directory": str(self.directory),
            "data": f"{self.base_name}_figure_data.csv",
        }
