"""Native-case grouping and independent-replica data for standard figures."""

from copy import deepcopy
from dataclasses import dataclass

import numpy as np
import pandas as pd

from electron_swarm.core.identity import content_digest
from electron_swarm.core.statistics import (
    summarize_replicates,
    conservative_rebin_probability_mass,
)
from electron_swarm.grids.energy import reported_cell_edges

FIGURE_DATA_KINDS = {
    "figure_id": "string",
    "panel": "string",
    "curve": "string",
    "x": "number",
    "x_lower": "number",
    "x_upper": "number",
    "y": "number",
    "ci95_low": "number",
    "ci95_high": "number",
    "replicates": "integer",
    "estimate_status": "string",
    "numerical_status": "string",
}


def optional(value):
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


@dataclass(frozen=True)
class FigureCase:
    solver: str
    case_id: str
    field: float
    context: dict
    observables: dict
    metadata: dict
    replicate: int | None

    @property
    def key(self):
        return [self.solver, self.case_id, self.field]

    @property
    def label(self):
        angular = (
            self.metadata.get("angular_model")
            or self.context["angular_scattering"]["model"]
        )
        fidelity = (
            self.metadata.get("angular_scattering_fidelity") or "fidelity unverified"
        )
        lmax = self.metadata.get("lmax")
        return (
            f"{self.solver}"
            + (f" L{lmax}" if lmax is not None else "")
            + f"; {angular}; {fidelity}"
        )


def source_cases(artifact):
    summary = artifact.tables.get("summary_csv")
    if summary is None:
        return []
    cases, keys = [], set()
    for item in artifact.metadata.get("cases", []):
        key = (item["solver"], item["case_id"], item["E_over_N_Td"])
        if key in keys:
            raise ValueError("figures require unique solver/case/E/N identities")
        keys.add(key)
        context = item.get("physical_context")
        if not context or context.get("e_over_n_Td") != key[2]:
            raise ValueError(
                "figures require a matching physical context for every case"
            )
        row = case_rows(summary, key)
        if len(row) != 1:
            raise ValueError("figure case has no unique summary row")
        meta = {
            name[5:]: optional(value)
            for name, value in row.iloc[0].items()
            if name.startswith("meta_")
        }
        meta["solver_method"] = optional(row.iloc[0].get("solver_method"))
        cases.append(
            FigureCase(
                *key,
                deepcopy(context),
                deepcopy(item.get("transport_observables", {})),
                meta,
                item.get("replicate"),
            )
        )
    return cases


def case_rows(frame, key):
    solver, case_id, field = key
    axis = "E_over_N_Td" if "E_over_N_Td" in frame else "E_over_N_rms_Td"
    return frame.loc[
        (frame.solver == solver) & (frame.case_id == case_id) & (frame[axis] == field)
    ]


def physical_groups(cases, *, include_field):
    groups = {}
    for case in cases:
        context = deepcopy(case.context)
        if not include_field:
            context.pop("e_over_n_Td")
        digest = content_digest(context)
        groups.setdefault(digest, (context, []))[1].append(case)
    return [groups[key] for key in sorted(groups)]


def independent(cases):
    return (
        len(cases) > 1
        and cases[0].solver == "monte_carlo"
        and all(case.replicate is not None for case in cases)
        and len({case.replicate for case in cases}) == len(cases)
        and all(
            case.metadata.get("monte_carlo_case_seed") is not None for case in cases
        )
        and len({case.metadata["monte_carlo_case_seed"] for case in cases})
        == len(cases)
    )


def status(cases, *, transport=False):
    states = {case.metadata.get("numerical_status") or "unknown" for case in cases}
    if transport and any(
        case.metadata.get("electron_electron_transport_stale") for case in cases
    ):
        states.add("transport_stale")
    return "; ".join(sorted(states))


def scalar_points(cases, values, *, transport=False, censored=False):
    """Keep unproved replicas separate; never report a false confidence band."""
    points = []
    for field in sorted({case.field for case in cases}):
        selected = [
            case
            for case in cases
            if case.field == field and values.get(tuple(case.key)) is not None
        ]
        if not selected:
            continue
        if independent(selected):
            stats = summarize_replicates(
                [values[tuple(case.key)] for case in selected],
                zero_is_censored=censored,
            )
            points.append(
                dict(
                    x=field,
                    y=stats.mean,
                    ci95_low=stats.ci95_low,
                    ci95_high=stats.ci95_high,
                    replicates=len(selected),
                    estimate_status=stats.estimate_status,
                    numerical_status=status(selected, transport=transport),
                )
            )
        else:
            for case in selected:
                points.append(
                    dict(
                        x=field,
                        y=values[tuple(case.key)],
                        ci95_low=None,
                        ci95_high=None,
                        replicates=1,
                        estimate_status="single_replica"
                        if case.solver == "monte_carlo"
                        else "native_case",
                        numerical_status=status([case], transport=transport),
                    )
                )
    return points


def eedf_curves(cases, frame):
    """Average only independent replicas on their common covered energy domain."""
    groups = [cases] if independent(cases) else [[case] for case in cases]
    curves = []
    for group in groups:
        data = [case_rows(frame, case.key).sort_values("energy_eV") for case in group]
        if any(item.empty for item in data):
            raise ValueError("requested EEDF figure has missing native cells")
        edges = [
            reported_cell_edges(item.energy_eV, item.energy_width_eV) for item in data
        ]
        same = all(
            np.array_equal(edges[0], edge)
            and np.array_equal(data[0].energy_eV, item.energy_eV)
            for edge, item in zip(edges, data)
        )
        outside = []
        if same:
            grid = edges[0]
            centers = np.asarray(data[0].energy_eV, dtype=float)
            values = np.array([item.eedf.to_numpy(dtype=float) for item in data])
            outside = [0.0] * len(data)
        else:
            end = min(edge[-1] for edge in edges)
            grid = np.unique(
                np.concatenate(
                    [np.array([0.0, end])]
                    + [edge[(edge > 0) & (edge < end)] for edge in edges]
                )
            )
            centers = 0.5 * (grid[1:] + grid[:-1])
            values = []
            for item, edge in zip(data, edges):
                original = item.eedf.to_numpy(
                    dtype=float
                ) * item.energy_width_eV.to_numpy(dtype=float)
                mass = conservative_rebin_probability_mass(edge, original, grid)
                outside.append(float(original.sum() - mass.sum()))
                values.append(mass / np.diff(grid))
            values = np.array(values)
        if len(group) > 1:
            stats = [
                summarize_replicates(values[:, index], zero_is_censored=True)
                for index in range(len(centers))
            ]
            density = np.array([s.mean for s in stats])
            low = [s.ci95_low for s in stats]
            high = [s.ci95_high for s in stats]
        else:
            density = values[0]
            low = high = [None] * len(centers)
        curves.append(
            (
                group,
                centers,
                np.diff(grid),
                density,
                low,
                high,
                {
                    "grid_operation": "native"
                    if same
                    else "conservative_common_overlap",
                    "source_mass_outside_common_domain": outside,
                    "normalization": float(np.sum(density * np.diff(grid))),
                    "negative_mass": float(
                        np.sum(np.maximum(-density, 0) * np.diff(grid))
                    ),
                },
            )
        )
    return curves


def energy_display_limit(energy, widths, values):
    """Choose a native cell edge containing 99.9% of the absolute cell mass."""
    edges = reported_cell_edges(energy, widths)
    mass = np.abs(np.asarray(values, dtype=float)) * widths
    total = float(mass.sum())
    index = min(int(np.searchsorted(np.cumsum(mass), 0.999 * total)), len(mass) - 1)
    return float(edges[index + 1])


def context_label(context):
    mixture = "/".join(
        f"{name}:{fraction:g}" for name, fraction, _ in context["species"]
    )
    text = f"{mixture}; Tg={context['gas_temperature_K']:g} K; N={context['gas_number_density_m3']:.4g} m^-3"
    if "e_over_n_Td" in context:
        text += f"; E/N{' RMS' if context['field_type'] != 'dc' else ''}={context['e_over_n_Td']:g} Td"
    if context.get("field_frequency_Hz"):
        text += f"; f={context['field_frequency_Hz']:g} Hz"
    return text


def curve_metadata(cases, label=None):
    return {
        "label": label or cases[0].label,
        "solver": cases[0].solver,
        "source_cases": [case.key for case in cases],
        "source_status": status(cases),
    }


def channel_label(key):
    species, process, process_type, threshold = key
    return f"{species}: {process}\n{process_type}" + (
        f"; threshold {threshold:g} eV" if threshold is not None else ""
    )
