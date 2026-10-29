"""Standard Swarm figures from verified native values, without solver execution."""

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from textwrap import fill

import numpy as np
import pandas as pd

from electron_swarm import __version__
from electron_swarm.core.numerics import eepf_from_eedf
from electron_swarm.grids.energy import reported_cell_edges
from .figure_data import (
    FIGURE_DATA_KINDS,
    source_cases,
    case_rows,
    physical_groups,
    eedf_curves,
    curve_metadata,
    scalar_points,
    energy_display_limit,
    context_label,
    channel_label,
    optional,
    status,
)
from .reaction_rates import rates_from_product
from .artifacts import FIGURE_ARTIFACT_SCHEMA


@dataclass(frozen=True)
class RenderedFigures:
    data: pd.DataFrame
    figures: dict
    paths: dict[str, Path]


class _Figures:
    def __init__(self, artifact):
        self.artifact = artifact
        self.cases = source_cases(artifact)
        self.specs, self.rows = {}, []
        self.counts = defaultdict(int)

    def new(self, kind, context, *, note=""):
        ident = f"{kind}_{self.counts[kind]:04d}"
        self.counts[kind] += 1
        self.specs[ident] = {
            "kind": kind,
            "physical_context": context,
            "panels": {},
            "curves": {},
            "notes": [note] if note else [],
        }
        return ident

    def panel(self, ident, key, label, x_unit, y_unit, **extra):
        self.specs[ident]["panels"][key] = {
            "label": label,
            "x_unit": x_unit,
            "y_unit": y_unit,
            **extra,
        }

    def curve(self, ident, panel, metadata, points):
        curve = f"curve_{len(self.specs[ident]['curves']):04d}"
        self.specs[ident]["curves"][curve] = metadata
        self.rows.extend(
            {
                "figure_id": ident,
                "panel": panel,
                "curve": curve,
                "x_lower": None,
                "x_upper": None,
                **point,
            }
            for point in points
        )


def _require(builder, table, kind):
    frame = builder.artifact.tables.get(table)
    if frame is None or frame.empty:
        raise ValueError(f"{kind} figure requires nonempty {table}")
    return frame


def _distribution_points(group, energy, widths, density, low, high):
    edges = reported_cell_edges(energy, widths)
    return [
        dict(
            x=float(x),
            x_lower=float(edges[i]),
            x_upper=float(edges[i + 1]),
            y=float(density[i]),
            ci95_low=low[i],
            ci95_high=high[i],
            replicates=len(group),
            estimate_status=(
                "censored_all_zero"
                if density[i] == 0 and high[i] is None
                else "independent_replica_estimate"
                if len(group) > 1
                else "single_replica"
            )
            if group[0].solver == "monte_carlo"
            else "native_case",
            numerical_status=status(group),
        )
        for i, x in enumerate(energy)
    ]


def _eedf(builder):
    frame = _require(builder, "eedf_csv", "eedf")
    for context, cases in physical_groups(builder.cases, include_field=True):
        ident = builder.new(
            "eedf",
            context,
            note="Native cell values; no renormalization. Bands are pointwise 95% Student-t intervals when independent MC replicas exist.",
        )
        prepared = []
        for solver in sorted({case.solver for case in cases}):
            prepared.extend(
                eedf_curves([case for case in cases if case.solver == solver], frame)
            )
        limit = max(
            energy_display_limit(energy, widths, values)
            for _, energy, widths, values, *_ in prepared
        )
        builder.panel(
            ident, "eedf", "EEDF", "eV", "1/eV", x_limit=[0, limit], logarithmic=True
        )
        builder.panel(
            ident,
            "eepf",
            "EEPF (source convention)",
            "eV",
            "1/eV^1.5",
            x_limit=[0, limit],
            logarithmic=True,
        )
        builder.panel(
            ident,
            "tail",
            "EEDF full energy support",
            "eV",
            "1/eV",
            logarithmic=True,
            log_x=True,
        )
        omitted = []
        for group, energy, widths, values, low, high, info in prepared:
            edges = reported_cell_edges(energy, widths)
            outside = float(
                np.sum(
                    np.abs(values[edges[:-1] >= limit]) * widths[edges[:-1] >= limit]
                )
            )
            omitted.append(outside)
            meta = {
                **curve_metadata(group),
                **info,
                "outside_display_absolute_mass": outside,
            }
            points = _distribution_points(group, energy, widths, values, low, high)
            builder.curve(ident, "eedf", meta, points)
            builder.curve(ident, "tail", meta, points)
            factor = 1 / np.sqrt(np.maximum(energy, 1e-30))
            builder.curve(
                ident,
                "eepf",
                meta,
                _distribution_points(
                    group,
                    energy,
                    widths,
                    eepf_from_eedf(energy, values),
                    [
                        None if value is None else value * factor[i]
                        for i, value in enumerate(low)
                    ],
                    [
                        None if value is None else value * factor[i]
                        for i, value in enumerate(high)
                    ],
                ),
            )
        builder.specs[ident]["notes"].append(
            f"Zoom: 0-{limit:g} eV; outside-view absolute EEDF mass <= {max(omitted):.4g}. Full support is retained in the third panel and CSV."
        )


_TRANSPORT = (
    ("mean_energy_eV", None, "Mean electron energy", "eV"),
    ("drift_velocity_m_s", "drift_velocity_m_s", "Flux drift", "m/s"),
    (
        "reduced_mobility_m2_V_s_m3",
        "mobility_m2_V_s",
        "Reduced mobility N mu",
        "1/(V m s)",
    ),
    (
        "reduced_diffusion_L_m2_s_m3",
        "diffusion_L_m2_s",
        "Reduced longitudinal diffusion",
        "1/(m s)",
    ),
    (
        "reduced_diffusion_T_m2_s_m3",
        "diffusion_T_m2_s",
        "Reduced transverse diffusion",
        "1/(m s)",
    ),
    (
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "electron_energy_mobility_m2_V_s",
        "Reduced electron-energy mobility",
        "1/(V m s)",
    ),
    (
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "electron_energy_diffusion_L_m2_s",
        "Reduced longitudinal energy diffusion",
        "1/(m s)",
    ),
    (
        "reduced_electron_energy_diffusion_T_m2_s_m3",
        "electron_energy_diffusion_T_m2_s",
        "Reduced transverse energy diffusion",
        "1/(m s)",
    ),
)


def _transport(builder):
    frame = _require(builder, "summary_csv", "transport")
    for context, cases in physical_groups(builder.cases, include_field=False):
        panels, omitted = [], []
        for column, observable, title, unit in _TRANSPORT:
            meanings = defaultdict(list)
            for case in cases:
                row = case_rows(frame, case.key).iloc[0]
                definition = (
                    case.observables.get(observable)
                    if observable
                    else {
                        "quantity": "mean_energy",
                        "method": "source_mean_energy",
                        "approximate": False,
                    }
                )
                if column not in row or optional(row[column]) is None or not definition:
                    omitted.append(f"{case.solver}: {title}")
                    continue
                meanings[definition["quantity"]].append(
                    (case, float(row[column]), definition)
                )
            for quantity, members in meanings.items():
                panels.append((column, title, unit, quantity, members))
        if not panels:
            raise ValueError("transport figure has no defined observable values")
        for start in range(0, len(panels), 6):
            ident = builder.new(
                "transport",
                context,
                note="Different observable quantities occupy different panels. Markers are native E/N points; lines are visual guides.",
            )
            builder.specs[ident]["unavailable"] = sorted(set(omitted))
            if omitted:
                builder.specs[ident]["notes"].append(
                    "Not returned / definition unknown: "
                    + "; ".join(sorted(set(omitted)))
                )
            for index, (column, title, unit, quantity, members) in enumerate(
                panels[start : start + 6]
            ):
                panel = f"panel_{index:04d}"
                axis = "E/N RMS" if context["field_type"] != "dc" else "E/N"
                builder.panel(
                    ident,
                    panel,
                    f"{title}\n{quantity}",
                    "Td",
                    unit,
                    x_label=axis,
                    quantity=quantity,
                )
                for solver in sorted({item[0].solver for item in members}):
                    selected = [item for item in members if item[0].solver == solver]
                    refs = [item[0] for item in selected]
                    values = {tuple(case.key): value for case, value, _ in selected}
                    methods = sorted({item[2]["method"] for item in selected})
                    approximate = any(item[2].get("approximate") for item in selected)
                    label = (
                        refs[0].label
                        + "; "
                        + "/".join(methods)
                        + (" [approximate]" if approximate else "")
                    )
                    if column != "mean_energy_eV" and any(
                        case.metadata.get("electron_electron_transport_stale")
                        for case in refs
                    ):
                        label += " [transport stale]"
                    builder.curve(
                        ident,
                        panel,
                        {**curve_metadata(refs, label), "definition": selected[0][2]},
                        scalar_points(
                            refs, values, transport=column != "mean_energy_eV"
                        ),
                    )


def _rates(builder):
    if (
        not builder.cases
        and "reaction_rates" in builder.artifact.metadata
        and "reaction_rates_csv" in builder.artifact.tables
    ):
        datasets = [
            (
                builder.artifact.metadata["reaction_rates"],
                builder.artifact.tables["reaction_rates_csv"],
            )
        ]
    else:
        _require(builder, "rates_csv", "reaction_rates")
        datasets = []
        for solver in sorted({case.solver for case in builder.cases}):
            data = rates_from_product(
                builder.artifact, solver=solver, source={"kind": "figure_source"}
            )
            datasets.append((data.metadata, data.table))
    groups = {}
    from electron_swarm.core.identity import content_digest

    for metadata, table in datasets:
        for series_id, series in metadata["series"].items():
            context = series["physical_context"]
            if context["field_type"] != "dc":
                raise NotImplementedError(
                    "reaction_rates figures require DC; use rf_phase for periodic electron response"
                )
            digest = content_digest(context)
            group = groups.setdefault(digest, (context, []))[1]
            group.append((metadata, series, table.loc[table.series_id == series_id]))
    for digest in sorted(groups):
        context, sources = groups[digest]
        channels = sorted(
            {
                (row.species, row.process, row.process_type, optional(row.threshold_eV))
                for _, _, frame in sources
                for row in frame.itertuples()
            },
            key=str,
        )
        for start in range(0, len(channels), 6):
            ident = builder.new(
                "reaction_rates",
                context,
                note="Per-target k, not x*k. Binary event rate = n_e*n_target*k. Species source terms require external stoichiometry. No input qualification is inferred from the plot.",
            )
            for index, channel in enumerate(channels[start : start + 6]):
                panel = f"panel_{index:04d}"
                builder.panel(
                    ident,
                    panel,
                    channel_label(channel),
                    "Td",
                    "m^3/s",
                    x_label="E/N",
                    logarithmic=True,
                    channel=list(channel),
                )
                for metadata, series, frame in sources:
                    selected = frame.loc[
                        (frame.species == channel[0])
                        & (frame.process == channel[1])
                        & (frame.process_type == channel[2])
                        & (
                            frame.threshold_eV.isna()
                            if channel[3] is None
                            else frame.threshold_eV == channel[3]
                        )
                    ].sort_values("E_over_N_Td")
                    if selected.empty:
                        continue
                    case_meta = series["cases"][0]["solver_metadata"]
                    angular = (
                        case_meta.get("angular_model")
                        or context["angular_scattering"]["model"]
                    )
                    label = f"{metadata['solver']}; {angular}; {case_meta.get('angular_scattering_fidelity', 'fidelity unverified')}"
                    rows = [
                        dict(
                            x=row.E_over_N_Td,
                            y=row.rate_coefficient_m3_s,
                            ci95_low=optional(row.ci95_low_m3_s),
                            ci95_high=optional(row.ci95_high_m3_s),
                            replicates=int(row.replicates),
                            estimate_status=row.estimate_status,
                            numerical_status=row.numerical_status,
                        )
                        for row in selected.itertuples()
                    ]
                    refs = [
                        [metadata["solver"], case["case_id"], case["E_over_N_Td"]]
                        for case in series["cases"]
                    ]
                    builder.curve(
                        ident,
                        panel,
                        {
                            "label": label,
                            "solver": metadata["solver"],
                            "source_cases": refs,
                            "target_species_fraction": sorted(
                                set(selected.target_species_fraction)
                            ),
                        },
                        rows,
                    )


def _energy_angle(builder):
    frame = _require(builder, "energy_angle_distribution_csv", "energy_angle")
    found = False
    for case in builder.cases:
        selected = case_rows(frame, case.key)
        if selected.empty:
            continue
        found = True
        ident = builder.new(
            "energy_angle",
            case.context,
            note="Native energy-angle cell density; no renormalization. Mu is the solver-native polar cosine; ordinary integral cross sections retain their declared angular closure.",
        )
        builder.panel(
            ident,
            "density",
            case.label,
            "eV",
            "1",
            grid="energy_angle",
            data_table="energy_angle_distribution_csv",
            source_case=case.key,
            numerical_status=status([case]),
        )
    if not found:
        raise ValueError("energy_angle figure has no native energy-angle result")


def _rf_phase(builder):
    frame = _require(builder, "rf_phase_csv", "rf_phase")
    _require(builder, "rf_phase_eedf_csv", "rf_phase")
    found = False
    for case in builder.cases:
        selected = case_rows(frame, case.key).sort_values("phase_index")
        if selected.empty:
            continue
        found = True
        ident = builder.new(
            "rf_phase",
            case.context,
            note="Saved periodic electron response. This is not chemical-species time evolution. Heatmap uses the nearest saved phase sample on periodic midpoint cells.",
        )
        for key, title, unit in (
            ("mean_energy_eV", "Phase mean energy", "eV"),
            ("instantaneous_E_over_N_Td", "Instantaneous E/N", "Td"),
            (
                "ionization_rate_coefficient_m3_s",
                "Reported phase ionization coefficient",
                "m^3/s",
            ),
        ):
            builder.panel(ident, key, title, "1", unit, x_label="RF phase / cycle")
            builder.curve(
                ident,
                key,
                curve_metadata([case]),
                [
                    dict(
                        x=row.phase_fraction,
                        y=getattr(row, key),
                        ci95_low=None,
                        ci95_high=None,
                        replicates=1,
                        estimate_status="native_case",
                        numerical_status=status([case]),
                    )
                    for row in selected.itertuples()
                ],
            )
        builder.panel(
            ident,
            "density",
            "Phase EEDF",
            "1",
            "eV",
            grid="rf_phase",
            data_table="rf_phase_eedf_csv",
            source_case=case.key,
            numerical_status=status([case]),
        )
    if not found:
        raise ValueError("rf_phase figure has no native RF phase result")


_PREPARERS = {
    "eedf": _eedf,
    "transport": _transport,
    "reaction_rates": _rates,
    "energy_angle": _energy_angle,
    "rf_phase": _rf_phase,
}


def _grid(axis, panel, artifact):
    from matplotlib.colors import LogNorm

    data = case_rows(artifact.tables[panel["data_table"]], panel["source_case"])
    energy_data = data.drop_duplicates("energy_eV").sort_values("energy_eV")
    energy = energy_data.energy_eV.to_numpy(dtype=float)
    widths = energy_data.energy_width_eV.to_numpy(dtype=float)
    energy_edges = reported_cell_edges(energy, widths)
    if panel["grid"] == "energy_angle":
        coordinate = data.drop_duplicates("mu").sort_values("mu")
        coordinate_edges = np.r_[
            -1.0, -1.0 + np.cumsum(coordinate.mu_width.to_numpy(dtype=float))
        ]
        if not np.isclose(coordinate_edges[-1], 1.0):
            raise ValueError(
                "energy-angle figure requires native mu widths spanning [-1,1]"
            )
        values = (
            data.pivot(
                index="energy_eV", columns="mu", values="energy_angle_density_eV_inv"
            )
            .sort_index()
            .sort_index(axis=1)
            .to_numpy(dtype=float)
        )
        xlabel, ylabel, cunit = (
            "Electron energy [eV]",
            "Native polar cosine mu [1]",
            "density [1/eV]",
        )
        marginal = values @ coordinate.mu_width.to_numpy(dtype=float)
        limit = energy_display_limit(energy, widths, marginal)
        absolute_mass = (
            np.abs(values) @ coordinate.mu_width.to_numpy(dtype=float) * widths
        )
        x, y, colors = energy_edges, coordinate_edges, values.T
        axis.set_xlim(0, limit)
    else:
        coordinate = data.drop_duplicates("phase_index").sort_values("phase_index")
        fractions = coordinate.phase_fraction.to_numpy(dtype=float)
        if fractions[0] < 0 or np.any(np.diff(fractions) <= 0) or fractions[-1] >= 1:
            raise ValueError(
                "RF figure requires saved phase fractions increasing within [0,1)"
            )
        first_edge = 0.5 * (fractions[0] + fractions[-1] - 1)
        coordinate_edges = np.r_[
            first_edge, 0.5 * (fractions[1:] + fractions[:-1]), first_edge + 1
        ]
        values = (
            data.pivot(index="energy_eV", columns="phase_index", values="eedf")
            .sort_index()
            .sort_index(axis=1)
            .to_numpy(dtype=float)
        )
        xlabel, ylabel, cunit = (
            "RF phase / cycle [1]",
            "Electron energy [eV]",
            "EEDF [1/eV]",
        )
        limit = max(
            energy_display_limit(energy, widths, values[:, index])
            for index in range(values.shape[1])
        )
        absolute_mass = np.max(np.abs(values), axis=1) * widths
        colors = values
        if coordinate_edges[0] > 0:
            coordinate_edges = np.r_[coordinate_edges[-2] - 1, coordinate_edges]
            colors = np.column_stack([values[:, -1], colors])
        if coordinate_edges[-1] < 1:
            coordinate_edges = np.r_[coordinate_edges, coordinate_edges[1] + 1]
            colors = np.column_stack([colors, values[:, 0]])
        x, y = coordinate_edges, energy_edges
        axis.set_xlim(0, 1)
        axis.set_ylim(0, limit)
    if np.any(~np.isfinite(values)):
        raise ValueError("figure grid has missing or nonfinite native cells")
    visible = values[energy_edges[:-1] < limit]
    positive = visible[visible > 0]
    logarithmic = (
        len(positive) > 0 and np.all(colors >= 0) and positive.max() > positive.min()
    )
    image = axis.pcolormesh(
        x,
        y,
        np.ma.masked_less_equal(colors, 0) if logarithmic else colors,
        shading="flat",
        cmap="viridis",
        norm=LogNorm(vmin=positive.min(), vmax=positive.max()) if logarithmic else None,
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.figure.colorbar(image, ax=axis, label=cunit)
    panel["display_energy_max_eV"] = limit
    panel["outside_display_absolute_mass_bound"] = float(
        absolute_mass[energy_edges[:-1] >= limit].sum()
    )
    panel["colour_scale"] = (
        "log; nonpositive cells omitted" if logarithmic else "linear"
    )
    panel["colour_limits"] = list(image.get_clim())


def _draw(ident, spec, points, artifact, directory, base_name):
    # FigureCanvasAgg avoids changing a user's pyplot backend or open figures.
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib import rc_context

    panels = list(spec["panels"])
    columns = 2 if len(panels) == 4 else min(3, len(panels))
    rows = int(np.ceil(len(panels) / columns))
    with rc_context({"font.size": 9, "svg.fonttype": "none", "text.usetex": False}):
        figure = Figure(figsize=(5.2 * columns, 3.8 * rows + 2.1))
        FigureCanvasAgg(figure)
        axes = np.atleast_1d(figure.subplots(rows, columns)).ravel()
        notes = list(spec["notes"])
        for axis, key in zip(axes, panels):
            panel = spec["panels"][key]
            axis.set_title(panel["label"].replace("$", r"\$"), fontsize=10)
            if "grid" in panel:
                _grid(axis, panel, artifact)
                axis.text(
                    0.01,
                    0.01,
                    fill(panel["numerical_status"], 50),
                    transform=axis.transAxes,
                    fontsize=6.5,
                    va="bottom",
                    bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
                )
                notes.append(
                    f"{key}: no rescaling; display energy <= {panel['display_energy_max_eV']:g} eV; outside-view absolute mass <= {panel['outside_display_absolute_mass_bound']:.4g}."
                )
                continue
            selected = points.loc[(points.figure_id == ident) & (points.panel == key)]
            for curve_id in selected.curve.unique():
                frame = selected.loc[selected.curve == curve_id].sort_values("x")
                label = fill(spec["curves"][curve_id]["label"].replace("$", r"\$"), 65)
                x, y = frame.x.to_numpy(dtype=float), frame.y.to_numpy(dtype=float)
                style = (
                    "--" if frame.numerical_status.str.contains("stale").any() else "-"
                )
                if frame.x_lower.notna().all() and frame.x_upper.notna().all():
                    edges = np.r_[
                        frame.x_lower.to_numpy(dtype=float),
                        float(frame.x_upper.iloc[-1]),
                    ]
                    line = axis.stairs(
                        y,
                        edges,
                        baseline=None,
                        linestyle=style,
                        linewidth=1.2,
                        label=label,
                    )
                    color = line.get_edgecolor()
                else:
                    line = axis.plot(
                        x,
                        y,
                        marker="o",
                        markersize=3,
                        linestyle=style if len(set(x)) == len(x) else "None",
                        linewidth=1.2,
                        label=label,
                    )[0]
                    color = line.get_color()
                low, high = (
                    frame.ci95_low.to_numpy(dtype=float),
                    frame.ci95_high.to_numpy(dtype=float),
                )
                available = np.isfinite(low) & np.isfinite(high)
                if panel.get("logarithmic"):
                    available &= (low > 0) & (high > 0)
                if available.any():
                    if frame.x_lower.notna().all() and frame.x_upper.notna().all():
                        axis.fill_between(
                            edges,
                            np.r_[low, low[-1]],
                            np.r_[high, high[-1]],
                            where=np.r_[available, available[-1]],
                            step="post",
                            color=color,
                            alpha=0.16,
                        )
                    else:
                        axis.fill_between(
                            x, low, high, where=available, color=color, alpha=0.16
                        )
            if panel.get("logarithmic") and not selected.empty:
                values = selected.y.to_numpy(dtype=float)
                if np.any(values < 0):
                    axis.set_yscale(
                        "symlog",
                        linthresh=max(float(np.max(np.abs(values))) * 1e-6, 1e-300),
                    )
                    panel["scale_y"] = "symlog; signed source values retained"
                elif np.any(values > 0):
                    axis.set_yscale("log")
                    panel["scale_y"] = "log; zero values not displayed"
            if panel.get("log_x"):
                axis.set_xscale("log")
            if "x_limit" in panel:
                axis.set_xlim(*panel["x_limit"])
                lower, upper = panel["x_limit"]
                visible_points = selected.loc[
                    (selected.x_upper > lower) & (selected.x_lower < upper), "y"
                ]
                visible = visible_points.to_numpy(dtype=float)
                if len(visible) and np.any(visible != 0):
                    low, high = float(visible.min()), float(visible.max())
                    if axis.get_yscale() == "log":
                        intervals = (
                            selected.loc[
                                visible_points.index, ["ci95_low", "ci95_high"]
                            ]
                            .to_numpy(dtype=float)
                            .ravel()
                        )
                        positive = np.r_[
                            visible[visible > 0],
                            intervals[np.isfinite(intervals) & (intervals > 0)],
                        ]
                        axis.set_ylim(
                            float(positive.min()) / 2, float(positive.max()) * 2
                        )
                    elif low < high:
                        axis.set_ylim(
                            low - 0.05 * (high - low), high + 0.05 * (high - low)
                        )
                panel["display_y_limits"] = list(axis.get_ylim())
            axis.set_xlabel(
                f"{panel.get('x_label', 'Electron energy')} [{panel['x_unit']}]"
            )
            axis.set_ylabel(f"{panel['label'].splitlines()[0]} [{panel['y_unit']}]")
            axis.grid(True, alpha=0.22)
            if not selected.empty:
                axis.legend(fontsize=7, loc="best")
                states = "; ".join(sorted(set(selected.numerical_status)))
                axis.text(
                    0.01,
                    0.01,
                    fill(states, 65),
                    transform=axis.transAxes,
                    fontsize=6.5,
                    va="bottom",
                    alpha=0.8,
                )
        for axis in axes[len(panels) :]:
            axis.set_visible(False)
        title = context_label(spec["physical_context"])
        heading = fill(f"{spec['kind']} | {title}", 50 * columns)
        figure.suptitle(heading, fontsize=11, y=0.985)
        notes.append(
            "Calculated electron data; numerical status and statistical intervals are not independent physical validation."
        )
        if any(panel.get("logarithmic") for panel in spec["panels"].values()):
            notes.append(
                "Log axes omit zero values and interval bounds <= 0; signed data use symlog. Values remain unchanged in CSV."
            )
        figure.text(
            0.02,
            0.015,
            fill(" ".join(notes), max(80, 50 * columns)),
            fontsize=7.5,
            va="bottom",
        )
        top = 1 - (0.18 * len(heading.splitlines()) + 0.15) / figure.get_figheight()
        figure.tight_layout(rect=(0, 0.18, 1, top))
        files, paths = {}, {}
        for extension in ("png", "svg"):
            path = directory / f"{base_name}_figure_{ident}.{extension}"
            figure.savefig(path, dpi=150, metadata={"Creator": "electron-swarm"})
            files[extension] = {
                "name": path.name,
                "sha256": sha256(path.read_bytes()).hexdigest(),
            }
            paths[f"figure_{ident}_{extension}"] = path
        figure.clear()
    return files, paths


def render_figures(artifact, requested, directory, base_name) -> RenderedFigures:
    builder = _Figures(artifact)
    for kind in requested:
        if kind not in _PREPARERS:
            raise ValueError(f"unsupported figure kind: {kind}")
        before = len(builder.specs)
        _PREPARERS[kind](builder)
        if len(builder.specs) == before:
            raise ValueError(f"requested {kind} figure has no usable source cases")
    points = pd.DataFrame(builder.rows, columns=FIGURE_DATA_KINDS, dtype=object)
    paths = {}
    for ident, spec in builder.specs.items():
        files, generated = _draw(ident, spec, points, artifact, directory, base_name)
        spec["files"] = files
        spec["generator_version"] = __version__
        spec["curve_data"] = "figure_data_csv"
        paths.update(generated)
    return RenderedFigures(points, builder.specs, paths)


def read_figures(manifest_path):
    from .publication import read_generation
    from .artifacts import PRODUCT_ARTIFACT_SCHEMA
    from .reaction_rates import RATE_ARTIFACT_SCHEMA

    path = Path(manifest_path).resolve()
    if not path.name.endswith("_manifest.json"):
        raise ValueError("expected a figure _manifest.json")
    result = read_generation(
        path,
        path.name[: -len("_manifest.json")],
        {
            PRODUCT_ARTIFACT_SCHEMA: "swarm_product_run",
            RATE_ARTIFACT_SCHEMA: "reaction_rates",
            FIGURE_ARTIFACT_SCHEMA: "swarm_figures",
        },
    )
    if not result.metadata.get("figures"):
        raise ValueError("manifest contains no figures")
    return result
