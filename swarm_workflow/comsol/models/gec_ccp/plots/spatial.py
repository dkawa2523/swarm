"""Spatial, domain, radial, phase, and waveform GEC CCP rendering."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from swarm_workflow.comsol.models.gec_ccp.plots import contracts as _contracts
from swarm_workflow.comsol.models.gec_ccp.plots import data as _data


def _plot_spatial_cut_overview(
    axis_paths: tuple[Path, Path],
    radial_paths: tuple[Path, Path],
    output_stem: Path,
    plt: Any,
    *,
    external_label: str = "External closure",
) -> list[Path]:
    """Plot four period-averaged fields along axial and radial cuts.

    Chart contract: compare absolute spatial profiles from the built-in and
    Swarm-table closures at fixed 1 W.  Rows are physical quantities and
    columns are the axial centerline and radial midplane cuts.  All four fields
    use linear axes so absolute profile magnitudes remain directly comparable.
    """
    axis_series = tuple(_data._comsol_spatial_data(path) for path in axis_paths)
    radial_series = tuple(_data._comsol_spatial_data(path) for path in radial_paths)
    count = min(
        *(fields.shape[1] for _, fields in (*axis_series, *radial_series)),
        4,
    )
    if count < 1:
        raise _contracts.GecCcpPlotError("COMSOL spatial cuts contain no result fields")

    labels = (
        "Electron density",
        "Electron temperature",
        "Electric potential",
        "Ionization source",
    )
    units = (r"m$^{-3}$", "eV", "V", r"m$^{-3}$ s$^{-1}$")
    cuts = (
        ("Axial centerline (r = 0 cm)", axis_series, "Axial position z [cm]"),
        ("Radial midplane (z = 1.27 cm)", radial_series, "Radius r [cm]"),
    )
    fig, axes = plt.subplots(
        count,
        2,
        figsize=(13, 2.75 * count),
        sharey="row",
        constrained_layout=True,
    )
    if count == 1:
        axes = [axes]
    for column, (title, series, x_label) in enumerate(cuts):
        axes[0][column].set_title(title)
        for index in range(count):
            chart = axes[index][column]
            for (coordinate, fields), color, style, name in (
                (series[0], "#555555", "--", "built-in Druyvesteyn"),
                (series[1], "#1769aa", "-", external_label),
            ):
                chart.plot(
                    coordinate * 100.0,
                    fields[:, index],
                    color=color,
                    linestyle=style,
                    linewidth=2,
                    label=name,
                )
            chart.grid(True, which="both", color="#777777", alpha=0.18)
            if column == 0:
                chart.set_ylabel(f"{labels[index]}\n[{units[index]}]")
            if index == count - 1:
                chart.set_xlabel(x_label)
    axes[0][0].legend(frameon=False, ncol=2, fontsize=9)
    fig.suptitle(
        "Argon GEC CCP period-averaged spatial distributions\n"
        f"fixed 1 W; dashed: built-in Druyvesteyn, solid: {external_label}",
        fontsize=14,
    )
    return _data._save(fig, output_stem, plt)


def _plot_gec_geometry_2d(
    domain_path: Path,
    output_stem: Path,
    plt: Any,
) -> list[Path]:
    """Draw the axisymmetric GEC-CCP gas domain and electrical boundaries."""
    from matplotlib.patches import Polygon, Rectangle

    coordinates, _ = _data._comsol_domain_data(domain_path)
    geometry = _gec_domain_geometry(coordinates)
    polygon = _gec_domain_polygon(geometry)
    fig, axis = plt.subplots(figsize=(9, 8), constrained_layout=True)
    axis.add_patch(
        Polygon(
            polygon,
            closed=True,
            facecolor="#dbeaf5",
            edgecolor="#222222",
            linewidth=1.8,
            label="Argon gas / plasma domain",
        )
    )
    gap_low = geometry["gap_z_min_cm"]
    gap_high = geometry["gap_z_max_cm"]
    powered_radius = geometry["powered_radius_cm"]
    throat = geometry["throat_radius_cm"]
    chamber = geometry["chamber_radius_cm"]
    axis.add_patch(
        Rectangle(
            (0.0, gap_low - 0.55),
            powered_radius,
            0.55,
            facecolor="#dddddd",
            edgecolor="#555555",
            hatch="///",
            linewidth=1.0,
        )
    )
    axis.add_patch(
        Rectangle(
            (0.0, gap_high),
            throat,
            0.55,
            facecolor="#dddddd",
            edgecolor="#555555",
            hatch="///",
            linewidth=1.0,
        )
    )
    _draw_gec_boundaries(axis, geometry)
    axis.text(
        2.5,
        0.5 * (gap_low + gap_high),
        "inter-electrode\nplasma",
        ha="center",
        va="center",
        fontsize=11,
        color="#0b3558",
    )
    axis.text(
        7.7,
        0.5 * (gap_low + gap_high),
        "outer chamber\nplasma",
        ha="center",
        va="center",
        fontsize=11,
        color="#0b3558",
    )
    axis.annotate(
        "powered electrode",
        xy=(2.5, gap_low),
        xytext=(2.5, gap_low - 1.0),
        ha="center",
        arrowprops={"arrowstyle": "->", "color": _contracts.POWERED_BOUNDARY_COLOR},
        color=_contracts.POWERED_BOUNDARY_COLOR,
    )
    axis.annotate(
        "grounded electrode",
        xy=(2.7, gap_high),
        xytext=(2.7, gap_high + 1.0),
        ha="center",
        arrowprops={"arrowstyle": "->", "color": _contracts.GROUNDED_BOUNDARY_COLOR},
        color=_contracts.GROUNDED_BOUNDARY_COLOR,
    )
    axis.annotate(
        "grounded chamber wall",
        xy=(chamber, 1.2),
        xytext=(8.1, 4.9),
        ha="center",
        arrowprops={"arrowstyle": "->", "color": _contracts.GROUNDED_BOUNDARY_COLOR},
        color=_contracts.GROUNDED_BOUNDARY_COLOR,
    )
    axis.annotate(
        "axis of symmetry",
        xy=(0.0, 1.2),
        xytext=(1.0, 4.4),
        arrowprops={"arrowstyle": "->", "color": "#777777"},
        color="#555555",
    )
    axis.annotate(
        "",
        xy=(4.6, gap_low),
        xytext=(4.6, gap_high),
        arrowprops={"arrowstyle": "<->", "color": "#333333"},
    )
    axis.text(
        4.72,
        0.5 * (gap_low + gap_high),
        f"gap {gap_high - gap_low:.2f} cm",
        rotation=90,
        va="center",
        fontsize=9,
    )
    axis.annotate(
        "",
        xy=(0.0, -1.0),
        xytext=(powered_radius, -1.0),
        arrowprops={"arrowstyle": "<->", "color": "#333333"},
    )
    axis.text(
        0.5 * powered_radius,
        -1.18,
        f"powered radius {powered_radius:.2f} cm",
        ha="center",
        va="top",
        fontsize=9,
    )
    axis.set_xlim(-0.5, chamber + 0.7)
    axis.set_ylim(
        geometry["outer_z_min_cm"] - 0.5,
        geometry["outer_z_max_cm"] + 0.5,
    )
    axis.set_aspect("equal")
    axis.set_xlabel("Radius r [cm]")
    axis.set_ylabel("Axial coordinate z [cm]")
    axis.set_title(
        "Axisymmetric GEC-CCP computational domain\n"
        "boundaries inferred from the exported COMSOL plasma mesh"
    )
    axis.grid(True, color="#777777", alpha=0.12)
    return _data._save(fig, output_stem, plt)


def _plot_domain_2d_comparison(
    baseline_path: Path,
    external_path: Path,
    output_stem: Path,
    plt: Any,
    *,
    baseline_label: str = "Original COMSOL\nBuilt-in Druyvesteyn",
    external_label: str = "External closure",
    additional_path: Path | None = None,
    additional_label: str | None = None,
    identify_electrodes: bool = False,
    include_ionization_source: bool = True,
    include_svg: bool = False,
) -> tuple[dict[str, Any], list[Path]]:
    """Plot COMSOL-style full-domain period-average comparisons.

    Each physical quantity is exported as one dedicated multi-panel figure,
    matching the single-surface focus of COMSOL's Graphics pane.  Both
    variants share the same linear normalization and RainbowLight
    approximation.
    """
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    baseline_coordinates, baseline = _data._comsol_domain_data(baseline_path)
    external_coordinates, external = _data._comsol_domain_data(external_path)
    if (additional_path is None) != (additional_label is None):
        raise _contracts.GecCcpPlotError(
            "additional COMSOL domain path and label must be supplied together"
        )
    coordinates = [baseline_coordinates, external_coordinates]
    field_tables = [baseline, external]
    panel_titles = [
        baseline_label,
        external_label,
    ]
    if additional_path is not None:
        additional_coordinates, additional = _data._comsol_domain_data(additional_path)
        coordinates.append(additional_coordinates)
        field_tables.append(additional)
        panel_titles.append(str(additional_label))
    if any(fields.shape[1] < 4 for fields in field_tables):
        raise _contracts.GecCcpPlotError("COMSOL 2D domain export needs four fields")
    geometry = _gec_domain_geometry(external_coordinates)
    triangulations = [_gec_domain_triangulation(item, geometry) for item in coordinates]
    rainbow_light = LinearSegmentedColormap.from_list(
        "comsol_rainbow_light",
        _contracts.COMSOL_RAINBOW_LIGHT_COLORS,
        N=1024,
    )
    all_specifications = (
        (0, "electron_density", "Electron density", r"n$_e$ [m$^{-3}$]", True),
        (
            1,
            "electron_temperature",
            "Electron temperature",
            "T$_e$ [eV]",
            False,
        ),
        (2, "electric_potential", "Electric potential", "V [V]", False),
        (
            3,
            "ionization_source",
            "Ionization source",
            r"R$_i$ [m$^{-3}$ s$^{-1}$]",
            True,
        ),
    )
    all_field_names = (
        "electron_density_m3",
        "electron_temperature_eV",
        "potential_V",
        "ionization_source_m3_s",
    )
    specifications = (
        all_specifications if include_ionization_source else all_specifications[:3]
    )
    field_names = all_field_names if include_ionization_source else all_field_names[:3]
    summary: dict[str, Any] = {
        "baseline_nodes": int(baseline_coordinates.shape[0]),
        "external_nodes": int(external_coordinates.shape[0]),
        "coordinate_match": bool(
            baseline_coordinates.shape == external_coordinates.shape
            and np.allclose(baseline_coordinates, external_coordinates)
        ),
        "geometry_cm": geometry,
        "presentation": {
            "style": "COMSOL 6.x Graphics-pane approximation",
            "absolute_color_table": "RainbowLight approximation",
            "color_levels": 1024,
            "layout": f"one {len(field_tables)}-panel figure per quantity",
            "shared_limits_between_closures": True,
            "color_scale": "linear",
            "electrical_boundaries_identified": identify_electrodes,
        },
        "reference_label": baseline_label,
        "comparison_label": external_label,
        "fields": {},
    }
    if additional_path is not None:
        summary["additional_nodes"] = int(coordinates[2].shape[0])
        summary["additional_coordinate_match"] = bool(
            baseline_coordinates.shape == coordinates[2].shape
            and np.allclose(baseline_coordinates, coordinates[2])
        )
        summary["additional_label"] = additional_label
    individual_paths: list[Path] = []
    for index, slug, title, colorbar_label, zero_based in specifications:
        raw_values = tuple(fields[:, index] for fields in field_tables)
        finite = np.concatenate([values[np.isfinite(values)] for values in raw_values])
        if finite.size == 0:
            raise _contracts.GecCcpPlotError(f"COMSOL 2D field is empty: {title}")
        low = float(np.min(finite))
        high = float(np.max(finite))
        if zero_based and low >= 0.0:
            low = 0.0
        if high <= low:
            high = low + max(abs(low) * 1.0e-12, 1.0e-12)
        norm = Normalize(vmin=low, vmax=high)
        field_fig, field_axes = plt.subplots(
            1,
            len(field_tables),
            figsize=(5.6 * len(field_tables), 5.3),
            constrained_layout=True,
        )
        field_meshes: list[Any] = []
        for field_axis, triangulation, values, panel_title in zip(
            field_axes,
            triangulations,
            raw_values,
            panel_titles,
        ):
            field_meshes.append(
                field_axis.tripcolor(
                    triangulation,
                    values,
                    shading="gouraud",
                    cmap=rainbow_light,
                    norm=norm,
                )
            )
            _draw_gec_boundaries(
                field_axis,
                geometry,
                comsol_style=True,
                identify_electrodes=identify_electrodes,
            )
            _style_comsol_surface_axis(field_axis, geometry, show_ylabel=True)
            field_axis.set_title(panel_title, fontsize=10.5, pad=7)
        if identify_electrodes:
            from matplotlib.lines import Line2D

            field_axes[0].legend(
                handles=(
                    Line2D(
                        [0],
                        [0],
                        color=_contracts.POWERED_BOUNDARY_COLOR,
                        linewidth=2.0,
                        label="powered electrode",
                    ),
                    Line2D(
                        [0],
                        [0],
                        color=_contracts.GROUNDED_BOUNDARY_COLOR,
                        linewidth=2.0,
                        label="grounded electrode / wall",
                    ),
                    Line2D(
                        [0],
                        [0],
                        color="#666666",
                        linewidth=1.0,
                        linestyle=":",
                        label="symmetry axis",
                    ),
                ),
                loc="lower left",
                bbox_to_anchor=(0.0, -0.27),
                ncol=3,
                frameon=False,
                fontsize=7.5,
            )
        colorbar = field_fig.colorbar(
            field_meshes[-1],
            ax=field_axes,
            shrink=0.88,
            pad=0.02,
        )
        colorbar.set_label(colorbar_label)
        colorbar.outline.set_linewidth(0.7)
        colorbar.ax.tick_params(labelsize=8.5, width=0.7, length=3)
        field_fig.suptitle(
            f"Surface: {title}\nPeriod average, P0 = 1 W — shared linear legend range",
            x=0.04,
            ha="left",
            fontsize=13,
            color="#222222",
        )
        individual_paths.extend(
            _data._save(
                field_fig,
                output_stem.with_name(f"comsol_style_{slug}_comparison"),
                plt,
                include_svg=include_svg,
            )
        )
        baseline_values = raw_values[0]
        external_values = raw_values[1]
        field_summary: dict[str, Any] = {
            "baseline_min": float(np.nanmin(baseline_values)),
            "baseline_max": float(np.nanmax(baseline_values)),
            "external_min": float(np.nanmin(external_values)),
            "external_max": float(np.nanmax(external_values)),
            "shared_display_range": [low, high],
            "display_scale": "linear",
            "relative_l2_difference": (
                float(np.linalg.norm(external_values - baseline_values))
                / float(np.linalg.norm(baseline_values))
                if summary["coordinate_match"]
                and float(np.linalg.norm(baseline_values)) > 0.0
                else None
            ),
        }
        if additional_path is not None:
            additional_values = raw_values[2]
            field_summary.update(
                {
                    "additional_min": float(np.nanmin(additional_values)),
                    "additional_max": float(np.nanmax(additional_values)),
                    "additional_relative_l2_difference": (
                        float(np.linalg.norm(additional_values - baseline_values))
                        / float(np.linalg.norm(baseline_values))
                        if summary["additional_coordinate_match"]
                        and float(np.linalg.norm(baseline_values)) > 0.0
                        else None
                    ),
                }
            )
        summary["fields"][field_names[index]] = field_summary
    return summary, individual_paths


def _plot_radial_midplane_solver_comparison(
    paths: tuple[Path, ...],
    labels: tuple[str, ...],
    output_stem: Path,
    plt: Any,
    *,
    geometry: dict[str, float],
    case_ids: tuple[str, ...] | None = None,
    include_ionization_source: bool = True,
    include_svg: bool = False,
) -> tuple[dict[str, Any], list[Path]]:
    """Plot accepted closures along the GEC gap midplane.

    Each quantity receives a dedicated COMSOL-like line figure.  Every y-axis
    is linear; nonnegative source and density fields start at zero so the
    absolute peak magnitudes can be compared without logarithmic expansion.
    """
    import numpy as np

    if len(paths) < 2 or len(paths) != len(labels):
        raise _contracts.GecCcpPlotError(
            "radial solver comparison needs matching paths and labels for "
            "at least two cases"
        )
    if case_ids is None:
        case_ids = tuple(f"case_{index + 1}" for index in range(len(paths)))
    if len(case_ids) != len(paths) or len(set(case_ids)) != len(case_ids):
        raise _contracts.GecCcpPlotError(
            "radial solver comparison case ids must be unique and match paths"
        )

    center_height_cm = 0.5 * (geometry["gap_z_min_cm"] + geometry["gap_z_max_cm"])
    series = tuple(
        _data._comsol_radial_midplane_data(path, center_height_cm) for path in paths
    )
    if any(fields.shape[1] < 4 for _, fields, _ in series):
        raise _contracts.GecCcpPlotError(
            "COMSOL radial-midplane export needs four fields"
        )

    available_styles = (
        ("#3f3f3f", "--"),
        ("#1769aa", "-"),
        ("#d97706", "-."),
    )
    styles = tuple(
        available_styles[index % len(available_styles)] for index in range(len(paths))
    )
    all_specifications = (
        (
            0,
            "electron_density",
            "Electron density",
            r"n$_e$ [m$^{-3}$]",
            "zero_based",
        ),
        (
            1,
            "electron_temperature",
            "Electron temperature",
            r"T$_e$ [eV]",
            "focused",
        ),
        (
            2,
            "electric_potential",
            "Electric potential",
            "V [V]",
            "include_zero",
        ),
        (
            3,
            "ionization_source",
            "Ionization source (ptp.Re_av)",
            r"R$_i$ [m$^{-3}$ s$^{-1}$]",
            "zero_based",
        ),
    )
    specifications = (
        all_specifications if include_ionization_source else all_specifications[:3]
    )
    baseline_radius = series[0][0]
    coordinate_match = all(
        radius.shape == baseline_radius.shape
        and np.allclose(radius, baseline_radius, rtol=0.0, atol=1.0e-10)
        for radius, _, _ in series[1:]
    )
    summary: dict[str, Any] = {
        "center_height_cm": center_height_cm,
        "actual_center_heights_cm": {
            case_id: actual_height
            for case_id, (_, _, actual_height) in zip(case_ids, series)
        },
        "coordinate_match": coordinate_match,
        "sample_counts": {
            case_id: int(radius.size)
            for case_id, (radius, _, _) in zip(case_ids, series)
        },
        "radius_range_cm": [
            float(min(np.min(radius) for radius, _, _ in series)),
            float(max(np.max(radius) for radius, _, _ in series)),
        ],
        "presentation": {
            "style": "COMSOL 6.x 1D Plot Group approximation",
            "layout": f"one {len(paths)}-series figure per quantity",
            "y_scale": "linear",
            "shared_axis_between_closures": True,
            "density_and_source_zero_based": True,
            "geometry_markers_cm": {
                "powered_electrode_edge": geometry["powered_radius_cm"],
                "chamber_throat": geometry["throat_radius_cm"],
            },
        },
        "electron_temperature_unit_note": (
            "COMSOL exports ptp.Teav in V; it is displayed as the "
            "numerically equivalent electron-temperature energy in eV"
        ),
        "fields": {},
    }
    figures: list[Path] = []
    for index, slug, title, y_label, range_policy in specifications:
        values_by_case = tuple(fields[:, index] for _, fields, _ in series)
        finite = np.concatenate(
            [values[np.isfinite(values)] for values in values_by_case]
        )
        if finite.size == 0:
            raise _contracts.GecCcpPlotError(f"COMSOL radial field is empty: {title}")
        data_low = float(np.min(finite))
        data_high = float(np.max(finite))
        span = data_high - data_low
        padding = 0.05 * span if span > 0.0 else max(abs(data_high) * 0.05, 1.0)
        if range_policy == "zero_based" and data_low >= 0.0:
            display_low = 0.0
            display_high = data_high + max(0.05 * data_high, padding)
        else:
            display_low = data_low - padding
            display_high = data_high + padding
            if range_policy == "include_zero":
                display_low = min(0.0, display_low)
                display_high = max(0.0, display_high)
        if display_high <= display_low:
            display_high = display_low + max(abs(display_low) * 0.05, 1.0)

        fig, axis = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
        for (radius, fields, _), label, (color, linestyle) in zip(
            series,
            labels,
            styles,
        ):
            values = fields[:, index]
            valid = np.isfinite(radius) & np.isfinite(values)
            axis.plot(
                radius[valid],
                values[valid],
                color=color,
                linestyle=linestyle,
                linewidth=2.2,
                label=label,
            )
        axis.axvline(
            geometry["powered_radius_cm"],
            color=_contracts.POWERED_BOUNDARY_COLOR,
            linestyle=":",
            linewidth=1.2,
        )
        axis.axvline(
            geometry["throat_radius_cm"],
            color=_contracts.GROUNDED_BOUNDARY_COLOR,
            linestyle=":",
            linewidth=1.2,
        )
        axis.text(
            geometry["powered_radius_cm"],
            0.02,
            "powered edge",
            color=_contracts.POWERED_BOUNDARY_COLOR,
            fontsize=7.5,
            rotation=90,
            ha="right",
            va="bottom",
            transform=axis.get_xaxis_transform(),
        )
        axis.text(
            geometry["throat_radius_cm"],
            0.02,
            "chamber throat",
            color=_contracts.GROUNDED_BOUNDARY_COLOR,
            fontsize=7.5,
            rotation=90,
            ha="left",
            va="bottom",
            transform=axis.get_xaxis_transform(),
        )
        axis.set_xlim(*summary["radius_range_cm"])
        axis.set_ylim(display_low, display_high)
        axis.set_xlabel("Radius r [cm]")
        axis.set_ylabel(y_label)
        axis.set_title(
            f"Line: {title}\n"
            f"Period average at gap midplane z = {center_height_cm:.2f} cm, "
            "P0 = 1 W — linear y-axis",
            loc="left",
            fontsize=12.5,
            color="#222222",
            pad=10,
        )
        axis.grid(True, color="#777777", alpha=0.18, linewidth=0.7)
        axis.ticklabel_format(
            axis="y",
            style="sci",
            scilimits=(-3, 4),
            useMathText=True,
        )
        axis.legend(loc="upper right", frameon=False, fontsize=9)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        figures.extend(
            _data._save(
                fig,
                output_stem.with_name(
                    f"comsol_style_radial_midplane_{slug}_comparison"
                ),
                plt,
                include_svg=include_svg,
            )
        )
        summary["fields"][slug] = {
            "display_scale": "linear",
            "display_range": [display_low, display_high],
            "series": {
                case_id: {
                    "minimum": float(np.nanmin(values)),
                    "maximum": float(np.nanmax(values)),
                }
                for case_id, values in zip(case_ids, values_by_case)
            },
        }
    return summary, figures


def _style_comsol_surface_axis(
    axis: Any,
    geometry: dict[str, float],
    *,
    show_ylabel: bool,
) -> None:
    """Apply a restrained COMSOL Graphics-pane-like 2D axis style."""
    axis.set_aspect("equal")
    axis.set_xlim(-0.12, geometry["chamber_radius_cm"] + 0.12)
    axis.set_ylim(
        geometry["outer_z_min_cm"] - 0.12,
        geometry["outer_z_max_cm"] + 0.12,
    )
    axis.set_xlabel("r (cm)", fontsize=9)
    axis.set_ylabel("z (cm)" if show_ylabel else "", fontsize=9)
    axis.set_facecolor("white")
    axis.grid(False)
    axis.tick_params(
        axis="both",
        direction="out",
        colors="#333333",
        labelsize=8.5,
        width=0.7,
        length=3.5,
    )
    for spine in axis.spines.values():
        spine.set_color("#4a4a4a")
        spine.set_linewidth(0.7)


def _plot_domain_2d_difference(
    baseline_path: Path,
    external_path: Path,
    output_stem: Path,
    plt: Any,
    *,
    external_label: str = "External closure",
) -> tuple[dict[str, Any], list[Path]]:
    """Plot signed spatial differences between Swarm and original COMSOL."""
    import numpy as np
    import matplotlib.tri as mtri
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    baseline_coordinates, baseline = _data._comsol_domain_data(baseline_path)
    external_coordinates, external = _data._comsol_domain_data(external_path)
    if baseline_coordinates.shape != external_coordinates.shape or not np.allclose(
        baseline_coordinates, external_coordinates
    ):
        raise _contracts.GecCcpPlotError(
            "2D difference maps require matching original and Swarm meshes"
        )
    geometry = _gec_domain_geometry(external_coordinates)
    triangulation = _gec_domain_triangulation(external_coordinates, geometry)
    diverging = LinearSegmentedColormap.from_list(
        "comsol_difference_light",
        _contracts.COMSOL_DIFFERENCE_COLORS,
        N=1024,
    )
    density_threshold = 1.0e-3 * float(np.nanmax(baseline[:, 0]))
    source_threshold = 1.0e-3 * float(np.nanmax(baseline[:, 3]))
    density_ratio = np.divide(
        external[:, 0],
        baseline[:, 0],
        out=np.full_like(external[:, 0], np.nan),
        where=baseline[:, 0] > 0.0,
    )
    source_ratio = np.divide(
        external[:, 3],
        baseline[:, 3],
        out=np.full_like(external[:, 3], np.nan),
        where=baseline[:, 3] > 0.0,
    )
    specifications = (
        (
            "Electron density change",
            100.0 * (density_ratio - 1.0),
            baseline[:, 0] >= density_threshold,
            r"100(n$_{e,external}$/n$_{e,original}$ - 1) [%]",
        ),
        (
            "Electron temperature change",
            external[:, 1] - baseline[:, 1],
            np.isfinite(baseline[:, 1]),
            r"T$_{e,external}$ - T$_{e,original}$ [eV]",
        ),
        (
            "Electric potential change",
            external[:, 2] - baseline[:, 2],
            np.isfinite(baseline[:, 2]),
            r"V$_{external}$ - V$_{original}$ [V]",
        ),
        (
            "Ionization-source ratio",
            np.log10(source_ratio),
            baseline[:, 3] >= source_threshold,
            r"log$_{10}$(R$_{i,external}$/R$_{i,original}$)",
        ),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    summary: dict[str, Any] = {
        "definition": "external closure minus original COMSOL on the common mesh",
        "external_label": external_label,
        "ratio_mask": (
            "electron-density and ionization ratios exclude nodes below "
            "0.1% of the respective original peak"
        ),
        "fields": {},
    }
    for axis, (title, raw_difference, valid, colorbar_label) in zip(
        axes.flat, specifications
    ):
        difference = np.where(
            valid & np.isfinite(raw_difference), raw_difference, np.nan
        )
        finite = difference[np.isfinite(difference)]
        if finite.size == 0:
            raise _contracts.GecCcpPlotError(f"COMSOL 2D difference is empty: {title}")
        extent = max(float(np.max(np.abs(finite))), 1.0e-12)
        norm = TwoSlopeNorm(vmin=-extent, vcenter=0.0, vmax=extent)
        metric_tri = mtri.Triangulation(
            triangulation.x,
            triangulation.y,
            triangles=triangulation.triangles,
        )
        base_mask = (
            np.zeros(metric_tri.triangles.shape[0], dtype=bool)
            if triangulation.mask is None
            else np.asarray(triangulation.mask, dtype=bool)
        )
        invalid_triangles = ~np.all(
            np.isfinite(difference)[metric_tri.triangles], axis=1
        )
        metric_tri.set_mask(base_mask | invalid_triangles)
        plot_values = np.where(np.isfinite(difference), difference, 0.0)
        mesh = axis.tripcolor(
            metric_tri,
            plot_values,
            shading="gouraud",
            cmap=diverging,
            norm=norm,
        )
        _draw_gec_boundaries(axis, geometry, comsol_style=True)
        _style_comsol_surface_axis(axis, geometry, show_ylabel=True)
        axis.set_title(title)
        colorbar = fig.colorbar(
            mesh,
            ax=axis,
            shrink=0.88,
            pad=0.015,
            label=colorbar_label,
        )
        colorbar.outline.set_linewidth(0.7)
        summary["fields"][title] = {
            "minimum": float(np.min(finite)),
            "maximum": float(np.max(finite)),
            "valid_nodes": int(finite.size),
        }
    fig.suptitle(
        f"Argon GEC CCP — {external_label} minus original COMSOL\n"
        "red: external closure is higher; blue: external closure is lower",
        x=0.04,
        ha="left",
        fontsize=14,
        color="#222222",
    )
    return summary, _data._save(fig, output_stem, plt, include_svg=False)


def _gec_domain_geometry(coordinates: Any) -> dict[str, float]:
    """Infer the stepped GEC gas-domain dimensions from exported mesh nodes."""
    import numpy as np

    radius_cm = coordinates[:, 0] * 100.0
    axial_cm = coordinates[:, 1] * 100.0
    chamber_radius = float(np.max(radius_cm))
    axis_radius = float(np.min(radius_cm))
    axis_tolerance = max(1.0e-8, 1.0e-6 * chamber_radius)
    on_axis = np.isclose(radius_cm, axis_radius, atol=axis_tolerance)
    if int(np.sum(on_axis)) < 2:
        raise _contracts.GecCcpPlotError(
            "COMSOL domain export does not resolve the symmetry axis"
        )
    gap_low = float(np.min(axial_cm[on_axis]))
    gap_high = float(np.max(axial_cm[on_axis]))
    outside_gap = (axial_cm < gap_low - 1.0e-6) | (axial_cm > gap_high + 1.0e-6)
    if not np.any(outside_gap):
        raise _contracts.GecCcpPlotError(
            "COMSOL domain export does not resolve the outer chamber"
        )
    throat = float(np.min(radius_cm[outside_gap]))
    powered_radius = float(radius_cm[np.argmin(np.abs(radius_cm - 5.08))])
    return {
        "axis_radius_cm": axis_radius,
        "powered_radius_cm": powered_radius,
        "throat_radius_cm": throat,
        "chamber_radius_cm": chamber_radius,
        "gap_z_min_cm": gap_low,
        "gap_z_max_cm": gap_high,
        "outer_z_min_cm": float(np.min(axial_cm)),
        "outer_z_max_cm": float(np.max(axial_cm)),
    }


def _gec_domain_polygon(geometry: dict[str, float]) -> list[tuple[float, float]]:
    axis = geometry["axis_radius_cm"]
    throat = geometry["throat_radius_cm"]
    chamber = geometry["chamber_radius_cm"]
    gap_low = geometry["gap_z_min_cm"]
    gap_high = geometry["gap_z_max_cm"]
    outer_low = geometry["outer_z_min_cm"]
    outer_high = geometry["outer_z_max_cm"]
    return [
        (axis, gap_low),
        (throat, gap_low),
        (throat, outer_low),
        (chamber, outer_low),
        (chamber, outer_high),
        (throat, outer_high),
        (throat, gap_high),
        (axis, gap_high),
    ]


def _gec_domain_triangulation(
    coordinates: Any,
    geometry: dict[str, float],
) -> Any:
    import matplotlib.tri as mtri
    import numpy as np

    radius = coordinates[:, 0] * 100.0
    axial = coordinates[:, 1] * 100.0
    triangulation = mtri.Triangulation(radius, axial)
    triangles = triangulation.triangles
    center_r = np.mean(radius[triangles], axis=1)
    center_z = np.mean(axial[triangles], axis=1)
    throat = geometry["throat_radius_cm"]
    in_gap = (
        (center_r <= throat + 1.0e-8)
        & (center_z >= geometry["gap_z_min_cm"] - 1.0e-8)
        & (center_z <= geometry["gap_z_max_cm"] + 1.0e-8)
    )
    in_outer = (
        (center_r >= throat - 1.0e-8)
        & (center_r <= geometry["chamber_radius_cm"] + 1.0e-8)
        & (center_z >= geometry["outer_z_min_cm"] - 1.0e-8)
        & (center_z <= geometry["outer_z_max_cm"] + 1.0e-8)
    )
    triangulation.set_mask(~(in_gap | in_outer))
    return triangulation


def _draw_gec_boundaries(
    axis: Any,
    geometry: dict[str, float],
    *,
    comsol_style: bool = False,
    identify_electrodes: bool = False,
) -> None:
    polygon = _gec_domain_polygon(geometry)
    closed = [*polygon, polygon[0]]
    edge_color = "#303030" if comsol_style else "#222222"
    use_boundary_colors = identify_electrodes or not comsol_style
    powered_color = (
        _contracts.POWERED_BOUNDARY_COLOR if use_boundary_colors else edge_color
    )
    grounded_color = (
        _contracts.GROUNDED_BOUNDARY_COLOR if use_boundary_colors else edge_color
    )
    axis.plot(
        [point[0] for point in closed],
        [point[1] for point in closed],
        color=edge_color,
        linewidth=0.9 if comsol_style else 1.5,
        zorder=5,
    )
    gap_low = geometry["gap_z_min_cm"]
    gap_high = geometry["gap_z_max_cm"]
    powered_radius = geometry["powered_radius_cm"]
    throat = geometry["throat_radius_cm"]
    chamber = geometry["chamber_radius_cm"]
    axis.plot(
        [0.0, powered_radius],
        [gap_low, gap_low],
        color=powered_color,
        linewidth=1.5 if comsol_style else 4.0,
        solid_capstyle="butt",
        zorder=7,
    )
    axis.plot(
        [0.0, throat],
        [gap_high, gap_high],
        color=grounded_color,
        linewidth=1.2 if comsol_style else 3.0,
        solid_capstyle="butt",
        zorder=7,
    )
    grounded_segments = (
        ([throat, chamber], [geometry["outer_z_min_cm"]] * 2),
        ([throat, chamber], [geometry["outer_z_max_cm"]] * 2),
        ([chamber, chamber], [geometry["outer_z_min_cm"], geometry["outer_z_max_cm"]]),
    )
    for x_values, y_values in grounded_segments:
        axis.plot(
            x_values,
            y_values,
            color=grounded_color,
            linewidth=1.0 if comsol_style else 2.4,
            zorder=6,
        )
    axis.plot(
        [geometry["axis_radius_cm"]] * 2,
        [gap_low, gap_high],
        color="#666666" if comsol_style else "#777777",
        linewidth=0.9 if comsol_style else 1.5,
        linestyle=":",
        zorder=8,
    )


def _plot_phase_axial_distributions(
    phase_path: Path,
    output_stem: Path,
    plt: Any,
    *,
    external_label: str = "External closure",
) -> tuple[dict[str, Any], list[Path]]:
    """Plot accepted Swarm-table fields versus axial position and RF phase.

    Chart contract: expose the space-time structure hidden by period averages
    for three directly solved COMSOL fields.  The grain is 20 centerline nodes
    by 51 RF phases.  Density is shown as log10 magnitude, temperature uses a
    sequential scale, and signed potential uses a zero-centered diverging
    scale when both signs occur.
    """
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    sequential = LinearSegmentedColormap.from_list(
        "neutral_blue", ("#f7f7f7", "#a8c7e6", "#1769aa", "#0b3558")
    )
    diverging = LinearSegmentedColormap.from_list(
        "blue_neutral_orange", ("#1769aa", "#f7f7f7", "#e69500")
    )
    specifications = (
        ("ptp.ne", "Electron density", r"log$_{10}$(n$_e$ [m$^{-3}$])", "log"),
        ("ptp.Te", "Electron temperature", "T$_e$ [eV]", "linear"),
        ("V", "Electric potential", "V [V]", "diverging"),
    )
    fig, axes = plt.subplots(
        len(specifications),
        1,
        figsize=(10, 9),
        sharex=True,
        constrained_layout=True,
    )
    summary: dict[str, Any] = {
        "source": external_label,
        "sample_scope": "center-axis nodes x RF phase",
        "fields": {},
    }
    for axis, (expression, title, unit, scale) in zip(axes, specifications):
        coordinate, values = _data._comsol_sorted_phase_field(phase_path, expression)
        phase = np.linspace(0.0, 1.0, values.shape[1])
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            raise _contracts.GecCcpPlotError(
                f"COMSOL phase field is empty: {expression}"
            )
        summary["fields"][expression] = {
            "minimum": float(np.min(finite)),
            "maximum": float(np.max(finite)),
            "nodes": int(values.shape[0]),
            "phase_samples": int(values.shape[1]),
        }
        if scale == "log":
            displayed = np.where(values > 0.0, np.log10(values), np.nan)
            mesh = axis.pcolormesh(
                phase,
                coordinate * 100.0,
                displayed,
                shading="auto",
                cmap=sequential,
            )
        elif scale == "diverging" and np.min(finite) < 0.0 < np.max(finite):
            mesh = axis.pcolormesh(
                phase,
                coordinate * 100.0,
                values,
                shading="auto",
                cmap=diverging,
                norm=TwoSlopeNorm(
                    vmin=float(np.min(finite)),
                    vcenter=0.0,
                    vmax=float(np.max(finite)),
                ),
            )
        else:
            mesh = axis.pcolormesh(
                phase,
                coordinate * 100.0,
                values,
                shading="auto",
                cmap=sequential,
            )
        axis.set_title(title)
        axis.set_ylabel("Axial position z [cm]")
        axis.set_xlim(0.0, 1.0)
        fig.colorbar(mesh, ax=axis, pad=0.015, label=unit)
    axes[-1].set_xlabel("RF period fraction")
    fig.suptitle(
        "Argon GEC CCP phase-resolved axial distributions\n"
        f"{external_label}, fixed 1 W; 20 nodes x 51 RF phases",
        fontsize=14,
    )
    return summary, _data._save(fig, output_stem, plt)


def _plot_waveform_comparison(
    baseline_path: Path,
    external_path: Path,
    output_stem: Path,
    plt: Any,
    *,
    external_label: str = "External closure",
) -> list[Path]:
    baseline_phase, baseline_voltage, baseline_current = _data._comsol_waveform_data(
        baseline_path
    )
    external_phase, external_voltage, external_current = _data._comsol_waveform_data(
        external_path
    )
    fig, axes = plt.subplots(
        2, 1, figsize=(11, 7), sharex=True, constrained_layout=True
    )
    for phase, voltage, current, color, style, label in (
        (
            baseline_phase,
            baseline_voltage,
            baseline_current,
            "#555555",
            "--",
            "built-in Druyvesteyn",
        ),
        (
            external_phase,
            external_voltage,
            external_current,
            "#1769aa",
            "-",
            external_label,
        ),
    ):
        axes[0].plot(phase, voltage, color=color, linestyle=style, label=label)
        axes[1].plot(phase, current, color=color, linestyle=style, label=label)
    axes[0].set_ylabel("Electrode voltage [V]")
    axes[1].set_ylabel("Electrode current [A]")
    axes[1].set_xlabel("RF period fraction")
    axes[0].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.grid(True, alpha=0.2)
    fig.suptitle("Argon GEC CCP powered-electrode waveform")
    return _data._save(fig, output_stem, plt)
