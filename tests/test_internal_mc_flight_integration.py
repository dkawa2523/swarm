from __future__ import annotations

import math

import numpy as np
import pytest

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import CrossSectionProcess, ProcessType
from electron_swarm.solvers.monte_carlo.cross_section_table import (
    NUMBA_AVAILABLE,
    PreparedCrossSectionTable,
)
from electron_swarm.solvers.monte_carlo.flight_integration import (
    accumulate_dc_flight_histogram,
    integrate_dc_flight_rates,
    new_rate_integration_workspace,
)


_SPEED_PER_SQRT_EV = math.sqrt(2.0 * EV_TO_J / ELECTRON_MASS_KG)


def _velocity(sqrt_energy_eV: float) -> np.ndarray:
    return np.asarray([_SPEED_PER_SQRT_EV * sqrt_energy_eV, 0.0, 0.0])


def _histogram_flight(
    before: np.ndarray,
    after: np.ndarray,
    *,
    duration_s: float,
    weight: float,
    edges_eV: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bins = edges_eV.size - 1
    counts = np.zeros(bins, dtype=np.int64)
    weighted = np.zeros(bins, dtype=float)
    squares = np.zeros(bins, dtype=float)
    scratch = np.zeros(bins, dtype=float)
    touched = np.empty(bins, dtype=np.int64)
    accumulate_dc_flight_histogram(
        before,
        after,
        duration_s,
        weight,
        edges_eV,
        counts,
        weighted,
        squares,
        scratch,
        touched,
    )
    return counts, weighted, squares


def test_exact_energy_crossings_conserve_residence_duration() -> None:
    edges = np.asarray([0.0, 1.0, 4.0, 9.0, 16.0])
    duration = 2.0
    weight = 3.0

    counts, weighted, squares = _histogram_flight(
        _velocity(0.5),
        _velocity(3.0),
        duration_s=duration,
        weight=weight,
        edges_eV=edges,
    )
    expected = weight * duration * np.asarray([0.2, 0.4, 0.4, 0.0])
    assert weighted == pytest.approx(expected, rel=2.0e-14, abs=1.0e-15)
    assert float(np.sum(weighted)) == pytest.approx(weight * duration, rel=1.0e-15)
    assert counts.tolist() == [1, 1, 1, 0]
    assert squares == pytest.approx(expected * expected)

    # A decelerating/reaccelerating flight visits each occupied bin twice.
    # Its contributions remain one correlated flight for effective counts.
    counts, weighted, squares = _histogram_flight(
        _velocity(2.0),
        _velocity(-2.0),
        duration_s=duration,
        weight=weight,
        edges_eV=edges,
    )
    expected = weight * duration * np.asarray([0.5, 0.5, 0.0, 0.0])
    assert weighted == pytest.approx(expected, rel=2.0e-14, abs=1.0e-15)
    assert float(np.sum(weighted)) == pytest.approx(weight * duration, rel=1.0e-15)
    assert counts.tolist() == [1, 1, 0, 0]
    assert squares == pytest.approx(expected * expected)

    counts, weighted, _ = _histogram_flight(
        _velocity(4.0),
        _velocity(4.0),
        duration_s=duration,
        weight=weight,
        edges_eV=edges,
    )
    assert counts.tolist() == [0, 0, 0, 1]
    assert weighted.tolist() == pytest.approx([0.0, 0.0, 0.0, weight * duration])


def _piecewise_linear_speed_integral(
    energy_eV: np.ndarray,
    sigma_m2: np.ndarray,
    *,
    start_sqrt_energy: float,
    stop_sqrt_energy: float,
    duration_s: float,
) -> float:
    total = 0.0
    for index in range(energy_eV.size - 1):
        left = max(start_sqrt_energy, math.sqrt(float(energy_eV[index])))
        right = min(stop_sqrt_energy, math.sqrt(float(energy_eV[index + 1])))
        if right <= left:
            continue
        slope = (sigma_m2[index + 1] - sigma_m2[index]) / (
            energy_eV[index + 1] - energy_eV[index]
        )
        intercept = sigma_m2[index] - slope * energy_eV[index]
        total += 0.5 * intercept * (right * right - left * left)
        total += 0.25 * slope * (right**4 - left**4)
    return (
        duration_s
        * _SPEED_PER_SQRT_EV
        * total
        / (stop_sqrt_energy - start_sqrt_energy)
    )


def test_piecewise_xs_and_threshold_partition_gives_analytic_hazard() -> None:
    energy = np.asarray([0.0, 1.0, 4.0, 9.0])
    sigma = 1.0e-20 * np.asarray([0.5, 2.0, 1.0, 3.0])
    threshold = 2.25
    process = CrossSectionProcess(
        species="X",
        process="piecewise",
        process_type=ProcessType.EXCITATION,
        energy_eV=energy,
        cross_section_m2=sigma,
        threshold_eV=threshold,
    )
    table = PreparedCrossSectionTable.build((process,))
    breakpoints = np.unique(np.concatenate((table.grid_eV, [threshold])))
    rate, loss, quadrature, intervals, depths = new_rate_integration_workspace(1)
    duration = 7.0e-9

    converged = integrate_dc_flight_rates(
        _velocity(0.5),
        _velocity(3.0),
        duration,
        breakpoints,
        table.grid_eV,
        table.values_m2,
        table.slopes_m2_eV,
        table.process_min_eV,
        table.process_max_eV,
        table.right_values_m2,
        np.asarray([threshold]),
        np.asarray([np.nan]),
        rate,
        loss,
        quadrature,
        intervals,
        depths,
    )
    expected = _piecewise_linear_speed_integral(
        energy,
        sigma,
        start_sqrt_energy=math.sqrt(threshold),
        stop_sqrt_energy=3.0,
        duration_s=duration,
    )

    assert converged
    assert table.evaluate(threshold - 1.0e-6)[0] == 0.0
    assert table.evaluate(threshold)[0] > 0.0
    assert rate[0] == pytest.approx(expected, rel=2.0e-13)
    assert loss[0] == pytest.approx(threshold * expected, rel=2.0e-13)


def test_zero_at_xs_knot_microinterval_converges_in_local_coordinates() -> None:
    energy_offset = 1.0e-8
    duration = 1.0e-9
    sigma_slope = 1.0e-20

    # These include the two knots exposed by the representative Ar/N2 run.
    # The construction is generic: a channel opens linearly at either knot.
    for knot in (13.273, 14.0):
        process = CrossSectionProcess(
            species="X",
            process="zero-at-knot",
            process_type=ProcessType.EXCITATION,
            energy_eV=np.asarray([knot - 1.0, knot, knot + 1.0]),
            cross_section_m2=sigma_slope * np.asarray([0.0, 0.0, 1.0]),
            threshold_eV=knot,
        )
        table = PreparedCrossSectionTable.build((process,))
        workspace = new_rate_integration_workspace(1)
        before_energy = knot - energy_offset
        after_energy = knot + energy_offset

        converged = integrate_dc_flight_rates(
            _velocity(math.sqrt(before_energy)),
            _velocity(math.sqrt(after_energy)),
            duration,
            table.grid_eV,
            table.grid_eV,
            table.values_m2,
            table.slopes_m2_eV,
            table.process_min_eV,
            table.process_max_eV,
            table.right_values_m2,
            np.asarray([knot]),
            np.asarray([np.nan]),
            *workspace,
        )
        sqrt_before = math.sqrt(before_energy)
        sqrt_after = math.sqrt(after_energy)
        sqrt_span = (after_energy - before_energy) / (
            sqrt_after + sqrt_before
        )
        expected_rate = (
            duration
            * _SPEED_PER_SQRT_EV
            * sigma_slope
            * (after_energy - knot) ** 2
            / (4.0 * sqrt_span)
        )

        assert converged
        assert workspace[0][0] == pytest.approx(expected_rate, rel=5.0e-7)
        assert workspace[1][0] == pytest.approx(knot * expected_rate, rel=5.0e-7)


@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="Numba is not installed")
def test_python_and_numba_flight_integrators_are_numerically_identical() -> None:
    energy = np.asarray([0.0, 0.7, 2.2, 5.0, 12.0])
    first = CrossSectionProcess(
        species="X",
        process="first",
        process_type=ProcessType.EXCITATION,
        energy_eV=energy,
        cross_section_m2=1.0e-20 * np.asarray([0.0, 0.4, 2.0, 1.2, 3.0]),
        threshold_eV=0.9,
    )
    second = CrossSectionProcess(
        species="X",
        process="second",
        process_type=ProcessType.IONIZATION,
        energy_eV=np.asarray([0.0, 1.3, 4.4, 12.0]),
        cross_section_m2=1.0e-21 * np.asarray([0.0, 0.0, 3.0, 5.0]),
        threshold_eV=1.3,
    )
    table = PreparedCrossSectionTable.build((first, second))
    breakpoints = np.unique(
        np.concatenate((table.grid_eV, [first.threshold_eV, second.threshold_eV]))
    )
    before = np.asarray([1.1e6, -0.7e6, 0.2e6])
    after = np.asarray([-0.6e6, 1.4e6, -0.3e6])
    arguments = (
        before,
        after,
        3.0e-9,
        breakpoints,
        table.grid_eV,
        table.values_m2,
        table.slopes_m2_eV,
        table.process_min_eV,
        table.process_max_eV,
        table.right_values_m2,
        np.asarray([first.threshold_eV, second.threshold_eV]),
        np.asarray([np.nan, np.nan]),
    )
    compiled_workspace = new_rate_integration_workspace(2)
    python_workspace = new_rate_integration_workspace(2)

    assert integrate_dc_flight_rates(*arguments, *compiled_workspace)
    assert integrate_dc_flight_rates.py_func(*arguments, *python_workspace)
    assert compiled_workspace[0] == pytest.approx(
        python_workspace[0], rel=2.0e-15, abs=0.0
    )
    assert compiled_workspace[1] == pytest.approx(
        python_workspace[1], rel=2.0e-15, abs=0.0
    )

    edges = np.linspace(0.0, 20.0, 81)
    compiled_histogram = _histogram_flight(
        before,
        after,
        duration_s=3.0e-9,
        weight=2.5,
        edges_eV=edges,
    )
    bins = edges.size - 1
    python_histogram = (
        np.zeros(bins, dtype=np.int64),
        np.zeros(bins),
        np.zeros(bins),
    )
    accumulate_dc_flight_histogram.py_func(
        before,
        after,
        3.0e-9,
        2.5,
        edges,
        *python_histogram,
        np.zeros(bins),
        np.empty(bins, dtype=np.int64),
    )
    for compiled, python in zip(compiled_histogram, python_histogram, strict=True):
        assert compiled == pytest.approx(python, rel=2.0e-15, abs=0.0)
