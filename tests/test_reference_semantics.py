from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from tools.benchmarks.benchmark_ar_external_references import _rate_metric
from tools.benchmarks.benchmark_ar_external_references import _mc_confidence_status
from tools.benchmarks.references.common import (
    reference_cases_from_frame,
    reference_comparison_case,
    reference_comparison_metrics,
    reference_comparison_status,
)


def _frame():
    return pd.DataFrame(
        {
            "case_id": ["NA", "NA"],
            "E_over_N_Td": [50.0, 50.0],
            "energy_eV": [0.5, 1.5],
            "eedf_eV_inv": [0.5, 0.5],
        }
    )


def _read(frame):
    return reference_cases_from_frame(
        frame, reference_id="synthetic", convention="eedf"
    )[0]


def test_reference_convention_cannot_be_overridden_by_a_column():
    with pytest.raises(ValueError, match="contradicts reference column"):
        reference_cases_from_frame(_frame(), reference_id="test", convention="eepf")


@pytest.mark.parametrize(
    "change",
    [
        "missing_id",
        "missing_field",
        "bad_scalar",
        "inconsistent_scalar",
        "negative_distribution",
        "negative_uncertainty",
    ],
)
def test_bad_reference_rows_are_not_silently_dropped_or_coerced(change):
    frame = _frame()
    if change == "missing_id":
        frame.loc[0, "case_id"] = None
    elif change == "missing_field":
        frame.loc[0, "E_over_N_Td"] = np.nan
    elif change == "bad_scalar":
        frame["mobility_m2_V_s"] = ["bad", ""]
    elif change == "inconsistent_scalar":
        frame["mobility_m2_V_s"] = [1.0, 2.0]
    elif change == "negative_distribution":
        frame["eedf_eV_inv"] = [-0.1, 1.1]
    else:
        frame["mean_energy_eV_ci95"] = -1.0
    with pytest.raises(ValueError):
        _read(frame)


def test_missing_transport_and_rates_remain_unavailable():
    reference = _read(_frame())
    candidate = reference_comparison_case(reference)
    metrics = reference_comparison_metrics(reference, candidate)
    assert candidate.mobility_m2_V_s is None
    assert "mobility_relative_difference" not in metrics
    assert "drift_velocity_relative_difference" not in metrics
    assert "major_rate_relative_difference" not in metrics
    verdict = reference_comparison_status(
        metrics, {"eedf_relative_l1": 0.1, "major_rate_relative_difference": 0.1}
    )
    assert verdict["status"] == "PARTIAL"


def test_rate_outputs_have_distinct_named_errors():
    frame = _frame()
    frame["rate__ionization"] = 1.0
    frame["rate__excitation"] = 2.0
    reference = _read(frame)
    candidate_reference = deepcopy(reference)
    candidate_reference.rates["ionization"] = 1.5
    metrics = reference_comparison_metrics(
        reference, reference_comparison_case(candidate_reference)
    )
    assert metrics["major_rate_relative_difference"] == 0.5
    assert _rate_metric(metrics, "ionization") == 0.5
    assert _rate_metric(metrics, "excitation") == 0.0


def test_scalar_confidence_cannot_override_failed_distribution_gate():
    verdict = reference_comparison_status(
        {"mean_energy_relative_difference": 0.2, "eedf_relative_l1": 0.5},
        {"mean_energy_relative_difference": 0.1, "eedf_relative_l1": 0.1},
        confidence_status="pass_within_mc_uncertainty",
    )
    assert verdict["status"] == "FAIL"
    assert "eedf_relative_l1" in verdict["failed_metrics"]


def test_declared_unavailable_uncertainty_is_not_used_for_acceptance():
    reference = _read(_frame())
    reference.metadata.update(
        {"uncertainty_unavailable": True, "scalar_ci95": {"mean_energy_eV": 1.0}}
    )
    status, _ = _mc_confidence_status(
        reference, reference_comparison_case(reference), {}
    )
    assert status == "unknown"


def test_benchmark_comparison_preserves_disjoint_probability_support():
    reference = _read(_frame())
    shifted = _frame()
    shifted["energy_eV"] = [2.5, 3.5]
    metrics = reference_comparison_metrics(
        reference, reference_comparison_case(_read(shifted))
    )
    assert metrics["eedf_l1"] == pytest.approx(2.0)
    assert metrics["eedf_mean_energy_relative_difference"] == pytest.approx(2.0)


def test_reference_preserves_reported_cells_instead_of_inferring_new_ones():
    frame = _frame()
    frame["energy_eV"] = [0.2, 1.8]
    frame["energy_width_eV"] = [1.0, 1.0]
    reference = _read(frame)
    assert reference.metadata["distribution_geometry"] == "reported_cells"
    assert reference_comparison_case(reference).energy_widths_eV.tolist() == [1.0, 1.0]
    assert np.dot(reference.eedf_eV_inv, reference.energy_widths_eV) == 1.0
