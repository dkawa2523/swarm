"""Shared result models for Monte Carlo and Boltzmann solvers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from electron_swarm.core.numerics import eepf_from_eedf
from electron_swarm.core.transport import ElectronTransport


def _finite_scalar(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(slots=True)
class RateResult:
    """Rate coefficient for one collision/reaction process."""

    solver: str
    case_id: str
    e_over_n_Td: float
    species: str
    process: str
    process_type: str
    threshold_eV: float | None
    rate_coefficient_m3_s: float
    target_species_fraction: float = 1.0
    energy_loss_eV: float | None = None
    gas_number_density_m3: float | None = None
    tail_fraction: float | None = None

    def __post_init__(self) -> None:
        if self.energy_loss_eV is None:
            process_type = self.process_type.lower()
            threshold = float(self.threshold_eV or 0.0)
            if process_type in {"excitation", "ionization"}:
                self.energy_loss_eV = threshold
            elif process_type == "superelastic":
                self.energy_loss_eV = -abs(threshold)
            else:
                self.energy_loss_eV = 0.0

    @property
    def mixture_weighted_rate_m3_s(self) -> float:
        return self.rate_coefficient_m3_s * self.target_species_fraction

    @property
    def energy_loss_rate_coefficient_eV_m3_s(self) -> float:
        return self.rate_coefficient_m3_s * float(self.energy_loss_eV or 0.0)

    @property
    def frequency_s_inv(self) -> float | None:
        if self.gas_number_density_m3 is None:
            return None
        return self.gas_number_density_m3 * self.mixture_weighted_rate_m3_s

    @property
    def power_loss_eV_s(self) -> float | None:
        if self.gas_number_density_m3 is None:
            return None
        return (
            self.gas_number_density_m3
            * self.target_species_fraction
            * self.energy_loss_rate_coefficient_eV_m3_s
        )


@dataclass(slots=True)
class RfPhaseResult:
    """Phase-resolved diagnostics for one periodic homogeneous EEDF case."""

    phase_fraction: np.ndarray = field(repr=False)
    phase_rad: np.ndarray = field(repr=False)
    instantaneous_e_over_n_Td: np.ndarray = field(repr=False)
    mean_energy_eV: np.ndarray = field(repr=False)
    ionization_rate_coefficient_m3_s: np.ndarray = field(repr=False)
    eedf: np.ndarray = field(repr=False)


@dataclass(slots=True)
class EnergyAngleDistributionResult:
    """Normalized energy-angle density integrated over azimuth."""

    energy_eV: np.ndarray = field(repr=False)
    energy_widths_eV: np.ndarray = field(repr=False)
    mu: np.ndarray = field(repr=False)
    mu_widths: np.ndarray = field(repr=False)
    density_eV_inv: np.ndarray = field(repr=False)


@dataclass(slots=True)
class SwarmCaseResult:
    """A single solver result for one E/N case."""

    solver: str
    case_id: str
    e_over_n_Td: float
    mean_energy_eV: float
    net_ionization_frequency_s: float
    effective_townsend_m2: float
    transport: ElectronTransport
    energy_eV: np.ndarray = field(repr=False)
    eedf: np.ndarray = field(repr=False)  # normalized energy distribution, 1/eV
    energy_widths_eV: np.ndarray | None = field(default=None, repr=False)
    eedf_counts: np.ndarray | None = field(default=None, repr=False)
    eedf_effective_counts: np.ndarray | None = field(default=None, repr=False)
    rates: list[RateResult] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    rf_phase: RfPhaseResult | None = field(default=None, repr=False)
    energy_angle_distribution: EnergyAngleDistributionResult | None = field(
        default=None,
        repr=False,
    )
    schema_version: str = "2"

    @property
    def drift_velocity_m_s(self) -> float:
        return self.transport.drift_velocity_m_s

    @property
    def mobility_m2_V_s(self) -> float:
        return self.transport.mobility_m2_V_s

    @property
    def reduced_mobility_m2_V_s_m3(self) -> float:
        return self.transport.reduced_mobility_m2_V_s_m3

    @property
    def diffusion_L_m2_s(self) -> float | None:
        return self.transport.diffusion_L_m2_s

    @property
    def diffusion_T_m2_s(self) -> float | None:
        return self.transport.diffusion_T_m2_s

    @property
    def reduced_diffusion_L_m2_s_m3(self) -> float | None:
        return self.transport.reduced_diffusion_L_m2_s_m3

    @property
    def reduced_diffusion_T_m2_s_m3(self) -> float | None:
        return self.transport.reduced_diffusion_T_m2_s_m3

    @property
    def reduced_electron_energy_mobility_m2_V_s_m3(self) -> float | None:
        return self.transport.reduced_electron_energy_mobility_m2_V_s_m3

    @property
    def reduced_electron_energy_diffusion_m2_s_m3(self) -> float | None:
        return self.transport.reduced_electron_energy_diffusion_m2_s_m3

    @property
    def reduced_electron_energy_diffusion_L_m2_s_m3(self) -> float | None:
        return self.transport.reduced_electron_energy_diffusion_L_m2_s_m3

    @property
    def reduced_electron_energy_diffusion_T_m2_s_m3(self) -> float | None:
        return self.transport.reduced_electron_energy_diffusion_T_m2_s_m3

    @property
    def electron_energy_mobility_m2_V_s(self) -> float | None:
        return self.transport.electron_energy_mobility_m2_V_s

    @property
    def electron_energy_diffusion_m2_s(self) -> float | None:
        return self.transport.electron_energy_diffusion_m2_s

    @property
    def electron_energy_diffusion_L_m2_s(self) -> float | None:
        return self.transport.electron_energy_diffusion_L_m2_s

    @property
    def electron_energy_diffusion_T_m2_s(self) -> float | None:
        return self.transport.electron_energy_diffusion_T_m2_s

    @property
    def eepf(self) -> np.ndarray:
        return eepf_from_eedf(self.energy_eV, self.eedf)


@dataclass(slots=True)
class SwarmRunResult:
    """All results from a YAML run."""

    cases: list[SwarmCaseResult]
    metadata: dict[str, Any] = field(default_factory=dict)

    def by_solver(self) -> dict[str, list[SwarmCaseResult]]:
        grouped: dict[str, list[SwarmCaseResult]] = {}
        for case in self.cases:
            grouped.setdefault(case.solver, []).append(case)
        return grouped


def validate_case_result(case: SwarmCaseResult) -> None:
    """Validate the canonical numerical result produced by any solver."""

    if not case.solver or not case.case_id:
        raise ValueError("solver results require solver and case_id")
    if _finite_scalar(case.e_over_n_Td, "E/N") <= 0.0:
        raise ValueError("solver result E/N must be positive")
    if _finite_scalar(case.mean_energy_eV, "mean energy") < 0.0:
        raise ValueError("solver result mean energy must be nonnegative")
    _finite_scalar(case.net_ionization_frequency_s, "net ionization frequency")
    _finite_scalar(case.effective_townsend_m2, "effective Townsend coefficient")

    energy = np.asarray(case.energy_eV, dtype=float)
    eedf = np.asarray(case.eedf, dtype=float)
    widths = (
        None
        if case.energy_widths_eV is None
        else np.asarray(case.energy_widths_eV, dtype=float)
    )
    if energy.ndim != 1 or energy.size == 0:
        raise ValueError("solver result energy grid must be a non-empty 1-D array")
    if eedf.shape != energy.shape:
        raise ValueError("solver result EEDF shape must match the energy grid")
    if widths is None or widths.shape != energy.shape:
        raise ValueError("solver result energy widths must match the energy grid")
    if (
        np.any(~np.isfinite(energy))
        or np.any(energy < 0.0)
        or np.any(np.diff(energy) <= 0.0)
    ):
        raise ValueError("solver result energy grid must be finite and increasing")
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
        raise ValueError("solver result energy widths must be finite and positive")
    if np.any(~np.isfinite(eedf)):
        raise ValueError("solver result EEDF must be finite")
    # Direct PN discretizations can leave a very small negative cell mass even
    # after convergence.  Judge that defect in the conserved measure rather
    # than by a grid-dependent point value; the multi-term core uses the same
    # physical acceptance bound and reports the value in its diagnostics.
    negative_mass = float(np.sum(np.clip(-eedf, 0.0, None) * widths))
    if negative_mass > 1.0e-7:
        raise ValueError("solver result EEDF negative mass exceeds tolerance")
    normalization = float(np.sum(eedf * widths))
    if not np.isclose(normalization, 1.0, rtol=1.0e-6, atol=1.0e-8):
        raise ValueError(
            "solver result EEDF must be normalized on its reported cell widths"
        )

    for name, values in (
        ("EEDF counts", case.eedf_counts),
        ("EEDF effective counts", case.eedf_effective_counts),
    ):
        if values is None:
            continue
        array = np.asarray(values)
        if array.shape != energy.shape:
            raise ValueError(f"solver result {name} shape must match the energy grid")
        if np.any(~np.isfinite(array)) or np.any(array < 0.0):
            raise ValueError(f"solver result {name} must be finite and nonnegative")

    transport = case.transport
    if _finite_scalar(transport.gas_number_density_m3, "gas number density") <= 0.0:
        raise ValueError("solver result gas number density must be positive")
    _finite_scalar(transport.drift_velocity_m_s, "drift velocity")
    _finite_scalar(transport.reduced_mobility_m2_V_s_m3, "reduced mobility")
    for name in (
        "reduced_diffusion_L_m2_s_m3",
        "reduced_diffusion_T_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    ):
        value = getattr(transport, name)
        if value is not None:
            _finite_scalar(value, name)

    for rate in case.rates:
        if rate.solver != case.solver or rate.case_id != case.case_id:
            raise ValueError("rate identity must match its solver result")
        if not np.isclose(float(rate.e_over_n_Td), float(case.e_over_n_Td)):
            raise ValueError("rate E/N must match its solver result")
        _finite_scalar(rate.rate_coefficient_m3_s, "rate coefficient")
        fraction = _finite_scalar(rate.target_species_fraction, "species fraction")
        if not 0.0 <= fraction <= 1.0:
            raise ValueError("rate species fraction must be in [0, 1]")

    phase = case.rf_phase
    if phase is not None:
        phase_count = len(phase.phase_fraction)
        if phase_count == 0 or any(
            len(values) != phase_count
            for values in (
                phase.phase_rad,
                phase.instantaneous_e_over_n_Td,
                phase.mean_energy_eV,
                phase.ionization_rate_coefficient_m3_s,
            )
        ):
            raise ValueError("RF phase arrays must have one common nonzero length")
        if phase.eedf.shape != (phase_count, energy.size):
            raise ValueError("RF phase EEDF shape must match phase and energy grids")

    angular = case.energy_angle_distribution
    if angular is not None:
        expected = (len(angular.energy_eV), len(angular.mu))
        if angular.density_eV_inv.shape != expected:
            raise ValueError("energy-angle density shape must match its coordinates")
        if len(angular.energy_widths_eV) != expected[0] or len(angular.mu_widths) != expected[1]:
            raise ValueError("energy-angle widths must match their coordinates")
