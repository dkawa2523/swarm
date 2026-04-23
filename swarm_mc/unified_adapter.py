"""Bridge the existing swarm_mc engine to the unified electron_swarm schema."""

from __future__ import annotations

from copy import deepcopy
from time import perf_counter

import numpy as np
import scipy.constants as csts

from electron_swarm.core.results import RateResult, SwarmCaseResult
from swarm_mc.simulation import Simulation


def _fraction_by_species(gases: list[str], fractions: list[float]) -> dict[str, float]:
    mapping: dict[str, float] = {}
    for gas, fraction in zip(gases, fractions, strict=False):
        mapping[gas] = mapping.get(gas, 0.0) + float(fraction)
    return mapping


def _as_float(value, default=np.nan) -> float:
    if value is None:
        return float(default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _infer_pressure_pa(config) -> float:
    if config.conditions.pressure_Pa is not None:
        return float(config.conditions.pressure_Pa)
    assert config.conditions.gas_number_density_m3 is not None
    return float(
        config.conditions.gas_number_density_m3
        * csts.Boltzmann
        * config.conditions.gas_temperature_K
    )


def _infer_max_cross_section_energy(cross_sections) -> float:
    return float(
        max(float(np.max(process.energy_eV)) for process in cross_sections.processes)
    )


def _effective_frequency(rate_coefficient_m3_s: float, gas_number_density: float) -> float:
    if np.isfinite(rate_coefficient_m3_s):
        return float(rate_coefficient_m3_s * gas_number_density)
    return float(np.nan)


def _effective_townsend(
    frequency_s: float, drift_velocity_m_s: float, gas_number_density: float
) -> tuple[float, float]:
    if not np.isfinite(frequency_s) or drift_velocity_m_s == 0.0:
        return float(np.nan), float(np.nan)
    return (
        float(frequency_s / (abs(drift_velocity_m_s) * gas_number_density)),
        float(frequency_s / abs(drift_velocity_m_s)),
    )


def _build_mc_config(config, cross_sections, *, case_id: str, e_over_n_Td: float, base_config: dict) -> dict:
    cfg = deepcopy(base_config)
    gases = cfg.setdefault("input_gases", {})
    gases.setdefault(
        "gases", [component.species for component in config.conditions.gas_mixture]
    )
    gases.setdefault(
        "paths_to_cross_section_files",
        [str(file_cfg.path) for file_cfg in config.cross_sections.files],
    )
    gases.setdefault(
        "fractions",
        [float(component.fraction) for component in config.conditions.gas_mixture],
    )
    gases.setdefault(
        "max_cross_section_energy", _infer_max_cross_section_energy(cross_sections)
    )

    output = cfg.setdefault("output", {})
    output_dir = config.output.directory / "_mc_runs" / case_id
    output.setdefault("output_directory", str(output_dir))
    output["base_name"] = case_id
    output.setdefault("save_simulation_pickle", False)
    output.setdefault("save_temporal_evolution", False)
    output.setdefault("save_swarm_parameters", False)
    output.setdefault("save_energy_distribution", False)
    output.setdefault("output_format", "csv")

    physical = cfg.setdefault("physical_conditions", {})
    physical["EN"] = float(e_over_n_Td)
    physical["pressure"] = _infer_pressure_pa(config)
    physical["temperature"] = float(config.conditions.gas_temperature_K)

    initial = cfg.setdefault("initial_state", {})
    initial.setdefault("num_e_initial", 5e3)
    initial.setdefault("initial_pos_electrons", [0.0, 0.0, 0.0])
    initial.setdefault("initial_std_electrons", [0.0, 0.0, 0.0])

    simulation = cfg.setdefault("simulation_settings", {})
    simulation.setdefault("num_energy_bins", 768)
    simulation.setdefault("energy_sharing_factor", 0.5)
    simulation.setdefault("isotropic_scattering", False)
    simulation.setdefault("use_aniso_scatter_lut", True)
    simulation.setdefault("conserve", False)
    simulation.setdefault("num_e_max", 1e6)
    simulation.setdefault("seed", "run_id")

    end = cfg.setdefault("end_conditions", {})
    end.setdefault("w_tol", 0.03)
    end.setdefault("DN_tol", 0.03)
    end.setdefault("num_col_max", 1e5)

    return cfg


def _build_rate_results(sim: Simulation, case_id: str, e_over_n_Td: float) -> list[RateResult]:
    per_reaction = list(getattr(sim.output.rates_conv, "per_reaction", []) or [])
    if not per_reaction:
        return []

    meta_lookup = {meta.label: meta for meta in sim.gas_mixture.cross_section_meta}
    fractions = _fraction_by_species(sim.config.gases, sim.config.fractions)
    thresholds = np.asarray(getattr(sim.gas_mixture, "thresholds", []), dtype=float)
    results: list[RateResult] = []
    for reaction in per_reaction:
        meta = meta_lookup.get(reaction.label)
        species = meta.species if meta is not None else "unknown"
        threshold = None
        if meta is not None and meta.index < thresholds.size:
            threshold = float(thresholds[meta.index])
        mixture_weighted = float(reaction.rate)
        fraction = float(fractions.get(species, 1.0))
        rate = mixture_weighted / fraction if fraction > 0.0 else mixture_weighted
        results.append(
            RateResult(
                solver="monte_carlo",
                case_id=case_id,
                e_over_n_Td=float(e_over_n_Td),
                species=species,
                process=reaction.label,
                process_type=str(reaction.ctype.name).lower(),
                threshold_eV=threshold,
                rate_coefficient_m3_s=rate,
                mixture_weighted_rate_m3_s=mixture_weighted,
            )
        )
    return results


def _fallback_distribution(sim: Simulation) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    energies = np.asarray(sim.electrons.energy, dtype=float)
    if energies.size == 0:
        energy = np.array([0.0, 1.0], dtype=float)
        eedf = np.array([1.0, 0.0], dtype=float)
        eepf = np.array([1.0e15, 0.0], dtype=float)
        return energy, eedf, eepf

    upper = max(float(np.max(energies)), 1.0)
    bins = np.linspace(0.0, upper, 128)
    hist, edges = np.histogram(energies, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist = hist / max(np.trapezoid(hist, centers), 1.0e-30)
    eepf = hist / np.sqrt(np.maximum(centers, 1.0e-30))
    return centers, hist, eepf


def _build_case_result(sim: Simulation, case_id: str, e_over_n_Td: float, run_time_s: float) -> SwarmCaseResult:
    summary = sim.output.swarm_parameters_dict() if sim.output else {}
    energy = getattr(sim.output.energy_distribution, "energy_bin_centers", None)
    eedf = getattr(sim.output.energy_distribution, "eedf", None)
    eepf = getattr(sim.output.energy_distribution, "eepf", None)
    if energy is None or eedf is None or eepf is None:
        energy, eedf, eepf = _fallback_distribution(sim)

    gas_number_density = float(sim.config.gas_number_density)
    electric_field = float(sim.electric_field)
    flux_drift = _as_float(summary.get("flux drift velocity (m.s-1)"), float(sim.output.flux.w[2]))
    bulk_drift = _as_float(summary.get("bulk drift velocity (m.s-1)"), float(sim.output.bulk.w[2]))
    reduced_flux_diff_l = _as_float(
        summary.get("flux L diffusion coeff. * N (m-1.s-1)"),
        float(sim.output.flux.DN[2]),
    )
    reduced_flux_diff_t = _as_float(
        summary.get("flux T diffusion coeff. * N (m-1.s-1)"),
        float(sim.output.flux.DN[0]),
    )
    reduced_bulk_diff_l = _as_float(
        summary.get("bulk L diffusion coeff. * N (m-1.s-1)"),
        float(sim.output.bulk.DN[2]),
    )
    reduced_bulk_diff_t = _as_float(
        summary.get("bulk T diffusion coeff. * N (m-1.s-1)"),
        float(sim.output.bulk.DN[0]),
    )
    mobility = flux_drift / electric_field if electric_field != 0.0 else np.nan
    reduced_mobility = mobility * gas_number_density if np.isfinite(mobility) else np.nan
    diff_l = reduced_flux_diff_l / gas_number_density
    diff_t = reduced_flux_diff_t / gas_number_density
    counted_effective_rate = _as_float(sim.output.rates_count.effective, np.nan)
    convolution_effective_rate = _as_float(sim.output.rates_conv.effective, np.nan)
    counted_nu_eff = _effective_frequency(counted_effective_rate, gas_number_density)
    convolution_nu_eff = _effective_frequency(
        convolution_effective_rate, gas_number_density
    )
    primary_nu_eff = (
        convolution_nu_eff
        if np.isfinite(convolution_nu_eff)
        else counted_nu_eff
    )
    reduced_townsend, townsend_1_m = _effective_townsend(
        primary_nu_eff, flux_drift, gas_number_density
    )

    metadata = {
        "run_label": case_id,
        "sweep_param": "E_over_N_Td",
        "sweep_value": float(e_over_n_Td),
        "run_time_s": float(run_time_s),
        "seed_used": sim.seed_used,
        "seed_source": sim.seed_source,
        "gas_number_density_m-3": gas_number_density,
        "electric_field_V_m": electric_field,
        "mean_energy_error_eV": _as_float(summary.get("mean energy error (eV)")),
        "bulk_drift_velocity_m_s": bulk_drift,
        "bulk_drift_velocity_error_m_s": _as_float(
            summary.get("bulk drift velocity error (m.s-1)")
        ),
        "bulk_reduced_diffusion_L_m-1_s-1": reduced_bulk_diff_l,
        "bulk_reduced_diffusion_L_error_m-1_s-1": _as_float(
            summary.get("bulk L diffusion coeff. error * N (m-1.s-1)")
        ),
        "bulk_reduced_diffusion_T_m-1_s-1": reduced_bulk_diff_t,
        "bulk_reduced_diffusion_T_error_m-1_s-1": _as_float(
            summary.get("bulk T diffusion coeff. error * N (m-1.s-1)")
        ),
        "flux_drift_velocity_error_m_s": _as_float(
            summary.get("flux drift velocity error (m.s-1)")
        ),
        "flux_reduced_diffusion_L_error_m-1_s-1": _as_float(
            summary.get("flux L diffusion coeff. error * N (m-1.s-1)")
        ),
        "flux_reduced_diffusion_T_error_m-1_s-1": _as_float(
            summary.get("flux T diffusion coeff. error * N (m-1.s-1)")
        ),
        "counted_effective_rate_coefficient_m3_s": counted_effective_rate,
        "counted_effective_rate_error_m3_s": _as_float(
            sim.output.rates_count.effective_err
        ),
        "counted_ionization_rate_coefficient_m3_s": _as_float(
            sim.output.rates_count.ionization
        ),
        "counted_ionization_rate_error_m3_s": _as_float(
            sim.output.rates_count.ionization_err
        ),
        "counted_attachment_rate_coefficient_m3_s": _as_float(
            sim.output.rates_count.attachment
        ),
        "counted_attachment_rate_error_m3_s": _as_float(
            sim.output.rates_count.attachment_err
        ),
        "convolution_effective_rate_coefficient_m3_s": convolution_effective_rate,
        "convolution_ionization_rate_coefficient_m3_s": _as_float(
            sim.output.rates_conv.ionization
        ),
        "convolution_attachment_rate_coefficient_m3_s": _as_float(
            sim.output.rates_conv.attachment
        ),
        "counted_net_ionization_frequency_s-1": counted_nu_eff,
        "convolution_net_ionization_frequency_s-1": convolution_nu_eff,
        "net_ionization_frequency_model": (
            "convolution_effective_rate"
            if np.isfinite(convolution_nu_eff)
            else "counted_effective_rate"
        ),
        "effective_townsend_1_m": townsend_1_m,
    }
    return SwarmCaseResult(
        solver="monte_carlo",
        case_id=case_id,
        e_over_n_Td=float(e_over_n_Td),
        mean_energy_eV=_as_float(summary.get("mean energy (eV)"), float(sim.electrons.mean_energy)),
        drift_velocity_m_s=flux_drift,
        mobility_m2_V_s=mobility,
        reduced_mobility_m2_V_s_m3=reduced_mobility,
        diffusion_L_m2_s=diff_l,
        diffusion_T_m2_s=diff_t,
        reduced_diffusion_L_m2_s_m3=reduced_flux_diff_l,
        reduced_diffusion_T_m2_s_m3=reduced_flux_diff_t,
        net_ionization_frequency_s=primary_nu_eff,
        effective_townsend_m2=reduced_townsend,
        energy_eV=np.asarray(energy, dtype=float),
        eedf=np.asarray(eedf, dtype=float),
        eepf=np.asarray(eepf, dtype=float),
        rates=_build_rate_results(sim, case_id, e_over_n_Td),
        metadata=metadata,
    )


def run_swarm(config, cross_sections, **kwargs) -> list[SwarmCaseResult]:
    """Run the current swarm_mc solver and return unified result objects."""

    base_config = deepcopy(kwargs.get("base_config", {}))
    results: list[SwarmCaseResult] = []
    for index, e_over_n_Td in enumerate(config.run.e_over_n_Td):
        case_id = f"{config.run.case_prefix}_{index:04d}"
        mc_config = _build_mc_config(
            config,
            cross_sections,
            case_id=case_id,
            e_over_n_Td=float(e_over_n_Td),
            base_config=base_config,
        )
        started = perf_counter()
        simulation = Simulation(mc_config, gas_mixture_cache=None, run_label=case_id)
        simulation.run()
        run_time_s = perf_counter() - started
        results.append(
            _build_case_result(simulation, case_id, float(e_over_n_Td), run_time_s)
        )
    return results
