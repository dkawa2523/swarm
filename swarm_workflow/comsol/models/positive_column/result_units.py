"""Positive-column field meanings and explicit conversions, without unit inference."""

from dataclasses import dataclass

from electron_swarm.core.constants import E_CHARGE_C


@dataclass(frozen=True, slots=True)
class ResultQuantity:
    unit: str
    label: str
    logarithmic: bool = False
    display_unit: str | None = None
    display_scale: float = 1.0


QUANTITIES = {
    "x": ResultQuantity("m", "Axial position"),
    "electron_density": ResultQuantity("1/m^3", "Electron density", True),
    "mean_electron_energy": ResultQuantity("eV", "Mean electron energy"),
    "electric_potential": ResultQuantity("V", "Electric potential"),
    "E_over_N": ResultQuantity("V*m^2", "Reduced electric field", True, "Td", 1e21),
    "electron_current_density": ResultQuantity("A/m^2", "Electron current density"),
    "ion_current_density": ResultQuantity("A/m^2", "Ion current density"),
    "total_current_density": ResultQuantity(
        "A/m^2", "Electron + Ar-ion conductive current density"
    ),
    "excitation_source": ResultQuantity("1/(m^3*s)", "Excitation source", True),
    "ionization_source": ResultQuantity("1/(m^3*s)", "Ionization source", True),
    "reduced_mobility": ResultQuantity("1/(V*m*s)", "Reduced mobility"),
    "reduced_diffusion_L": ResultQuantity("1/(m*s)", "Reduced longitudinal diffusion"),
    "excitation_townsend": ResultQuantity(
        "m^2", "Reduced excitation Townsend coefficient"
    ),
    "ionization_townsend": ResultQuantity(
        "m^2", "Reduced ionization Townsend coefficient"
    ),
    "applied_voltage": ResultQuantity("V", "Applied voltage"),
    "gas_pressure": ResultQuantity("Pa", "Gas pressure"),
}

_SCALES = {
    "m": {"m": 1.0, "cm": 0.01, "mm": 0.001, "um": 1e-6},
    "1/m^3": {"1/m^3": 1.0, "m^-3": 1.0, "1/cm^3": 1e6, "cm^-3": 1e6},
    "eV": {"eV": 1.0, "J": 1.0 / E_CHARGE_C},
    "V": {"V": 1.0, "kV": 1e3, "mV": 1e-3},
    "V*m^2": {"V*m^2": 1.0, "Td": 1e-21},
    "A/m^2": {"A/m^2": 1.0, "A/cm^2": 1e4, "mA/m^2": 1e-3, "mA/cm^2": 10.0},
    "1/(m^3*s)": {"1/(m^3*s)": 1.0, "1/(cm^3*s)": 1e6},
    "1/(V*m*s)": {"1/(V*m*s)": 1.0, "1/(V*cm*s)": 100.0},
    "1/(m*s)": {"1/(m*s)": 1.0, "1/(cm*s)": 100.0},
    "m^2": {"m^2": 1.0, "cm^2": 1e-4},
    "Pa": {"Pa": 1.0, "kPa": 1e3, "hPa": 100.0},
}


def canonical_unit(name: str, unit: str) -> tuple[str, float]:
    """Return an explicit scale; custom field names retain their declared unit."""
    if not isinstance(unit, str) or not unit.strip():
        raise ValueError(f"result field {name!r} requires a nonempty unit")
    if name not in QUANTITIES:
        return unit, 1.0
    target = QUANTITIES[name].unit
    normalized = "".join(unit.split())
    scales = _SCALES[target]
    if normalized not in scales:
        raise ValueError(
            f"result field {name!r}: unit {unit!r} is unsupported; choose {list(scales)}"
        )
    return target, scales[normalized]
