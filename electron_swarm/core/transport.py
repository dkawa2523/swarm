"""Canonical electron transport result container."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ElectronTransport:
    definition: str
    gas_number_density_m3: float
    drift_velocity_m_s: float
    reduced_mobility_m2_V_s_m3: float
    reduced_diffusion_L_m2_s_m3: float
    reduced_diffusion_T_m2_s_m3: float
    reduced_electron_energy_mobility_eV_m2_V_s_m3: float | None = None
    reduced_electron_energy_diffusion_eV_m2_s_m3: float | None = None

    @classmethod
    def from_actual(
        cls,
        *,
        definition: str,
        gas_number_density_m3: float,
        drift_velocity_m_s: float,
        mobility_m2_V_s: float,
        diffusion_L_m2_s: float,
        diffusion_T_m2_s: float,
        electron_energy_mobility_eV_m2_V_s: float | None = None,
        electron_energy_diffusion_eV_m2_s: float | None = None,
    ) -> "ElectronTransport":
        density = float(gas_number_density_m3)
        return cls(
            definition=definition,
            gas_number_density_m3=density,
            drift_velocity_m_s=float(drift_velocity_m_s),
            reduced_mobility_m2_V_s_m3=float(mobility_m2_V_s) * density,
            reduced_diffusion_L_m2_s_m3=float(diffusion_L_m2_s) * density,
            reduced_diffusion_T_m2_s_m3=float(diffusion_T_m2_s) * density,
            reduced_electron_energy_mobility_eV_m2_V_s_m3=(
                None
                if electron_energy_mobility_eV_m2_V_s is None
                else float(electron_energy_mobility_eV_m2_V_s) * density
            ),
            reduced_electron_energy_diffusion_eV_m2_s_m3=(
                None
                if electron_energy_diffusion_eV_m2_s is None
                else float(electron_energy_diffusion_eV_m2_s) * density
            ),
        )

    @property
    def mobility_m2_V_s(self) -> float:
        return self.reduced_mobility_m2_V_s_m3 / self.gas_number_density_m3

    @property
    def diffusion_L_m2_s(self) -> float:
        return self.reduced_diffusion_L_m2_s_m3 / self.gas_number_density_m3

    @property
    def diffusion_T_m2_s(self) -> float:
        return self.reduced_diffusion_T_m2_s_m3 / self.gas_number_density_m3

    @property
    def electron_energy_mobility_eV_m2_V_s(self) -> float | None:
        if self.reduced_electron_energy_mobility_eV_m2_V_s_m3 is None:
            return None
        return (
            self.reduced_electron_energy_mobility_eV_m2_V_s_m3
            / self.gas_number_density_m3
        )

    @property
    def electron_energy_diffusion_eV_m2_s(self) -> float | None:
        if self.reduced_electron_energy_diffusion_eV_m2_s_m3 is None:
            return None
        return (
            self.reduced_electron_energy_diffusion_eV_m2_s_m3
            / self.gas_number_density_m3
        )

    @property
    def characteristic_energy_L_eV(self) -> float | None:
        mobility = self.mobility_m2_V_s
        if mobility == 0.0:
            return None
        return self.diffusion_L_m2_s / abs(mobility)

    @property
    def characteristic_energy_T_eV(self) -> float | None:
        mobility = self.mobility_m2_V_s
        if mobility == 0.0:
            return None
        return self.diffusion_T_m2_s / abs(mobility)
