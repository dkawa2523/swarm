"""Canonical electron transport result container."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ElectronTransport:
    definition: str
    gas_number_density_m3: float
    drift_velocity_m_s: float
    reduced_mobility_m2_V_s_m3: float
    reduced_diffusion_L_m2_s_m3: float | None
    reduced_diffusion_T_m2_s_m3: float | None
    reduced_electron_energy_mobility_m2_V_s_m3: float | None = None
    reduced_electron_energy_diffusion_m2_s_m3: float | None = None
    reduced_electron_energy_diffusion_L_m2_s_m3: float | None = None
    reduced_electron_energy_diffusion_T_m2_s_m3: float | None = None

    def __post_init__(self) -> None:
        scalar = self.reduced_electron_energy_diffusion_m2_s_m3
        if scalar is None:
            return
        if self.reduced_electron_energy_diffusion_L_m2_s_m3 is None:
            object.__setattr__(
                self,
                "reduced_electron_energy_diffusion_L_m2_s_m3",
                scalar,
            )
        if self.reduced_electron_energy_diffusion_T_m2_s_m3 is None:
            object.__setattr__(
                self,
                "reduced_electron_energy_diffusion_T_m2_s_m3",
                scalar,
            )

    @classmethod
    def from_actual(
        cls,
        *,
        definition: str,
        gas_number_density_m3: float,
        drift_velocity_m_s: float,
        mobility_m2_V_s: float,
        diffusion_L_m2_s: float | None,
        diffusion_T_m2_s: float | None,
        electron_energy_mobility_m2_V_s: float | None = None,
        electron_energy_diffusion_m2_s: float | None = None,
        electron_energy_diffusion_L_m2_s: float | None = None,
        electron_energy_diffusion_T_m2_s: float | None = None,
    ) -> "ElectronTransport":
        density = float(gas_number_density_m3)
        return cls(
            definition=definition,
            gas_number_density_m3=density,
            drift_velocity_m_s=float(drift_velocity_m_s),
            reduced_mobility_m2_V_s_m3=float(mobility_m2_V_s) * density,
            reduced_diffusion_L_m2_s_m3=(
                None
                if diffusion_L_m2_s is None
                else float(diffusion_L_m2_s) * density
            ),
            reduced_diffusion_T_m2_s_m3=(
                None
                if diffusion_T_m2_s is None
                else float(diffusion_T_m2_s) * density
            ),
            reduced_electron_energy_mobility_m2_V_s_m3=(
                None
                if electron_energy_mobility_m2_V_s is None
                else float(electron_energy_mobility_m2_V_s) * density
            ),
            reduced_electron_energy_diffusion_m2_s_m3=(
                None
                if electron_energy_diffusion_m2_s is None
                else float(electron_energy_diffusion_m2_s) * density
            ),
            reduced_electron_energy_diffusion_L_m2_s_m3=(
                float(
                    electron_energy_diffusion_m2_s
                    if electron_energy_diffusion_L_m2_s is None
                    else electron_energy_diffusion_L_m2_s
                )
                * density
                if (
                    electron_energy_diffusion_L_m2_s is not None
                    or electron_energy_diffusion_m2_s is not None
                )
                else None
            ),
            reduced_electron_energy_diffusion_T_m2_s_m3=(
                float(
                    electron_energy_diffusion_m2_s
                    if electron_energy_diffusion_T_m2_s is None
                    else electron_energy_diffusion_T_m2_s
                )
                * density
                if (
                    electron_energy_diffusion_T_m2_s is not None
                    or electron_energy_diffusion_m2_s is not None
                )
                else None
            ),
        )

    @property
    def mobility_m2_V_s(self) -> float:
        return self.reduced_mobility_m2_V_s_m3 / self.gas_number_density_m3

    @property
    def diffusion_L_m2_s(self) -> float | None:
        if self.reduced_diffusion_L_m2_s_m3 is None:
            return None
        return self.reduced_diffusion_L_m2_s_m3 / self.gas_number_density_m3

    @property
    def diffusion_T_m2_s(self) -> float | None:
        if self.reduced_diffusion_T_m2_s_m3 is None:
            return None
        return self.reduced_diffusion_T_m2_s_m3 / self.gas_number_density_m3

    @property
    def electron_energy_mobility_m2_V_s(self) -> float | None:
        if self.reduced_electron_energy_mobility_m2_V_s_m3 is None:
            return None
        return (
            self.reduced_electron_energy_mobility_m2_V_s_m3
            / self.gas_number_density_m3
        )

    @property
    def electron_energy_diffusion_m2_s(self) -> float | None:
        if self.reduced_electron_energy_diffusion_m2_s_m3 is None:
            return None
        return (
            self.reduced_electron_energy_diffusion_m2_s_m3
            / self.gas_number_density_m3
        )

    @property
    def electron_energy_diffusion_L_m2_s(self) -> float | None:
        value = self.reduced_electron_energy_diffusion_L_m2_s_m3
        if value is None:
            value = self.reduced_electron_energy_diffusion_m2_s_m3
        return None if value is None else value / self.gas_number_density_m3

    @property
    def electron_energy_diffusion_T_m2_s(self) -> float | None:
        value = self.reduced_electron_energy_diffusion_T_m2_s_m3
        if value is None:
            value = self.reduced_electron_energy_diffusion_m2_s_m3
        return None if value is None else value / self.gas_number_density_m3

    @property
    def characteristic_energy_L_eV(self) -> float | None:
        mobility = self.mobility_m2_V_s
        diffusion = self.diffusion_L_m2_s
        if mobility == 0.0 or diffusion is None:
            return None
        return diffusion / abs(mobility)

    @property
    def characteristic_energy_T_eV(self) -> float | None:
        mobility = self.mobility_m2_V_s
        diffusion = self.diffusion_T_m2_s
        if mobility == 0.0 or diffusion is None:
            return None
        return diffusion / abs(mobility)
