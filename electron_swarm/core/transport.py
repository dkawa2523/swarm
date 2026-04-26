"""Explicit flux, bulk, and source-gradient transport containers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


def _mobility_from_drift_and_field(
    drift_velocity_m_s: float, electric_field_V_m: float
) -> float:
    return (
        drift_velocity_m_s / electric_field_V_m
        if electric_field_V_m != 0.0
        else float("nan")
    )


def _characteristic_energy(
    diffusion_m2_s: float | None, mobility_m2_V_s: float
) -> float | None:
    if diffusion_m2_s is None or mobility_m2_V_s == 0.0:
        return None
    return float(diffusion_m2_s) / abs(float(mobility_m2_V_s))


def _transport_values_with_characteristic_energies(
    drift_velocity_m_s: float,
    mobility_m2_V_s: float,
    diffusion_longitudinal_m2_s: float | None,
    diffusion_transverse_m2_s: float | None,
) -> tuple[float, float, float | None, float | None, float | None, float | None]:
    return (
        float(drift_velocity_m_s),
        float(mobility_m2_V_s),
        diffusion_longitudinal_m2_s,
        diffusion_transverse_m2_s,
        _characteristic_energy(diffusion_longitudinal_m2_s, mobility_m2_V_s),
        _characteristic_energy(diffusion_transverse_m2_s, mobility_m2_V_s),
    )


@dataclass(frozen=True, slots=True)
class FluxTransport:
    drift_velocity_m_s: float
    mobility_m2_V_s: float
    diffusion_longitudinal_m2_s: float | None = None
    diffusion_transverse_m2_s: float | None = None
    characteristic_energy_longitudinal_eV: float | None = None
    characteristic_energy_transverse_eV: float | None = None

    @classmethod
    def from_drift_and_field(
        cls,
        drift_velocity_m_s: float,
        electric_field_V_m: float,
        diffusion_longitudinal_m2_s: float | None = None,
        diffusion_transverse_m2_s: float | None = None,
    ) -> "FluxTransport":
        mobility = _mobility_from_drift_and_field(
            drift_velocity_m_s, electric_field_V_m
        )
        return cls.with_characteristic_energies(
            drift_velocity_m_s,
            mobility,
            diffusion_longitudinal_m2_s,
            diffusion_transverse_m2_s,
        )

    @classmethod
    def with_characteristic_energies(
        cls,
        drift_velocity_m_s: float,
        mobility_m2_V_s: float,
        diffusion_longitudinal_m2_s: float | None = None,
        diffusion_transverse_m2_s: float | None = None,
    ) -> "FluxTransport":
        return cls(
            *_transport_values_with_characteristic_energies(
                drift_velocity_m_s,
                mobility_m2_V_s,
                diffusion_longitudinal_m2_s,
                diffusion_transverse_m2_s,
            )
        )


@dataclass(frozen=True, slots=True)
class BulkTransport:
    drift_velocity_m_s: float
    mobility_m2_V_s: float
    diffusion_longitudinal_m2_s: float | None = None
    diffusion_transverse_m2_s: float | None = None
    characteristic_energy_longitudinal_eV: float | None = None
    characteristic_energy_transverse_eV: float | None = None

    @classmethod
    def from_drift_and_field(
        cls,
        drift_velocity_m_s: float,
        electric_field_V_m: float,
        diffusion_longitudinal_m2_s: float | None = None,
        diffusion_transverse_m2_s: float | None = None,
    ) -> "BulkTransport":
        mobility = _mobility_from_drift_and_field(
            drift_velocity_m_s, electric_field_V_m
        )
        return cls.with_characteristic_energies(
            drift_velocity_m_s,
            mobility,
            diffusion_longitudinal_m2_s,
            diffusion_transverse_m2_s,
        )

    @classmethod
    def with_characteristic_energies(
        cls,
        drift_velocity_m_s: float,
        mobility_m2_V_s: float,
        diffusion_longitudinal_m2_s: float | None = None,
        diffusion_transverse_m2_s: float | None = None,
    ) -> "BulkTransport":
        return cls(
            *_transport_values_with_characteristic_energies(
                drift_velocity_m_s,
                mobility_m2_V_s,
                diffusion_longitudinal_m2_s,
                diffusion_transverse_m2_s,
            )
        )


@dataclass(frozen=True, slots=True)
class SourceGradientTransport:
    ionization_frequency_s_inv: float = 0.0
    attachment_frequency_s_inv: float = 0.0
    effective_growth_frequency_s_inv: float = 0.0
    gradient_velocity_m_s: float | None = None
    curvature_diffusion_longitudinal_m2_s: float | None = None
    curvature_diffusion_transverse_m2_s: float | None = None

    @classmethod
    def from_flux_bulk(
        cls,
        flux: FluxTransport,
        bulk: BulkTransport,
        ionization_frequency_s_inv: float = 0.0,
        attachment_frequency_s_inv: float = 0.0,
    ) -> "SourceGradientTransport":
        dl = None
        if (
            bulk.diffusion_longitudinal_m2_s is not None
            and flux.diffusion_longitudinal_m2_s is not None
        ):
            dl = bulk.diffusion_longitudinal_m2_s - flux.diffusion_longitudinal_m2_s
        dt = None
        if (
            bulk.diffusion_transverse_m2_s is not None
            and flux.diffusion_transverse_m2_s is not None
        ):
            dt = bulk.diffusion_transverse_m2_s - flux.diffusion_transverse_m2_s
        return cls(
            ionization_frequency_s_inv=float(ionization_frequency_s_inv),
            attachment_frequency_s_inv=float(attachment_frequency_s_inv),
            effective_growth_frequency_s_inv=float(
                ionization_frequency_s_inv - attachment_frequency_s_inv
            ),
            gradient_velocity_m_s=bulk.drift_velocity_m_s
            - flux.drift_velocity_m_s,
            curvature_diffusion_longitudinal_m2_s=dl,
            curvature_diffusion_transverse_m2_s=dt,
        )


@dataclass(frozen=True, slots=True)
class TransportMetadata:
    solver: str
    coefficient_definition: Literal["flux", "bulk", "flux_bulk_source"] = (
        "flux_bulk_source"
    )
    swarm_condition: str = "hydrodynamic"
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class TransportSet:
    flux: FluxTransport
    bulk: BulkTransport | None
    source: SourceGradientTransport
    metadata: TransportMetadata | None = None

    def require_bulk(self) -> BulkTransport:
        if self.bulk is None:
            raise RuntimeError(
                "Bulk transport coefficients were not computed. Use a "
                "hydrodynamic or Monte Carlo bulk mode."
            )
        return self.bulk

    @classmethod
    def from_flux_only(
        cls,
        flux: FluxTransport,
        ionization_frequency_s_inv: float = 0.0,
        attachment_frequency_s_inv: float = 0.0,
        metadata: TransportMetadata | None = None,
    ) -> "TransportSet":
        return cls(
            flux=flux,
            bulk=None,
            source=SourceGradientTransport(
                ionization_frequency_s_inv=ionization_frequency_s_inv,
                attachment_frequency_s_inv=attachment_frequency_s_inv,
                effective_growth_frequency_s_inv=ionization_frequency_s_inv
                - attachment_frequency_s_inv,
            ),
            metadata=metadata,
        )

    @classmethod
    def from_flux_bulk(
        cls,
        flux: FluxTransport,
        bulk: BulkTransport,
        ionization_frequency_s_inv: float = 0.0,
        attachment_frequency_s_inv: float = 0.0,
        metadata: TransportMetadata | None = None,
    ) -> "TransportSet":
        return cls(
            flux=flux,
            bulk=bulk,
            source=SourceGradientTransport.from_flux_bulk(
                flux,
                bulk,
                ionization_frequency_s_inv,
                attachment_frequency_s_inv,
            ),
            metadata=metadata,
        )
