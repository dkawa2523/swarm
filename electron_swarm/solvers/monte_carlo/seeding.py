"""Deterministic random-stream identities for Monte Carlo field cases."""

from __future__ import annotations

from hashlib import sha256

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence


def derive_case_seed(*, base_seed: int, e_over_n_Td: float) -> int:
    """Derive a stable uint32 seed without depending on batch order."""

    material = (
        f"{_mc_evidence.MC_CASE_SEED_DERIVATION}\0monte_carlo\0{int(base_seed)}\0"
        f"{float(e_over_n_Td).hex()}"
    ).encode("ascii")
    return int.from_bytes(sha256(material).digest()[:4], "big", signed=False)
