from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import glob

_SOLVER_TAGS = {
    "monte_carlo": "mc",
    "boltzmann_two_term": "boltzmann",
    "multiterm_boltzmann": "multiterm",
}


@dataclass
class DiscoveredInputs:
    transport_csv: Optional[Path] = None
    rates_csv: Optional[Path] = None
    eedf_stacked_csv: Optional[Path] = None
    eedf_files_glob: Optional[str] = None


def resolve_rel(run_dir: Path, maybe_rel: Optional[str]) -> Optional[Path]:
    if maybe_rel in (None, "", "null"):
        return None
    p = Path(maybe_rel)
    if not p.is_absolute():
        p = run_dir / p
    return p


def _is_default_name(value: Optional[str], defaults: Iterable[str]) -> bool:
    if value in (None, "", "null"):
        return True
    try:
        name = Path(value).name
    except TypeError:
        return False
    return name in set(defaults)


def resolve_with_fallbacks(
    run_dir: Path,
    maybe_rel: Optional[str],
    fallbacks: Sequence[str],
    default_names: Sequence[str],
) -> Optional[Path]:
    p = resolve_rel(run_dir, maybe_rel)
    if p is not None and p.exists():
        return p
    if not _is_default_name(maybe_rel, default_names):
        return p
    for name in fallbacks:
        candidate = run_dir / name
        if candidate.exists():
            return candidate
    return p


def _primary_solver_name(input_cfg: Dict) -> Optional[str]:
    value = input_cfg.get("primary_solver", input_cfg.get("solver"))
    if value in (None, "", "null"):
        return None
    return str(value)


def _solver_first_fallbacks(primary_solver: Optional[str], prefix: str) -> tuple[str, ...]:
    names = (
        f"{prefix}.csv",
        f"{prefix}_mc.csv",
        f"{prefix}_boltzmann.csv",
        f"{prefix}_multiterm.csv",
    )
    if primary_solver is None:
        return names
    tag = _SOLVER_TAGS.get(primary_solver)
    if tag is None:
        raise ValueError(f"Unknown COMSOL input primary_solver: {primary_solver}")
    preferred = f"{prefix}_{tag}.csv"
    return (preferred,) + tuple(name for name in names if name != preferred)


def _solver_first_table_fallbacks(
    primary_solver: Optional[str], table_prefix: str
) -> tuple[str, ...]:
    names = (
        f"{table_prefix}.csv",
        f"{table_prefix}_mc.csv",
        f"{table_prefix}_boltzmann.csv",
        f"{table_prefix}_multiterm.csv",
    )
    if primary_solver is None:
        return names
    tag = _SOLVER_TAGS.get(primary_solver)
    if tag is None:
        raise ValueError(f"Unknown COMSOL input primary_solver: {primary_solver}")
    preferred = f"{table_prefix}_{tag}.csv"
    return (preferred,) + tuple(name for name in names if name != preferred)


def _ensure_unambiguous_solver_fallback(
    run_dir: Path,
    maybe_rel: Optional[str],
    default_names: Sequence[str],
    primary_solver: Optional[str],
) -> None:
    if primary_solver is not None or not _is_default_name(maybe_rel, default_names):
        return
    p = resolve_rel(run_dir, maybe_rel)
    if p is not None and p.exists():
        return
    if (run_dir / "summary.csv").exists():
        return
    existing = [
        name
        for name in ("summary_mc.csv", "summary_boltzmann.csv", "summary_multiterm.csv")
        if (run_dir / name).exists()
    ]
    if len(existing) > 1:
        raise ValueError(
            "Multiple solver-specific summary files were found "
            f"({existing}). Set comsol_export.input.primary_solver explicitly."
        )


def discover_inputs(run_dir: Path, input_cfg: Dict) -> DiscoveredInputs:
    primary_solver = _primary_solver_name(input_cfg)
    _ensure_unambiguous_solver_fallback(
        run_dir,
        input_cfg.get("transport_csv"),
        ("transport.csv",),
        primary_solver,
    )
    transport_csv = resolve_with_fallbacks(
        run_dir,
        input_cfg.get("transport_csv"),
        fallbacks=_solver_first_fallbacks(primary_solver, "summary")
        + ("transport_table.csv",),
        default_names=("transport.csv",),
    )
    _ensure_unambiguous_solver_fallback(
        run_dir,
        input_cfg.get("rates_csv"),
        ("rates.csv",),
        primary_solver,
    )
    rates_csv = resolve_with_fallbacks(
        run_dir,
        input_cfg.get("rates_csv"),
        fallbacks=_solver_first_fallbacks(primary_solver, "summary")
        + ("rates_table.csv",),
        default_names=("rates.csv",),
    )

    eedf_cfg = dict(input_cfg.get("eedf", {}))
    mode = str(eedf_cfg.get("mode", "auto")).lower()

    eedf_stacked_csv = resolve_with_fallbacks(
        run_dir,
        eedf_cfg.get("stacked_csv"),
        fallbacks=_solver_first_table_fallbacks(primary_solver, "eedf_table")
        + _solver_first_table_fallbacks(primary_solver, "energy_table"),
        default_names=("eedf.csv",),
    )
    eedf_files_glob = eedf_cfg.get("files_glob", None)
    if eedf_files_glob and not Path(eedf_files_glob).is_absolute():
        eedf_files_glob = str(run_dir / eedf_files_glob)

    # Auto mode: prefer stacked if exists, else per_file glob
    if mode == "auto":
        if eedf_stacked_csv and eedf_stacked_csv.exists():
            mode = "stacked_csv"
        elif eedf_files_glob:
            mode = "per_file"

    if mode == "stacked_csv":
        if eedf_stacked_csv is None:
            raise ValueError("eedf.mode=stacked_csv but eedf.stacked_csv is not set")
        if not eedf_stacked_csv.exists():
            raise FileNotFoundError(str(eedf_stacked_csv))
    elif mode == "per_file":
        if not eedf_files_glob:
            raise ValueError("eedf.mode=per_file but eedf.files_glob is not set")
        files = glob.glob(eedf_files_glob)
        if not files:
            raise FileNotFoundError(f"No EEDF files matched: {eedf_files_glob}")
    else:
        raise ValueError(f"Unknown eedf.mode: {mode}")

    return DiscoveredInputs(
        transport_csv=transport_csv,
        rates_csv=rates_csv,
        eedf_stacked_csv=eedf_stacked_csv if mode == "stacked_csv" else None,
        eedf_files_glob=eedf_files_glob if mode == "per_file" else None,
    )
