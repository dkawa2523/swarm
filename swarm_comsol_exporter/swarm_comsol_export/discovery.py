from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import glob


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


def discover_inputs(run_dir: Path, input_cfg: Dict) -> DiscoveredInputs:
    transport_csv = resolve_with_fallbacks(
        run_dir,
        input_cfg.get("transport_csv"),
        fallbacks=("summary.csv", "transport_table.csv"),
        default_names=("transport.csv",),
    )
    rates_csv = resolve_with_fallbacks(
        run_dir,
        input_cfg.get("rates_csv"),
        fallbacks=("summary.csv", "rates_table.csv"),
        default_names=("rates.csv",),
    )

    eedf_cfg = dict(input_cfg.get("eedf", {}))
    mode = str(eedf_cfg.get("mode", "auto")).lower()

    eedf_stacked_csv = resolve_with_fallbacks(
        run_dir,
        eedf_cfg.get("stacked_csv"),
        fallbacks=("eedf_table.csv", "energy_table.csv"),
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
