from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import glob
import re

import numpy as np
import pandas as pd

from .columns import require_column
from .io import read_table


@dataclass
class EEDFLoadResult:
    df: pd.DataFrame  # columns: eps_eV, epsbar_eV, f_eV_m32
    meta: Dict


def load_eedf_stacked(
    path: Path,
    delimiter: Optional[str],
    columns_cfg: Dict[str, List[str]],
) -> EEDFLoadResult:
    df = read_table(path, delimiter=delimiter)
    c_eps = require_column(df, columns_cfg["eps_eV"], "eedf.eps_eV")
    c_epsbar = require_column(df, columns_cfg["epsbar_eV"], "eedf.epsbar_eV")
    c_f = require_column(df, columns_cfg["f"], "eedf.f")

    out = df[[c_eps, c_epsbar, c_f]].copy()
    out.columns = ["eps_eV", "epsbar_eV", "f"]
    out = out.dropna()

    meta = {
        "mode": "stacked_csv",
        "source": str(path),
        "rows": int(len(out)),
    }
    return EEDFLoadResult(df=out, meta=meta)


def load_eedf_per_file(
    files_glob: str,
    delimiter: Optional[str],
    columns_cfg: Dict[str, List[str]],
    file_mean_energy_regex: Optional[str] = None,
) -> EEDFLoadResult:
    paths = [Path(p) for p in glob.glob(files_glob)]
    if not paths:
        raise FileNotFoundError(f"No files matched: {files_glob}")

    regex = re.compile(file_mean_energy_regex) if file_mean_energy_regex else None

    frames: List[pd.DataFrame] = []
    for p in sorted(paths):
        df = read_table(p, delimiter=delimiter)
        c_eps = require_column(df, columns_cfg["eps_eV"], "eedf.eps_eV")
        c_f = require_column(df, columns_cfg["f"], "eedf.f")

        if any(c.lower() in [c_.lower() for c_ in df.columns] for c in columns_cfg.get("epsbar_eV", [])):
            c_epsbar = require_column(df, columns_cfg["epsbar_eV"], "eedf.epsbar_eV")
            epsbar_series = df[c_epsbar]
        else:
            if regex is None:
                raise ValueError(
                    f"EEDF file {p} does not contain epsbar column and file_mean_energy_regex is not provided."
                )
            m = regex.search(p.name)
            if not m:
                raise ValueError(f"Regex did not match filename for epsbar: {p.name} / {regex.pattern}")
            epsbar_val = float(m.group("epsbar"))
            epsbar_series = pd.Series([epsbar_val] * len(df))

        out = pd.DataFrame({
            "eps_eV": df[c_eps].astype(float),
            "epsbar_eV": epsbar_series.astype(float),
            "f": df[c_f].astype(float),
        })
        out = out.dropna()
        frames.append(out)

    out_all = pd.concat(frames, ignore_index=True)
    meta = {
        "mode": "per_file",
        "files_glob": files_glob,
        "files": [str(p) for p in paths],
        "rows": int(len(out_all)),
    }
    return EEDFLoadResult(df=out_all, meta=meta)


def normalize_eedf_pdf_to_comsol(df: pd.DataFrame, eps_min: float = 1e-12) -> pd.DataFrame:
    """Convert a PDF defined as p(ε) [eV^-1] into COMSOL-like f(ε) [eV^-3/2] via f = p / sqrt(ε).

    This is a *best-effort* conversion. Many Boltzmann solvers already output f(ε) in eV^-3/2, in which case
    you should keep eedf_normalization='none'.
    """
    out = df.copy()
    eps = np.maximum(out["eps_eV"].to_numpy(dtype=float), eps_min)
    out["f"] = out["f"].to_numpy(dtype=float) / np.sqrt(eps)
    return out
