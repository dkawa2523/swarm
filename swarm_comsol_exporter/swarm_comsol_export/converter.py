from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from .columns import require_column
from .discovery import discover_inputs
from .eedf import (
    EEDFLoadResult,
    load_eedf_per_file,
    load_eedf_stacked,
    normalize_eedf_pdf_to_comsol,
)
from .io import ensure_dir, read_table, write_comsol_spreadsheet_csv


@dataclass
class ExportPaths:
    out_dir: Path
    transport_csv: Path
    rates_csv: Path
    eedf_csv: Path
    report_json: Path


def _deduplicate_epsbar(df: pd.DataFrame, mode: str, tol: float) -> pd.DataFrame:
    # Group by epsbar within tolerance by rounding after scaling
    if df.empty:
        return df

    epsbar = df["epsbar_eV"].to_numpy(dtype=float)
    scale = 1.0 / max(tol, 1e-30)
    key = np.round(epsbar * scale).astype(np.int64)
    df = df.copy()
    df["_epsbar_key"] = key

    if mode == "error":
        if df["_epsbar_key"].duplicated().any():
            dups = df[df["_epsbar_key"].duplicated(keep=False)].sort_values("_epsbar_key")
            raise ValueError(f"Duplicate epsbar detected within tol={tol}.\n{dups.head(20)}")
        return df.drop(columns=["_epsbar_key"])
    elif mode == "first":
        return df.drop_duplicates("_epsbar_key", keep="first").drop(columns=["_epsbar_key"])
    elif mode == "average":
        numeric_cols = [c for c in df.columns if c not in ["_epsbar_key"]]
        # Average all numeric columns; for epsbar, take mean as well
        grouped = df.groupby("_epsbar_key", as_index=False)[numeric_cols].mean(numeric_only=True)
        return grouped.drop(columns=["_epsbar_key"])
    else:
        raise ValueError(f"Unknown deduplicate mode: {mode}")


def _sort_epsbar(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values("epsbar_eV").reset_index(drop=True)


def _extend_candidates(cfg: Dict[str, List[str]], key: str, extra: Iterable[str]) -> None:
    existing = cfg.get(key)
    if existing is None:
        cfg[key] = list(extra)
        return
    if not isinstance(existing, list):
        existing = [existing]
    for cand in extra:
        if cand not in existing:
            existing.append(cand)
    cfg[key] = existing


def _infer_rate_columns(columns: Iterable[str]) -> List[str]:
    inferred: List[str] = []
    for c in columns:
        name = str(c)
        lower = name.lower()
        if "rate coeff" in lower and "error" not in lower:
            inferred.append(c)
    return inferred


def export_to_comsol(
    run_dir: Path,
    output_dir: Path,
    input_cfg: Dict,
    processing: Dict,
    output_cfg: Dict,
) -> ExportPaths:
    out_dir = ensure_dir(output_dir)

    transport_filename = output_cfg.get("transport_filename", "comsol_transport.csv")
    rates_filename = output_cfg.get("rates_filename", "comsol_rates.csv")
    eedf_filename = output_cfg.get("eedf_filename", "comsol_eedf.csv")
    report_filename = output_cfg.get("report_filename", "comsol_export_report.json")
    delimiter = output_cfg.get("delimiter", ",")
    float_format = output_cfg.get("float_format", "%.8e")

    paths = ExportPaths(
        out_dir=out_dir,
        transport_csv=out_dir / transport_filename,
        rates_csv=out_dir / rates_filename,
        eedf_csv=out_dir / eedf_filename,
        report_json=out_dir / report_filename,
    )

    discovered = discover_inputs(run_dir, input_cfg)

    # Load transport
    transport_meta = {}
    if discovered.transport_csv is None:
        raise ValueError("input.transport_csv is not set (or not found). Please set comsol_export.input.transport_csv.")
    transport_df = read_table(discovered.transport_csv, delimiter=None if delimiter == "auto" else delimiter)
    col_cfg = dict(input_cfg.get("columns", {}))
    _extend_candidates(
        col_cfg,
        "epsbar_eV",
        [
            "epsbar_eV",
            "mean_energy_eV",
            "Ebar_eV",
            "mean_eps_eV",
            "mean energy (eV)",
        ],
    )
    _extend_candidates(
        col_cfg,
        "mu",
        [
            "mu_m2_V_s",
            "mobility",
            "mu",
            "bulk drift velocity (m.s-1)",
            "flux drift velocity (m.s-1)",
        ],
    )
    _extend_candidates(
        col_cfg,
        "De",
        [
            "De_m2_s",
            "diffusion",
            "D",
            "bulk L diffusion coeff. * N (m-1.s-1)",
            "flux L diffusion coeff. * N (m-1.s-1)",
        ],
    )
    c_epsbar = require_column(transport_df, col_cfg["epsbar_eV"], "transport.epsbar_eV")
    c_mu = require_column(transport_df, col_cfg["mu"], "transport.mu")
    c_De = require_column(transport_df, col_cfg["De"], "transport.De")

    transport_out = transport_df[[c_epsbar, c_mu, c_De]].copy()
    transport_out.columns = ["epsbar_eV", "mu_m2_V_s", "De_m2_s"]
    transport_out = transport_out.dropna()
    transport_out = transport_out.astype(float)

    # Optional processing
    if bool(processing.get("sort_by_epsbar", True)):
        transport_out = _sort_epsbar(transport_out)

    dedup_mode = str(processing.get("deduplicate_epsbar", "average")).lower()
    tol = float(processing.get("epsbar_tol", 1e-8))
    transport_out = _deduplicate_epsbar(transport_out, dedup_mode, tol)

    transport_comment = [
        "COMSOL Interpolation (Spreadsheet) | args: epsbar_eV | cols: epsbar_eV, mu_m2_V_s, De_m2_s",
        f"Source: {discovered.transport_csv}",
    ]
    transport_comment.append("Column map: 1=epsbar_eV, 2=mu_m2_V_s, 3=De_m2_s")
    transport_comment.append(f"Source columns: epsbar='{c_epsbar}', mu='{c_mu}', De='{c_De}'")
    if c_mu in {"bulk drift velocity (m.s-1)", "flux drift velocity (m.s-1)"}:
        transport_comment.append(f"Note: mu source column='{c_mu}' (velocity). Check units before COMSOL import.")
    if c_De in {"bulk L diffusion coeff. * N (m-1.s-1)", "flux L diffusion coeff. * N (m-1.s-1)"}:
        transport_comment.append(f"Note: De source column='{c_De}' (D*N). Check units before COMSOL import.")
    write_comsol_spreadsheet_csv(
        transport_out[["epsbar_eV", "mu_m2_V_s", "De_m2_s"]],
        paths.transport_csv,
        comment_lines=transport_comment,
        delimiter=delimiter if delimiter != "auto" else ",",
        float_format=float_format,
    )
    transport_meta.update({
        "source": str(discovered.transport_csv),
        "rows": int(len(transport_out)),
        "epsbar_min": float(transport_out["epsbar_eV"].min()) if len(transport_out) else None,
        "epsbar_max": float(transport_out["epsbar_eV"].max()) if len(transport_out) else None,
        "column_map": ["1=epsbar_eV", "2=mu_m2_V_s", "3=De_m2_s"],
        "source_columns": {"epsbar": c_epsbar, "mu": c_mu, "De": c_De},
    })

    # Load rates
    rates_meta = {}
    if discovered.rates_csv is None:
        raise ValueError("input.rates_csv is not set (or not found). Please set comsol_export.input.rates_csv.")
    rates_df = read_table(discovered.rates_csv, delimiter=None if delimiter == "auto" else delimiter)
    c_epsbar_r = require_column(rates_df, col_cfg["epsbar_eV"], "rates.epsbar_eV")

    rate_prefix = str(col_cfg.get("rate_prefix", "k_"))
    rate_prefix_used = rate_prefix
    rate_cols = [c for c in rates_df.columns if str(c).startswith(rate_prefix)]
    if not rate_cols:
        rate_cols = _infer_rate_columns(rates_df.columns)
        if rate_cols:
            rate_prefix_used = "<auto>"
    if not rate_cols:
        raise ValueError(f"No rate columns found with prefix '{rate_prefix}'. Available={list(rates_df.columns)}")

    rates_out = rates_df[[c_epsbar_r] + rate_cols].copy()
    rates_out = rates_out.dropna()
    rates_out = rates_out.astype(float)
    rates_out = rates_out.rename(columns={c_epsbar_r: "epsbar_eV"})

    if bool(processing.get("sort_by_epsbar", True)):
        rates_out = _sort_epsbar(rates_out)

    rates_out = _deduplicate_epsbar(rates_out, dedup_mode, tol)

    # Ensure column order: epsbar first
    rates_out = rates_out[["epsbar_eV"] + rate_cols]

    rates_comment = [
        "COMSOL Interpolation (Spreadsheet) | args: epsbar_eV | cols: epsbar_eV, k_* ...",
        f"Source: {discovered.rates_csv}",
        f"rate_prefix: {rate_prefix_used}",
    ]
    # Add a column map so COMSOL users can identify each reaction rate column.
    col_map = ["1=epsbar_eV"] + [f"{i+2}={c}" for i, c in enumerate(rate_cols)]
    rates_comment.append("Column map: " + ", ".join(col_map))
    write_comsol_spreadsheet_csv(
        rates_out,
        paths.rates_csv,
        comment_lines=rates_comment,
        delimiter=delimiter if delimiter != "auto" else ",",
        float_format=float_format,
    )
    rates_meta.update({
        "source": str(discovered.rates_csv),
        "rows": int(len(rates_out)),
        "epsbar_min": float(rates_out["epsbar_eV"].min()) if len(rates_out) else None,
        "epsbar_max": float(rates_out["epsbar_eV"].max()) if len(rates_out) else None,
        "rate_columns": rate_cols,
        "rate_prefix_used": rate_prefix_used,
        "column_map": ["1=epsbar_eV"] + [f"{i+2}={c}" for i, c in enumerate(rate_cols)],
    })

    # Load EEDF
    eedf_cfg = dict(input_cfg.get("eedf", {}))
    eedf_cols_cfg = dict(eedf_cfg.get("columns", {}))
    eedf_cols_cfg.setdefault("eps_eV", ["eps_eV", "energy_eV", "E_eV", "epsilon"])
    eedf_cols_cfg.setdefault("epsbar_eV", ["epsbar_eV", "mean_energy_eV", "Ebar_eV", "mean_eps_eV"])
    eedf_cols_cfg.setdefault("f", ["f_eV_m32", "eedf", "f", "F"])
    _extend_candidates(
        eedf_cols_cfg,
        "eps_eV",
        [
            "energy (eV)",
        ],
    )
    _extend_candidates(
        eedf_cols_cfg,
        "epsbar_eV",
        [
            "mean energy (eV)",
        ],
    )
    _extend_candidates(
        eedf_cols_cfg,
        "f",
        [
            "eedf (eV-1)",
        ],
    )

    if discovered.eedf_stacked_csv is not None:
        eedf_res = load_eedf_stacked(
            discovered.eedf_stacked_csv,
            delimiter=None if delimiter == "auto" else delimiter,
            columns_cfg=eedf_cols_cfg,
        )
    elif discovered.eedf_files_glob is not None:
        eedf_res = load_eedf_per_file(
            discovered.eedf_files_glob,
            delimiter=None if delimiter == "auto" else delimiter,
            columns_cfg=eedf_cols_cfg,
            file_mean_energy_regex=eedf_cfg.get("file_mean_energy_regex"),
        )
    else:
        raise RuntimeError("EEDF discovery failed unexpectedly.")

    eedf_norm = str(processing.get("eedf_normalization", "none")).lower()
    eedf_df = eedf_res.df.copy()
    # standardize types
    eedf_df = eedf_df.astype(float)
    if eedf_norm == "pdf_to_comsol":
        eedf_df = normalize_eedf_pdf_to_comsol(eedf_df)
    elif eedf_norm == "none":
        pass
    else:
        raise ValueError(f"Unknown eedf_normalization: {eedf_norm}")

    # Sort for nicer file (not required)
    eedf_df = eedf_df.sort_values(["epsbar_eV", "eps_eV"]).reset_index(drop=True)

    eedf_comment = [
        "COMSOL Interpolation (Spreadsheet) | args: eps_eV, epsbar_eV | col: f (EEDF)",
        f"Source mode: {eedf_res.meta.get('mode')}",
        f"Normalization: {eedf_norm}",
    ]
    if "source" in eedf_res.meta:
        eedf_comment.append(f"Source: {eedf_res.meta['source']}")
    if "files_glob" in eedf_res.meta:
        eedf_comment.append(f"files_glob: {eedf_res.meta['files_glob']}")
    eedf_comment.append("Column map: 1=eps_eV, 2=epsbar_eV, 3=f_eV_m32")

    # Output columns in strict order
    eedf_out = eedf_df[["eps_eV", "epsbar_eV", "f"]].copy()
    eedf_out.columns = ["eps_eV", "epsbar_eV", "f_eV_m32"]
    write_comsol_spreadsheet_csv(
        eedf_out,
        paths.eedf_csv,
        comment_lines=eedf_comment,
        delimiter=delimiter if delimiter != "auto" else ",",
        float_format=float_format,
    )

    eedf_meta = dict(eedf_res.meta)
    eedf_meta.update({
        "rows": int(len(eedf_out)),
        "eps_min": float(eedf_out["eps_eV"].min()) if len(eedf_out) else None,
        "eps_max": float(eedf_out["eps_eV"].max()) if len(eedf_out) else None,
        "epsbar_min": float(eedf_out["epsbar_eV"].min()) if len(eedf_out) else None,
        "epsbar_max": float(eedf_out["epsbar_eV"].max()) if len(eedf_out) else None,
        "output_unit_note": "f_eV_m32 is intended as eV^(-3/2) (COMSOL blog/user guide convention).",
        "column_map": ["1=eps_eV", "2=epsbar_eV", "3=f_eV_m32"],
    })

    # Export report
    report = {
        "run_dir": str(run_dir),
        "output_dir": str(out_dir),
        "transport": transport_meta,
        "rates": rates_meta,
        "eedf": eedf_meta,
        "timestamp_utc": __import__("datetime").datetime.utcnow().isoformat() + "Z",
    }
    paths.report_json.write_text(__import__("json").dumps(report, indent=2), encoding="utf-8")

    return paths
