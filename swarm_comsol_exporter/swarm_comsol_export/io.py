from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import pandas as pd


def guess_delimiter(sample_line: str) -> str:
    # Prefer comma if present, else tab, else whitespace.
    if "," in sample_line:
        return ","
    if "\t" in sample_line:
        return "\t"
    if ";" in sample_line:
        return ";"
    return r"\s+"


def read_table(path: Union[str, Path], delimiter: Optional[str] = None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(str(path))

    # Try to infer delimiter if not provided
    if delimiter is None:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("%"):
                    continue
                delimiter = guess_delimiter(line)
                break
        if delimiter is None:
            delimiter = ","

    # If delimiter is regex for whitespace, use sep with regex engine
    if delimiter == r"\s+":
        return pd.read_csv(path, sep=delimiter, comment="%", engine="python")
    return pd.read_csv(path, sep=delimiter, comment="%")


def ensure_dir(path: Union[str, Path]) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_comsol_spreadsheet_csv(
    df: pd.DataFrame,
    path: Union[str, Path],
    comment_lines: Sequence[str] = (),
    delimiter: str = ",",
    float_format: str = "%.8e",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for line in comment_lines:
            f.write(f"% {line}\n")
        # No header, numeric only
        df.to_csv(
            f,
            index=False,
            header=False,
            sep=delimiter,
            float_format=float_format,
        )
