from __future__ import annotations

from typing import Iterable, List, Optional

import pandas as pd


def find_column(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand is None:
            continue
        key = str(cand).lower()
        if key in cols:
            return cols[key]
    return None


def require_column(df: pd.DataFrame, candidates: List[str], what: str) -> str:
    c = find_column(df, candidates)
    if c is None:
        raise KeyError(f"Required column for '{what}' not found. Candidates={candidates}. Available={list(df.columns)}")
    return c
