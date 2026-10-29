"""Content identity for scientific arrays and small structured records."""

from dataclasses import fields, is_dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def content_digest(value: object) -> str:
    digest = sha256()

    def update(item):
        if is_dataclass(item):
            update({field.name: getattr(item, field.name) for field in fields(item)})
        elif isinstance(item, np.ndarray):
            update(("array", item.dtype.str, item.shape))
            if item.dtype.hasobject:
                update(item.tolist())
            else:
                digest.update(np.ascontiguousarray(item).tobytes())
        elif isinstance(item, dict):
            digest.update(b"{")
            for key in sorted(item, key=lambda key: (type(key).__name__, str(key))):
                update(key)
                update(item[key])
            digest.update(b"}")
        elif isinstance(item, (tuple, list)):
            digest.update(b"[")
            for element in item:
                update(element)
            digest.update(b"]")
        else:
            if isinstance(item, np.generic):
                item = item.item()
            if isinstance(item, Path):
                item = str(item)
            digest.update(
                json.dumps(item, sort_keys=True, allow_nan=True).encode("utf-8") + b"\0"
            )

    update(value)
    return digest.hexdigest()
