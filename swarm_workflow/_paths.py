"""Small path helpers shared by workflow modules."""

from __future__ import annotations

from pathlib import Path


def discover_repo_root(path: str | Path) -> Path:
    """Return the nearest parent containing project markers.

    The workflow accepts mapping files from examples, model directories, or
    temporary test folders.  Project-relative paths should resolve from the
    closest repository root when one is visible, otherwise from the file's
    parent directory.
    """

    resolved = Path(path).resolve()
    anchor = resolved if resolved.is_dir() else resolved.parent
    for parent in (anchor, *anchor.parents):
        if (parent / "pyproject.toml").exists() or (parent / ".git").exists():
            return parent.resolve()
    return anchor.resolve()


def resolve_path(value: str | Path, root: str | Path) -> Path:
    """Resolve ``value`` against ``root`` when it is relative."""

    path = Path(value)
    if not path.is_absolute():
        path = Path(root) / path
    return path.resolve()
