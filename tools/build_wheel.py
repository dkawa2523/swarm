"""Build current runtime sources in a clean staging directory."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="wheel output directory")
    args = parser.parse_args()
    output = args.output.resolve()
    root = Path(__file__).resolve().parents[1]
    output.mkdir(parents=True, exist_ok=True)
    # Only this owned temporary directory is cleaned up by the context manager.
    with tempfile.TemporaryDirectory(prefix="swarm-wheel-") as temporary:
        staging = Path(temporary).resolve()
        for package in ("electron_swarm", "swarm_workflow"):
            shutil.copytree(
                root / package,
                staging / package,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
            )
        for name in ("pyproject.toml", "README.md"):
            shutil.copy2(root / name, staging / name)
        subprocess.run(
            [sys.executable, "-m", "pip", "wheel", str(staging), "--no-deps",
             "--no-build-isolation", "--wheel-dir", str(output)],
            check=True,
        )


if __name__ == "__main__":
    main()
