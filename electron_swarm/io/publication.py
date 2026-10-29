"""One locked, atomic publication mechanism for canonical table generations."""

from contextlib import contextmanager
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import tempfile
import re
from copy import deepcopy

from electron_swarm.core.output import FIGURE_KINDS

from .artifacts import (
    ProductRunArtifact,
    read_scalar_table,
    write_scalar_table,
    json_scalar,
)

OUTPUT_SUFFIXES = (
    "summary.csv",
    "eedf.csv",
    "energy_angle_distribution.csv",
    "rates.csv",
    "solver_plan.csv",
    "rf_phase.csv",
    "rf_phase_eedf.csv",
    "comparison_summary.csv",
    "manifest.json",
    "reaction_rates.csv",
    "figure_data.csv",
)


@contextmanager
def publication_lock(directory: Path, base_name: str):
    with (directory / f".{base_name}.lock").open("a+b") as handle:
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle, fcntl.LOCK_UN)


@contextmanager
def staged_publication(output):
    directory = Path(output.directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(prefix=f".{output.base_name}.staging-", dir=directory)
    ).resolve()
    if stage.parent != directory:
        raise ValueError("output staging escaped its target directory")
    paths = {}
    keep_recovery = False
    try:
        yield stage, paths
        manifest = paths["manifest_json"]
        backup = stage / "previous"
        backup.mkdir()
        targets = {path.name for path in paths.values()}
        managed = {f"{output.base_name}_{suffix}" for suffix in OUTPUT_SUFFIXES}
        with publication_lock(directory, output.base_name):
            pattern = re.compile(
                re.escape(output.base_name)
                + r"_figure_(?:"
                + "|".join(FIGURE_KINDS)
                + r")_\d+\.(png|svg)$"
            )
            managed.update(
                path.name
                for path in directory.iterdir()
                if path.is_file() and pattern.fullmatch(path.name)
            )
            for name in managed:
                existing = directory / name
                if existing.exists():
                    shutil.copy2(existing, backup / name)
            changed = []
            try:
                for path in paths.values():
                    if path != manifest:
                        os.replace(path, directory / path.name)
                        changed.append(path.name)
                for name in sorted(managed - targets):
                    stale = directory / name
                    if stale.exists():
                        stale.unlink()
                        changed.append(name)
                os.replace(manifest, directory / manifest.name)
                changed.append(manifest.name)
            except Exception:
                try:
                    for name in reversed(changed):
                        prior, target = backup / name, directory / name
                        if prior.exists():
                            os.replace(prior, target)
                        elif target.exists():
                            target.unlink()
                except OSError as error:
                    keep_recovery = True
                    raise RuntimeError(
                        f"output rollback failed; recovery copies retained at {backup}"
                    ) from error
                raise
        paths.update({key: directory / path.name for key, path in paths.items()})
    finally:
        if not keep_recovery:
            shutil.rmtree(stage)


def read_generation(
    path: Path, base_name: str, schemas: dict[str, str], *, verify_figures: bool = True
) -> ProductRunArtifact:
    with publication_lock(path.parent, base_name):
        metadata = json.loads(path.read_text(encoding="utf-8"))
        expected_role = schemas.get(metadata.get("artifact_schema"))
        if expected_role is None or expected_role != metadata.get("artifact_role"):
            raise ValueError("unsupported output manifest")
        tables = {}
        for key, entry in metadata["files"].items():
            name = entry["name"]
            if not isinstance(name, str) or Path(name).name != name:
                raise ValueError("invalid product artifact path")
            content = (path.parent / name).read_bytes()
            if sha256(content).hexdigest() != entry["sha256"]:
                raise ValueError(
                    f"product output generation is incomplete or altered: {name}"
                )
            tables[key] = read_scalar_table(content, entry["schema"])
        if verify_figures:
            for figure in metadata.get("figures", {}).values():
                for entry in figure["files"].values():
                    name = entry["name"]
                    if Path(name).name != name:
                        raise ValueError("invalid figure artifact path")
                    if (
                        sha256((path.parent / name).read_bytes()).hexdigest()
                        != entry["sha256"]
                    ):
                        raise ValueError(
                            f"figure output generation is incomplete or altered: {name}"
                        )
        return ProductRunArtifact(metadata, tables)


def publish_tables(metadata, tables, kinds, output):
    """Commit native tables, requested figures and one manifest together."""
    output = deepcopy(output)
    output.validate()
    output.directory = Path(output.directory)
    metadata, tables, kinds = deepcopy((metadata, tables, kinds))
    canonical_tables = {
        suffix[:-4] + "_csv" for suffix in OUTPUT_SUFFIXES if suffix.endswith(".csv")
    }
    if not set(tables) <= canonical_tables:
        raise ValueError("publication tables must use canonical output roles")
    tables.pop("figure_data_csv", None)
    kinds.pop("figure_data_csv", None)
    metadata.pop("figures", None)
    with staged_publication(output) as (stage, paths):
        if output.figures:
            from .figures import render_figures, FIGURE_DATA_KINDS

            rendered = render_figures(
                ProductRunArtifact(metadata, tables),
                output.figures,
                stage,
                output.base_name,
            )
            tables["figure_data_csv"] = rendered.data
            kinds["figure_data_csv"] = FIGURE_DATA_KINDS
            metadata["figures"] = rendered.figures
            paths.update(rendered.paths)
        files = {}
        for key, frame in tables.items():
            path = stage / f"{output.base_name}_{key[:-4]}.csv"
            schema = write_scalar_table(frame, path, kinds[key])
            files[key] = {
                "name": path.name,
                "sha256": sha256(path.read_bytes()).hexdigest(),
                "schema": schema,
            }
            paths[key] = path
        metadata["files"] = files
        metadata.setdefault("publication", {}).update(
            {
                "base_name": output.base_name,
                "created_directory": str(output.directory.resolve()),
                "figures": list(output.figures),
            }
        )
        manifest = stage / f"{output.base_name}_manifest.json"
        manifest.write_text(
            json.dumps(metadata, indent=2, default=json_scalar, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        paths["manifest_json"] = manifest
    return paths
