"""Simulator-neutral DC reaction-coefficient points and their physical context."""

from copy import deepcopy
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path

import pandas as pd
import numpy as np
from electron_swarm import __version__

from electron_swarm.core.identity import content_digest
from electron_swarm.core.output import OutputConfig
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS
from electron_swarm.core.statistics import summarize_replicates
from .artifacts import PRODUCT_ARTIFACT_SCHEMA, ProductRunArtifact, json_scalar
from .publication import read_generation, publish_tables


RATE_ARTIFACT_SCHEMA = "swarm.reaction_rates.v1"
RATE_COLUMNS = {
    "series_id": "string",
    "solver": "string",
    "E_over_N_Td": "number",
    "species": "string",
    "process": "string",
    "process_type": "string",
    "threshold_eV": "number",
    "target_species_fraction": "number",
    "rate_coefficient_m3_s": "number",
    "mixture_weighted_rate_m3_s": "number",
    "replicates": "integer",
    "standard_error_m3_s": "number",
    "ci95_low_m3_s": "number",
    "ci95_high_m3_s": "number",
    "estimate_status": "string",
    "numerical_status": "string",
}
RATE_UNITS = {
    "E_over_N_Td": "Td",
    "threshold_eV": "eV",
    "target_species_fraction": "1",
    **{
        name: "m^3/s"
        for name in (
            "rate_coefficient_m3_s",
            "mixture_weighted_rate_m3_s",
            "standard_error_m3_s",
            "ci95_low_m3_s",
            "ci95_high_m3_s",
        )
    },
}


@dataclass(frozen=True, slots=True)
class RateSample:
    solver: str
    case_id: str
    e_over_n_Td: float
    physical_context: dict
    rates: list[dict]
    metadata: dict
    replicate: int = 0
    seed: int | None = None


@dataclass(frozen=True, slots=True)
class ReactionRatesArtifact:
    metadata: dict
    table: pd.DataFrame
    _digest: str | None = field(default=None, init=False, repr=False)

    def _seal(self):
        object.__setattr__(self, "_digest", self._content_digest())
        return self

    def _content_digest(self):
        return content_digest((self.metadata, self.table.to_dict("records")))


def case_rate_sample(case, *, replicate=0, seed=None) -> RateSample:
    return RateSample(
        case.solver,
        case.case_id,
        case.e_over_n_Td,
        asdict(case.physical_context) if case.physical_context else {},
        [asdict(rate) for rate in case.rates],
        deepcopy(case.metadata),
        replicate,
        seed,
    )


def _number(value, name):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise ValueError(f"reaction_rates {name} must be finite numeric data")
    return float(value)


def _channels(sample):
    result = {}
    species = {s[0]: s[1] for s in sample.physical_context["species"]}
    for row in sample.rates:
        key = tuple(
            row[name] for name in ("species", "process", "process_type", "threshold_eV")
        )
        if any(not isinstance(v, str) or not v for v in key[:3]):
            raise ValueError("reaction_rates requires explicit species/process/type")
        if key[3] is not None:
            _number(key[3], "threshold_eV")
        if key in result:
            raise ValueError(
                "reaction_rates has duplicate channel identities; cannot guess a sum"
            )
        fraction = _number(row["target_species_fraction"], "target_species_fraction")
        if fraction != species.get(row["species"]):
            raise ValueError(
                "reaction_rates target fraction differs from its physical context"
            )
        value = _number(row["rate_coefficient_m3_s"], "rate_coefficient_m3_s")
        if value < 0:
            raise ValueError("reaction rate coefficient must be nonnegative")
        result[key] = (value, fraction)
    if not result:
        raise ValueError("selected solver has no reaction coefficients")
    return result


def build_reaction_rates(
    samples: list[RateSample], *, solver: str, source: dict
) -> ReactionRatesArtifact:
    """Keep native points; never infer a reaction network or a transport closure."""
    if solver not in CANONICAL_SOLVER_IDS:
        raise ValueError("reaction_rates solver must be canonical")
    samples = deepcopy(samples)
    grouped = {}
    contexts = {}
    for sample in samples:
        if sample.solver != solver:
            continue
        context = deepcopy(sample.physical_context)
        if not context or context.get("e_over_n_Td") != sample.e_over_n_Td:
            raise ValueError(
                "reaction_rates requires matching physical context for every case"
            )
        if context.get("field_type") != "dc":
            raise NotImplementedError(
                "reaction_rates export currently supports DC E/N only"
            )
        field_value = _number(sample.e_over_n_Td, "E_over_N_Td")
        if field_value <= 0:
            raise ValueError("reaction_rates requires positive DC E/N")
        context.pop("e_over_n_Td")
        digest = content_digest(context)
        contexts[digest] = context
        grouped.setdefault(digest, {}).setdefault(field_value, []).append(sample)
    if not grouped:
        raise ValueError(f"no {solver} cases are available for reaction_rates export")
    rows, series = [], {}
    for index, digest in enumerate(sorted(grouped)):
        series_id = f"series_{index:04d}"
        points = grouped[digest]
        series[series_id] = {
            "physical_context": contexts[digest],
            "E_over_N_Td": sorted(points),
            "range_Td": [min(points), max(points)],
            "cases": [],
        }
        inventory = None
        for e_over_n, replicas in sorted(points.items()):
            if solver != "monte_carlo" and len(replicas) != 1:
                raise ValueError(
                    "deterministic reaction_rates requires one case per physical point"
                )
            if solver == "monte_carlo" and len(replicas) > 1:
                if (
                    len({r.replicate for r in replicas}) != len(replicas)
                    or any(r.seed is None for r in replicas)
                    or len({r.seed for r in replicas}) != len(replicas)
                ):
                    raise ValueError(
                        "MC reaction_rates requires distinct independent replica seeds"
                    )
            channels = [_channels(sample) for sample in replicas]
            if inventory is None:
                inventory = set(channels[0])
            if any(set(c) != inventory for c in channels):
                raise ValueError(
                    "reaction channel coverage differs across points or replicas"
                )
            statuses = {r.metadata.get("numerical_status", "unknown") for r in replicas}
            numerical_status = next(iter(statuses)) if len(statuses) == 1 else "mixed"
            numerical_status = {
                "postprocessed_transport_stale": "postprocessed_rates",
                "requires_independent_ensemble": "statistical_precision_unassessed",
            }.get(numerical_status, numerical_status)
            for sample in replicas:
                series[series_id]["cases"].append(
                    {
                        "E_over_N_Td": e_over_n,
                        "case_id": sample.case_id,
                        "replicate": sample.replicate,
                        "seed": sample.seed,
                        "solver_metadata": deepcopy(sample.metadata),
                    }
                )
            for key in sorted(inventory, key=lambda k: json.dumps(k)):
                values = [c[key][0] for c in channels]
                fraction = channels[0][key][1]
                stats = summarize_replicates(
                    values, zero_is_censored=solver == "monte_carlo"
                )
                monte_carlo = solver == "monte_carlo"
                estimate = (
                    "deterministic"
                    if not monte_carlo
                    else "censored_all_zero"
                    if stats.estimate_status == "censored_all_zero"
                    else "independent_replica_estimate"
                    if len(replicas) > 1
                    else "single_replica"
                )
                rows.append(
                    {
                        "series_id": series_id,
                        "solver": solver,
                        "E_over_N_Td": e_over_n,
                        **dict(
                            zip(
                                ("species", "process", "process_type", "threshold_eV"),
                                key,
                            )
                        ),
                        "target_species_fraction": fraction,
                        "rate_coefficient_m3_s": stats.mean,
                        "mixture_weighted_rate_m3_s": stats.mean * fraction,
                        "replicates": len(replicas),
                        "standard_error_m3_s": stats.standard_error
                        if monte_carlo
                        else None,
                        "ci95_low_m3_s": stats.ci95_low if monte_carlo else None,
                        "ci95_high_m3_s": stats.ci95_high if monte_carlo else None,
                        "estimate_status": estimate,
                        "numerical_status": numerical_status,
                    }
                )
    metadata = {
        "schema": RATE_ARTIFACT_SCHEMA,
        "package_version": __version__,
        "solver": solver,
        "quantity": "reaction_rates",
        "independent_variable": "E_over_N_Td",
        "units": deepcopy(RATE_UNITS),
        "coefficient_basis": "k per target particle; mixture_weighted_rate=x*k; binary source=n_e*n_target*k",
        "data_role": "coefficient_points",
        "qualification": "numerical status and independent-replica uncertainty are reported; physical and downstream validation are not established",
        "interpolation": "none",
        "extrapolation": "unsupported",
        "uncertainty_method": "95% Student-t interval across independent MC replicas; unavailable for one replica or censored all-zero estimates",
        "series": series,
        "source": deepcopy(source),
    }
    metadata = json.loads(json.dumps(metadata, default=json_scalar, allow_nan=False))
    return ReactionRatesArtifact(
        metadata, pd.DataFrame(rows, columns=RATE_COLUMNS, dtype=object)
    )._seal()


def rates_from_product(
    artifact: ProductRunArtifact, *, solver: str, source: dict
) -> ReactionRatesArtifact:
    """Use the verified original rates, not the transport-dependent COMSOL tables."""
    if "rates_csv" not in artifact.tables:
        raise ValueError("product run contains no rates")
    samples = []
    rates = artifact.tables["rates_csv"]
    summary = artifact.tables["summary_csv"]
    for case in artifact.metadata["cases"]:
        if case["solver"] != solver:
            continue
        key = (
            (summary.solver == solver)
            & (summary.case_id == case["case_id"])
            & (summary.E_over_N_Td == case["E_over_N_Td"])
        )
        selected = summary.loc[key]
        if len(selected) != 1:
            raise ValueError("product rates do not have a unique case context")
        metadata = {
            k[5:]: _optional(v)
            for k, v in selected.iloc[0].items()
            if k.startswith("meta_")
        }
        mask = (
            (rates.solver == solver)
            & (rates.case_id == case["case_id"])
            & (rates.E_over_N_Td == case["E_over_N_Td"])
        )
        records = [
            {k: _optional(v) for k, v in row.items()}
            for row in rates.loc[mask].to_dict("records")
        ]
        samples.append(
            RateSample(
                solver,
                case["case_id"],
                case["E_over_N_Td"],
                case["physical_context"],
                records,
                metadata,
                case.get("replicate") or 0,
                metadata.get("monte_carlo_case_seed"),
            )
        )
    return build_reaction_rates(samples, solver=solver, source=source)


def _optional(value):
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


def write_reaction_rates(
    artifact: ReactionRatesArtifact, output: OutputConfig
) -> dict[str, Path]:
    """Publish a neutral coefficient artifact, preserving source and alteration scope."""
    output = deepcopy(output)
    output.validate()
    artifact = deepcopy(artifact)
    if output.export is not None and output.export.solver != artifact.metadata.get(
        "solver"
    ):
        raise ValueError("output.export.solver differs from the reaction_rates data")
    if artifact._digest is None or artifact._content_digest() != artifact._digest:
        status = "provided" if artifact._digest is None else "derived"
        artifact.metadata["source"] = {
            "kind": status,
            "source": artifact.metadata.get("source"),
        }
        artifact.table["numerical_status"] = status
    _validate_artifact(artifact)
    metadata = {
        "artifact_schema": RATE_ARTIFACT_SCHEMA,
        "artifact_role": "reaction_rates",
        "reaction_rates": artifact.metadata,
    }
    return publish_tables(
        metadata,
        {"reaction_rates_csv": artifact.table},
        {"reaction_rates_csv": RATE_COLUMNS},
        output,
    )


def read_reaction_rates(manifest_path: str | Path) -> ReactionRatesArtifact:
    from .artifacts import FIGURE_ARTIFACT_SCHEMA

    path = Path(manifest_path).resolve()
    suffix = "_manifest.json"
    if not path.name.endswith(suffix):
        raise ValueError("expected a reaction rates or product _manifest.json")
    artifact = read_generation(
        path,
        path.name[: -len(suffix)],
        {
            RATE_ARTIFACT_SCHEMA: "reaction_rates",
            PRODUCT_ARTIFACT_SCHEMA: "swarm_product_run",
            FIGURE_ARTIFACT_SCHEMA: "swarm_figures",
        },
    )
    metadata = artifact.metadata.get("reaction_rates")
    if not metadata or metadata.get("schema") != RATE_ARTIFACT_SCHEMA:
        raise ValueError("manifest contains no reaction_rates export")
    result = ReactionRatesArtifact(metadata, artifact.tables["reaction_rates_csv"])
    _validate_artifact(result)
    return result._seal()


def _validate_artifact(artifact: ReactionRatesArtifact) -> None:
    metadata, table = artifact.metadata, artifact.table
    if metadata.get("units") != RATE_UNITS or list(table) != list(RATE_COLUMNS):
        raise ValueError(
            "reaction_rates columns and units must use the canonical schema"
        )
    if metadata.get("solver") not in CANONICAL_SOLVER_IDS or set(table.solver) != {
        metadata["solver"]
    }:
        raise ValueError("reaction_rates requires one explicit canonical solver")
    if table.duplicated(
        [
            "series_id",
            "E_over_N_Td",
            "species",
            "process",
            "process_type",
            "threshold_eV",
        ]
    ).any():
        raise ValueError("reaction_rates has duplicate physical points/channels")
    for row in table.to_dict("records"):
        coefficient = _number(row["rate_coefficient_m3_s"], "rate_coefficient_m3_s")
        fraction = _number(row["target_species_fraction"], "target_species_fraction")
        if coefficient < 0 or not 0 <= fraction <= 1:
            raise ValueError(
                "reaction_rates requires nonnegative k and fractions in [0, 1]"
            )
        if row["mixture_weighted_rate_m3_s"] != coefficient * fraction:
            raise ValueError("reaction_rates mixture coefficient must equal x*k")
        series = metadata["series"].get(row["series_id"])
        if not series or row["E_over_N_Td"] not in series["E_over_N_Td"]:
            raise ValueError(
                "reaction_rates point has no matching physical series context"
            )
