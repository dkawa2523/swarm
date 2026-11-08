"""Screen direct electron-impact pathways using native Swarm rate points."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict
import json
import math

import pandas as pd

from electron_swarm.chemistry.catalog import (
    CHEMICAL_RATE_TYPES,
    PathwayScreen,
    Reaction,
    Stoichiometry,
    read_reaction_catalog,
)
from .reaction_rates import CHANNEL_COLUMNS, ReactionRatesArtifact


PATHWAY_COLUMNS = {
    "series_id": "string",
    "solver": "string",
    "E_over_N_Td": "number",
    **{
        name: "number" if name == "threshold_eV" else "string"
        for name in CHANNEL_COLUMNS
    },
    "reaction_id": "string",
    "category": "string",
    "resolution": "string",
    "equation": "string",
    "reactants": "string",
    "products": "string",
    "reference": "string",
    "event_frequency_s_inv": "number",
    "focus_net_stoichiometry": "integer",
    "focus_source_per_electron_s_inv": "number",
    "decision": "string",
    "ranking_group": "string",
    "ranking_basis": "string",
    "ranking_score_s_inv": "number",
    "rank": "integer",
}
_SCATTERING_TYPES = {"elastic", "effective", "momentum"}


def _equation_side(terms: Stoichiometry) -> str:
    return " + ".join(
        name if count == 1 else f"{count} {name}" for name, count in terms
    )


def _ranking_coefficient(rate: dict) -> tuple[float | None, str]:
    status = rate["estimate_status"]
    numerical = rate["numerical_status"]
    if numerical == "numerically_converged" and status == "deterministic":
        return rate["rate_coefficient_m3_s"], "converged_point"
    if numerical != "statistical_precision_unassessed":
        return None, "unranked_numerical_status"
    if status == "censored_all_zero":
        return None, "mc_censored"
    if status != "independent_replica_estimate" or rate["replicates"] < 2:
        return None, "mc_requires_independent_replicas"
    low, high = rate["ci95_low_m3_s"], rate["ci95_high_m3_s"]
    if (
        low is None
        or high is None
        or not math.isfinite(low)
        or not math.isfinite(high)
        or not low <= rate["rate_coefficient_m3_s"] <= high
    ):
        return None, "mc_interval_unavailable"
    if low <= 0:
        return None, "mc_interval_includes_zero"
    return low, "independent_mc_ci95_lower"


def _screen_channel(
    rate: dict, reaction: Reaction | None, screen: PathwayScreen, target_density: float
) -> dict:
    row = {
        column: rate[column]
        for column in ("series_id", "solver", "E_over_N_Td", *CHANNEL_COLUMNS)
    }
    if rate["process_type"] in _SCATTERING_TYPES:
        return {**row, "decision": "scattering_only"}
    if rate["process_type"] not in CHEMICAL_RATE_TYPES:
        return {**row, "decision": "unsupported_collision_family"}

    event = target_density * rate["rate_coefficient_m3_s"]
    if not math.isfinite(event):
        raise ValueError("reaction pathway event frequency exceeds the numeric range")
    row["event_frequency_s_inv"] = event
    if reaction is None:
        return {**row, "decision": "unmapped_channel"}
    row.update(
        reaction_id=reaction.id,
        category=reaction.category,
        resolution=reaction.resolution,
        equation=f"{_equation_side(reaction.reactants)} -> {_equation_side(reaction.products)}",
        reactants=json.dumps(dict(reaction.reactants), ensure_ascii=False),
        products=json.dumps(dict(reaction.products), ensure_ascii=False),
        reference=reaction.reference,
    )
    multiplier = 1
    group = reaction.category
    if screen.focus_species is not None:
        net = reaction.net(screen.focus_species)
        row.update(
            focus_net_stoichiometry=net, focus_source_per_electron_s_inv=net * event
        )
        if not net:
            return {**row, "decision": "does_not_change_focus"}
        multiplier = abs(net)
        group = "production" if net > 0 else "loss"
    row["ranking_group"] = group

    coefficient, basis = _ranking_coefficient(rate)
    if coefficient is None:
        return {**row, "decision": basis}
    score = multiplier * (target_density * coefficient)
    if not math.isfinite(score):
        raise ValueError("reaction pathway ranking score exceeds the numeric range")
    if score == 0:
        decision = "zero_activity"
    elif score < screen.minimum_activity_s_inv:
        decision = "below_minimum_activity"
    else:
        decision = "ranked"
    return {
        **row,
        "ranking_score_s_inv": score,
        "ranking_basis": basis,
        "decision": decision,
    }


def screen_reaction_pathways(
    rates: ReactionRatesArtifact, screen: PathwayScreen
) -> tuple[dict, pd.DataFrame]:
    """Return auditable decisions at each saved point with no cross-point ranking."""
    catalog = read_reaction_catalog(screen.catalog_path)
    if catalog.sha256 != screen.catalog_sha256:
        raise ValueError(
            "reaction catalog changed after planning; reload the configuration"
        )
    source = rates.metadata["source"]
    input_files = source.get("input_files")
    input_evidence = source.get("input_evidence")
    if (
        not isinstance(input_files, dict)
        or not isinstance(input_evidence, dict)
        or not input_evidence
        or any(path not in input_files for path in input_evidence)
    ):
        raise ValueError("reaction pathways require recorded collision-file provenance")
    screen.verify_collision_files(input_files[path] for path in input_evidence)

    rows = []
    for series_id, series_rows in rates.table.groupby("series_id", sort=True):
        context = rates.metadata["series"][series_id]["physical_context"]
        density = context["gas_number_density_m3"]
        if (
            not isinstance(density, (float, int))
            or not math.isfinite(density)
            or density <= 0
        ):
            raise ValueError(
                "reaction pathways require a finite positive background gas density"
            )
        fractions = {species: fraction for species, fraction, _ in context["species"]}
        observed = set(
            series_rows[list(CHANNEL_COLUMNS)].itertuples(index=False, name=None)
        )
        missing = [
            reaction.id
            for key, reaction in catalog.reactions.items()
            if key[0] in fractions and key not in observed
        ]
        if missing:
            raise ValueError(
                f"reaction catalog has no matching channel for active targets in {series_id}: "
                f"{', '.join(missing)}; check process labels, categories, and thresholds"
            )

        for _, point in series_rows.groupby("E_over_N_Td", sort=True):
            candidates: dict[str, list[dict]] = defaultdict(list)
            for rate in point.to_dict("records"):
                key = tuple(rate[column] for column in CHANNEL_COLUMNS)
                fraction = rate["target_species_fraction"]
                if fraction <= 0 or fraction != fractions.get(key[0]):
                    raise ValueError(
                        "reaction pathway target fraction differs from its saved gas state"
                    )
                row = _screen_channel(
                    rate, catalog.reactions.get(key), screen, density * fraction
                )
                rows.append(row)
                if row["decision"] == "ranked":
                    candidates[row["ranking_group"]].append(row)
            for group in candidates.values():
                levels = sorted(
                    {row["ranking_score_s_inv"] for row in group}, reverse=True
                )
                ranks = {score: rank for rank, score in enumerate(levels, 1)}
                for row in group:
                    row["rank"] = ranks[row["ranking_score_s_inv"]]

    decisions = Counter(row["decision"] for row in rows)
    metadata = {
        "schema": "swarm.reaction_pathways.v1",
        **screen.plan(),
        "solver": rates.metadata["solver"],
        "status": "ranked_channels_available"
        if decisions["ranked"]
        else "no_ranked_channels",
        "decision_counts": dict(sorted(decisions.items())),
        "datasets": [asdict(dataset) for dataset in catalog.datasets],
        "species": {
            name: {"atoms": dict(item.atoms), "charge": item.charge}
            for name, item in catalog.species.items()
        },
        "event_frequency": "n_gas * x_target * k per electron",
        "ranking": "mapped eligible channels only; descending within point and group; ties share a rank",
        "ranking_score": "n_gas * x_target * k_basis, times abs(net stoichiometry) when focused",
        "uncertainty": "MC ranking uses a positive lower endpoint of the independent-replica 95% interval",
        "units": {
            "E_over_N_Td": "Td",
            "threshold_eV": "eV",
            "event_frequency_s_inv": "1/s",
            "focus_source_per_electron_s_inv": "1/s",
            "ranking_score_s_inv": "1/s",
        },
    }
    return metadata, pd.DataFrame(rows, columns=PATHWAY_COLUMNS, dtype=object)
