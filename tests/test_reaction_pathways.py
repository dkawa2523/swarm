"""Exercise catalog binding, fixed-gas screening, and published decisions."""

from hashlib import sha256
import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from electron_swarm import load_config
from electron_swarm.chemistry.catalog import CATALOG_SCHEMA, PathwayScreen
from electron_swarm.io import read_reaction_rates
from electron_swarm.io.pathways import screen_reaction_pathways
from electron_swarm.io.reaction_rates import RateSample, build_reaction_rates
from product_helpers import base_product_config, write_config
from swarm_workflow.cli import main


def _species(atoms, charge=0):
    return {"atoms": atoms, "charge": charge}


def _channel(target, process, threshold, category, products):
    return {
        "target": target,
        "process": process,
        "threshold_eV": threshold,
        "category": category,
        "products": products,
        "resolution": "resolved",
        "reference": "test dataset, known channel",
    }


def _catalog(path, digest, species, channels):
    return write_config(
        path,
        {
            "schema": CATALOG_SCHEMA,
            "datasets": {"test": {"sha256": digest, "reference": "test dataset"}},
            "species": species,
            "reactions": channels,
        },
        "catalog.yaml",
    )


@pytest.mark.parametrize("campaign", [False, True], ids=["direct", "campaign"])
def test_gas_to_pathways_and_saved_rescreen(tmp_path, capsys, campaign):
    config = base_product_config(tmp_path)
    data_path = Path(config["cross_sections"]["files"][0]["path"])
    digest = sha256(data_path.read_bytes()).hexdigest()
    catalog = _catalog(
        tmp_path,
        digest,
        {
            "Ar": _species({"Ar": 1}),
            "Ar+": _species({"Ar": 1}, 1),
        },
        {
            "ar_ionization": _channel(
                "Ar", "Ar ionization 15.76 eV", 15.76, "ionization", {"e": 2, "Ar+": 1}
            )
        },
    )
    config["output"]["export"] = {
        "solver": "two_term",
        "quantities": ["reaction_pathways"],
        "pathways": {"catalog": catalog.name, "focus_species": "Ar+"},
    }
    if campaign:
        config["campaign"] = {"database": "swarm.sqlite"}
    path = write_config(tmp_path, config)
    main(["run", str(path), "--plan", "--json"])
    plan = json.loads(capsys.readouterr().out)
    assert plan["export"]["pathways"]["focus_species"] == "Ar+"
    assert plan["export"]["pathways"]["collision_file_sha256"] == [digest]
    assert plan["export"]["pathways_csv"] == "prod_reaction_pathways.csv"
    assert not (tmp_path / "prod_manifest.json").exists()

    main(["run", str(path)])
    assert "reaction_pathways: ranked_channels_available" in capsys.readouterr().out
    table = pd.read_csv(tmp_path / "prod_reaction_pathways.csv")
    assert set(table.decision) == {"ranked", "unmapped_channel", "scattering_only"}
    ionization = table[table.decision == "ranked"].iloc[0]
    assert ionization["rank"] == 1
    assert ionization["ranking_group"] == "production"
    assert json.loads(ionization["products"]) == {"e": 2, "Ar+": 1}

    manifest = tmp_path / "prod_manifest.json"
    native = read_reaction_rates(manifest)
    context = native.metadata["series"][ionization.series_id]["physical_context"]
    rate = native.table.loc[native.table.process_type == "ionization"].iloc[0]
    expected = (
        context["gas_number_density_m3"]
        * rate.target_species_fraction
        * rate.rate_coefficient_m3_s
    )
    assert ionization.focus_source_per_electron_s_inv == pytest.approx(
        expected, rel=1e-12, abs=0
    )
    source = tmp_path / "swarm.sqlite" if campaign else manifest
    main(
        [
            "export-reaction-rates",
            str(source),
            "--solver",
            "two_term",
            "--output",
            str(tmp_path / "saved"),
            "--pathway-catalog",
            str(catalog),
            "--focus-species",
            "Ar+",
            "--minimum-activity-s-inv",
            "1e20",
        ]
    )
    assert "reaction_pathways: no_ranked_channels" in capsys.readouterr().out
    reread = pd.read_csv(tmp_path / "saved/swarm_reaction_pathways.csv")
    assert "below_minimum_activity" in set(reread.decision)
    assert reread["rank"].isna().all()

    # This catalog would still parse and have the same names. The bytes it
    # claims to describe are different, so a new calculation fails preflight.
    catalog.write_text(
        catalog.read_text(encoding="utf-8").replace(digest, "f" * 64), encoding="utf-8"
    )
    if not campaign:
        with pytest.raises(ValueError, match="SHA-256 set differs"):
            load_config(path)


@pytest.fixture
def oxygen_catalog(tmp_path):
    digest = "1" * 64
    channels = {
        "dissociation": _channel(
            "O2", "dissociation label", 5, "dissociation", {"e": 1, "O": 2}
        ),
        "attachment": _channel(
            "O2", "attachment label", None, "dissociative_attachment", {"O-": 1, "O": 1}
        ),
        "excitation": _channel(
            "O2", "excitation label", 2, "excitation", {"e": 1, "O2*": 1}
        ),
    }
    catalog = _catalog(
        tmp_path,
        digest,
        {
            "O2": _species({"O": 2}),
            "O2*": _species({"O": 2}),
            "O": _species({"O": 1}),
            "O-": _species({"O": 1}, -1),
        },
        channels,
    )
    source = {"input_evidence": {"input": {}}, "input_files": {"input": digest}}
    return catalog, source


def _sample(solver, replica, status, values, field=50):
    channels = (
        ("dissociation label", "excitation", 5.0),
        ("attachment label", "attachment", None),
        ("excitation label", "excitation", 2.0),
    )
    return RateSample(
        solver,
        f"r{replica}",
        float(field),
        {
            "species": [["O2", 0.1, 31.998], ["Ar", 0.9, 39.948]],
            "field_type": "dc",
            "gas_number_density_m3": 1e21,
            "e_over_n_Td": float(field),
        },
        [
            {
                "species": "O2",
                "process": label,
                "process_type": kind,
                "threshold_eV": threshold,
                "rate_coefficient_m3_s": value,
                "target_species_fraction": 0.1,
            }
            for (label, kind, threshold), value in zip(channels, values, strict=True)
        ],
        {"numerical_status": status},
        replica,
        17 + replica,
    )


def test_focus_accounts_for_products_and_does_not_rank_invalid_points(oxygen_catalog):
    path, source = oxygen_catalog
    samples = [
        _sample("two_term", 0, "numerically_converged", [2e-16, 3e-16, 1e-15]),
        _sample("two_term", 0, "tail_limited", [2e-16, 3e-16, 1e-15], field=60),
    ]
    rates = build_reaction_rates(samples, solver="two_term", source=source)
    metadata, table = screen_reaction_pathways(rates, PathwayScreen(path, "O"))
    assert metadata["status"] == "ranked_channels_available"
    at_50 = table[table.E_over_N_Td == 50].set_index("reaction_id")
    assert at_50.loc[
        "dissociation", "focus_source_per_electron_s_inv"
    ] == pytest.approx(4e4)
    assert at_50.loc["dissociation", "rank"] == 1
    assert at_50.loc["attachment", "rank"] == 2
    assert at_50.loc["excitation", "decision"] == "does_not_change_focus"
    assert set(table.loc[table.E_over_N_Td == 60, "decision"]) == {
        "unranked_numerical_status",
        "does_not_change_focus",
    }
    _, electron = screen_reaction_pathways(rates, PathwayScreen(path, "e"))
    attachment = electron.loc[
        (electron.E_over_N_Td == 50) & (electron.reaction_id == "attachment")
    ].iloc[0]
    assert attachment.ranking_group == "loss"
    assert attachment.focus_source_per_electron_s_inv == pytest.approx(-3e4)

    # A threshold edit produces a distinct collision channel, even if its
    # formula is identical. An active unmatched catalog entry fails closed.
    path.write_text(
        path.read_text(encoding="utf-8").replace(
            "threshold_eV: 5", "threshold_eV: 5.1"
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no matching channel"):
        screen_reaction_pathways(rates, PathwayScreen(path))


def test_mc_priority_uses_interval_and_preserves_censoring(oxygen_catalog):
    path, source = oxygen_catalog
    samples = [
        _sample(
            "monte_carlo",
            index,
            "requires_independent_ensemble",
            [value, 0.0, uncertain],
        )
        for index, (value, uncertain) in enumerate(
            zip([2e-16, 2.01e-16, 1.99e-16], [0.0, 0.0, 1e-15])
        )
    ]
    rates = build_reaction_rates(samples, solver="monte_carlo", source=source)
    _, table = screen_reaction_pathways(rates, PathwayScreen(path))
    rows = table.set_index("reaction_id")
    assert rows.loc["dissociation", "ranking_basis"] == "independent_mc_ci95_lower"
    assert (
        0
        < rows.loc["dissociation", "ranking_score_s_inv"]
        < rows.loc["dissociation", "event_frequency_s_inv"]
    )
    assert rows.loc["attachment", "decision"] == "mc_censored"
    assert rows.loc["excitation", "decision"] == "mc_interval_includes_zero"
    assert rows.loc[["attachment", "excitation"], "rank"].isna().all()


@pytest.mark.parametrize(
    ("products", "reason"),
    [
        ({"e": 1, "O": 1}, "product counts"),
        ({"e": 1, "O": 1, "O-": 1}, "charge balance"),
        ({"e": 1, "O": 1, "O2": 1}, "element or charge balance"),
    ],
)
def test_catalog_rejects_category_and_stoichiometry_errors(
    oxygen_catalog, products, reason
):
    path, _ = oxygen_catalog
    catalog = yaml.safe_load(path.read_text(encoding="utf-8"))
    catalog["reactions"]["dissociation"]["products"] = products
    path.write_text(yaml.safe_dump(catalog), encoding="utf-8")
    with pytest.raises(ValueError, match=reason):
        PathwayScreen(path)
