"""The simulator-neutral reader uses the exact qualified table values."""

import csv
import hashlib
import json

import pytest

from electron_swarm.core.constants import BOLTZMANN_J_K
from swarm_workflow.campaign.sweep import run_sweep
from swarm_workflow.tables import (
    CoefficientLookupError,
    build_tables,
    load_coefficients,
)

from product_helpers import ROOT, base_product_config, write_config


GAS = {
    "gas_temperature_K": 300.0,
    "gas_number_density_m3": 100.0 / (BOLTZMANN_J_K * 300.0),
    "mole_fractions": {"Ar": 1.0},
}


@pytest.fixture(scope="module")
def source_tables(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("coefficient-source")
    data = base_product_config(tmp)
    data["run"]["e_over_n_Td"] = [10, 20, 30]
    data["cross_sections"]["files"][0]["path"] = str(
        ROOT / "examples" / "cross_sections" / "argon_application_library.csv"
    )
    data["campaign"] = {"database": "source.sqlite"}
    run_sweep(write_config(tmp, data))
    build_tables(tmp / "source.sqlite", tmp / "tables", source="two_term")
    return tmp / "tables/mixture_0000"


def test_lookup_preserves_units_anchor_values_and_reaction_identity(source_tables):
    tables = load_coefficients(source_tables, **GAS)
    reaction = next(
        item for item in tables.reactions if item.process_type == "ionization"
    )
    series = tables.series("rate_coefficient_m3_s", reaction=reaction)
    with (source_tables / "rates_vs_en.csv").open(newline="") as handle:
        anchors = {
            float(row["E_over_N_Td"]): float(row["rate_coefficient_m3_s"])
            for row in csv.DictReader(handle)
            if row["process_type"] == "ionization"
        }
    assert series.unit == "m^3/s"
    assert series.at(20) == anchors[20]
    assert series.at(15) == pytest.approx((anchors[10] + anchors[20]) / 2)
    assert tables.series("reduced_mobility_m2_V_s_m3").at(25) > 0
    with pytest.raises(CoefficientLookupError, match="no qualified"):
        series.at(30.1)
    with pytest.raises(CoefficientLookupError, match="differs"):
        load_coefficients(source_tables, **{**GAS, "gas_temperature_K": 301})
    with pytest.raises(CoefficientLookupError, match="mole fractions differ"):
        load_coefficients(
            source_tables, **{**GAS, "mole_fractions": {"Ar": 0.9, "He": 0.1}}
        )


def test_lookup_checks_selected_csv_bytes_and_units(tmp_path, source_tables):
    name = "transport_vs_en.csv"
    manifest = json.loads((source_tables / "manifest.json").read_text())
    record = manifest["tables"][name]
    (tmp_path / name).write_bytes((source_tables / name).read_bytes() + b"\n")
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(CoefficientLookupError, match="SHA-256"):
        load_coefficients(tmp_path, **GAS).series("mean_energy_eV")
    record["units"]["mean_energy_eV"] = "K"
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(CoefficientLookupError, match="units"):
        load_coefficients(tmp_path, **GAS).series("mean_energy_eV")


def test_lookup_refuses_a_gap_and_restricted_mc_evidence(tmp_path, source_tables):
    name = "transport_vs_en.csv"
    manifest = json.loads((source_tables / "manifest.json").read_text())
    manifest["source"] = "monte_carlo"
    manifest["source_policy"].update(
        source="monte_carlo",
        qualification_profile="function_eedf_restricted_lmea",
        transport_eligible_E_over_N_Td=[10, 30],
        transport_excluded_E_over_N_Td={"20": ["failed"]},
    )
    lines = (source_tables / name).read_bytes().splitlines(keepends=True)
    content = b"".join((lines[0], lines[1], lines[3]))
    (tmp_path / name).write_bytes(content)
    manifest["tables"][name]["sha256"] = hashlib.sha256(content).hexdigest()
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    tables = load_coefficients(tmp_path, **GAS)
    series = tables.series("reduced_mobility_m2_V_s_m3")
    assert series.at(30) > 0
    with pytest.raises(CoefficientLookupError, match="no qualified"):
        series.at(25)
    with pytest.raises(CoefficientLookupError, match="restricted MC qualification"):
        tables.series("reduced_diffusion_L_m2_s_m3")
