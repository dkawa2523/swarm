"""Build and execute the reader-facing benchmark analysis notebook."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Sequence

_LOCAL_JUPYTER_ROOT = Path(__file__).resolve().parents[3] / ".tmp" / "report_jupyter"
_LOCAL_JUPYTER_ROOT.mkdir(parents=True, exist_ok=True)
for _name in ("ipython", "data", "config", "runtime"):
    (_LOCAL_JUPYTER_ROOT / _name).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("IPYTHONDIR", str(_LOCAL_JUPYTER_ROOT / "ipython"))
os.environ.setdefault("JUPYTER_DATA_DIR", str(_LOCAL_JUPYTER_ROOT / "data"))
os.environ.setdefault("JUPYTER_CONFIG_DIR", str(_LOCAL_JUPYTER_ROOT / "config"))
os.environ.setdefault("JUPYTER_RUNTIME_DIR", str(_LOCAL_JUPYTER_ROOT / "runtime"))
os.environ.setdefault(
    "JUPYTER_PATH", str(_LOCAL_JUPYTER_ROOT / "prefix" / "share" / "jupyter")
)

import nbformat
from nbclient import NotebookClient
from ipykernel.kernelspec import install as install_kernel

from benchmark_analysis import DEFAULT_CONFIG, PROJECT_ROOT, resolve_path, run_analysis


def build_notebook(config_path: Path, output_path: Path | None = None) -> Path:
    summary = run_analysis(config_path)
    config = json.loads(config_path.resolve().read_text(encoding="utf-8"))
    kernel_name = "swarm-report-python"
    kernel_prefix = _LOCAL_JUPYTER_ROOT / "prefix"
    install_kernel(
        kernel_name=kernel_name,
        display_name="Swarm report Python",
        prefix=str(kernel_prefix),
    )
    notebook_path = (
        output_path.resolve()
        if output_path is not None
        else resolve_path(config["outputs"]["notebook"])
    )
    data_dir = Path(config["outputs"]["data_dir"]).as_posix()
    evidence_status = summary["evidence_status"]
    profile_evidence_status = summary["data_vintages"]["profile_comparison"]
    lookup_evidence_status = summary["data_vintages"]["lookup_tables"]
    qa_status = summary["qa"]["data_quality_status"]
    gates = summary["formal_gates"]
    derived_table_count = len(summary["derived_tables"])
    figure_count = len(summary["figures"])
    tail_path = resolve_path(Path(data_dir) / "tail_convergence_audit.csv")
    with tail_path.open("r", encoding="utf-8", newline="") as handle:
        tail_rows = list(csv.DictReader(handle))
    tail_case_count = len(tail_rows)
    tail_grid_pass_count = sum(
        row["formal_tail_grid_gate_pass"].strip().lower() == "true"
        for row in tail_rows
    )
    solver_pass_count = sum(
        row["formal_solver_convergence_gate_pass"].strip().lower() == "true"
        for row in tail_rows
    )
    all_case_pass_count = sum(
        row["formal_swarm_case_gate_pass"].strip().lower() == "true"
        for row in tail_rows
    )
    notebook = nbformat.v4.new_notebook()
    notebook.metadata.update(
        {
            "kernelspec": {
                "display_name": "Swarm report Python",
                "language": "python",
                "name": kernel_name,
            },
            "language_info": {"name": "python", "version": "3"},
            "benchmark": {
                "config": config_path.resolve().relative_to(PROJECT_ROOT).as_posix(),
                "evidence_status": evidence_status,
            },
        }
    )
    notebook.cells = [
        nbformat.v4.new_markdown_cell(
            "# External Swarm tables in a COMSOL 1D argon DC glow discharge\n\n"
            "## tl;dr\n\n"
            f"- Evidence status: **{evidence_status}**; data-quality decision: "
            f"**{qa_status}**.\n"
            f"- Profile comparison/runtime evidence: **{profile_evidence_status}**; "
            f"lookup curves: **{lookup_evidence_status}**.\n"
            "- The application-library basis is the full electrode-gap **1D DC glow "
            "discharge** (cathode fall, dark-space/column structure, and anode side), "
            "not an isolated uniform positive-column calculation.\n"
            "- The external path is a **selected-channel hybrid**: Swarm supplies four "
            "transport quantities and direct eir2/eir4 Townsend data, while COMSOL "
            "retains the LEA fluid equations and residual chemistry/energy treatment.\n"
            "- The bundle carries **7 lookup quantities = 7 X-Y pairs = 14 numerical "
            "arrays**. In archived LEA, the E/N-to-mean-energy pair is stored but "
            "inactive; written properties and equation-level activation are distinct.\n"
            f"- Tail/grid criteria pass **{tail_grid_pass_count}/{tail_case_count}** "
            f"points, but the iterative solver convergence flag passes only "
            f"**{solver_pass_count}/{tail_case_count}**; therefore the full Swarm "
            f"case gate passes **{all_case_pass_count}/{tail_case_count}**.\n"
            f"- Runtime repetition gate: **{gates['runtime_repetition_gate_met']}**; "
            f"raw 200/400 mesh gate: **{gates['raw_mesh_200_400_gate_met']}**.\n"
            "- Exact piecewise-linear L2 and spatial integral ratio are reported "
            "separately; the built-in result is a reference, not ground truth."
        ),
        nbformat.v4.new_markdown_cell(
            "## Context & Methods\n\n"
            f"This companion notebook regenerates **{derived_table_count} derived "
            f"CSVs and {figure_count} static figures** without modifying raw COMSOL "
            "or Swarm outputs. The archived COMSOL profiles are final-time BDF "
            "profiles; a stationarity residual was not exported, so they are not "
            "promoted to independently verified steady states.\n\n"
            "### Key Assumptions\n\n"
            "- Profiles overlap over a common 1D interval and are linearly "
            "interpolated onto the union of their x coordinates.\n"
            "- Integrals of the reconstructed profiles, their squares, and their "
            "products are evaluated exactly on every piecewise-linear interval. "
            "Trapezoidal integration of squared nodal samples is retained only as "
            "a sensitivity diagnostic.\n"
            "- E/N in exported profiles is V m² and is converted using 1 Td = "
            "10⁻²¹ V m².\n"
            "- Central 80% means 10% of physical domain length removed at each end; "
            "it is a geometric window, not a physical bulk/sheath definition.\n"
            "- `excitation_source` and `ionization_source` are eir2 direct "
            "excitation and eir4 direct ionization, not net reaction sources. "
            "`total_current_density` is electron + Ar⁺ conductive current.\n"
            "- E/N > 500 Td is COMSOL's typical drift-diffusion regime guidance, "
            "not a hard validity cutoff.\n"
            "- The Townsend source representation is "
            r"\(S_T=(\alpha/N)N|J_e|/e\). The notebook's "
            r"\(Nn_e k\) reconstruction is a counterfactual postprocess, not a "
            "matched COMSOL rerun and not evidence that one representation is better.\n"
            "- Timing rows are independent only when the run process that produced "
            "them was independently launched; the notebook cannot infer this from "
            "numbers alone."
        ),
        nbformat.v4.new_code_cell(
            "from pathlib import Path\n"
            "import json\n"
            "import pandas as pd\n"
            "from IPython.display import Image, display\n\n"
            "PROJECT_ROOT = Path.cwd().resolve()\n"
            f"CONFIG = PROJECT_ROOT / {config_path.resolve().relative_to(PROJECT_ROOT).as_posix()!r}\n"
            "REPORT_ROOT = PROJECT_ROOT / 'reports/comsol_swarm_benchmark_2026'\n"
            f"DATA_DIR = PROJECT_ROOT / {data_dir!r}\n"
            "assert CONFIG.exists() and DATA_DIR.exists()\n"
            "print(f'Config: {CONFIG.relative_to(PROJECT_ROOT)}')"
        ),
        nbformat.v4.new_markdown_cell(
            "## Data\n\n"
            "### 1. Regenerate analysis products\n\n"
            "The checked-in analysis module is the single implementation of metric "
            "definitions used by this notebook and the report."
        ),
        nbformat.v4.new_code_cell(
            "import importlib.util\n"
            "import sys\n\n"
            "module_path = REPORT_ROOT / 'repro/benchmark_analysis.py'\n"
            "spec = importlib.util.spec_from_file_location('benchmark_analysis_runtime', module_path)\n"
            "analysis = importlib.util.module_from_spec(spec)\n"
            "sys.modules[spec.name] = analysis\n"
            "spec.loader.exec_module(analysis)\n"
            "summary = analysis.run_analysis(CONFIG)\n"
            "{key: summary[key] for key in ['evidence_status', 'data_vintages', 'formal_gates']}"
        ),
        nbformat.v4.new_markdown_cell(
            "### 2. Verify the 19-table contract, source hashes, and data-quality checks"
        ),
        nbformat.v4.new_code_cell(
            "DERIVED_CSV_CONTRACT = [\n"
            "    'aligned_profiles_long.csv', 'clamp_impact.csv', 'clamp_metrics.csv',\n"
            "    'comparison_metrics.csv', 'current_rsd_sensitivity.csv',\n"
            "    'current_uniformity.csv', 'fluid_applicability_diagnostics.csv',\n"
            "    'historical_activation_audit.csv', 'operating_point_diagnostics.csv',\n"
            "    'partial_electron_energy_audit.csv', 'qa_summary.csv',\n"
            "    'regime_impact.csv', 'regime_metrics.csv',\n"
            "    'regional_error_attribution.csv', 'runtime_speedup.csv',\n"
            "    'runtime_summary.csv', 'source_hashes.csv',\n"
            "    'tail_convergence_audit.csv',\n"
            "    'townsend_source_representation_audit.csv',\n"
            "]\n"
            "actual_csvs = sorted(path.name for path in DATA_DIR.glob('*.csv'))\n"
            "assert actual_csvs == sorted(DERIVED_CSV_CONTRACT), (actual_csvs, DERIVED_CSV_CONTRACT)\n"
            "print(f'Derived CSV contract: {len(actual_csvs)}/19 files present')\n\n"
            "source_hashes = pd.read_csv(DATA_DIR / 'source_hashes.csv')\n"
            "qa = pd.read_csv(DATA_DIR / 'qa_summary.csv')\n"
            "display(source_hashes)\n"
            "display(qa.loc[~qa['passed'], ['source_id', 'check', 'observed', 'expected', 'severity_if_failed']])"
        ),
        nbformat.v4.new_markdown_cell(
            "## Results\n\n"
            "### 3. Exact piecewise-linear spatial metrics\n\n"
            "Relative L2 integrates the square of the continuous piecewise-linear "
            "reconstruction exactly. The adjacent signed integral ratio tests net "
            "magnitude, and the correlation tests co-variation after exact spatial "
            "centering. The trapezoidal-of-squared-samples column quantifies numerical "
            "sensitivity, especially for edge-local current spikes."
        ),
        nbformat.v4.new_code_cell(
            "metrics = pd.read_csv(DATA_DIR / 'comparison_metrics.csv')\n"
            "display(metrics[['quantity_label', 'relative_L2_spatial_weighted', "
            "'relative_L2_trapezoidal_of_squared_samples_sensitivity', "
            "'integral_ratio_external_over_reference', "
            "'shape_correlation_spatial_weighted', "
            "'legacy_L2_unweighted_discrete']].round(6))"
        ),
        nbformat.v4.new_markdown_cell(
            "### 4. Lookup coverage and Swarm tail/solver gates\n\n"
            "The LEA equations query active transport and Townsend data by mean "
            "electron energy. Therefore the mean-energy boundary-policy rows are "
            "primary; the low-E/N row is secondary because the E/N-to-mean-energy "
            "pair is inactive in the archived LEA formulation. Tail/grid adequacy "
            "and iterative solver convergence are audited independently."
        ),
        nbformat.v4.new_code_cell(
            "clamp = pd.read_csv(DATA_DIR / 'clamp_metrics.csv')\n"
            "clamp_impact = pd.read_csv(DATA_DIR / 'clamp_impact.csv')\n"
            "tail = pd.read_csv(DATA_DIR / 'tail_convergence_audit.csv')\n"
            "tail_summary = pd.Series({\n"
            "    'formal points': len(tail),\n"
            "    'tail/grid gate pass': int(tail['formal_tail_grid_gate_pass'].sum()),\n"
            "    'solver convergence gate pass': int(tail['formal_solver_convergence_gate_pass'].sum()),\n"
            "    'all-case gate pass': int(tail['formal_swarm_case_gate_pass'].sum()),\n"
            "}, name='count')\n"
            "display(clamp)\n"
            "display(clamp_impact[['region', 'quantity_label', "
            "'absolute_integral_share']].round(6))\n"
            "display(tail_summary.to_frame())\n"
            "display(tail.loc[~tail['formal_swarm_case_gate_pass'], [\n"
            "    'E_over_N_Td', 'solver_converged', 'iterations', 'residual_L1',\n"
            "    'grid_max_eV', 'tail_probability', 'edge_to_peak',\n"
            "    'formal_tail_grid_gate_pass',\n"
            "    'formal_solver_convergence_gate_pass',\n"
            "]])"
        ),
        nbformat.v4.new_markdown_cell(
            "### 5. Operating point, source representation, partial energy, and "
            "fluid-applicability diagnostics\n\n"
            "These are diagnostic decompositions of the archived profiles. The "
            "Townsend/rate comparison is a representation sensitivity only; the "
            "energy table is intentionally partial; and the net-flux-speed ratio is "
            "an operational boundary/locality flag rather than a validity theorem."
        ),
        nbformat.v4.new_code_cell(
            "operating = pd.read_csv(DATA_DIR / 'operating_point_diagnostics.csv')\n"
            "townsend = pd.read_csv(DATA_DIR / 'townsend_source_representation_audit.csv')\n"
            "partial_energy = pd.read_csv(DATA_DIR / 'partial_electron_energy_audit.csv')\n"
            "fluid = pd.read_csv(DATA_DIR / 'fluid_applicability_diagnostics.csv')\n"
            "display(operating[[\n"
            "    'path_label', 'source_voltage_V', 'gap_endpoint_voltage_V',\n"
            "    'ballast_current_A', 'internal_max_potential_V',\n"
            "    'anode_side_potential_overshoot_V',\n"
            "    'high_field_absolute_potential_variation_share',\n"
            "]].round(6))\n"
            "display(townsend[[\n"
            "    'channel_label', 'observed_townsend_source_integral_m2_s',\n"
            "    'townsend_reconstruction_relative_difference',\n"
            "    'counterfactual_rate_source_integral_m2_s',\n"
            "    'counterfactual_rate_over_observed_townsend',\n"
            "]])\n"
            "display(partial_energy[[\n"
            "    'path_label', 'electron_field_power_signed_W_m2',\n"
            "    'eir2_threshold_weighted_loss_W_m2',\n"
            "    'eir4_threshold_weighted_loss_W_m2',\n"
            "    'direct_loss_to_field_power_ratio',\n"
            "]].round(6))\n"
            "display(fluid[[\n"
            "    'path_label', 'max_electron_to_neutral_ratio',\n"
            "    'max_net_flux_speed_ratio',\n"
            "    'spatial_fraction_net_flux_speed_ratio_above_guide',\n"
            "    'spatial_fraction_net_flux_speed_ratio_above_unity',\n"
            "]])"
        ),
        nbformat.v4.new_markdown_cell(
            "### 6. Current, spatial robustness, regime, activation, and runtime"
        ),
        nbformat.v4.new_code_cell(
            "current = pd.read_csv(DATA_DIR / 'current_uniformity.csv')\n"
            "current_sensitivity = pd.read_csv(DATA_DIR / 'current_rsd_sensitivity.csv')\n"
            "regional = pd.read_csv(DATA_DIR / 'regional_error_attribution.csv')\n"
            "regime_metrics = pd.read_csv(DATA_DIR / 'regime_metrics.csv')\n"
            "regime_impact = pd.read_csv(DATA_DIR / 'regime_impact.csv')\n"
            "activation = pd.read_csv(DATA_DIR / 'historical_activation_audit.csv')\n"
            "runtime = pd.read_csv(DATA_DIR / 'runtime_summary.csv')\n"
            "display(current)\n"
            "display(current_sensitivity.head(12))\n"
            "display(regional)\n"
            "display(regime_metrics)\n"
            "display(regime_impact[regime_impact['quantity'].isin("
            "['excitation_source', 'ionization_source'])])\n"
            "display(activation)\n"
            "display(runtime)"
        ),
        nbformat.v4.new_markdown_cell(
            "### 7. Publication figures\n\n"
            "Each PNG also has a vector PDF peer. Figure contracts, exact source "
            "hashes, output hashes, and evidence labels are in "
            "`data/derived/figure_metadata.json`."
        ),
        nbformat.v4.new_code_cell(
            "for figure_id in [\n"
            "    'fig01_workflow_boundary', 'fig02_lookup_domain',\n"
            "    'fig03_profile_state', 'fig04_sources_currents_mesh',\n"
            "    'fig05_metric_summary', 'fig06_runtime_breakdown',\n"
            "    'fig07_spatial_robustness', 'fig08_activation_regime',\n"
            "    'fig09_fluid_applicability']:\n"
            "    print(figure_id)\n"
            "    display(Image(filename=str(REPORT_ROOT / 'figures' / f'{figure_id}.png'), width=850))"
        ),
        nbformat.v4.new_markdown_cell(
            "## Takeaways\n\n"
            f"- The present dataset remains **{evidence_status}**; failed formal "
            "gates are visible above and must not be hidden in publication claims.\n"
            f"- Tail/grid adequacy is **{tail_grid_pass_count}/{tail_case_count}**, "
            f"but solver convergence is only **{solver_pass_count}/{tail_case_count}**; "
            "the lookup bundle is tail-audited, not a fully converged formal input.\n"
            "- Shape correlation, weighted relative L2, and integral ratio are "
            "not interchangeable; conclusions must identify which metric supports them.\n"
            "- In LEA, active lookup coverage must be assessed on the mean-energy "
            "axis; the explicit low-energy floor is a boundary policy, not a direct "
            "Swarm result.\n"
            "- The >500 Td analysis separates lookup coverage from the physical "
            "regime risk of drift-diffusion/local closures.\n"
            "- The source-form and partial-energy audits identify sensitivity and "
            "missing balance terms; neither is a matched rerun or a complete balance.\n"
            "- The activation audit is historical: unresolved Java/class provenance "
            "prevents it from becoming a formal functional-activation result.\n"
            "- Replace paths and evidence status in `repro/analysis_inputs.json` "
            "after clean runs; no analysis code change is required."
        ),
    ]
    notebook_path.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, notebook_path)
    client = NotebookClient(
        notebook,
        timeout=600,
        kernel_name=kernel_name,
        resources={"metadata": {"path": str(PROJECT_ROOT)}},
    )
    client.execute()
    nbformat.write(notebook, notebook_path)
    return notebook_path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    path = build_notebook(args.config, args.output)
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
