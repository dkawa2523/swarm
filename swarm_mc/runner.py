# Copyright (c) 2020-2021 ETH Zurich
"""
Batch-oriented runner that drives swarm_mc simulations from YAML configs.

The goal is to make multi-point E/N sweeps easy to configure and reproducible while
keeping the low-level Monte Carlo kernels unchanged.
"""

from __future__ import annotations

import os
import time
from copy import deepcopy
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import yaml

from swarm_mc.caching import GasMixtureCache
from swarm_mc.config import Config
from swarm_mc.simulation import Simulation


def _arange_inclusive(start: float, stop: float, step: float) -> List[float]:
    """
    Inclusive arange with a small tolerance to capture the stop value.
    """
    if step == 0:
        raise ValueError("step cannot be zero for sweep definitions")
    count = int(np.floor((stop - start) / step)) + 1
    return [start + i * step for i in range(count)]


@dataclass
class SweepDefinition:
    param: str
    values: List[float]
    label_template: str = "{param}_{value}"

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SweepDefinition":
        """
        Parse sweep definitions.

        Supported shapes:
        sweep:
          param: EN
          start: 20
          stop: 200
          step: 20

        sweep:
          EN:
            start: 20
            stop: 200
            step: 20
        """
        if not data:
            raise ValueError("sweep definition is empty")

        label_template = data.get("label_template", "{param}_{value}")
        if "param" in data:
            param = data["param"]
            values = cls._extract_values(data)
        else:
            if len(data) != 1:
                raise ValueError("sweep must define exactly one parameter")
            param, payload = next(iter(data.items()))
            label_template = payload.get("label_template", label_template)
            values = cls._extract_values(payload)
        return cls(param=param, values=values, label_template=label_template)

    @staticmethod
    def _extract_values(data: Dict[str, Any]) -> List[float]:
        if "values" in data and data["values"] is not None:
            return list(data["values"])
        if {"start", "stop", "step"} <= data.keys():
            return _arange_inclusive(float(data["start"]),
                                     float(data["stop"]),
                                     float(data["step"]))
        raise ValueError("sweep requires either values or start/stop/step fields")


@dataclass
class PlannedRun:
    config: Dict[str, Any]
    label: str
    sweep_param: Optional[str] = None
    sweep_value: Optional[float] = None


@dataclass
class ExperimentConfig:
    base_config: Dict[str, Any]
    name: str
    output_root: Path
    sweep: Optional[SweepDefinition] = None
    base_name: Optional[str] = None
    output_format: Optional[str] = None
    plotting: Optional[Dict[str, Any]] = None
    parallel_workers: int = 1

    def build_runs(self) -> List[PlannedRun]:
        runs: List[PlannedRun] = []
        base_cfg = deepcopy(self.base_config)
        base_name = self.base_name or base_cfg["output"]["base_name"]
        if self.output_format:
            base_cfg["output"]["output_format"] = self.output_format

        if not self.sweep:
            cfg = self._prepare_output(deepcopy(base_cfg), base_name)
            runs.append(PlannedRun(cfg, label=base_name))
            return runs

        for value in self.sweep.values:
            cfg = deepcopy(base_cfg)
            self._apply_param(cfg, self.sweep.param, value)
            label = self.sweep.label_template.format(param=self.sweep.param,
                                                     value=value)
            run_name = f"{base_name}_{label}"
            cfg = self._prepare_output(cfg, run_name)
            runs.append(PlannedRun(cfg, label=run_name,
                                   sweep_param=self.sweep.param,
                                   sweep_value=value))
        return runs

    def _prepare_output(self, cfg: Dict[str, Any], run_name: str) -> Dict[str, Any]:
        root = self.output_root / self.name
        run_dir = root / run_name
        run_dir.mkdir(parents=True, exist_ok=True)
        cfg["output"]["base_name"] = run_name
        cfg["output"]["output_directory"] = str(run_dir) + os.sep
        if self.output_format:
            cfg["output"]["output_format"] = self.output_format
        return cfg

    @staticmethod
    def _apply_param(cfg: Dict[str, Any], param: str, value: Any) -> None:
        """
        Apply a parameter change. Supports dot-notation (category.key) or plain keys.
        """
        if "." in param:
            category, key = param.split(".", 1)
            if category not in cfg or key not in cfg[category]:
                raise ValueError(f"cannot apply sweep for parameter '{param}'")
            cfg[category][key] = value
            return

        for category, payload in cfg.items():
            if isinstance(payload, dict) and param in payload:
                payload[param] = value
                return
        raise ValueError(f"{param} is not a valid configuration parameter name.")


class SimulationRunner:
    """
    Drive one or more simulations from a YAML (or json/json5) experiment description.
    """

    def __init__(self, experiment: ExperimentConfig,
                 gas_mixture_cache: Optional[GasMixtureCache] = None):
        self.experiment = experiment
        self.gas_mixture_cache = gas_mixture_cache or GasMixtureCache()
        self.run_stats: List[Dict[str, Any]] = []
        self.total_runtime_s: Optional[float] = None
        self.log_path: Optional[Path] = None

    @classmethod
    def from_yaml(cls, path: str) -> "SimulationRunner":
        with open(path, "r") as fh:
            data = yaml.safe_load(fh)

        if "base_config" in data:
            base_config = data["base_config"]
        else:
            base_config = {k: v for k, v in data.items()
                           if k not in ("experiment", "sweep")}
        experiment_meta = data.get("experiment", {})
        sweep = SweepDefinition.from_dict(data["sweep"]) if "sweep" in data else None
        base_name = experiment_meta.get("base_name",
                                        base_config["output"]["base_name"])
        output_root = Path(experiment_meta.get("output_root",
                                               base_config["output"]
                                               .get("output_directory", "results")))
        output_format = experiment_meta.get("output_format",
                                            base_config["output"].get("output_format"))
        plotting = experiment_meta.get("plotting", {})
        parallel_workers = int(experiment_meta.get("parallel_workers", 1))

        exp = ExperimentConfig(
            base_config=base_config,
            name=experiment_meta.get("name", base_name),
            output_root=output_root,
            sweep=sweep,
            base_name=base_name,
            output_format=output_format,
            plotting=plotting,
            parallel_workers=parallel_workers,
        )
        return cls(exp)

    def run_all(self) -> List[Simulation]:
        runs = self.experiment.build_runs()
        self.run_stats = []
        self.total_runtime_s = None
        log_dir = self.experiment.output_root / self.experiment.name
        if log_dir.exists():
            shutil.rmtree(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        self.log_path = log_dir / "execution.log"
        self.log_path.write_text(f"Execution log for {self.experiment.name}\n")

        def log(message: str) -> None:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line = f"[{timestamp}] {message}"
            print(line)
            with open(self.log_path, "a") as log_file:
                log_file.write(line + "\n")

        log(f"Starting sweep with {len(runs)} run(s)"
            f"{f' in parallel (workers={self.experiment.parallel_workers})' if self.experiment.parallel_workers > 1 else ''}")
        total_start = time.perf_counter()
        simulations: List[Simulation] = []
        summary_rows: List[Dict[str, Any]] = []
        eedf_rows: List[Dict[str, Any]] = []
        eepf_rows: List[Dict[str, Any]] = []
        workers = max(1, int(self.experiment.parallel_workers))

        if workers == 1:
            for planned in runs:
                run_en = planned.config.get("physical_conditions", {}).get("EN")
                if run_en is None and "EN" in planned.config:
                    run_en = planned.config["EN"]
                log(f"Starting run '{planned.label}'"
                    f"{f' (E/N={run_en} Td)' if run_en is not None else ''}")
                run_start = time.perf_counter()
                sim = Simulation(planned.config,
                                 gas_mixture_cache=self.gas_mixture_cache,
                                 run_label=planned.label)
                if sim.seed_used is not None:
                    log(f"Using RNG seed {sim.seed_used} (source={sim.seed_source})")
                sim.run()
                duration = time.perf_counter() - run_start
                self.run_stats.append({
                    "run_label": planned.label,
                    "sweep_param": planned.sweep_param,
                    "sweep_value": planned.sweep_value,
                    "E/N (Td)": run_en,
                    "run_time (s)": duration,
                })
                log(f"Finished run '{planned.label}' in {duration:.3f} s")
                simulations.append(sim)

                summary = None
                if sim.output and hasattr(sim.output, "summary_row"):
                    summary = sim.output.summary_row()
                    if summary:
                        summary["run_label"] = planned.label
                        summary["sweep_param"] = planned.sweep_param
                        summary["sweep_value"] = planned.sweep_value
                        summary["run_time (s)"] = duration
                        summary_rows.append(summary)

                run_dir = Path(planned.config["output"]["output_directory"])
                base = planned.config["output"]["base_name"]
                ed_path = run_dir / f"{base}_energy_distribution.csv"
                if ed_path.exists():
                    ed_df = pd.read_csv(ed_path)
                    mean_energy_val = summary.get("mean energy (eV)", np.nan) if summary else np.nan
                    for _, row in ed_df.iterrows():
                        common = {
                            "run_label": planned.label,
                            "sweep_param": planned.sweep_param,
                            "sweep_value": planned.sweep_value,
                            "E/N (Td)": run_en,
                            "mean energy (eV)": mean_energy_val,
                            "energy (eV)": row["energy (eV)"],
                        }
                        eedf_rows.append({
                            **common,
                            "eedf (eV-1)": row["eedf (eV-1)"],
                        })
                        eepf_rows.append({
                            **common,
                            "eepf (eV-3/2)": row["eepf (eV-3/2)"],
                        })
        else:
            futures = {}
            with ProcessPoolExecutor(max_workers=workers) as executor:
                for planned in runs:
                    run_en = planned.config.get("physical_conditions", {}).get("EN")
                    if run_en is None and "EN" in planned.config:
                        run_en = planned.config["EN"]
                    log(f"Queue run '{planned.label}'"
                        f"{f' (E/N={run_en} Td)' if run_en is not None else ''}")
                    fut = executor.submit(_run_planned, planned)
                    futures[fut] = planned

                for fut in as_completed(futures):
                    planned = futures[fut]
                    result = fut.result()
                    duration = result["run_time (s)"]
                    run_en = result.get("E/N (Td)")
                    log(f"Finished run '{planned.label}' in {duration:.3f} s")
                    self.run_stats.append({
                        "run_label": planned.label,
                        "sweep_param": planned.sweep_param,
                        "sweep_value": planned.sweep_value,
                        "E/N (Td)": run_en,
                        "run_time (s)": duration,
                    })
                    summary = result.get("summary")
                    if summary:
                        summary_rows.append(summary)
                    eedf_rows.extend(result.get("eedf_rows", []))
                    eepf_rows.extend(result.get("eepf_rows", []))

        self.total_runtime_s = time.perf_counter() - total_start

        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            if "E/N (Td)" in summary_df.columns:
                summary_df = summary_df.sort_values(by="E/N (Td)")
            summary_dir = self.experiment.output_root / self.experiment.name
            summary_dir.mkdir(parents=True, exist_ok=True)
            summary_df.to_csv(summary_dir / "summary.csv", index=False)
        if eedf_rows:
            eedf_df = pd.DataFrame(eedf_rows)
            eedf_df.to_csv(log_dir / "energy_table.csv", index=False)
            eedf_df[["mean energy (eV)", "energy (eV)", "eedf (eV-1)"]].to_csv(
                log_dir / "eedf_table.csv", index=False
            )
        if eepf_rows:
            eepf_df = pd.DataFrame(eepf_rows)
            eepf_df[["mean energy (eV)", "energy (eV)", "eepf (eV-3/2)"]].to_csv(
                log_dir / "eepf_table.csv", index=False
            )
        if self.run_stats:
            timing_df = pd.DataFrame(self.run_stats)
            if self.total_runtime_s is not None:
                timing_df = pd.concat(
                    [timing_df, pd.DataFrame([{
                        "run_label": "TOTAL",
                        "sweep_param": None,
                        "sweep_value": None,
                        "E/N (Td)": None,
                        "run_time (s)": self.total_runtime_s,
                    }])],
                    ignore_index=True,
                )
            timing_df.to_csv(log_dir / "execution_times.csv", index=False)
        log(f"Completed sweep in {self.total_runtime_s:.3f} s")

        return simulations


def _run_planned(planned: PlannedRun) -> Dict[str, Any]:
    """
    Execute a single planned run in a separate process.
    Returns summary, timing, and aggregated EEDF/EEPF rows.
    """
    run_en = planned.config.get("physical_conditions", {}).get("EN")
    if run_en is None and "EN" in planned.config:
        run_en = planned.config["EN"]

    run_start = time.perf_counter()
    sim = Simulation(planned.config, gas_mixture_cache=None, run_label=planned.label)
    sim.run()
    duration = time.perf_counter() - run_start

    summary = None
    if sim.output and hasattr(sim.output, "summary_row"):
        summary = sim.output.summary_row()
        if summary:
            summary["run_label"] = planned.label
            summary["sweep_param"] = planned.sweep_param
            summary["sweep_value"] = planned.sweep_value
            summary["run_time (s)"] = duration
            summary["E/N (Td)"] = run_en

    eedf_rows: List[Dict[str, Any]] = []
    eepf_rows: List[Dict[str, Any]] = []
    run_dir = Path(planned.config["output"]["output_directory"])
    base = planned.config["output"]["base_name"]
    ed_path = run_dir / f"{base}_energy_distribution.csv"
    mean_energy_val = summary.get("mean energy (eV)", np.nan) if summary else np.nan
    if ed_path.exists():
        ed_df = pd.read_csv(ed_path)
        for _, row in ed_df.iterrows():
            common = {
                "run_label": planned.label,
                "sweep_param": planned.sweep_param,
                "sweep_value": planned.sweep_value,
                "E/N (Td)": run_en,
                "mean energy (eV)": mean_energy_val,
                "energy (eV)": row["energy (eV)"],
            }
            eedf_rows.append({
                **common,
                "eedf (eV-1)": row["eedf (eV-1)"],
            })
            eepf_rows.append({
                **common,
                "eepf (eV-3/2)": row["eepf (eV-3/2)"],
            })

    return {
        "run_time (s)": duration,
        "E/N (Td)": run_en,
        "summary": summary,
        "eedf_rows": eedf_rows,
        "eepf_rows": eepf_rows,
    }
