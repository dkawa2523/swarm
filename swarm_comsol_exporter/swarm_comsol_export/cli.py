from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict

from .config import ComsolExportConfig, load_yaml
from .converter import export_to_comsol


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="swarm-comsol-export",
        description="Convert swarm/Boltzmann outputs to COMSOL Plasma Module interpolation tables.",
    )
    p.add_argument("--config", required=True, help="Path to YAML config (examples/comsol_export.yaml)")
    p.add_argument("--run-dir", default=None, help="Override run_dir in YAML")
    p.add_argument("--dry-run", action="store_true", help="Only validate inputs; do not write outputs")
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    cfg_dict: Dict[str, Any] = load_yaml(args.config)
    cfg = ComsolExportConfig.from_dict(cfg_dict, run_dir_override=args.run_dir)

    if not cfg.enabled:
        print("[swarm-comsol-export] comsol_export.enabled=false; nothing to do.")
        return 0

    out_dir = cfg.output_dir or (cfg.run_dir / "comsol")

    if args.dry_run:
        # Try to load and validate without writing
        print("[swarm-comsol-export] Dry run: validating inputs only.")
        export_to_comsol(
            run_dir=cfg.run_dir,
            output_dir=out_dir,
            input_cfg=cfg.input_cfg,
            processing=cfg.processing,
            output_cfg={**cfg.output, "transport_filename": "_dry_transport.csv", "rates_filename": "_dry_rates.csv", "eedf_filename": "_dry_eedf.csv", "report_filename": "_dry_report.json"},
        )
        # remove created outputs
        for f in out_dir.glob("_dry_*"):
            try:
                f.unlink()
            except Exception:
                pass
        print("[swarm-comsol-export] OK.")
        return 0

    paths = export_to_comsol(
        run_dir=cfg.run_dir,
        output_dir=out_dir,
        input_cfg=cfg.input_cfg,
        processing=cfg.processing,
        output_cfg=cfg.output,
    )
    print("[swarm-comsol-export] Export completed:")
    print(f"  - transport: {paths.transport_csv}")
    print(f"  - rates:     {paths.rates_csv}")
    print(f"  - eedf:      {paths.eedf_csv}")
    print(f"  - report:    {paths.report_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
