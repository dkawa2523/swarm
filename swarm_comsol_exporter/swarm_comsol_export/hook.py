from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Union

from .config import ComsolExportConfig, load_yaml
from .converter import export_to_comsol


def maybe_export_comsol(cfg: Union[str, Path, Dict[str, Any]], run_dir: Optional[Union[str, Path]] = None) -> Optional[Path]:
    """Swarm側から呼ぶための最小フック。

    Parameters
    ----------
    cfg:
        - YAMLパス (str/Path) または dict
    run_dir:
        - Swarm出力ディレクトリ（未指定なら cfg から）
    Returns
    -------
    output_dir (Path) or None
    """
    if isinstance(cfg, (str, Path)):
        cfg_dict = load_yaml(cfg)
    else:
        cfg_dict = cfg

    c = ComsolExportConfig.from_dict(cfg_dict, run_dir_override=run_dir)

    if not c.enabled:
        return None

    out_dir = c.output_dir or (c.run_dir / "comsol")
    export_to_comsol(
        run_dir=c.run_dir,
        output_dir=out_dir,
        input_cfg=c.input_cfg,
        processing=c.processing,
        output_cfg=c.output,
    )
    return out_dir
