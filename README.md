# swarm

この repo には 2 つの実行系があります。

- 既存 Monte Carlo workflow: `run_sweep.py` + `configs/*.yaml`
- unified workflow: `python -m electron_swarm` + `configs/unified/*.yaml`

`electron_swarm` は既存 `swarm_mc` を残したまま、Monte Carlo と Boltzmann two-term を同じ結果 schema で扱うための追加レイヤです。

## Setup

`uv` を使う場合:

```powershell
uv venv
uv sync
```

`uv` を使わない場合は `py -3 -m ...` でも実行できます。

## Workflow 1: Existing Monte Carlo

```powershell
uv run python run_sweep.py --config configs\Ar_N2_sweep.yaml
```

主な出力は `outputs\<experiment.name>\` 以下の `summary.csv`、各 run の CSV、`execution.log`、`execution_times.csv` です。

## Workflow 2: Unified electron_swarm

Boltzmann only:

```powershell
uv run python -m electron_swarm configs\unified\boltzmann_only.yaml
```

Monte Carlo + Boltzmann:

```powershell
uv run python -m electron_swarm configs\unified\both_template.yaml
```

write を切る場合:

```powershell
py -3 -m electron_swarm configs\unified\boltzmann_only.yaml --no-write
```

### Unified configs

- `configs/unified/boltzmann_only.yaml`
- `configs/unified/both_template.yaml`
- `configs/unified/ar_n2_both.yaml`

`both_template.yaml` と `ar_n2_both.yaml` では Monte Carlo 側を `swarm_mc.unified_adapter:run_swarm` に接続しています。既存 `swarm_mc` 用の詳細設定は `monte_carlo.passthrough.base_config` にそのまま渡せます。

### Unified outputs

`electron_swarm` は unified 出力と legacy 互換出力を同時に書き出します。

Unified files:

- `<base>_summary.csv`
- `<base>_eedf.csv`
- `<base>_rates.csv`
- comparison plots

Legacy compatibility files:

- `summary_mc.csv`
- `summary_boltzmann.csv`
- `eedf_table_mc.csv`
- `eedf_table_boltzmann.csv`
- `energy_table_mc.csv`
- `energy_table_boltzmann.csv`

Legacy aliases:

- `summary.csv`
- `eedf_table.csv`
- `energy_table.csv`

`summary.csv` / `eedf_table.csv` / `energy_table.csv` は `output.compatibility.primary_solver` が指す solver を alias します。既定は `monte_carlo` です。

`<base>_summary.csv` の `net_ionization_frequency_s` は、Monte Carlo と Boltzmann の両方で比較できるように convolution/integral effective rate を基準にしています。Monte Carlo の counted rate と Boltzmann の growth frequency は `meta_*` 列に残します。

### Cross sections

`electron_swarm` は canonical CSV に加えて、既存 repo が使っていた `txt|lxcat|bolsig` 形式も読めます。これにより `cross_sections/Ar_Biagi.txt` などを変換せずに unified runner から使えます。

## COMSOL export

既存 `run_sweep.py` の `comsol_export.enabled: true` はそのまま使えます。

```powershell
uv run python run_sweep.py --config configs\example_visual_jit_comsol.yaml --plot --save-plots
```

`swarm_comsol_exporter` は unified 出力の `summary_mc.csv` / `summary_boltzmann.csv` / `eedf_table_mc.csv` / `eedf_table_boltzmann.csv` も fallback として認識します。

## EEDF fitting

`EEDF_optimze.py` は従来の `eedf_table.csv` に加えて、`eedf_table_mc.csv` / `eedf_table_boltzmann.csv` にも対応しています。combined output を使う場合は `--solver` で対象 solver を絞れます。

```powershell
uv run python EEDF_optimze.py --table outputs\unified\argon_both\eedf_table.csv --solver monte_carlo
```

## Docs

- `docs/boltzmann_two_term_design.md`
- `docs/comsol.md`
- `docs/README_optimize.md`
- `docs/plan_fast2.md`
- `plan_Bolz.md`
