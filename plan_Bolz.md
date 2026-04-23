既存 repo への統合は、**最初は既存 MC コードをほぼ触らず、Boltzmann solver を新しい独立 package として横に追加する**のが安全です。その後、MC の入出力だけを adapter に接続し、最後に共通 core へ徐々に寄せる流れがよいです。

実装 ZIP はこれです。

[swarm_boltzmann_production_bolsig.zip](sandbox:/mnt/data/swarm_boltzmann_production_bolsig.zip)

## 1. 推奨統合方針

今回の実装は、既存 repo に対して以下を追加する構成です。

```text
electron_swarm/                      # 新しい共通 package
  core/                              # config / cross section / result model
  solvers/
    boltzmann_two_term.py            # BOLSIG+ 相当 backend
    monte_carlo_adapter.py           # 既存 MC との接続口
  io/
  plotting.py
  runner.py

configs/
  boltzmann_only.yaml
  both_template.yaml

examples/cross_sections/
  argon_minimal.csv

docs/
  boltzmann_two_term_design.md
```

既存の粒子 Monte Carlo 実装は、最初は直接書き換えず、`monte_carlo_adapter.py` から

```yaml
monte_carlo:
  command: ...
```

または

```yaml
monte_carlo:
  python_api: "module:function"
```

で呼ぶ形にします。
この方針なら、既存 MC solver の計算ロジックを壊さずに、Boltzmann solver と出力統合だけを先に検証できます。

## 2. まず統合作業用ブランチを切る

既存 repo で作業します。

```bash
cd /path/to/swarm

git status
git switch update3-1
git pull
git switch -c feature/unified-mc-boltzmann
```

`git pull` が使えない環境なら、現在の `update3-1` のローカル状態からそのままブランチを切って構いません。

## 3. ZIP を展開する

```bash
mkdir -p /tmp/swarm_boltzmann_overlay
unzip /path/to/swarm_boltzmann_production_bolsig.zip -d /tmp/swarm_boltzmann_overlay
cd /tmp/swarm_boltzmann_overlay/swarm_boltzmann_impl
```

ZIP の中の top directory は `swarm_boltzmann_impl/` です。

## 4. いきなり `--force` で上書きしない

ZIP 内には統合スクリプトがあります。

```bash
python tools/integrate_overlay.py /path/to/swarm
```

ただし、このスクリプトは保守的に「既存ファイルがあれば止まる」作りです。最初はこれで問題ありません。

```bash
python tools/integrate_overlay.py /path/to/swarm
```

もし

```text
Refusing to overwrite ...
```

と出た場合、すぐに `--force` は使わず、次の手動コピー方式にしてください。

## 5. 安全な手動コピー方式

既存 repo 側に `electron_swarm/` や `configs/` が既にある場合は、こちらを推奨します。

```bash
OVERLAY=/tmp/swarm_boltzmann_overlay/swarm_boltzmann_impl
REPO=/path/to/swarm

rsync -aivn \
  --exclude='__pycache__' \
  --exclude='*.pyc' \
  "$OVERLAY/electron_swarm/" \
  "$REPO/electron_swarm/"
```

`-n` は dry-run です。表示された差分を確認して問題なければ実行します。

```bash
rsync -aiv \
  --exclude='__pycache__' \
  --exclude='*.pyc' \
  "$OVERLAY/electron_swarm/" \
  "$REPO/electron_swarm/"
```

設定・例・ドキュメントもコピーします。

```bash
mkdir -p "$REPO/configs" "$REPO/examples/cross_sections" "$REPO/docs" "$REPO/tests" "$REPO/tools"

cp "$OVERLAY/configs/boltzmann_only.yaml" "$REPO/configs/boltzmann_only.yaml"
cp "$OVERLAY/configs/both_template.yaml" "$REPO/configs/both_template.yaml"
cp "$OVERLAY/examples/cross_sections/argon_minimal.csv" "$REPO/examples/cross_sections/argon_minimal.csv"
cp "$OVERLAY/docs/boltzmann_two_term_design.md" "$REPO/docs/boltzmann_two_term_design.md"

cp "$OVERLAY/tests/test_boltzmann_two_term.py" "$REPO/tests/test_boltzmann_two_term.py"
cp "$OVERLAY/tools/validate_against_bolos.py" "$REPO/tools/validate_against_bolos.py"
```

不要な `__pycache__` が混ざった場合は削除します。

```bash
cd "$REPO"
find electron_swarm -type d -name '__pycache__' -prune -exec rm -rf {} +
find electron_swarm -name '*.pyc' -delete
```

## 6. 依存関係を既存 repo に追加する

今回の追加 package が必要とする依存関係は以下です。

```text
numpy>=1.23
scipy>=1.10
pandas>=1.5
PyYAML>=6.0
matplotlib>=3.6
```

既存 repo が `requirements.txt` を使っているなら追加します。

```text
numpy>=1.23
scipy>=1.10
pandas>=1.5
PyYAML>=6.0
matplotlib>=3.6
```

既存 repo が `pyproject.toml` を使っているなら、`[project] dependencies` に追加します。

```toml
dependencies = [
  "numpy>=1.23",
  "scipy>=1.10",
  "pandas>=1.5",
  "PyYAML>=6.0",
  "matplotlib>=3.6",
]
```

BOLOS との optional 比較検証も行うなら、任意で追加します。

```text
bolos>=0.2
pytest>=7
```

開発環境では以下で入れます。

```bash
cd /path/to/swarm
python -m pip install -e .
python -m pip install scipy pandas PyYAML matplotlib pytest
```

既存 repo が package 化されていない場合は、ひとまず以下でも動作確認できます。

```bash
cd /path/to/swarm
python -m pip install numpy scipy pandas PyYAML matplotlib pytest
```

## 7. まず Boltzmann 単独で動作確認する

最初に MC とは接続せず、Boltzmann solver だけを動かしてください。

```bash
cd /path/to/swarm
python -m electron_swarm configs/boltzmann_only.yaml
```

成功すると、設定上は以下に出力されます。

```text
outputs/argon_boltzmann/argon_summary.csv
outputs/argon_boltzmann/argon_eedf.csv
outputs/argon_boltzmann/argon_rates.csv
outputs/argon_boltzmann/argon_mean_energy.png
outputs/argon_boltzmann/argon_drift_velocity.png
outputs/argon_boltzmann/argon_reduced_mobility.png
outputs/argon_boltzmann/argon_reduced_diffusion.png
outputs/argon_boltzmann/argon_eedf.png
```

この段階で見るべき列は `argon_summary.csv` の以下です。

```text
solver
E_over_N_Td
mean_energy_eV
drift_velocity_m_s
reduced_mobility_m2_V_s_m3
diffusion_L_m2_s
diffusion_T_m2_s
net_ionization_frequency_s
effective_townsend_m2
meta_converged
meta_residual_L1
meta_tail_probability
meta_edge_to_peak
meta_grid_max_eV
```

特に本番用途では、まず以下を確認してください。

```text
meta_converged == True
meta_tail_probability が十分小さい
meta_edge_to_peak が十分小さい
mean_energy_eV が異常に高すぎない
EEDF tail が energy grid 上端で切れていない
```

## 8. 単体テストを通す

```bash
cd /path/to/swarm
python -m pytest tests/test_boltzmann_two_term.py -q
```

このテストは、Boltzmann solver が最低限以下を満たすことを確認します。

```text
YAML config が読める
断面積が読める
BOLSIG-like backend が収束する
summary / EEDF / rates が生成される
```

## 9. 既存 MC と接続する

ここが最も重要です。

`configs/both_template.yaml` のこの部分を、既存 repo の MC 実行方法に合わせて変更します。

```yaml
monte_carlo:
  enabled: true
  command: "python main.py --config {config}"
  working_directory: ".."
  timeout_s: 3600
  output_summary_csv: "../outputs/mc_summary.csv"
  output_eedf_csv: "../outputs/mc_eedf.csv"
```

`{config}` には現在の YAML path が渡されます。
`{output_dir}` には `output.directory` が渡されます。

例えば、既存 MC の実行が

```bash
python run.py input.yaml
```

なら、

```yaml
monte_carlo:
  enabled: true
  command: "python run.py {config}"
  working_directory: ".."
  timeout_s: 3600
  output_summary_csv: "../outputs/mc_summary.csv"
  output_eedf_csv: "../outputs/mc_eedf.csv"
```

のようにします。

既存 MC が既存形式の設定ファイルしか受け付けない場合は、`command` に wrapper script を挟むのが安全です。

```yaml
monte_carlo:
  enabled: true
  command: "python tools/run_existing_mc_from_unified_yaml.py {config}"
  working_directory: ".."
  timeout_s: 3600
  output_summary_csv: "../outputs/mc_summary.csv"
  output_eedf_csv: "../outputs/mc_eedf.csv"
```

この wrapper で、統一 YAML から既存 MC 用 input に変換し、既存 MC を実行し、最後に `mc_summary.csv` と `mc_eedf.csv` を出す形にします。

## 10. MC 側が出すべき summary CSV

`monte_carlo_adapter.py` が最低限要求する summary CSV は、`E_over_N_Td` 列です。

推奨 schema は以下です。

```csv
case_id,E_over_N_Td,mean_energy_eV,drift_velocity_m_s,mobility_m2_V_s,reduced_mobility_m2_V_s_m3,diffusion_L_m2_s,diffusion_T_m2_s,reduced_diffusion_L_m2_s_m3,reduced_diffusion_T_m2_s_m3,net_ionization_frequency_s,effective_townsend_m2
Ar_0000,20,1.23,10000,0.1,1.0e-21,0.2,0.2,2.0e-21,2.0e-21,0.0,0.0
Ar_0001,50,2.34,30000,0.2,2.0e-21,0.4,0.4,4.0e-21,4.0e-21,1.0e5,3.0e-20
```

最低限はこれでも動きます。

```csv
E_over_N_Td,mean_energy_eV,drift_velocity_m_s
20,1.23,10000
50,2.34,30000
```

ただし、Boltzmann と比較するには、できるだけ同じ列を出すべきです。

## 11. MC 側が出すべき EEDF CSV

EEDF も統合したい場合、MC 側は以下のような CSV を出してください。

```csv
case_id,E_over_N_Td,energy_eV,eedf
Ar_0000,20,0.001,0.12
Ar_0000,20,0.01,0.34
Ar_0000,20,0.1,0.56
Ar_0001,50,0.001,0.09
Ar_0001,50,0.01,0.21
Ar_0001,50,0.1,0.44
```

`case_id` がない場合は、`E_over_N_Td` で対応付けます。
`eedf` は

```text
∫ f(ε) dε = 1
```

になるように規格化しておくのが望ましいです。adapter 側でも一応再規格化します。

## 12. `both` mode を動かす

`configs/both_template.yaml` を編集したら実行します。

```bash
cd /path/to/swarm
python -m electron_swarm configs/both_template.yaml
```

成功すると、1つの出力ファイルに MC と Boltzmann の結果がまとまります。

```text
outputs/argon_both/argon_summary.csv
outputs/argon_both/argon_eedf.csv
outputs/argon_both/argon_rates.csv
```

`argon_summary.csv` は例えば次のような構造になります。

```csv
solver,case_id,E_over_N_Td,mean_energy_eV,drift_velocity_m_s,...
boltzmann_two_term,Ar_0000,20,...
boltzmann_two_term,Ar_0001,50,...
monte_carlo,Ar_0000,20,...
monte_carlo,Ar_0001,50,...
```

つまり、既存 MC と Boltzmann の違いは `solver` 列で識別します。
グラフも同じ `solver` 列で分けて描画されます。

## 13. より堅い接続は `python_api` 方式

長期的には、`command` で外部プロセスとして呼ぶより、既存 MC solver に Python API を1つ生やす方が保守しやすいです。

YAML はこうします。

```yaml
monte_carlo:
  enabled: true
  python_api: "existing_mc.unified_adapter:run_swarm"
  passthrough:
    n_particles: 100000
    random_seed: 1
```

既存 repo 側に、例えば以下を追加します。

```python
# existing_mc/unified_adapter.py

from electron_swarm.core.results import SwarmCaseResult

def run_swarm(config, cross_sections, **kwargs):
    """
    config: electron_swarm.core.config.SwarmConfig
    cross_sections: electron_swarm.core.cross_sections.CrossSectionSet
    kwargs: YAML の monte_carlo.passthrough
    """

    # 1. config から既存 MC 用入力を作る
    # 2. 既存 MC solver を呼ぶ
    # 3. 結果を SwarmCaseResult の list に変換して返す

    results: list[SwarmCaseResult] = []

    for case in config.run.e_over_n_Td:
        # 既存 MC 結果から値を詰める
        result = SwarmCaseResult(
            solver="monte_carlo",
            case_id=f"{config.run.case_prefix}_{case:g}Td",
            e_over_n_Td=case,
            mean_energy_eV=...,
            drift_velocity_m_s=...,
            mobility_m2_V_s=...,
            reduced_mobility_m2_V_s_m3=...,
            diffusion_L_m2_s=...,
            diffusion_T_m2_s=...,
            reduced_diffusion_L_m2_s_m3=...,
            reduced_diffusion_T_m2_s_m3=...,
            net_ionization_frequency_s=...,
            effective_townsend_m2=...,
            energy_eV=...,
            eedf=...,
            eepf=...,
            rates=[],
            metadata={
                "n_particles": kwargs.get("n_particles"),
                "random_seed": kwargs.get("random_seed"),
            },
        )
        results.append(result)

    return results
```

この方式の利点は、MC 出力ファイルの列名変換や一時ファイルに依存しなくなることです。

## 14. 既存 MC input と統一 YAML の対応

統一 YAML では、物理条件は solver 共通部に置きます。

```yaml
run:
  mode: both
  e_over_n_Td: [20, 50, 100, 200]

conditions:
  gas_temperature_K: 300.0
  pressure_Pa: 13.3322
  gas_mixture:
    - species: Ar
      fraction: 1.0
      mass_amu: 39.948

cross_sections:
  format: csv
  files:
    - path: ../examples/cross_sections/argon_minimal.csv
      species: Ar
      format: csv
```

MC 固有設定は `monte_carlo:` に閉じ込めます。

```yaml
monte_carlo:
  enabled: true
  python_api: "existing_mc.unified_adapter:run_swarm"
  passthrough:
    n_particles: 100000
    max_collisions: 1000000
    time_step_s: 1.0e-12
    random_seed: 1
```

Boltzmann 固有設定は `boltzmann_two_term:` に閉じ込めます。

```yaml
boltzmann_two_term:
  enabled: true
  backend: native_bolsig
  energy_grid:
    min_eV: 0.0001
    max_eV: 80.0
    n: 600
    spacing: quadratic
  adaptive_grid:
    enabled: true
    max_cycles: 4
    mean_energy_multiplier: 15.0
    tail_probability: 1.0e-8
  convergence:
    max_iterations: 120
    tolerance: 1.0e-8
    residual_tolerance: 1.0e-7
  nonconservative_model: growth
  ionization_energy_sharing: equal
```

既存 MC 側に粒子数、乱数 seed、サンプリング時間、初期電子温度などがある場合、それらは `monte_carlo.passthrough` に入れるのがよいです。
共通物理条件と solver 固有数値条件を混ぜないことが、後々の保守性に効きます。

## 15. 断面積ファイルの統合

今回の Boltzmann solver は、標準では CSV を読みます。

長形式 CSV は以下です。

```csv
species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2
Ar,Ar momentum transfer,momentum,0.0,39.948,0.1,2.0e-20
Ar,Ar excitation 11.55 eV,excitation,11.55,39.948,12.0,1.0e-21
Ar,Ar ionization,ionization,15.76,39.948,20.0,2.0e-21
```

対応 type は以下です。

```text
momentum
elastic
effective
excitation
ionization
attachment
superelastic
```

既存 MC が別形式の断面積を使っている場合、最初は MC 側の読み込みはそのままにして、Boltzmann 用に CSV 変換したファイルを用意する方が安全です。
その後、共通 `electron_swarm.core.cross_sections` を MC からも読む形に寄せるとよいです。

## 16. repo 構成としてはこうするのが理想

最終的には、既存 repo を次のような構成に寄せるのがよいです。

```text
swarm/
  electron_swarm/
    core/
      config.py
      constants.py
      cross_sections.py
      results.py

    solvers/
      base.py
      boltzmann_two_term.py
      monte_carlo_adapter.py
      monte_carlo_native.py      # 将来的に既存 MC をここへ移す、または薄い wrapper にする

    io/
      writers.py

    plotting.py
    runner.py

  existing_mc/                  # 既存 MC 実装。最初はそのまま
    ...

  configs/
    boltzmann_only.yaml
    both_template.yaml
    mc_only.yaml

  examples/
    cross_sections/

  tests/
    test_boltzmann_two_term.py
    test_mc_adapter.py
    test_both_mode.py

  docs/
    boltzmann_two_term_design.md
```

重要なのは、MC と Boltzmann を同じ階層に無理やり混ぜるのではなく、

```text
core        : 物理条件、断面積、結果 schema
solvers     : solver 固有処理
io          : 出力
plotting    : 可視化
runner      : orchestration
```

に分けることです。

## 17. 統合後に確認すべき `git diff`

コピー後、必ず差分を確認してください。

```bash
cd /path/to/swarm
git status
git diff --stat
git diff -- electron_swarm/core/config.py
git diff -- configs/both_template.yaml
```

期待される追加は主に以下です。

```text
new file: electron_swarm/...
new file: configs/boltzmann_only.yaml
new file: configs/both_template.yaml
new file: examples/cross_sections/argon_minimal.csv
new file: docs/boltzmann_two_term_design.md
new file: tests/test_boltzmann_two_term.py
new file: tools/validate_against_bolos.py
```

既存 MC の主要ファイルに大きな変更が入っていたら、いったん戻した方が安全です。

```bash
git checkout -- path/to/existing_mc_file.py
```

## 18. CI に追加するテスト

GitHub Actions や手元 CI があるなら、まず以下だけ追加します。

```bash
python -m pytest tests/test_boltzmann_two_term.py -q
python -m electron_swarm configs/boltzmann_only.yaml --no-write
```

`--no-write` は solver 実行だけを確認し、CSV/PNG を書かないモードです。

さらに MC 接続後に、以下のような統合テストを追加するとよいです。

```bash
python -m electron_swarm configs/both_template.yaml
python - <<'PY'
import pandas as pd
df = pd.read_csv("outputs/argon_both/argon_summary.csv")
assert {"boltzmann_two_term", "monte_carlo"} <= set(df["solver"])
assert df["E_over_N_Td"].nunique() >= 1
print(df.head())
PY
```

## 19. 最小コミット単位

破綻を避けるため、1回の巨大 commit にしない方がよいです。

推奨 commit 分割は以下です。

```text
commit 1: add electron_swarm core/result/config/cross_section modules
commit 2: add native_bolsig two-term Boltzmann solver
commit 3: add unified CSV/plot output
commit 4: add YAML examples and docs
commit 5: add MC adapter with command/python_api bridge
commit 6: connect existing MC output to adapter
commit 7: add both-mode regression tests
```

この分け方にすると、Boltzmann solver の問題、YAML schema の問題、既存 MC 接続の問題を切り分けやすくなります。

## 20. 実際に最初にやるべきコマンドまとめ

一番安全な流れだけまとめると、次です。

```bash
# 既存 repo
cd /path/to/swarm
git switch update3-1
git switch -c feature/unified-mc-boltzmann

# ZIP 展開
mkdir -p /tmp/swarm_boltzmann_overlay
unzip /path/to/swarm_boltzmann_production_bolsig.zip -d /tmp/swarm_boltzmann_overlay

# コピー
OVERLAY=/tmp/swarm_boltzmann_overlay/swarm_boltzmann_impl
REPO=/path/to/swarm

rsync -aiv \
  --exclude='__pycache__' \
  --exclude='*.pyc' \
  "$OVERLAY/electron_swarm/" \
  "$REPO/electron_swarm/"

mkdir -p "$REPO/configs" "$REPO/examples/cross_sections" "$REPO/docs" "$REPO/tests" "$REPO/tools"

cp "$OVERLAY/configs/boltzmann_only.yaml" "$REPO/configs/boltzmann_only.yaml"
cp "$OVERLAY/configs/both_template.yaml" "$REPO/configs/both_template.yaml"
cp "$OVERLAY/examples/cross_sections/argon_minimal.csv" "$REPO/examples/cross_sections/argon_minimal.csv"
cp "$OVERLAY/docs/boltzmann_two_term_design.md" "$REPO/docs/boltzmann_two_term_design.md"
cp "$OVERLAY/tests/test_boltzmann_two_term.py" "$REPO/tests/test_boltzmann_two_term.py"
cp "$OVERLAY/tools/validate_against_bolos.py" "$REPO/tools/validate_against_bolos.py"

# 依存関係
cd "$REPO"
python -m pip install numpy scipy pandas PyYAML matplotlib pytest

# Boltzmann 単独確認
python -m pytest tests/test_boltzmann_two_term.py -q
python -m electron_swarm configs/boltzmann_only.yaml

# 差分確認
git status
git diff --stat
```

その後に `configs/both_template.yaml` の `monte_carlo:` を既存 MC の実行方法に合わせて修正し、

```bash
python -m electron_swarm configs/both_template.yaml
```

を実行します。

この順番で進めれば、既存 MC を壊さずに、まず BOLSIG+ 相当の Boltzmann 二項近似 solver を repo に追加し、その後に MC と Boltzmann の同時実行・同一形式出力へ接続できます。
