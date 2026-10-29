# Swarmの設定と実行

Swarm用の設定は一つのYAMLで完結する。断面積データはそのYAMLから参照する。
計算をCOMSOLへ接続するときだけ、別のCOMSOLモデル設定を使う。
基本設定ファイルへの参照、scenario、preset、YAML継承は使わない。

```powershell
py -3 -m swarm_workflow.cli run examples\two_term.yaml --plan
py -3 -m swarm_workflow.cli run examples\two_term.yaml
```

`--plan`は計算も出力ディレクトリの作成も行わず、物理要求の扱い、ケース数、保存先を表示する。
`--plan --json`で実効ソルバー設定、物理条件、sampling単位、worker数とメモリ予算を確認できる。
ケース数は再開による既存ケースの省略前の要求数である。

## 通常計算

既存の`schema_version: 2`と四つの正規solver idを使う。
`run.solvers`は計算手法、`physics`は衝突・場・散乱などの物理要求を指定する。

```yaml
schema_version: 2
run:
  solvers:
    - id: two_term
  e_over_n_Td: [10, 30, 100]
conditions:
  gas_temperature_K: 300
  pressure_Pa: 13.3
  gas_mixture:
    - {species: Ar, fraction: 1.0, mass_amu: 39.948}
cross_sections:
  files:
    - {path: cross_sections/argon_minimal.csv, species: Ar}
output:
  directory: ../outputs/argon
  base_name: swarm
```

この例の断面積パスは設定YAMLの場所を基準にする。
実際の最小例は[こちら](../examples/two_term.yaml)。外部断面積は`files`に複数指定できる。
CSVとLXCat/BOLSIG形式の断面積テキストを受け付ける。外部ソルバーの結果ファイルを
断面積として指定する入口ではない。

通常計算は正規のsummary・rates・EEDF・plan CSVとmanifestを出力する。
`electron-swarm <config.yaml>`と`python -m electron_swarm <config.yaml>`も、単一混合比の
直接計算用として使える。campaignを含む設定は共通入口の`swarm-workflow run`を使う。

## 保存・再開と独立反復

同じYAMLに`campaign`を加える。E/N・気体条件・seedを書き直す必要はない。

```yaml
campaign:
  database: ../outputs/argon/cases.sqlite
```

同じコマンドを再実行すると、同じ計算・入力データ・実装に属する既存ケースを再利用し、
不足分だけを実行する。条件や計算実装が変わった場合は、新しいDB名を指定する。
古いDBの由来を書き換えて再利用しない。

MCでは`solvers.monte_carlo.seed`を直接計算とcampaignの共通のmaster seedとして指定する。
`campaign.mc`は反復・worker・必要なら点ごとのsampling予算を設定する。

```yaml
solvers:
  monte_carlo:
    seed: 20261025
    particles: 256
    production_trial_periods: 1024
campaign:
  database: ../outputs/argon/mc.sqlite
  mc:
    replicas: 4
    workers: 2
```

この部分例では、`run.solvers`にも`monte_carlo`を選ぶ。試行数は受理衝突数ではなく
null-collision trial period数である。seed未指定は計画時に失敗する。
各反復・E/N・混合比の独立seedと、各solverケースの実効seedを記録する。
同じmaster seedでも、直接計算のstreamとcampaignの反復0が同一という意味ではない。
例の予算は入力方法を示すもので、統計的な資格を保証しない。

必要な場合だけ`campaign.mc.sampling_plan`、`campaign.mc.convergence`、
`campaign.quality`、`campaign.execution.deterministic`を追加する。
通常は`run.e_over_n_Td`が全solverの計算点になる。
MCだけ別の点を計算する場合に限り`campaign.mc.e_over_n_Td`を明示できる。
その範囲とケース数も計画に表示される。
sampling予算の数値は一混合比あたりで、MCを実行する混合比数を掛けたものが総量になる。
平均エネルギー範囲の継続を設定した場合、初期計画点に加え、指定した最大step数まで追加点があり得る。

campaignではsolver比較CSVの生成は未実装である。`comparison.enabled: true`は計画時に拒否する。
比較出力が必要な直接計算では[比較例](../examples/compare_three_solvers.yaml)を使う。
campaign結果はまずSQLiteへ保存し、必要なら`aggregate`、`build-tables`へ渡す。

## 複数の混合比

`conditions.gas_mixture`のfractionを、同じ長さの配列にする。
同じ添字が一つの混合比であり、成分ごとの直積ではない。

```yaml
conditions:
  gas_temperature_K: 300
  pressure_Pa: 13.3
  gas_mixture:
    - {species: Ar, mass_amu: 39.948, fraction: [1.0, 0.99]}
    - {species: O2, mass_amu: 31.998, fraction: [0.0, 0.01]}
```

この例は純ArとAr/O2=99/1の二つを表す。全成分で配列を使い、各混合比の合計を1にする。
不揃い・空・重複の混合比、暗黙の重み正規化は受け付けない。
複数混合比では`campaign.database`を指定する。E/Nとの組合せを同じ計画器で列挙し、保存・再開する。
実行できる完全な例は[workflow_ar_o2.yaml](../examples/workflow_ar_o2.yaml)。

## 出力と図

0.10のSwarm用`output`は`directory`、`base_name`と任意の`export`、`figures`を受け付ける。
`export`の最初の対応は、指定した一つのsolverのDC反応速度係数である。
Swarmの図は同じYAMLで選ぶ。空リストまたは省略時は作図しない。
任意空間観測点・補間時刻は、プラズマモデル側の設定になる。
通常計算はsummary/rates/EEDFをCSVへ出し、対応する結果には角度分布とRF位相CSVも加える。
manifestと一緒に読むと、ケースの条件、近似、欠損、単位・量の定義を保持できる。

```python
from electron_swarm.io import read_product_run

# examples/two_term.yamlを実行した後
result = read_product_run("outputs/examples/two_term_manifest.json")
summary = result.tables["summary_csv"]
rates = result.tables["rates_csv"]
eedf = result.tables["eedf_csv"]
```

campaignも、solverが返すRF位相・位相別EEDFと角度分布を保存する。
RFの必要な子データが欠けたケースは完了扱いにせず、再実行でそのケースを計算する。
DBから一つの原結果を型付きで読む入口は次の通り。

```python
from swarm_workflow.campaign.store import WorkflowStore

with WorkflowStore("outputs/argon/cases.sqlite", readonly=True) as store:
    case = store.read_case(
        mixture_id=0, solver="two_term", e_over_n_Td=30.0, replicate=0,
    )
# case.rf_phase、case.energy_angle_distributionなどは、solverが返した場合だけ存在する。
```

読込結果は保存データの復元であり、新しいsolver実行の証拠ではない。
旧保存形式は黙って変換せず、新しいDBへ再計算する。
EEDF比較の資格評価、GEC-CCPの空間・RF図、GEC-ICP保存解の時間履歴には専用ツールがある。
positive-columnは接続先mapの`results.fields`に出力式と単位を指定し、
`dataset`と`solutions: all`または`last`を明示する。
CSVはdataset・保存解番号・点番号・時刻・解パラメータを保持する。
任意座標・補間時刻の指定は未対応として拒否する。
組成スキャンを組成の時間発展として扱わない。化学種の時系列は利用側プラズマモデルの結果になる。
Swarmの設定に化学反応網・濃度初期値・化学種時間積分を加える予定はない。
既存のRF位相別EEDFは電子応答として扱い、プラズマ結果の時系列図は接続先の後処理へ置く。
現状の対応範囲・再現した不足・改訂案は[出力仕様レビュー](dev/modeling_output_review.md)を参照する。

図を使う場合は`py -3 -m pip install -e ".[plot]"`を実行し、同じ設定に次を加える。

```yaml
output:
  directory: ../outputs/argon
  base_name: swarm
  figures: [eedf, transport, reaction_rates]
```

| 名前 | 図の内容 |
| --- | --- |
| `eedf` | 同一条件のsolverを比較するEEDF・EEPFと、全エネルギー域のtail |
| `transport` | E/Nに対する平均エネルギー、drift、利用可能な還元移動度・拡散・エネルギー輸送 |
| `reaction_rates` | DCの反応別`k [m^3/s]`。係数exportのsolver選択とは独立して全計算solverを表示 |
| `energy_angle` | propagatorが返したエネルギー・角度セル密度 |
| `rf_phase` | time-dependent two-termが返した位相別平均エネルギー・場・報告済み電離係数・EEDF |

図は`<base>_figure_<kind>_<番号>.png/svg`、曲線の数値は`<base>_figure_data.csv`となる。
単位、元ケース、表示範囲、近似・数値状態は同じmanifestの`figures`に記録する。
2次元の密度は元の角度/RF CSVを参照し、同じ配列を複製しない。
組成・温度・密度・物理条件が違う系列は分ける。組成配列を時間軸にしない。
輸送量はquantityの定義が違えば別パネルにし、未計算量を0で埋めない。
e-e後処理の古い輸送量は`transport stale`、近似は`approximate`と明記する。

MCの信頼区間は独立seed・反復番号を確認できる保存済み反復から求める。
一反復・全ゼロの区間は空欄。異なるEEDF格子は共通範囲へ確率質量を保存して再配置し、
範囲外の質量と負の微小値を保持する。表示のための再正規化・正値への置換は行わない。
図の帯は点ごとのStudent-t 95%区間であり、同時区間や入力資格を意味しない。
図に必要なデータがない明示要求はエラーとし、P1の拡散欠損を理由にEEDF/rates図を止めない。

保存済みのmanifestまたは互換schemaのSQLiteから、再計算せず描き直せる。

```powershell
swarm-workflow plot outputs/argon/cases.sqlite --figures eedf transport reaction_rates --output outputs/argon/figures
swarm-workflow plot outputs/argon/figures/swarm_manifest.json --figures eedf --output outputs/argon/figures
```

```python
from electron_swarm.core.output import OutputConfig
from electron_swarm.io import read_figures
from swarm_workflow.figures import plot_results

paths = plot_results("outputs/argon/cases.sqlite", OutputConfig(
    "outputs/argon/figures", figures=("eedf", "transport"),
))
data = read_figures(paths["manifest_json"])
curves = data.tables["figure_data_csv"]
```

同じ出力先・base名への再描画もでき、削除した画像を復元する。
数値CSVのhashは常に検証する。通常のreaderは画像hashも検証し、
数値だけを消費する`read_product_run(..., verify_figures=False)`では画像の検証を省ける。
旧実装のDBを解析しても、その由来を現在の実装に更新しない。
DB起点の出力は`swarm.figures.v1`/`saved_results_analysis`であり、新しいsolver実行ではない。
直接計算の再描画では元の数値と実行identityを維持する。図の生成版は別に記録する。
campaignに`figures`と`export`の両方がある場合は、保存した同じsnapshotから一世代として出力する。

## 結果を他のモデルへ渡す

必要な係数は同じSwarm YAMLで指定できる。

```yaml
output:
  directory: ../outputs/argon
  base_name: swarm
  export:
    solver: multi_term
    quantities: [reaction_rates]
```

このsolverは`run.solvers`でも有効にする。四solverすべてを同じ形式で扱い、
別solverへの代替は行わない。未対応量、RF RMS/位相、実行できない選択は計画時に拒否する。
元のsummary/rates/EEDFは維持し、選択した反応係数を`swarm_reaction_rates.csv`へ追加する。
直接実行とcampaignのどちらも`swarm_manifest.json`に結び付く。

```python
from electron_swarm.io import read_reaction_rates

data = read_reaction_rates("outputs/argon/swarm_manifest.json")
rates = data.table
contexts = data.metadata["series"]
```

保存済み結果から、計算をやり直さずに作成することもできる。

```powershell
swarm-workflow export-reaction-rates outputs/argon/swarm_manifest.json --solver multi_term --output outputs/argon/coefficients
swarm-workflow export-reaction-rates outputs/argon/cases.sqlite --solver monte_carlo --output outputs/argon/mc_coefficients
```

反応のspecies/process/type/thresholdと、固定組成・温度・密度・場・角度仮定・e-e条件を保持する。
系列ごとにE/Nの実計算点と範囲を示し、補間・外挿・平均エネルギー軸への変換は行わない。
`series_id`はそのartifact内の識別子であり、別artifactとの条件照合には物理contextを使う。
`k [m^3/s]`は対象粒子あたり、`mixture_weighted_rate_m3_s`は`x k`。
二体反応源`n_e n_target k`に対し、組成を二重に掛けない。
化学量論や化学種ごとの正味生成・消滅源は、利用側の反応網から与える。
輸送量、driftによるTownsend換算、平均エネルギーの単調性は、この出力の条件にしない。
e-e後処理の更新済みratesも出力し、元のtransport-stale情報は系列のcase metadataに保持する。

MCは同じ条件の独立反復を等重みで平均し、反復数、標準誤差、Student-tの95%区間を出す。
一反復や全ゼロ観測の区間は空欄で、全ゼロは`censored_all_zero`とする。
計画したMC反復の欠落、同じseedの重複、反応inventoryの不一致は拒否する。
数値状態と不確かさを付けた係数点データであり、MCの精度判定・独立物理検証・
接続先モデルの受入が完了した入力という意味ではない。
同じ出力先・base名は一つの世代として置換する。通常実行からcampaignへ切り替えたときも
古いmanifest・CSVを残さない。保存済みデータからのexportでは入力元の世代を上書きしない。
利用者が編集・提供したデータは`derived`/`provided`とし、元の収束主張を引き継がない。

Swarm計算とCOMSOLの形式・MPHへの対応・空間モデル実行は別工程である。
COMSOL側の温度・組成は入力との適合性を確認するもので、Swarm計算条件を上書きしない。
モデル設定内のbundleを使うときはCLIの`--bundle`再指定は不要。明示した差替えは実効mappingへ記録する。
手書き設定の相対パスはそのYAMLの場所基準、bundle内参照はbundleの場所基準である。

このCSV/APIの読込はCOMSOLを必要としない。COMSOL形式への変換と、
非COMSOLモデルでの範囲内係数取得の利用例は次の段階に残る。
現在のCOMSOLテーブルをそのまま汎用の化学反応入力とは扱わない。
速度係数`k`と組成重み付き`x k`、E/Nと平均エネルギー、未計算量とゼロを区別する。
既存COMSOLテーブルの輸送・統計・補間要件は維持している。

計算由来は`calculation_sha256`を使い、出力先やcampaignの品質評価条件と分けて記録する。
table/bundleはformat v2、Propagatorのtarget資格はv2である。古い成果物を新しい意味へ読み替えない。
数値モデルの資格、数値収束、独立した物理検証、接続先モデルの受入は別の判定である。
