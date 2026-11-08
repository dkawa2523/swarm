# COMSOL モデルと Swarm 入力

## 設定済みMPH

現行の `two_term` で計算した純Arの入力を組み込み、COMSOL **6.4.0.429** で下記の条件の受入検査を通過しました。
元の入力CSVを一時的に移動し、別フォルダに置いたMPHをCOMSOL Desktopから開いて再計算しました。

| モデル | 開くMPH | 受入条件 | GUIで選ぶスタディ |
| --- | --- | --- | --- |
| 正柱 | [positive_column_two_term.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/10e05e76287ce403cd749b9356fbbd341d03b118/comsol_models/positive_column/positive_column_two_term.mph) | 200 V、13.3322 Pa、293.15 K | **Study 1**（8秒） |
| GEC-CCP | [argon_gec_ccp_two_term.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/10e05e76287ce403cd749b9356fbbd341d03b118/comsol_models/gec_ccp/argon_gec_ccp_two_term.mph) | 1 W、300 K | **Time Periodic Study2**（235秒）、次に **Time Periodic to Time Dependent**（3秒） |
| GEC-ICP | [argon_gec_icp_two_term.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/10e05e76287ce403cd749b9356fbbd341d03b118/comsol_models/gec_icp/argon_gec_icp_two_term.mph) | 1500 W、13.56 MHz、2.66644 Pa、300 K、0.01 s | **Study 1**（56秒） |

MPHを開き、表のスタディをモデルビルダーで選んで **Compute** を押します。
CCPの外部Swarm用周期スタディのタグは `std3` です。変換 `std2` で電極波形を物理時刻で表示できます。
CCPの補間関数はEEDFのCSV見出し1行を読み飛ばす警告を出しました。

各 `release.json` に配布ファイルのバイト数とSHA-256、運転条件、GUI結果、受入記録をまとめています。
正柱は局所エネルギー近似でCOMSOLが平均エネルギーを解き、Swarmの輸送係数とTownsend係数を使います。
導電電流の一様性の受入範囲は正柱の中央80%です。電極近傍には400要素に増やしても
電流密度の振動が残りました。[正柱の比較結果](../docs/dev/improvement_plan.md#i4と独立照合で確認した結果)
CCPとICPはSwarmのFunction EEDFと移動度を使い、COMSOLが拡散、エネルギー、RF場と重粒子を扱います。
受入は表の条件と各モデルの閉包に対する判定です。入力表の範囲と物理的な利用範囲はそれぞれ確認してください。

## ダウンロード

MPH、SQLite、大きなEEDFと位相データはGit LFSで管理します。

```powershell
git lfs install
git clone --branch codex/swarm-0.12-configured-mph https://github.com/dkawa2523/swarm.git
cd swarm
git lfs pull
git lfs fsck
```

個別に取得する場合はMPHのダウンロードリンクを使います。`raw.githubusercontent.com` が
`version https://git-lfs.github.com/spec/v1` で始まる短いファイルを返した場合はLFSの参照テキストです。
サイズとSHA-256を `release.json`、元モデルは `sources.json`、過去のICPは `archive.json` で照合できます。
GUIが停止する場合は、取得ファイルのサイズ、COMSOLの完全なビルド番号、エラー表示を確認します。
描画で落ちる場合はソフトウェア描画で切り分けられます。[COMSOLの描画トラブル対策](https://www.comsol.com/support/knowledgebase/933)

## 入力と解析用出力

モデルごとのフォルダに、Swarmの設定 `swarm.yaml`、計算済みSQLite、COMSOL用 `bundle/`、
モデルの設定 `mapping.yaml`、保存済み解と解析データを収録しました。

| フォルダ | SwarmのE/N点 | 解析用データ |
| --- | --- | --- |
| [positive_column](positive_column/) | 45点、0.0002–5000 Td | 最終保存解の801点。電子密度、平均エネルギー、電位、電流、生成項。EEDFと輸送の図 |
| [gec_ccp](gec_ccp/) | 54点、0.05–4000 Td | 軸方向・半径方向・領域の周期平均、位相分布、電極波形、EEDFと空間分布の図 |
| [gec_icp](gec_icp/) | 14点、1–3000 Td | 0.001 sと0.01 sの空間分布、26時刻の電子・Ar*・Ar+粒子数、平均エネルギー・電力の図 |

`swarm-workflow run` に `swarm.yaml` を渡すと、出力先の新しいSQLiteに計算します。
`build-tables` と `export-comsol` で新しい表を作れます。配布のSQLiteは元の計算を確認するための読み取り用です。
各 `run-*` の `--dry-run` は実行前の入力検査とJava生成を行います。

```powershell
py -3 -m swarm_workflow.cli run-positive-column comsol_models/positive_column/mapping.yaml --pressure-Pa 13.3322 --dry-run
py -3 -m swarm_workflow.cli run-gec-ccp comsol_models/gec_ccp/mapping.yaml --dry-run
py -3 -m swarm_workflow.cli run-gec-icp comsol_models/gec_icp/mapping.yaml --dry-run
```

`export-positive-column` と `export-gec-icp` は保存解を読み、指定した結果を新しい出力フォルダに保存します。
`plot-gec-ccp` は `bundle/` と `results/` から図を作れます。描画には `plot` オプション依存関係を使います。
`compare-gec-icp` は移動先の `results/` から共通時刻の比較と時系列図を作れます。

```powershell
py -3 -m pip install ".[plot]"
py -3 -m swarm_workflow.cli plot-gec-ccp --bundle comsol_models/gec_ccp/bundle --comsol-results comsol_models/gec_ccp/results --output outputs/gec_ccp_plots
```

組成の時間発展はプラズマモデルが解きます。ここでのICPの時系列はCOMSOLの保存解です。
配布済みSQLiteは[外部モデルでの係数利用](../docs/coefficient_lookup.md)にも使えます。
具体的な設定とコマンドは[正柱](../docs/comsol_workflow.md)、[CCP](../docs/comsol_gec_ccp.md)、[ICP](../docs/comsol_gec_icp.md)を参照してください。

## 過去の計算済みICP

2026-09-13 UTCにCOMSOL 6.4.0.429で保存したモデルです。format-v1の元入力は当時の証跡として同梱しています。
元のCSVパスが表示されるため、表を再インポートするときは同梱の `bundle/` を指定してください。

- [ICP two-term](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_models/saved/gec_icp_two_term/argon_gec_icp_two_term_function_eedf.mph)
- [ICP Monte Carlo](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_models/saved/gec_icp_monte_carlo/argon_gec_icp_monte_carlo_function_eedf.mph)

## 元モデル

- [positive_column_1d.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/Model/positive_column_1d.mph) 正柱の元モデル
- [argon_gec_ccp.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_modes/argon_gec_ccp.mph) CCPの元モデル
- [argon_gec_icp.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_modes/argon_gec_icp.mph) ICPの元モデル
- [positive_column_1d_boltzmann.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/Model/positive_column_1d_boltzmann.mph) 正柱のCOMSOL組込Boltzmann参照
- [boltzmann_argon.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/Model/boltzmann_argon.mph) 空間一様のCOMSOL組込Boltzmann参照
