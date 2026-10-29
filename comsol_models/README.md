# COMSOL models and Swarm inputs

## ダウンロード

MPH、保存済みSwarm SQLiteと大きなEEDF表はGit LFSで管理する。Git LFSをインストールした環境で:

```powershell
git lfs install
git clone --branch codex/swarm-0.12-configured-mph https://github.com/dkawa2523/swarm.git
cd swarm
git lfs pull
git lfs fsck
```

MPHの実体と入力CSVを取得するため、GitHubのソースZIPではなく上の手順を使う。
CSV/JSON/YAMLのバイト列は改行変換から保護し、記録済みSHA-256を保つ。


1件ずつ取得する場合は、このページのMPHリンクを使う。リンクはGit LFSの実体を返す
ダウンロード先に固定している。`raw.githubusercontent.com` のURLから取得したものはMPHの拡張子でも
数百バイトの参照テキストの場合がある。先頭が
`version https://git-lfs.github.com/spec/v1`ならモデル本体ではない。
例えばGEC-CCPのMPH実体は **24,506,900バイト**。Rawから返る参照は133バイト。

### GUIで開くと停止・クラッシュする場合

まず取得したファイルのサイズと`SHA-256`を、元モデルは`sources.json`、保存済みICPは
各`archive.json`と照合する。配布元と別フォルダへGitHubから取得した7件について、
バイト列の一致とMPH内部全メンバーのCRCを確認した。これはGUIでの読込・再計算確認とは別である。

次にCOMSOLの完全なバージョン番号と、MPH内に記録された保存時のビルドを確認する。
設定済みICPは **6.4.0.429** で保存されている。こちらのPCには6.4.0.378と6.4.0.429の
インストールが混在しているため、ファイル関連付けと実際に起動する実行ファイルの確認も必要。
ビルド差だけでクラッシュ原因と断定しない。

COMSOLが起動してからモデル表示時に落ちる場合、通常のDesktop起動コマンドに
`-3drend sw`を付けてソフトウェア描画で切り分けられる。
これは[COMSOL公式の描画トラブル対策](https://www.comsol.com/support/knowledgebase/933)であり、
MPHや物理設定の変更ではない。起動時のライセンス期限切れは、MPH読込時のクラッシュと区別する。
問題のファイル名、サイズ、COMSOLのビルド、停止時のメッセージまたは該当ログを残して調査する。

GitHubのZIPにLFS実体が含まれるかはリポジトリ設定に依存する。
[GitHub公式説明](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/managing-repository-settings/managing-git-lfs-objects-in-archives-of-your-repository)を参照。

## 保存済みの設定済みMPH

次の2件は2026-09-13 UTCの計算済みICPモデルを変更せず保存したもの。
最終solve出力のSHA-256、内包EEDFと元CSVの一致を確認した。
`archive.json`にMPH・元入力・当時の計算証跡のハッシュをまとめている。

| Swarm入力 | 開くMPH | 大きさ | 保存時のCOMSOL | 状態 |
| --- | --- | --- | --- | --- |
| two_term | [ICP two-term](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_models/saved/gec_icp_two_term/argon_gec_icp_two_term_function_eedf.mph) | 29.9 MB | 6.4.0.429 | 過去の計算済みモデル |
| monte_carlo | [ICP Monte Carlo](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_models/saved/gec_icp_monte_carlo/argon_gec_icp_monte_carlo_function_eedf.mph) | 31.8 MB | 6.4.0.429 | 過去の計算済みモデル |

**これらは0.12.0で再設定・再計算したモデルではない。** 当時の受入結果と、現在の互換性・物理受入は分けて記録している。
MPHにはインポート済みEEDFと解が保存されている。元の絶対パス表示は残しているため、
外部ファイルを再インポートする際は同梱の`bundle/`へ参照先を指定し直す。
保存済みモデルの別環境での起動・再計算は、今回のライセンス期限切れのため未確認。
同梱のformat-v1 bundleは当時の証跡であり、現行アダプタに入力しない。

正柱の旧`Model/work/positive_column_1d_swarm_template.mph`は参照CSVが存在しないため、
配布する設定済みMPHには含めていない。CCPの設定済みMPHは見つからなかった。

## 0.12.0で新しく計算した入力

次のディレクトリには、現行コードで収束した点だけを保存したSQLite、
format-v2 COMSOL bundle、単一Swarm YAMLとモデル別マッピングがある。
`preparation.json`がファイルのハッシュ・入力範囲・停止理由を記録する。

| ディレクトリ | 点数 | E/N [Td] | 現在の状態 |
| --- | --- | --- | --- |
| [positive_column](positive_column/) | 26 | 0.0444–500 | オフライン準備通過、MPHへの反映未完了 |
| [gec_ccp](gec_ccp/) | 47 | 0.05–500 | bundle生成済み、アダプタの入力ハッシュ連鎖検査で停止 |
| [gec_icp](gec_icp/) | 10 | 1–300 | 元の作業場所で準備通過。別フォルダでは断面積のパス依存ハッシュで停止 |

**これらのディレクトリには設定済みMPHはない。** COMSOLの起動が
`The license has expired / License error: -10`で停止した。
生成Javaのコンパイル成功はMPHへの反映成功を意味しない。
保存済みモデルをZIP/XML編集で書き換える処理は行っていない。

これらは範囲を限定した入力例であり、従来のモデルの全運転範囲を資格付けするデータではない。
元の広範囲計算は未収束で停止し、追跡したICPでは1000 Tdで600反復の上限に達した。
未収束値の採用や古いデータへの置換は行っていない。範囲外の物理利用を認めるものではない。
粒子・組成の時間発展、RF場と空間方程式は引き続きCOMSOLが所有する。

有効なライセンスと必要な入力検査の修正後、モデル別に実行する。
出力は`outputs/configured_models/<model>/`に書かれ、上の過去のMPHを上書きしない。
生成Javaの絶対パスは元の準備時点の記録なので、別の場所では`--dry-run`で再生成する。

```powershell
py -3 -m swarm_workflow.cli run-positive-column comsol_models/positive_column/mapping.yaml --dry-run
py -3 -m swarm_workflow.cli run-gec-ccp comsol_models/gec_ccp/mapping.yaml --dry-run
py -3 -m swarm_workflow.cli run-gec-icp comsol_models/gec_icp/mapping.yaml --dry-run
```

CCPの上記コマンドは現時点で明示的に停止する。ICPもGitHubから別フォルダへ取得すると停止する。
ICPの断面積ファイル自体のSHA-256は元と同じだが、アダプタが比較する複合ハッシュに
絶対パスが含まれるためである。`gec_icp/evidence/relocation_preflight.json`に検証結果を記録した。
検査を外したり、bundleの由来を書き換えて進めたりしない。
`--dry-run`を外すと設定だけでなくモデルの計算も行うため、入力範囲・運転条件を確認して実行する。
Swarmを再計算する場合は各ディレクトリの`swarm.yaml`を`swarm-workflow run`に渡す。
保存済み`swarm_results.sqlite`は元結果の読み取り用コピーであり、別の場所への移動に合わせて改変していない。
bundleは既存の`build-tables`、`export-comsol`で新しい出力ディレクトリへ再生成する。

## 元モデル

次の5件もGit LFSで配布する。既存マッピングが参照するパスを維持している。

- [positive_column_1d.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/Model/positive_column_1d.mph): 正柱アダプタの元モデル。
- [argon_gec_ccp.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_modes/argon_gec_ccp.mph): CCPアダプタの元モデル。
- [argon_gec_icp.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/comsol_modes/argon_gec_icp.mph): ICPアダプタの元モデル。
- [positive_column_1d_boltzmann.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/Model/positive_column_1d_boltzmann.mph): COMSOL組込Boltzmann参照。
- [boltzmann_argon.mph](https://media.githubusercontent.com/media/dkawa2523/swarm/ecfd83fd3032da8cb6be6ef46933b246edef9d10/Model/boltzmann_argon.mph): COMSOLの空間一様Boltzmann参照。

元モデルの存在は外部Swarm入力の設定済みを意味しない。
詳細な責任分担と受入条件は[モデル別アダプタ](../docs/comsol_model_adapters.md)を参照。
