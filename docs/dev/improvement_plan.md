# Solver comparison improvement plan

2026-10-29、0.12の実装・モデル配布位置。schema v2と四つの正規solverを維持する。
各段階は狭い範囲を設定・計画・実行・出力・文書・振る舞いテストまで完成させる。

## 現在位置

Git LFSで元MPH 5件と過去の設定済みICP 2件を配布する準備を行った。
現行コードの範囲限定入力83点とbundleも[配布モデル一覧](../../comsol_models/README.md)にまとめた。
新しいMPHの保存はCOMSOLライセンス期限切れで未完了。過去のモデルを現行版の成功に扱っていない。
広範囲two-term計算の未収束（追跡したICPは1000 Tdで600反復上限）と、
現行two-term bundleに対するCCPの入力SHA-256連鎖検査不整合も見つかった。
次のI3b2では既存の残件に加え、この入力連鎖と収束範囲を確認し、ライセンス復旧後に
現行入力をMPHへ反映・再読込する。今回のオフライン準備はモデル受入完了ではない。

| 段階 | 範囲 | 状態 |
| --- | --- | --- |
| Phase 1--7h | 物理モデル、結果の意味、入出力・由来・配布基盤 | 0.6までに完了 |
| I1 | 入力YAML一元化、seed検査、保存・反復・再開 | 0.7で完了 |
| O1 | RF位相の保存・読込、positive-columnの保存解選択 | 0.8で完了 |
| I2 | 四solver共通のDC反応別`k(E/N)`点データ | 0.9で完了 |
| O2 | 同じYAMLと保存済みデータからのSwarm標準図 | 0.10で完了 |
| I3a | positive-columnの単位付き結果・監査・保存解比較 | 0.11で完了 |
| I3b1 | 資格付き元表とCOMSOL形式変換の分離 | **0.12で完了。全体回帰・配布後実行を確認済み** |
| I3b2 | 保存MPHの単独出力、CCP/ICPの結果入口 | **次に進める範囲。I3全体は未完了** |
| I4 | 非COMSOL向け係数の読込・利用例 | 係数取得までに限定 |
| Phase 8 | 独立物理参照、ライセンス下COMSOL受入 | 未完了。ソフトウェア検査とは別判定 |

利用手順は[Swarm入力・実行](../swarm_inputs.md)、
[COMSOL入出力](../comsol_workflow.md)、
対応範囲は[出力仕様レビュー](modeling_output_review.md)を参照する。

## 責任の境界

| 責任 | 所有者 |
| --- | --- |
| 固定した気体条件のEEDF、反応係数、輸送量、対応するRF電子応答 | Swarmの四solver |
| 資格付きE/N点データ、元の値・符号・不確かさ・由来 | `swarm_workflow/tables` |
| COMSOL用の独立変数・関数形式、入力bundle | `swarm_workflow/comsol/input` |
| 化学種濃度・組成の時間発展、粒子・エネルギー収支、壁・流入出 | COMSOLや外部Globalモデルなどのプラズマモデル |
| 既存プラズマ解からのデータ取得と図 | 対応する既存adapterの後処理 |

Swarmに反応網、濃度初期値、化学種の時間積分、自動的な組成結合ループを追加しない。
Swarm YAML一つと接続先mapを入口とし、図専用YAMLや汎用観測frameworkを増やさない。

## I3b1の対象と実装

共通層が作っていた平均エネルギー軸のCOMSOL用テーブルと、そのrole付与を
`comsol/input/export/coefficient_tables.py`へ移した。schema-v2 YAMLとCLI入口は維持する。
MPHのタグ、方程式、closure所有者、Swarmの数値モデル、COMSOL受入閾値は変更しない。

| 元の観測量 | 資格付き元表・独立変数 | COMSOL用の表・独立変数 |
| --- | --- | --- |
| 輸送量 | `transport_vs_en.csv` / E/N | `transport_vs_mean_energy.csv` / 平均エネルギー |
| 反応別係数 | `rates_vs_en.csv` / E/N。各行に平均エネルギーも保持 | `rates_vs_mean_energy.csv` / 平均エネルギー |
| 符号付き弾性損失 | `elastic_energy_loss.csv` / E/N | `elastic_energy_loss_vs_mean_energy.csv` / 平均エネルギー |

- 元表は形式v3、COMSOL bundleは独立した形式v2。古い元表は記録済みの現行DBから再生成する。
  数値実装とworkflowが一致するDBを使い、過去の資格hashや形式番号を付け替えない。
- 元表はE/Nの点データとして保持する。平均エネルギーの単調性を共通層の出力条件にしない。
- COMSOL変換時に量ごとの独立変数・単位・平均エネルギーの一価性を検査する。
  非単調な元データは保存し、変換だけを明示的に拒否する。点を修正・補間して通さない。
- 各反応行の実際の平均エネルギーを使う。MCの輸送とratesの採用範囲が異なっても、
  狭い輸送格子へ結合して点を失ったり、別のE/Nから値を補ったりしない。
- 値・符号・不確かさ・欠損を保持し、平均エネルギー軸へ並べ替える。
  派生表のsource filename/hashと投影方法を記録し、元表も同じbundleへ保持する。
- 既存COMSOLモデルの表名・消費量は維持する。受入資格はmodel adapterが別途判定する。

三つの資格付き`build-tables` sourceと、四solverの中立な反応係数exportは別の入口である。
`multi_term`をCOMSOL入力へ昇格させたり、欠けた輸送量を別solverから補ったりしない。

## 改良中に見つけて修正した事項

| 問題 | 0.12の動作 |
| --- | --- |
| 共通の元表生成がCOMSOL用の平均エネルギー軸表まで作る | 形式変換とCOMSOL roleをCOMSOL入力層へ集約 |
| ratesを輸送格子の平均エネルギーへ結合するとMCの独立な点を失う | 各行の観測値と量別のE/N範囲を保持 |
| 単位の宣言を確認せず、元データへSIラベルを付け直し得る | 宣言されたcanonical単位を検証。欠損・相違は拒否 |
| GEC-CCPがbundle版を共通テーブル層から取得する | 各形式の所有者を分け、COMSOL側の版を参照 |
| 非単調データの範囲を先頭・末尾だけで表し得る | 実際のmin/max、隣接差の大きさを記録 |
| 組成ごとの再出力途中で失敗しても古いroot成功が残る | 再出力開始時に古いroot manifestを無効化 |
| 汎用の単位名推測helperや共通層の重複した平均エネルギー表生成 | 未使用推測・重複生成・不要な単調性診断を削除 |

## 整理と移行

通常・compositeの両方から平均エネルギー軸表の重複生成を削除した。
再生成では古い派生表や未登録の残骸が元表ディレクトリへ残らないことを確認する。
旧CSV前提のテストは現行の値・単位・範囲・投影結果を検証する試験へ置き換える。
新しい単位ライブラリ、設定ファイル、化学種積分器は追加しない。

作業前snapshotはリポジトリ外へ保存した。元の断面積、MPH、利用者の資料は保持する。
旧v09 DB、手作業の図例、個別描画scriptの退避・削除は、前段階の自動承認レビューが
「blocked by policy」として拒否した。これらと同様に拒否された`build/`は再試行せず残している。
現行参照・配布対象からは外し、`tools/build_wheel.py`で隔離ビルドする。
今回の一時検証script・ログ17ファイルは、hash照合後にリポジトリ外の検証記録へ退避し、
作業ディレクトリから削除した。作業前snapshotと最新wheelへの参照は保持した。

## 検証と数値資格

- 最終`py -3 -m pytest -q`: **1232 passed, 2 skipped in 980.76 s**。
  Java動作試験にはCOMSOLのJDKを指定した。スキップはローカルCOMSOL実行設定を要する2件。
- 形式変換・単位・符号・範囲・移行の重点試験、GEC-CCPの112件も通過した。
- リポジトリ外に0.12 wheelをインストールし、299実行ファイルが最終ソースと一致した。
  四solverの数値モデルとworkflow identityは0.11から不変である。
- CCP/ICPの現行v10 DBのコピーから、各3点の元表9表・bundle16表を再生成した。
  3種類の派生係数表について、行・値・符号・欠損・元表hashの一致を確認した。
  旧派生表の除去、単位/旧形式の拒否、過去PN artifactからの中立rates出力も確認した。
  元のDB・map・MPH・資格JSONは変更していない。新しい数値計算の資格として再発行しない。
- Ruff E9/F、文書11件のローカルリンク、`git diff --check`を確認した。

[23-case P1 core資格](results/propagator_p1_core_qualification.json)、
[CCP](results/propagator_gec_ccp_target_qualification.json)・
[ICP](results/propagator_gec_icp_target_qualification.json)の数値証拠は維持する。
現行v10 DBの記録を再利用する際も、元DB・資格JSONのhashを変更しない。
ライセンス下COMSOL solveと独立した物理受入は今回の実装確認には含めない。

追加の実map preflightでは、元mapとMPHを変更せず、出力先だけを隔離した。
ICPはcheckoutからのprepareを通過したが、wheelでは`examples/cross_sections/argon_application_library.csv`
をインストール先から探して停止した。`Model/Ar_xsecs.txt`も同じcheckout依存を持つ。
CCPは保存資格のoperating bracketが3000--4000 Td、adapter要求が3000--3500 Tdであるため停止した。
元の資格はその宣言範囲で有効だが、今回の検証をCCPモデル受入済みとは扱わない。
両方とも最初の失敗条件を記録し、資格JSONの変更、map要求の緩和、ライセンス実行は行っていない。

## 次に進める内容と残件

1. **I3b2の先頭:** ICPの外部資産参照をcheckout固定から分離し、wheelでも宣言した入力へ到達できるようにする。
   CCPは要求範囲と資格範囲の差の意味を確認し、必要なscopeで正規の手順から再評価する。
   古い資格の書換えや、要求を緩めるだけの整合は行わない。
2. positive-columnの保存MPH単独exportを結果入口へ接続する。
   元モデル・解・単位の証拠と外部closure受入を区別し、再solveや成功ラベルの追加を要求しない。
3. CCP/ICPの既存結果について、保存解・必要な観測位置/時刻・図への入口を接続する。
   CCP overviewの全輸送量への依存と、ICP export/比較入口の分散を整理する。
4. I4は係数読込・単位/固定条件・範囲内の値取得例までとし、Globalモデル本体は作らない。

RF反応係数図、追加のエネルギー損失quantity、組成軸の比較図は未実装であり、必要な利用経路と別途扱う。
Phase 8には、旧argon断面積の原出典確認、同じ角度仮定・輸送定義の独立参照、
真の縦拡散応答とF0-gradient/weighted MCの区別、適格入力によるライセンス下COMSOL受入が残る。
