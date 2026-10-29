# モデリング利用者のための出力仕様レビュー

2026-10-27。0.8のO1、0.9のI2に続き、0.10でO2（Swarm標準図と再描画）を実装した。
O01--O03、O05の中立なrates出力に対応し、O04のSwarm図を全輸送量への依存から分離した。
0.11のI3aでpositive-columnの結果単位・保存解と比較入口を接続した。
0.12のI3b1で、元の資格付きE/N表とCOMSOL用の平均エネルギー軸変換を分離した。
元表は形式v3、bundleは独立した形式v2。反応行の平均エネルギーと量別の範囲を保持する。
COMSOL固有overviewの整理とO06のICP接続はI3b2に残る。

2026-10-26の範囲再確認で、I4に含めていた利用側の時系列・組成図をSwarm基盤の必須範囲から外した。
I4は外部向け係数の読込例に限定する。組成の時間発展はプラズマモデル側が解き、
その結果の描画は既存adapterの後処理として扱う。

## 判断

電子データの計算結果はCSV・SQLite・Pythonで取得できる。一方、利用者が一つの入口から
必要な量・図・観測点・時刻を選び、異なるシミュレータでも同じ意味で使う仕様には未到達である。
最優先は出力項目を増やすことより、返されたデータの保存漏れと曖昧な解・観測点の選択をなくすこと。

Swarmは基本的に一様な気体条件に対する電子データを計算する。
気体の組成配列やE/Nスキャンは独立条件の集合であり、組成の時間発展ではない。
RF位相別EEDFも、重粒子化学種の過渡反応・空間分布を解く機能とは区別する。
化学種の時間発展にはCOMSOLやGlobalモデルなど、粒子・エネルギー収支を解く利用側の結果が必要となる。

範囲再確認の時点では0.7の検証済みwheelと実行Pythonファイル289件が一致した。
その後の0.8改修は保存・読込と既存COMSOL結果の選択に限定する。
新しい化学反応網、化学種濃度の時間積分、Globalモデル本体は追加していない。
既存のRF位相別EEDFとMC電子軌道は電子計算の一部であり、組成の時間積分とは区別して維持する。

## 現行の入出力

| 利用目的 | 実装済みの内容 | 設定・利用上の限界 |
| --- | --- | --- |
| 各Swarmケースの物理量 | `summary.csv`。平均エネルギー、flux drift、移動度、利用可能な拡散・エネルギー輸送、反応由来の正味生成頻度など | 全solverが全量を返すわけではない。未計算・近似・staleの意味をmetadataと合わせて読む |
| 反応別係数 | `rates.csv`。species/process/type/threshold、対象粒子あたりの`k`、組成重み付き`x k`、損失係数など | 化学反応網の化学量論や反応源の全体を定義するファイルではない |
| EEDF/EEPF | `eedf.csv`。エネルギー座標、セル幅、分布、MCのcount情報 | EEDFは確率密度。countやweight effective countは独立反復の信頼区間ではない |
| エネルギー・角度分布 | solverが返す場合の`energy_angle_distribution.csv` | 全solver共通の角度分布出力ではない。座標・セル幅を保持する |
| RF位相 | 対応するtwo-termの直接CSVと、campaignの`rf_phases`/`rf_phase_eedf_bins` | 周期F0・瞬時F1の近似。必要な位相セルの欠落は再開時に再計算する |
| campaignの原結果・統計 | SQLiteのcases/rates/EEDF/energy-angle/RF、`WorkflowStore.read_case`、`aggregate`の係数・EEDF・quality CSV | 読込は原結果の復元であり、新しい実行証拠ではない。集計対象は原結果全体とは一致しない |
| Swarm出力のYAML指定 | `output.directory`、`base_name`、任意の`export`、`figures` | exportはDC反応係数点。標準図5種類に対応。任意空間観測点・補間時刻は利用側の責任。未知の項目は拒否 |
| positive-columnの結果 | `results.fields`に名前・COMSOL式・単位、`dataset`と`solutions: all|last`を明示。単位とCSV hashを持つrun summaryから読込・比較する | 比較時に複数保存解を明示選択。既知の単位を変換。任意座標や補間時刻は拒否する。外側解は1件、内側parameterは最大1軸 |
| GEC-CCPの結果 | 周期平均のdomain/axis/radial、位相別軸・domain/closure、電極波形など | 主な式と出力構成がadapter実装内で固定。任意の化学種・観測点をYAMLで選べない |
| GEC-ICPの結果 | 通常実行の体積積分・平均エネルギー範囲・coil powerの時間系列。別途、保存解の共通時刻・終端・全保存時刻の比較export | Ar/Ar*/Ar+モデル向けに固定。観測点、任意化学種、任意の組成基準を選ぶ共通入口はない |
| COMSOL入力 | E/N・平均エネルギーに対する係数、反応別rates/Townsend、EEDF、Function EEDF、solver-native弾性損失など | `build-tables`のsourceはtwo_term/propagator/monte_carlo。multi_termを含む中立exportは別入口であり、COMSOLへの受入を保証しない |
| 非COMSOL向け係数 | 四solver共通の`reaction_rates.csv`と`read_reaction_rates`。固定条件・単位・実点/範囲・数値状態・MC反復区間を保持 | 補間・外挿と係数取得の利用例は未実装。物理/接続先受入済みの入力とは扱わない |

読み出しはCSV単体より、`electron_swarm.io.read_product_run`でmanifestと一緒に行う。
ケースの物理条件、observableの定義、実装・入力の由来、近似・欠損を保持できる。
実行の完了、数値収束、独立した物理検証、利用先モデルの受入は異なる判定である。

## 図の現状

Swarm標準図は`output.figures`と`swarm-workflow plot`へ統一した。
直接計算・campaign・保存済みmanifest/DBで同じ描画処理を使う。
PNG/SVG、曲線CSV、元の2次元配列を同じmanifestへ結び付ける。
詳細な[設定・再描画例](../swarm_inputs.md#出力と図)を参照する。
次の専用コードは資格評価またはプラズマモデルの結果を扱うため維持する。

- `plots/eedf_comparison.py`と`tools/validation/run_eedf_comparison.py`は、two-term/P1/MCの
  EEDF・tail・差分をPNG/SVGにする。四solver任意の組合せを描く共通図ではない。
- `plot-gec-ccp`はbundleの係数・EEDFと、対応するCOMSOL結果の空間分布・RF位相・波形を描く。
  現行のoverviewには全輸送量への依存があり、P1のように一部の量を返さないsourceで止まる。
- GEC-ICPの`comparison_analysis`には、電子・Ar*・Ar+のinventory、平均エネルギー、電力の
  時間履歴PNG/SVGがある。保存解比較に限定され、一般の化学種組成プロットではない。
- 品質を通過しないデータの考察図を、適格なモデル入力・受入済み比較と分ける必要がある。
  MC EEDF比較にはraw値と失敗状態を残す経路がある一方、CCP結果図は受入済み結果を要求する。

inventoryは体積積分した粒子数、number densityは局所密度、fractionは指定した分母に対する割合。
Ar*のinventory時間履歴があることを、任意のガス組成比を時間追跡できることとは扱わない。

## 0.7で再現した課題と現在位置

| 優先度・ID | 0.7で再現した挙動 | 0.8の対応・残件 |
| --- | --- | --- |
| P0 / O01 | RF配列を持つケースをWorkflowStoreが保存しない | 修正済み。位相・EEDFを同一transactionで保存。欠落時の再計算、型付き読込、直接CSVとの一致を検証 |
| P0 / O02 | probeへ`point_m`、`time_s`、`dataset`を追加しても黙って捨てる | 修正済み。`fields`と明示的datasetへ変更。未知の項目と旧`probes`を拒否 |
| P0 / O03 | datasetを指定せず末尾解だけ出し、解・時刻を記録しない | 修正済み。`all`/`last`を必須にし、解・点・時刻・parameterを保持。外側sweep等は明示的に未対応 |
| P1 / O04 | P1のEEDF・rates・mobilityが存在しても、CCP overviewは電子エネルギー移動度を必須として停止する | O2でSwarm標準図を分離し、図ごとに必要量を判定。未計算の輸送量を0で埋めない。既存CCP overviewの整理はI3 |
| P1 / O05 | multi_termが連携テーブルのsourceに含まれず、ratesだけでも輸送・Townsend・平均エネルギー軸の条件を要求する | I2で中立な係数点の入口を追加。四solverのratesを輸送なしで出力。既存COMSOL表の条件は維持 |
| P1 / O06 | ICPの時間履歴・比較exportは別モジュールで、通常の出力設定・CLIと分散する | 通常利用の結果設定から観測・集計・時間選択を辿れるようにする。固定の受入出力と任意の解析出力を区別する |

O01は合成配列の損失なし往復に加え、two-termの実RF計算で保存・再開・欠落修復を確認した。
この検査はRF物理モデルの独立精度検証ではない。
O03の配列順は[COMSOL 6.3 NumericalFeature API](https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/api/com/comsol/model/NumericalFeature.html)
の`data[expr][solnum][vertex]`と照合した。実機での複数解・外側parameter sweepの挙動は、
対応するCOMSOLモデルで別途確認する必要がある。0.8では複数解×複数点のJava実行試験と
実COMSOL 6.4 APIへのcompileを行った。ライセンス下のsolveは実施していない。

改修中には、古いCSVをstdout経由で成功扱いできる経路と、異なる時刻の平均値が一致すれば
profile比較に通る経路も見つかった。前者を削除し、後者は一つの解と一定の運転条件・時刻を必須にした。
生EEDFの未使用列にあった`1/sqrt(ESS)`のRSEも削除した。独立MC反復の不確かさは維持する。

0.11ではrun summaryを比較入口にし、監査・比較で同じ単位変換を使う。
cm、密度のcm^-3、Td、J、kPaなど明示した対応単位を自動変換し、未知単位は拒否する。
詳細な[読込・比較手順](../comsol_workflow.md#profile-comparison)を参照する。
非一様格子の点ごとのノルムを空間積分へ改め、図と同じ曲線CSVを出力する。
電子とArイオンの伝導電流の空間一様性を、電荷保存やプラズマbulkの診断とは呼ばない。
単体CSVや旧summaryから単位・受入を推定する経路、重複した読込と常時ゼロの診断を削除した。
新しい形式へ接続されていないbuilt-in参照MPHの単独export、CCP/ICPとの接続は残件である。

## 改訂後の出力の責任

| 層 | 必要な内容 | 利用者の設定場所 |
| --- | --- | --- |
| 計算の原結果 | solver/case/混合比/E/N、返された分布と係数、単位、物理条件、欠損・近似、必要なMC統計、入力と実装の由来 | Swarm YAML。再開用DBへの保存と、選択したCSVの書出しを混同しない |
| 解析用の図 | 原結果に結び付く図、描画元の数値、軸・単位・系列、品質状態、表示範囲、補間した場合の方法 | 同じSwarm YAMLの任意の図設定。図のためだけに別の基本YAMLを要求しない |
| シミュレータ向け電子データ | 必要量・solver・独立変数・反応識別・係数の基準・補間/範囲外方針と量別資格 | 任意のexport設定。単純利用ではSwarm YAML内で完結させる |
| 既存プラズマ結果の取得・描画 | そのモデルが解いたdataset/solution、時刻/位相、位置、出力式・単位と必要な図 | 接続先モデルの結果設定。保存済み解の抽出と描画を扱う |

化学種の粒子・エネルギー収支、反応網、壁・流入出、濃度の初期値・時間積分はプラズマモデルが所有する。
Swarm YAMLへこれらの設定や、化学反応の時間積分を呼び出す入口を追加しない。
オンラインで組成を受け渡してSwarmを再実行する結合ループも今回の対象に含めない。

CSVと最小限のmanifestを基本とし、図はPNGとSVGを標準候補にする。
大きな分布への別形式追加は、実際のI/O量を確認して決める。複数の新形式を先行して必須化しない。
原結果・解析図・適格な投入表を同じ成功フラグでまとめない。

### 同じSwarm YAMLの任意設定

`export`は0.9、`figures`は0.10で対応した。図専用のYAMLは追加しない。

```yaml
output:
  directory: ../outputs/argon
  base_name: swarm
  export:
    solver: multi_term
    quantities: [reaction_rates]
```

この例はmulti_termのratesだけを利用側へ渡す。拡散・Townsend・COMSOL関数表の生成を要求しない。
計算するsolverの指定は引き続き`run.solvers`に置く。図やexportの指定で物理条件を再定義しない。
I2はDC E/Nを独立変数とし、`reaction_rates`を対象粒子あたりの`k [m^3/s]`と定義する。
`k`と`x k`、元の数値状態、独立MC反復のSE/CI95、物理条件と由来を保持する。
e-e後処理のratesは出力し、transport-staleは元caseのmetadataへ保持する。
反復不足・全ゼロの区間は空欄、欠けた反応の補完・別solverへの代替は行わない。
手順は[利用文書](../swarm_inputs.md#結果を他のモデルへ渡す)を参照する。
O2は任意の`output.figures: [eedf, transport, reaction_rates]`を同じYAMLへ追加する案とし、
PNG/SVGを既定にする。図用の別YAMLや有効化フラグを要求しない。
平均エネルギー軸などの追加は、それを必要とする利用経路の検証と合わせて別途実装する。

### 既存adapterの結果取得を整理する範囲

- O1でpositive-columnの出力式一覧を`fields`へ変更し、既存dataset、保存解選択を必須にした。
  点観測の意味は持たせず、未対応の座標・時刻指定は拒否する。
  汎用のpoint/line/domain観測機構は追加していない。
- 出力には選んだdataset・解番号・その解のパラメータ、時刻がある場合は実時刻を保持する。
  暗黙に末尾へ縮約しない。既存解を抽出する処理からプラズマの時間積分を起動しない。
- I3で実際に必要な点取得を実装する場合だけ、対応する既存モデル内で座標・単位・
  解・時刻の選択を追加する。範囲外の点や未対応の選択はエラーとする。
- 固定の受入用出力はadapterが所有する。利用者が図を省略しても受入証拠は失われない。

## Swarm標準図の仕様

| 図 | 表すもの | 誤読を防ぐ条件 |
| --- | --- | --- |
| EEDF/EEPF、tail | 分布のセル値と全エネルギー域 | 規約・単位・セル幅を保持。表示範囲外の質量を示し、表示のために再正規化しない |
| 輸送・平均エネルギー | 固定組成ごとのE/Nに対する移動度・drift・利用可能な拡散等 | quantityが違えば別パネル。近似、stale、未計算量を区別。線は実点間の視覚的補助 |
| 反応 | DCの反応別`k` | `x k`と混同しない。反応識別と係数の基準を保持。エネルギー損失の追加quantityは未実装 |
| エネルギー・角度分布 | 対応するsolverの分布と異方性 | 座標・向き・角度閉包を記載。ordinary XSをexact DCSと表示しない |
| RF位相 | 電子計算が返した位相別EEDF・平均エネルギー・rates | 周期内の電子応答として示し、気体組成の時間発展と呼ばない |

独立MC反復の不確かさは、その定義と水準を付けて図に示す。
未適格データも考察用には明示した状態で表示できるようにするが、その図を入力資格の証拠に昇格させない。
対数軸でゼロ・負値が表示できない場合も、データを正の値に置き換えない。

### プラズマモデル結果の図

こちらは既存adapterが取得した計算済み結果の後処理であり、Swarm標準図とは別に扱う。
新しい汎用化学種モデルや共通時間積分基盤の追加は要件としない。

| 図 | 入力に必要な既存結果 | 表示上の条件 |
| --- | --- | --- |
| 空間分布・点の履歴 | プラズマモデルが計算した位置・時刻と物理量 | 局所値・平均・体積積分、単位、選択した解を区別する |
| 化学種密度・inventoryの時間履歴 | プラズマモデルの化学種データ | 実時刻と集計方法を保持する。データのない化学種を計算して補わない |
| 組成比の時間変化 | プラズマモデルが出した割合、または基準を明示できる十分な密度データ | number/mass/moleの基準と分母を明示する。不足する化学種から全組成を推定しない |

## COMSOL・他のモデルへ渡す粒度

共通のrates出力は、solver・固定物理条件・独立変数の各点に対する反応別の行とする。
species/state/process/type/thresholdなどの識別を保持し、文字列の反応式から化学量論を推測しない。
`k`の単位は対象粒子密度を掛ける前の`m^3/s`。二体源`n_e n_target k`と`x k`の使用で
組成を二重に掛けない。電子エネルギー損失係数や符号付き弾性交換は別のquantityとして扱う。

E/N表に平均エネルギー軸の単調性を要求しない。平均エネルギーを独立変数にする場合だけ
対応関係の一価性を確認する。補間法・有効範囲・範囲外の拒否を保持し、元の点を単調化で書き換えない。
固定組成・温度・populationの表を、変化するGlobalモデルの任意状態に適用できるとは保証しない。

共通の量別判定の後に、COMSOL固有のFunction EEDF・関数表記・MPHへの対応を行う。
既存COMSOL受入は維持する。multi_termのratesや、e-e後処理後の更新済みEEDF/ratesを
利用するために、欠けた輸送量を捏造したりCOMSOL全体の受入条件を緩めたりしない。

## 実装順と完了条件

| 段階 | 完了させる範囲 | 主な確認 |
| --- | --- | --- |
| O1: 保存と選択の正しさ | 0.8で実装。RF位相保存・再開・読込、positive-columnの出力式と既存解選択 | 直接出力とDBのRF配列が一致。複数解を識別して保存し、未対応指定を拒否 |
| I2: 必要量別データ | 0.9で共通のDC反応係数点を実装 | PN ratesを輸送なしで出し、物理条件・由来・数値状態・反復区間を保持。投入表としての受入とは区別 |
| O2: Swarm標準図 | 0.10で実装。同じYAMLの任意設定と保存済み結果からの再描画 | EEDF・主要量・rates、対応する角度/RF図。PNG/SVGと元数値が一致し、未対応量を明示する |
| I3a: positive-column結果 | 0.11で単位付きsummary、保存解選択、監査・空間比較・PNG/SVGを接続 | 一つの正規化snapshotと空間積分を使う。実機の受入は別判定 |
| I3b1: 形式変換の分離 | 0.12で資格付き元表とCOMSOL形式を分離。値・単位・量別の範囲を保持 | 元の非単調E/N点を修正せず、COMSOL用の一価性は変換時に判定する |
| I3b2: 残る結果入口 | CCP/ICPの結果入口、保存MPHの単独export | 既存解の識別、必要な観測位置・時刻を実機で確認する |
| I4: 非COMSOL利用例 | 係数CSV/APIの読込、単位・固定条件の確認、範囲内の係数取得 | COMSOLなしで利用できる。新しいGlobalモデルや化学種時間積分を実装しない |

各段階は狭い範囲を設定・計画・実行・出力・文書・振る舞いテストまで完成させる。
次はI3b2。保存MPHの単独出力と、他の既存プラズマ解の取得・描画を接続する。
独立物理検証とライセンス下プラズマモデル受入は、引き続きPhase 8の別判定である。

## O2で見つけて修正した事項

| 問題 | 修正 |
| --- | --- |
| 係数exportのsolver選択で比較図まで絞られ得る | exportと図を分け、図は保存された各solverから作る |
| 異なる拡散定義を同じ量の曲線として重ね得る | typed quantityで別パネルにし、方法・近似・staleを示す |
| RF位相を0始まりの左端値と仮定すると元の座標とずれる | 保存位相を保持し、周期的な最近傍sampleを明記してheatmap化 |
| MCの格子差を補間・再正規化で見かけ上揃え得る | 共通範囲へ確率質量を保存して再配置し、除外された質量を記録 |
| 負の微小tailやゼロ観測を正のfloorへ置き換え得る | 元の符号とゼロをCSVへ保持し、symlog/logでの表示規約を記録 |
| 画像の削除で数値だけのexportも止まり得る | 数値検証と画像検証を分け、再描画・数値exportは欠損画像を許容 |
| 保存済みデータの再出力が現在の実行成功に見え得る | 元のidentityと解析scopeを保持し、図の生成版だけ別に記録 |
| PNG/SVGとCSVが異なる世代として残り得る | 一つのpublisherへ集約し、失敗時rollbackと旧図の削除を行う |

実データの標準図は`outputs/standard_figures_20261027/`へ保存する。
旧手作業の図例と個別描画scriptは標準APIの実行例で置き換えた。
旧生成物の退避・削除は自動承認レビューが拒否したため、現行文書からの参照を外して残す。
四solverの数値モデルは0.9から変更しない。実行基盤の変更に対し、CCP/ICPのmedium/fineを
v10 DBへ新規計算して資格を検証する。最終検証結果は[改良計画](improvement_plan.md)に記録する。

## 主な根拠

- [出力設定](../../electron_swarm/core/config_sections/output.py)、[OutputConfig](../../electron_swarm/core/output.py)、
  [product writer](../../electron_swarm/io/writers.py)、[result型](../../electron_swarm/core/results.py)。
- [campaign保存](../../swarm_workflow/campaign/store.py)、[集計](../../swarm_workflow/campaign/aggregate.py)。
- [positive-column設定](../../swarm_workflow/comsol/models/positive_column/config.py)、
  [field export生成Java](../../swarm_workflow/comsol/models/positive_column/java.py)。
- [CCP出力](../../swarm_workflow/comsol/models/gec_ccp/java_export_sections.py)、
  [CCP図](../../swarm_workflow/comsol/models/gec_ccp/plots/workflow.py)、
  [係数/EEDF overview](../../swarm_workflow/comsol/models/gec_ccp/plots/eedf.py)。
- [ICP通常出力](../../swarm_workflow/comsol/models/gec_icp/java.py)、
  [ICP保存解export](../../swarm_workflow/comsol/models/gec_icp/comparison_export.py)、
  [ICP時間履歴](../../swarm_workflow/comsol/models/gec_icp/comparison_analysis_render.py)。
- [共通テーブルの現行制約](../../swarm_workflow/tables/contracts.py)、
  [COMSOL出力の役割](../../swarm_workflow/comsol/input/export/contracts.py)。
