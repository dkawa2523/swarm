# モデリング利用者から見た入力・データ連携のレビュー

2026-10-27、0.10の現在位置。初回0.6レビューのうち、移行済みの古い設定手順は削除した。
現行の具体的な操作は[Swarm入力・実行](../swarm_inputs.md)へ一本化する。

## 入力の責任

通常のSwarm計算、E/N・混合比のスキャン、MC独立反復、SQLite保存・再開は一つの
schema-v2 YAMLと断面積データで設定する。`run.solvers`は四つの計算手法、
`physics`は衝突・場・散乱などの物理要求、`campaign`は保存とsamplingを所有する。
`output`は出力先、任意の係数exportと標準図を所有する。
scenario、preset、YAML継承、基本YAMLへの参照は使わない。

COMSOL形式への変換、MPHタグや保存解の選択、空間モデルの実行は別工程である。
Swarm側へMPH設定を加えず、通常は接続先mapのbundle参照を使う。
mapの相対パスはそのYAMLの場所基準、bundle内の参照はbundle基準となる。

## 調査項目の現在位置

| ID | 元の問題 | 現在の対応・残件 |
| --- | --- | --- |
| M01 | 基本設定とworkflowでE/N・組成を二重に指定する | I1で単一YAMLへ移行。旧基本設定9ファイルと重複した構造テストを削除 |
| M02 | seedなしMCが計画を通過し実行時に失敗する | 直接計算とcampaignの両方で計画時に必須検査 |
| M03 | multi_termのratesを他モデルへ渡す共通入口がない | I2で四solver共通のDC反応係数点をCSV/APIへ出力 |
| M04 | ratesだけでも拡散・Townsend換算を要求する | 中立なrates出力では不要。既存COMSOL一式の条件は維持 |
| M05 | E/N係数点にも平均エネルギー軸の単調性や表の補間条件が及ぶ | 係数点と補間表を分離。中立exportは実点と範囲のみを出す |
| M06 | e-e後処理で輸送が古いと更新済みratesも利用できない | 更新済みratesを出力し、元のtransport-stale情報を保持。図でも明記 |
| M07 | 共通のデータとCOMSOL固有の形式・受入が混在する | 中立な係数点・標準図・出版処理を分離。COMSOL固有の形式変換整理はI3 |
| M08 | モデルによって相対パスの基準が異なる | 手書き設定はYAMLの場所基準へ統一 |
| M09 | map内bundleをCLIでも重複指定する | 通常はmapの参照で準備。明示的差替えだけ実効mappingへ記録 |
| M10 | 計算前に総ケース数・予算を確認しにくい | 共通`run --plan`、詳細`--plan --json`へ統一 |

`output.export`は指定した一つのsolverのDC係数を出す。
`output.figures`は計算された各solverのEEDF・輸送量・反応係数と、対応する角度/RF図を選ぶ。
これらは独立した設定であり、同じYAML・同じ出力世代に共存する。
図専用の設定ファイルは作らず、保存済みmanifest/DBは`swarm-workflow plot`で再描画する。

## 利用時に区別する意味

| 論点 | 保持する意味 |
| --- | --- |
| 係数の基準 | `k [m^3/s]`は対象粒子あたり、`x k`は組成重み付き。二体源`n_e n_target k`に組成を二重に掛けない |
| 反応の対応 | species/process/type/thresholdを保持し、process文字列から化学量論を推定しない |
| 独立変数 | E/N、平均エネルギー、電子温度を区別。引数変換の一価性を単調化で偽装しない |
| 固定条件 | 気体組成・温度・密度・場・population・角度仮定・e-e条件を保持 |
| 輸送量 | flux、F0-gradient近似、restricted density-packetを異なるquantityとして扱う |
| 欠損・不確かさ | 未計算、stale、反復不足、全ゼロを0や精度保証へ変換しない |
| 図と由来 | 元の数値と条件を保持。過去のデータの再描画を現在の実行・入力資格としない |
| 適用範囲 | 係数点は補間・外挿・反応網や、利用先モデルの受入済み入力を意味しない |

化学種の濃度・組成の時間発展、反応網、壁・流入出、粒子・エネルギー収支は利用側の責任となる。
Swarmの組成配列は独立条件であり、時間軸ではない。RF位相図は既存の電子応答を表す。
新しいGlobalモデルや化学種時間積分を追加しない。

## 次の作業

O2の次はI3。COMSOL固有形式と必要量の対応を整理し、既存プラズマ解の識別・
観測位置/時刻・単位情報を比較・描画へ接続する。
続くI4は非COMSOLの係数読込・単位/固定条件確認・範囲内の値取得の例に限定する。
独立物理検証とライセンス下COMSOL受入はPhase 8として別に残る。
完了位置、実計算・配布後検証と残件は[改良計画](improvement_plan.md)を参照する。

## 主な実装根拠

- [統一入力](../../swarm_workflow/campaign/document.py)、[計画と実行](../../swarm_workflow/campaign/sweep.py)。
- [出力設定](../../electron_swarm/core/output.py)、[中立な係数点](../../electron_swarm/io/reaction_rates.py)。
- [標準図](../../electron_swarm/io/figures.py)、[保存済み再描画](../../swarm_workflow/figures.py)。
- [物理条件](../../electron_swarm/core/physical_context.py)、[輸送量の意味](../../electron_swarm/core/observables.py)。
- [COMSOL adapter](../comsol_model_adapters.md)、[物理的制約](../physics_limitations.md)。
