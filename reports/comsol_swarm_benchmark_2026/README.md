# 外部Swarm表–COMSOL 1D Ar DC glow-discharge benchmark

本directoryは、外部Swarm表をCOMSOL 6.4の1D Ar DC glow-discharge
Application Libraryモデルへ接続し、組込みBoltzmann経路と比較するための
手法・再現性benchmarkである。対象は電極間全域（cathode fall、dark-space／
positive-column構造、anode側を含む）であり、一様な正カラムだけを切り出した
モデルではない。2D GEC CCP、実験validation、solverの精度優越は対象外とする。

## 現在の判定

**`share_with_caveats`（暫定技術benchmark）**。

- COMSOL profile、mesh summary、時間は2026-07-14の履歴証拠である。
- 2026-07-31のclean COMSOL再実行はlicense error -10
  （`Product has expired`）でmodel apply前に停止した。
- 履歴runにはcompiler false-positiveがあり、Javaと実行classの対応を証明できない。
- 8.62倍／4.73倍は各経路1観測の候補比であり、正式な性能値ではない。
- 最新schema-v2 Swarm表は34/34点でtail probability、edge/peak、20 keV
  grid-limit条件を満たす一方、iterative solver convergenceは23/34点である。
  したがってtail監査済みではあるが、34点すべてが正式収束した入力ではない。

## 外部入力の責任境界

履歴COMSOL経路はLocal Energy Approximation（LEA）であり、平均電子エネルギー
方程式をCOMSOL側で解く。外部Swarmはすべての電子物理を置換するのではなく、
選択した係数だけを置換する**selected-channel hybrid**である。

bundleの契約は、次の**7 lookup quantities = 7 X-Y pairs = 14 numerical
arrays**である。

1. \(E/N\rightarrow\bar{\varepsilon}\)
2. 換算電子移動度
3. 縦方向換算電子拡散
4. 換算energy mobility
5. 換算energy diffusivity
6. eir2 direct-excitation reduced Townsend
7. eir4 direct-ionization reduced Townsend

LEAでは1番目のpairは保存されるが方程式上はinactiveで、2–7番目が候補active
集合となる。COMSOLは電子密度・平均エネルギー・Poisson／circuit／heavy-species
方程式、境界条件、残余反応・energy treatmentを保持する。したがってproperty
readbackは「書込み」の証拠であり、7量の同時functional activationの証明ではない。

## 成果物

- `report.html`: English abstractと印刷用CSSを含む自己完結HTML。
- `artifact.json`: HTML生成用canonical report artifact。
- `manuscript.md`: プラズマ物理・COMSOL定式化を明示した和文論文原稿。
- `benchmark_analysis_executed.ipynb`: 先頭から実行済みの比較解析notebook。
- `data/derived/`: 19 CSV。空間metric、lookup／tail監査、operating point、
  Townsend表現感度、partial electron-energy、fluid-applicability、runtime、
  QA、SHA-256を含む。
- `figures/`: 9図の投稿用PNGとvector PDF。
- `data/formal_swarm_bundle_solver_qa_20keV/`: 2026-07-31 schema-v2候補表。
- `data/formal_swarm_database_solver_qa_20keV.sqlite`: 34点のcase diagnosticsを
  含むworkflow database。
- `data/archive_2026-07-14/`: 暫定COMSOL profileと履歴log。
- `data/clean_run_attempts/`: 修復後runのlicense停止証跡。
- `repro/README.md`: Swarm表、外部／組込みCOMSOL、解析、検証の再現コマンド。
- `SOURCE_INVENTORY.md`: sourceごとの役割、証拠時点、採用可否。
- `CHART_MAP.md`: 9図の分析目的、表現契約、視覚QA。

正式再実行では、MPH rootの定式化readback、7 pairのproperty照合、
mode-aware activation、外部／組込み各3独立process、200／400要素raw profile、
同一条件、完全log・hashを必須gateとする。license更新後は
`repro/README.md`に従い、`repro/analysis_inputs.json`の暫定入力をclean runへ
差し替える。
