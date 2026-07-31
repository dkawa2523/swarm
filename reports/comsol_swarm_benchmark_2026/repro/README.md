# 外部Swarm–COMSOL benchmark再現手順

以下はrepository rootからPowerShellで実行する。Pythonは本報告作成時のlocal venvを例示する。COMSOL executableは環境に合わせて変更する。

```powershell
$PY = ".\.tmp\venv\Scripts\python.exe"
$COMSOL = "C:\Program Files\COMSOL\COMSOL64\Multiphysics_copy1\bin\win64\comsolbatch.exe"
```

本手順のCOMSOL対象はApplication Libraryの1D Ar DC glow dischargeであり、
電極間全域を解く。履歴外部経路は、COMSOLのLEA流体方程式と残余chemistryを
保持し、選択した輸送・Townsend係数をSwarm表で置換するselected-channel
hybridである。

## 1. schema v2 Swarm表

設定は`schema_version: 2`、`run.solvers: [{id: two_term}]`、純Ar、
293.15 K、13.3322 Pa、0.0444–5000 Tdの34点、adaptive energy-grid上限
20 keVである。workflow databaseのcase diagnosticsに対して、tail probability、
edge/peak、grid-limitとiterative solver convergenceを別々に判定する。

```powershell
& $PY -m swarm_workflow.cli sweep reports\comsol_swarm_benchmark_2026\repro\workflow_argon_comsol_benchmark_2026.yaml

& $PY -m swarm_workflow.cli build-tables `
  outputs\comsol_swarm_benchmark_2026\swarm_qualified_20keV.sqlite `
  --output outputs\comsol_swarm_benchmark_2026\tables_solver_qa_20keV `
  --source two_term

& $PY -m swarm_workflow.cli export-comsol `
  outputs\comsol_swarm_benchmark_2026\tables_solver_qa_20keV\mixture_0000 `
  --output outputs\comsol_swarm_benchmark_2026\bundle_solver_qa_20keV
```

同梱databaseの監査結果は、tail/grid gateが34/34点、solver convergence gateが
23/34点、両方を満たすall-case gateが23/34点である。したがって本bundleは
tail-audited候補であり、「全34点solver-converged」とは記述しない。

生成後、COMSOLで暗黙外挿させないため、mean energy 0、0.05、0.1、0.2、0.35 eVを持つ監査可能な派生floor bundleを作る。輸送と弾性行は最下計算点を保持し、非弾性rate／Townsendは0とする。この処理はSwarmの直接計算値ではない。

```powershell
& $PY reports\comsol_swarm_benchmark_2026\repro\prepare_low_energy_bundle.py `
  outputs\comsol_swarm_benchmark_2026\bundle_solver_qa_20keV `
  outputs\comsol_swarm_benchmark_2026\bundle_solver_qa_20keV_low_energy
```

scriptは既存output directoryを拒否する。再生成時は旧directoryを意図的に別名へ退避し、新しい空pathを指定する。

## 2. 外部Swarm表のCOMSOL適用

単発のdry runは、bundle、mapping、生成JavaをCOMSOL起動前に検査する。
mappingでは平均電子エネルギー定式化を明示する。本報告の履歴経路は
`LocalEnergyApproximationE`（LEA）であり、平均電子エネルギー方程式を解く。

```powershell
& $PY -m swarm_workflow.cli run-comsol `
  Model\maps\positive_column_external.yaml `
  --bundle outputs\comsol_swarm_benchmark_2026\bundle_solver_qa_20keV_low_energy `
  --comsol $COMSOL `
  --dry-run
```

正式runでは、1回のworkflow内でapply、表property verify、定式化に応じた
activation verify、run、exportを独立COMSOL batchとして実行する。bundleは
7 lookup quantities、7 X-Y pairs、14 numerical arraysを持つ。そのproperty
readbackは「書き込まれたこと」の監査であり、「7量が同時に方程式へ作用した
こと」の証明ではない。LEAでは
\(E/N\rightarrow\bar{\varepsilon}\) 表は非活性で、輸送4量とeir2/eir4
Townsendが候補活性集合となる。空間profileで直接照合するのは、安定したCOMSOL
変数を持つ換算移動度、換算拡散、eir2/eir4 Townsendの4量である。energy
mobility／diffusivityはproperty-levelで監査し、将来は個別摂動によるfunctional
activation testを行う。

Local Field Approximation（LFA）は別の方程式クラスである。LFAでは
\(E/N\rightarrow\bar{\varepsilon}\) 表が活性になる一方、平均エネルギー
方程式を解かないためenergy mobility／diffusivityは非活性となる。LEAとLFAを
同じactivation gateで評価しない。

次のrunnerはworkflow自体も独立Python processで3回起動し、各stage時間、
wall end-to-end、raw log、profile、hashを非破壊のrun archiveへ保存する。
`--output-root`は存在しないpathを指定する。

```powershell
& $PY reports\comsol_swarm_benchmark_2026\repro\run_external_benchmark.py `
  --mapping Model\maps\positive_column_external.yaml `
  --bundle outputs\comsol_swarm_benchmark_2026\bundle_solver_qa_20keV_low_energy `
  --comsol $COMSOL `
  --output-root outputs\comsol_swarm_benchmark_2026\external_200elem_runs_clean `
  --repetitions 3 `
  --pressure-Pa 13.3322 `
  --gas-temperature-K 293.15 `
  --mesh-elements 200
```

compilerの終了コードが0でも、標準出力に`Failed to compile java file`、
`Compilation failed`等があれば失敗する。各stageはrun専用class directoryを
使用し、事前class、空class、Javaより古いclass、変更されない既存output MPHを
拒否する。provenanceにはJava/class、mapping、input/output MPH、
bundle/configのSHA-256、COMSOL version/build、実行電圧列、Plasma rootの
`MeanElectronEnergyModel` readbackを記録する。期待定式化との不一致は明示的に
失敗させる。

400要素の感度runも別archiveへ保存する。

```powershell
& $PY reports\comsol_swarm_benchmark_2026\repro\run_external_benchmark.py `
  --mapping Model\maps\positive_column_external.yaml `
  --bundle outputs\comsol_swarm_benchmark_2026\bundle_solver_qa_20keV_low_energy `
  --comsol $COMSOL `
  --output-root outputs\comsol_swarm_benchmark_2026\external_400elem_runs_clean `
  --repetitions 1 `
  --pressure-Pa 13.3322 `
  --gas-temperature-K 293.15 `
  --mesh-elements 400
```

各`run_NN/`には`workflow_stdout.txt`、`workflow_stderr.txt`、command、
provenance、archive manifestに加え、`workflow/`以下へresult CSV、run summary、
effective mapping、生成Java、apply/verify/run/exportの全logとprovenance、
verification CSVをcopy-onlyで保存する。rootの`timing_runs.csv`と
`timing_summary.json`は中央値・最小・最大を持つ。`solve_s`という履歴列名が
残る場合も、その値はrun-stage COMSOL batch（class起動、model load、電圧
continuation、saveを含む）のwall timeであり、PDE kernelの純粋なsolution time
ではない。本文では「run-stage時間」と呼ぶ。`total_s`は
apply+verify+run+export、`workflow_wall_end_to_end_s`は独立subprocess全体で
ある。pressureはexport CSVの全点、temperatureは`pes1`の`T` property、meshは
再build後の`comp1/mesh1/edg1/dis1.elemcount`をread backする。タグ欠落や
不一致は明示失敗し、行数からmeshを推定しない。

本報告作成時（2026-07-31）は、clean class生成後の最初のCOMSOL applyがlicense error -10（`Product has expired`）で停止した。sandboxと通常ユーザー環境の双方の証跡は`../data/clean_run_attempts/`に保存した。したがって200／400要素の新raw profileを取得済みとは扱わない。

## 3. COMSOL組込みBoltzmann参照のクリーン反復

200要素を3回、独立batch processで実行する。

```powershell
& $PY reports\comsol_swarm_benchmark_2026\repro\run_builtin_benchmark.py `
  --input-mph Model\positive_column_1d_boltzmann.mph `
  --comsol $COMSOL `
  --output-root outputs\comsol_swarm_benchmark_2026\builtin_200elem_runs `
  --repetitions 3 `
  --mesh-elements 200 `
  --voltage-V 200 `
  --pressure-Pa 13.3322 `
  --gas-temperature-K 293.15
```

400要素の感度runを少なくとも1回実行する。正式なmesh間統計が必要なら`--repetitions 3`とする。

```powershell
& $PY reports\comsol_swarm_benchmark_2026\repro\run_builtin_benchmark.py `
  --input-mph Model\positive_column_1d_boltzmann.mph `
  --comsol $COMSOL `
  --output-root outputs\comsol_swarm_benchmark_2026\builtin_400elem_runs `
  --repetitions 1 `
  --mesh-elements 400 `
  --voltage-V 200 `
  --pressure-Pa 13.3322 `
  --gas-temperature-K 293.15
```

runnerは保存済み最終solver `sol3.runAll()`を使い、`comp1/mesh1/edg1/dis1`の`elemcount`を200または400へ明示設定する。各`run_NN/`には12列raw profile、Java/class、input/output MPHとlogのhash、COMSOL version/build、compile/batch/solve/end-to-end時間、条件一定性QA、provenanceを保存する。output root、run directory、事前classを再利用しない。

## 4. 比較入力の差し替え

`analysis_inputs.json`が解析の唯一のrouting fileである。正式run後は次を更新する。

- `comparison.external.path`: 外部200要素raw profile。
- `comparison.reference.path`: 組込み200要素raw profile。
- `mesh_profiles`: 200要素と400要素のraw profile。
- `mesh_summary_fallback`: 両raw profileが揃ったら`null`。
- `runtime`: `scope,path,seconds,repetition`列を持つCSV。`apply`、`verify`、`solve`、`export`、`solve_only`、`end_to_end`を区別する。
- `evidence_status`: 正式gateをすべて満たした場合に限り正式状態へ変更。

現在のdefaultは、2026-07-31の最新Swarm表と、2026-07-14の暫定COMSOL profile／単発時間を明示的に組み合わせる。

## 5. 解析、notebook、図

```powershell
$env:MPLBACKEND = "Agg"
& $PY reports\comsol_swarm_benchmark_2026\repro\benchmark_analysis.py
& $PY reports\comsol_swarm_benchmark_2026\repro\build_notebook.py
```

1つ目は`../data/derived/`へ19 CSVを出力する。内容は区分線形profileに対する
厳密空間積分のL2、積分比、相関、伝導電流RSD、幾何学的窓感度、領域別誤差、
高電界regime、履歴activation、tail/solver convergence、operating point、
Townsend source representation、partial electron-energy、fluid applicability、
runtime、QA、source hashである。`../figures/`へ9図の投稿用PNGとvector PDFを
出力する。旧「節点値二乗の台形積分」と点数等重みの指標は感度／追跡用として
別列に残す。2つ目は同じ解析を先頭から実行した
`../benchmark_analysis_executed.ipynb`を作る。raw入力はread-onlyで開き、
source、派生CSV、PNG、PDFをSHA-256で追跡する。

## 6. 検証

```powershell
$env:MPLBACKEND = "Agg"
& $PY -m pytest -q
```

正式採用条件は、clean compile、完全なraw log、最新class、root定式化readback、
定式化別のactive/inactive契約、7 lookup pairs（14 arrays）のwritten-property
照合、可能な空間照合、全34点のtail/grid gateとsolver convergence gate、
入力provenance、外部・参照各3反復、200/400要素raw profile、条件一致である。
また、\(E/N>500\ {\rm Td}\)の空間割合とdirect-source寄与を必ず報告し、
lookup上限がdrift-diffusion近似の妥当性を保証しないことを明記する。いずれかを
満たさない場合、2026-07-14値は「暫定技術ベンチマーク」のままとする。
