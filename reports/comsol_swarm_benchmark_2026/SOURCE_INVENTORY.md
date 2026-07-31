# 証拠ソース一覧と採用方針

## 報告範囲と判定

- **対象:** COMSOL 6.4の1D Ar DC glow-discharge Application Library modelを
  基礎とする外部Swarm表経路。cathode fall、dark-space／positive-column構造、
  anode側を含む全電極間を対象とし、未接続の2D GEC CCPは対象外。
- **数値vintage:** Swarm正式入力候補は2026-07-31、COMSOL profileと単発時間は
  2026-07-14。
- **参照:** COMSOL組込みBoltzmannは比較参照であり、実験ground truthではない。
- **判定:** `share_with_caveats`。2026-07-31 clean runはCOMSOL license
  error -10でmodel apply前に停止した。履歴runはcompiler false-positiveを含み、
  Java/class対応を証明できない。
- **Swarm候補表:** 34/34点でtail probability、edge/peak、20 keV grid-limit
  gateを満たすが、iterative solver convergenceは23/34点である。tail適合と
  solver収束を混同せず、all-case formal gateは未達とする。
- **formulation監査:** 基礎MPH、履歴外部表適用MPH、組込みBoltzmann MPHは
  いずれもroot property
  `ElectronProperties/MeanElectronEnergyModel=LocalEnergyApproximationE`。
  したがって履歴比較はLFAではなく、平均エネルギー方程式を保持したLEAである。

## 同梱した数値証拠

| ID | パス | 用途 | 採用状態と制約 |
|---|---|---|---|
| `formal_swarm_bundle_candidate` | `data/formal_swarm_bundle_solver_qa_20keV/` | schema v2、純Ar、293.15 K、13.3322 Pa、34点、0.0444–5000 Td、adaptive上限20 keV | tail/grid 34/34、solver convergence 23/34。COMSOL activation未完了で正式入力ではない |
| `formal_swarm_database_candidate` | `data/formal_swarm_database_solver_qa_20keV.sqlite` | 34点のcase diagnostics、solver・tail・energy-grid metadata | `tail_convergence_audit.csv`の一次数値source |
| `historical_swarm_conditions` | `data/archive_2026-07-14/historical_swarm_source_config.yaml`、`historical_swarm_provenance.json` | 履歴表の300 K、13.3 Pa条件 | COMSOL profileの293.15 K、13.3322 Paと不一致 |
| `historical_intended_bundle` | `outputs/validation/188_low_mean_energy_floor_bundle/mixture_0000/` | 2026-07-14 run summaryが指すintended table | activation residualの比較元。履歴compiler defectのため正式なclass-input証明ではない |
| `external_profile_history` | `data/archive_2026-07-14/positive_column_results.csv` | 外部表経路の801点profile | 暫定比較のみ |
| `builtin_profile_history` | `data/archive_2026-07-14/builtin_200V_fresh.csv` | 組込みBoltzmann参照の801点profile | 暫定比較のみ |
| `run_summary_history` | `data/archive_2026-07-14/positive_column_run_summary.json` | 条件、table range、単発stage時間 | timer scopeがPDE kernel単体とは限らない |
| `comparison_history` | `data/archive_2026-07-14/comparison_metrics.csv`、`comparison_summary.json` | 旧point-count-weighted比較 | 追跡用supplementのみ |
| `diagnostic_history` | `data/archive_2026-07-14/diagnostic_summary.json` | 旧200/400電流summary | 400 raw profile欠落。grid convergenceには不十分 |
| `activation_history` | `data/archive_2026-07-14/function_verify_summary.json` | `sw_meanE`、`sw_muN` control-point | property/function確認でありequation-level activation証明ではない |
| `compile_false_positive` | `data/archive_2026-07-14/compile_result.json`、`compile_stdout.txt` | return code 0とcompile失敗文字列の同居 | stale classを排除できない決定的制約 |
| `clean_attempts` | `data/clean_run_attempts/` | fresh Java/class、COMSOL起動、license停止 | model apply前に停止 |
| `derived_metrics` | `data/derived/` | 19 CSV: exact piecewise-linear metric、lookup/tail、operating point、Townsend表現、partial energy、fluid applicability、activation、runtime、QA、hash | 本文の派生値。formal gate未充足をmetadataに保持 |
| `executed_notebook` | `benchmark_analysis_executed.ipynb` | 派生CSV・図の再現 | 先頭から実行しerror outputがないことを最終QAする |
| `figures` | `figures/` | 本文9図 | PNG/PDF、source hash、evidence vintageをmetadata化 |

## formulationとactive-setの証拠

| Model | root mean-energy mode | solver/dependent-variable evidence | 解釈 |
|---|---|---|---|
| `Model/positive_column_1d.mph` | `LocalEnergyApproximationE` | vendor基礎モデル | 平均エネルギーPDEを解くLEA |
| `Model/work/positive_column_external_applied.mph` | `LocalEnergyApproximationE` | 履歴外部表適用MPH | mean-energy tableはstored/inactive。輸送4量とeir2/eir4 Townsendが候補active |
| `Model/positive_column_1d_boltzmann.mph` | `LocalEnergyApproximationE` | saved `sol3/std3`のsolve-forに`comp1.En`、`comp1.Ne`、`comp1.plas.F0`、heavy species、potential、circuit state | mean-energy PDEに加えてEEDF DOFを連成するbuilt-in参照 |

外部経路は、Swarmが選択した輸送・eir2/eir4 Townsendを供給し、COMSOLがLEA
流体方程式、境界、heavy species、残余chemistry／energy treatmentを保持する
selected-channel hybridである。bundleは7 lookup quantities、7 X-Y pairs、
14 numerical arraysを持つ。これらの書込みと、7量が同時に支配方程式へ作用する
ことは区別する。LEAでは \(E/N\rightarrow\bar{\varepsilon}\) pairがinactiveで
6量が候補active、LFAではenergy mobility／diffusivityがinactiveで5量が候補
activeとなる。いずれもfunctional activationにはroot readbackに加え、個別係数
摂動と解応答が必要である。

## 派生データ契約

| ファイル | 主な内容 | 注意 |
|---|---|---|
| `comparison_metrics.csv` | 区分線形profileに対するexact relative L2、積分比、相関、旧指標 | potential積分比は共通gauge下のdescriptive quantity |
| `current_uniformity.csv` | external/reference、full/中央窓のconductive-current RSD | conservation residualではない |
| `current_rsd_sensitivity.csv` | retained window幅に対するRSD | 物理的bulk/sheath区分ではない |
| `regional_error_attribution.csv` | cathode-side 10%、central 80%、anode-side 10%のsquared-error share | 幾何学窓 |
| `regime_metrics.csv`、`regime_impact.csv` | \(E/N>500\) Tdのnode/length割合とquantity積分寄与 | 500 TdはCOMSOL公式のtypical guidanceで、硬いinvalid閾値ではない |
| `fluid_applicability_diagnostics.csv` | 最大 \(n_e/N\)、net-flux-speed ratioの最大値と空間support | 弱電離proxyとboundary/locality heuristic。Knudsen数またはhard validity testではない |
| `operating_point_diagnostics.csv` | source voltage、gap endpoint voltage、ballast current、内部potential最大、絶対potential variation | 履歴final-time profileからの診断。source voltageをgap voltageと同一視しない |
| `townsend_source_representation_audit.csv` | \((\alpha/N)N|J_e|/e\) の再構成と \(Nn_e k\) counterfactual | matched COMSOL rerunではなく、rate form優越の証拠ではない |
| `partial_electron_energy_audit.csv` | signed \(J_eE\)、eir2/eir4 threshold-weighted loss | flux、elastic、残余反応、boundary、transientを欠くpartial audit |
| `tail_convergence_audit.csv` | 34点ごとのsolver、tail probability、edge/peak、grid-limit gate | tail/grid 34/34、solver 23/34、all-case 23/34 |
| `historical_activation_audit.csv` | intended table対exported profile residual | 履歴compiler defectを伴うdiagnostic |
| `runtime_summary.csv`、`runtime_speedup.csv` | run-stage／end-to-endの単発時間 | \(n=1\)、同一timer scope未立証 |
| `qa_summary.csv`、`source_hashes.csv` | formal gate、有限値、行数、条件、hash | failを隠さない |

## 実装・再現ソース

| ID | パス | 役割 |
|---|---|---|
| `schema_v2_config` | `data/formal_swarm_bundle_solver_qa_20keV/executed_argon_comsol_benchmark_2026.yaml` | `schema_version: 2`、`run.solvers: [{id: two_term}]`、34点、20 keV上限 |
| `low_energy_policy` | `repro/prepare_low_energy_bundle.py` | mean-energy下限の明示floor。直接Swarm計算点とは区別 |
| `mapping` | `Model/maps/positive_column_external.yaml`、`swarm_workflow/comsol_mapping.py` | mean-energy formulationとtable mapping |
| `java_apply` | `swarm_workflow/comsol_java.py` | property設定とroot formulation readback |
| `external_runner` | `repro/run_external_benchmark.py`、`swarm_workflow/comsol_positive_column.py`、`comsol_verify.py` | clean class、stale拒否、stage分離、provenance、反復 |
| `builtin_runner` | `repro/run_builtin_benchmark.py` | built-in saved solverの独立process反復 |
| `analysis` | `repro/benchmark_analysis.py`、`analysis_inputs.json` | exact spatial integration、regime／activation audit、図 |
| `notebook_builder` | `repro/build_notebook.py` | 実行済みnotebook生成 |

## 外部一次資料

| ID | 一次資料 | 使用目的 |
|---|---|---|
| `comsol_1d_model` | [COMSOL DC Glow Discharge, 1D](https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.positive_column_1d/positive_column_1d.html) | 方程式、反応、境界、mesh、cathode fall、1D制約 |
| `comsol_1d_boltzmann_model` | [COMSOL DC Glow Discharge, 1D, Coupled with a Boltzmann Equation](https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.positive_column_1d_boltzmann/positive_column_1d_boltzmann.html) | 組込みBoltzmann参照の連成構成とsolver context |
| `comsol_interface` | [The Drift Diffusion Interface](https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.02.html) | LEA/LFAの方程式差 |
| `comsol_lookup` | [Drift–Diffusion Model](https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.04.html) | reduced transport table、mean-energy specification |
| `comsol_theory` | [Introduction to Drift–Diffusion Theory](https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.17.html) | 10 mtorr、typically <500 Td、弱電離・collisional条件 |
| `comsol_guide` | [Plasma Module User’s Guide 6.4](https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/PlasmaModuleUsersGuide.pdf) | COMSOL property／unit／terminal-current契約 |
| `bolsig_paper` | [Hagelaar and Pitchford (2005)](https://doi.org/10.1088/0963-0252/14/4/011) | two-term Boltzmann、fluid係数 |
| `flux_bulk` | [Dujko and White (2008)](https://doi.org/10.1088/1742-6596/133/1/012005) | flux/bulk transport convention |
| `nonlocal` | [Dujko et al. (2008)](https://doi.org/10.1088/0022-3727/41/24/245205) | steady-state Townsend regimeにおけるnonhydrodynamic transport |
| `lxcat_paper` | [Pancheshnyi et al. (2012)](https://doi.org/10.1016/j.chemphys.2011.04.020) | cross-section／swarm provenance |

## 本文不採用または補足限定

- 出典・比較方向・run状態が本文contractと異なる旧summary画像。
- manual-edit MPH候補と出典不明slides。
- 旧8-panel図はsupplement閲覧だけ。
- 2D GEC CCP、実験精度、他gas、`multi_term`／`monte_carlo`のCOMSOL結果。
- 200/400 summaryだけからのGCIまたはformal grid-convergence主張。

## 正式採用に残るgap

1. 有効license下でapply、formulation readback、verify、run、exportを完走する。
2. 外部・built-inを同一timer contract、独立process、各3反復で測る。
3. 293.15 K、13.3322 Pa、同一mesh／反応在庫／solver stateへ可能な限り統制する。
4. 200/400（望ましくは800）要素のraw profileを両経路で保存する。
5. 7 lookup pairs（14 arrays）のwritten-property readbackとmode-aware active
   setを照合し、各tableを±1%摂動するfunctional activation testを行う。
6. \(n_{\rm Ar^+}\)、charge density、Debye length、terminal／displacement current、
   boundary fluxをexportし、物理的region、charge／particle balanceを定義する。
7. \(E/N>500\) Tdのlength/source shareとlocality diagnosticを再計算する。
8. 現候補でsolver未収束の11点を収束させ、tail/gridとsolverの両gateを34/34点で
   満たすbundleを再生成する。
9. 必要なら実験または独立文献profileでvalidationする。これは本稿の現範囲外。

これらを満たすまで、profile差は履歴・記述的比較、8.62倍／4.73倍は単発の
候補比としてのみ引用する。
