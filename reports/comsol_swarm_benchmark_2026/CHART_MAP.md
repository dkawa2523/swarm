# Chart map and visual QA contract

本報告はsingle-columnのportable HTMLを正本とする。青は外部Swarm表経路、
濃灰はCOMSOL組込み参照、goldは範囲・感度・regime警告に限定する。色だけに
依存せず、線種、marker、hatch、直接labelでも系列を区別する。解析artifactは
19 derived CSVと9 static figuresを契約とする。

## Technical-report structure map

| Required role | Visible section | Contract |
|---|---|---|
| Title | `External Swarm Tables in a COMSOL 1D Argon DC Glow Discharge: A Formulation-Aware Provisional Benchmark` | 電極間全域のDC glow discharge、LEA、暫定性をtitleで示す |
| Summary | `技術要約`、`English abstract` | 最初にformulation監査、次に数値結果、最後に採用制約 |
| Scope | `問題設定とモデル対照` | cathode fallからanode側までの1D Ar全電極間。2D GEC CCP、実験validation、solver優越を除外 |
| Method | `外部表契約`、`比較指標`、`再現手順` | selected-channel hybrid、7 quantities／7 X-Y pairs／14 arrays、written propertyとequation-level activationを分離 |
| Results | profile、direct sources、conductive current、operating point、source representation、partial energy、regional error、high-field／fluid regime、runtime | 全図・表にevidence date、母数、比較方向、単位を付ける |
| Limitations | `物理解釈の境界`、`正式採用gate` | high-field、履歴class、条件差、mesh、n=1を可視化 |
| Third-party input | `第三者Swarm adapter要件` | source provenanceとcanonical solver idを分離 |

## Figure contracts

| Figure | Analytical question | Form and source | Required interpretation |
|---|---|---|---|
| `fig01_workflow_boundary` | SwarmとCOMSOLの責任境界は何か | schema-v2 YAML → SQLite → tables → bundle → MPH apply/readback → COMSOL PDE/profile | 外部経路はselected-channel hybrid。7 lookup quantities＝7 X-Y pairs＝14 arraysの書込みと7量の同時活性化は同義でない |
| `fig02_lookup_domain` | LEAでactiveなlookup argumentは表域とboundary policyのどこを使うか | panel A: inactive \(E/N\rightarrow\bar{\varepsilon}\)。panel B–F: mean-energy軸のactive transport、energy transport、eir2/eir4 Townsend。最新表と履歴profile | 低mean-energyのgold bandは人工的floor policy、最初の直接Swarm計算点未満の範囲を明示する。500 Td超は別のregime警告で、table coverageは物理妥当性を保証しない |
| `fig03_profile_state` | state/field profileのshapeとamplitudeはどう違うか | density、mean energy、potential、\(E/N\) | potentialと\(E/N\)は微分関係で独立証拠ではない。potential積分比はgauge依存 |
| `fig04_sources_currents_mesh` | direct sourceとconductive currentの差はどこにあるか | eir2 direct excitation、eir4 direct ionization、\(J_e\)、\(J_e+J_{\mathrm{Ar}^+}\)、RSD、旧200/400 summary | 「全励起・全電離」「全電流」と呼ばない。中央80%は幾何学的窓で、bulk/sheath同定ではない |
| `fig05_metric_summary` | quantity別の局所差と総量差は何か | 区分線形厳密積分L2と積分比 | 分母はCOMSOL組込み参照で、ground truthではない。signed current比を明示 |
| `fig06_runtime_breakdown` | timing差はどの工程にあるか | apply/verify/run/exportとend-to-end | 各path \(n=1\)。42.914 sはPDE kernelのsolution timeでなくrun-stage wall time |
| `fig07_spatial_robustness` | 誤差はどの幾何学領域に集中し、current RSDは窓選択に頑健か | cathode-side 10%／central 80%／anode-side 10%のsquared-error share、retained-window RSD | 領域名は幾何学定義。物理的cathode fall/positive column/sheathの同定ではない |
| `fig08_activation_regime` | 履歴表は実際のprofileへ機能したか、高電界域は結果へどれだけ寄与するか | intended-table対exported-profile activation residual、\(E/N>500\) Tdの長さ・direct-source share | mean-energy table不一致はLEAで非活性というMPH readbackと整合。compiler provenance欠陥のため履歴監査であり正式activation証明ではない |
| `fig09_fluid_applicability` | 弱電離近似とboundary/localityの運用上の注意域はどこか | \(n_e/N\) profileと \(\lvert J_e\rvert/[e n_e\sqrt{2e\bar{\varepsilon}/m_e}]\)、0.1／1 guide | density ratioは弱電離proxy、net-flux-speed ratioはdrift・diffusion・boundary fluxを含むheuristicであり、Knudsen数やhard validity testではない |
| native metric charts | exact valuesをportable HTMLで読めるか | `comparison_metrics.csv` | static figureのsemantic table fallback |
| native runtime chart | scope別時間を読めるか | `runtime_summary.csv` | distributionを示すwhiskerを作らない |

## Numerical chart rules

- profileを節点間で線形とみなし、\(\int q\,dx\)、\(\int q^2dx\)、
  \(\int q_1q_2dx\)を区分ごとに厳密積分する。
- 「節点値二乗の台形積分」はcurrent spike感度を過大化し得るため、感度列として
  残し主値にしない。
- 旧point-count-weighted指標は追跡用supplementだけに置く。
- eir2/eir4はdirect reaction channelであり、net excitation／ionization sourceと
  labelしない。
- \(J_e+J_{\mathrm{Ar}^+}\)はconductive current。変位電流またはterminal currentを
  含む「全電流」とlabelしない。
- current RSDは空間一様性指標であり、charge-conservation residual、ambipolarity、
  grid convergenceの証明ではない。
- potentialの積分比は共通ground gauge下のdescriptive quantityと注記する。
- total/conductive-current correlationは参照がほぼ一定なら解釈しない。
- \(E/N>500\ {\rm Td}\)はCOMSOL公式の「typically below」guidanceに基づく
  regime flagで、硬いvalid/invalid閾値ではない。
- `tail_convergence_audit.csv`ではtail/grid gateとiterative solver convergence
  gateを分離する。現候補は34/34対23/34であり、前者だけを「全点収束」と
  labelしない。
- Townsend-formと\(Nn_e k\)の差はpostprocess representation sensitivityであり、
  matched COMSOL rerunまたはrate form優越の図として示さない。
- partial electron-energy auditはeir2/eir4 threshold lossだけを含み、完全な
  energy balanceとlabelしない。
- runtime \(n=1\)はpoint/barのみ。median、range、whiskerを暗示しない。

## Final-context QA checklist

- [ ] 全axisに量と単位がある。
- [ ] 外部／参照、分子／分母、signed／absoluteがcaptionだけで判別できる。
- [ ] 2026-07-31 Swarm表と2026-07-14 COMSOL profileのmixed vintageを明示する。
- [ ] provisional figureに`2026-07-14 archived`を表示する。
- [ ] active/inactiveはroot formulation readbackと整合する。
- [ ] central 80%をbulk、外側10%をsheathと呼ばない。
- [ ] eir2/eir4を全sourceと呼ばない。
- [ ] conductive currentをterminal/total currentと呼ばない。
- [ ] figure 8のactivation errorはlog scaleとfail thresholdを明記する。
- [ ] figure 9の0.1 guideを著者の運用値、1をscale comparisonとして明記する。
- [ ] tail/grid 34/34とsolver convergence 23/34を別々に表示する。
- [ ] narrow HTML、print CSS、English abstract、data URI図、source dialogを検証する。
- [ ] artifactのcanonical source path、query、metric definition、hashが揃う。

## Omitted-chart reasons

- 実験一致図は実験dataがないため作らない。
- timing uncertaintyは独立反復がないため描かない。
- quasi-neutrality、Debye-length、charge-conservation residual、global particle
  balanceは必要raw variablesが履歴exportにないため、将来の必須診断として残す。
- 200/400だけ、かつ400 raw欠落のためGCIやformal grid-convergence rateを出さない。
- `multi_term`、`monte_carlo`は今回COMSOLへ適用していないため空categoryを作らない。
