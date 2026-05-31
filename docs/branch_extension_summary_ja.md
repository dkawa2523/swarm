# `swarm-7_complete` ブランチ拡張内容

## 1. 目的

このブランチでは、既存の電子スウォーム計算コードを「個別ソルバーの実行ツール」から「複数ソルバーを同じ条件で比較する製品インターフェース」へ整理した。中心は schema v2 の導入、ソルバー選択と物理要求の分離、direct PN 型 `multi_term` の実行経路、内部 Monte Carlo、外部 BOLSIG+ / MCIG 参照ベンチである。

設計上の主語はソルバーではなく、次の 3 つの比較対象で固定した。

| 分類 | canonical id | 位置づけ |
|---|---:|---|
| 二項近似 Boltzmann | `two_term` | 基準となるエネルギー空間 SG ソルバー |
| 多項近似 Boltzmann | `multi_term` | direct PN による角度クロージャ比較ソルバー |
| 粒子 Monte Carlo | `monte_carlo` | 外部 adapter または内部粒子追跡 backend |

従来の単一文字列による実行モード指定や一括実行指定は、比較対象と物理モデルを混同しやすい。そのため、新しい設定では必ず

```yaml
schema_version: 2
run:
  solvers:
    - id: two_term
    - id: multi_term
    - id: monte_carlo
```

の形で、実行するソルバーを明示する。

---

## 2. 全体像

実行の流れは、設定を直接ソルバーへ渡すのではなく、一度「solve plan」を作ってから実行する構成になった。

```mermaid
flowchart LR
    A[schema v2 YAML] --> B[typed config validation]
    B --> C[capability matrix]
    C --> D[solve plan]
    D --> E[solver execution]
    E --> F[case metadata]
    F --> G[tail / EEDF diagnostics]
    G --> H[comparison summary]
    H --> I[canonical CSV outputs]
```

主な責務分担は次の通り。

| 層 | 主なファイル | 役割 |
|---|---|---|
| 設定 | `electron_swarm/core/config.py` | schema v2 の型定義、旧 schema 拒否、内部設定への変換 |
| 能力表 | `electron_swarm/core/capabilities.py` | ソルバーごとの対応範囲を `exact` / `approximate` / `unsupported` で定義 |
| 計画 | `electron_swarm/orchestration/plan.py` | 要求物理と能力表から runnable / skipped / degraded を決定 |
| 実行 | `electron_swarm/orchestration/executor.py` | solve plan に従ってソルバーを呼び出し、製品 metadata を付与 |
| 比較 | `electron_swarm/orchestration/comparison.py` | 同一 case 間のスカラー値、EEDF、角度モデル整合性を比較 |
| 出力 | `electron_swarm/io/writers.py` | canonical CSV のみを書き出す |

---

## 3. schema v2 への移行

schema v2 の狙いは、ソルバーの指定と物理機能の要求を分離することである。たとえば磁場、角度散乱、e-e、tail 診断は `physics` 側の要求であり、`two_term` / `multi_term` / `monte_carlo` のどれを走らせるかとは別に扱う。

| 項目 | 旧設計 | 新設計 |
|---|---|---|
| schema | v1 / 互換挙動あり | `schema_version: 2` 必須 |
| 実行指定 | 単一の mode 文字列 | `run.solvers` |
| solver id | 旧 two-term / multi-term 公開名 | `two_term`, `multi_term`, `monte_carlo` |
| 一括指定 | `both`, `all` | 廃止。配列で明示 |
| 出力互換 | compatibility block と互換別名 | 廃止 |
| 物理要求 | ソルバー設定と混在 | `physics.*` に集約 |
| 未対応物理 | 暗黙に無視される危険 | `feature_policy` に従って fail / skip / fallback |

旧 public id は migration error として拒否される。

| 拒否される入力 | 理由 |
|---|---|
| 旧 mode 指定 | solver-mode と比較目的が混ざるため |
| 旧一括 solver 指定 | 暗黙実行になり、metadata が曖昧になるため |
| 旧 two-term 公開 id | product schema の solver id ではないため |
| 旧 multi-term 公開 id | product schema の solver id ではないため |
| 互換目的の別名出力 | canonical output 以外を製品面に残さないため |

---

## 4. solver capability と feature policy

各ソルバーは、要求された物理 feature に対して対応レベルを返す。対応レベルは、solve plan と結果 metadata の両方に残る。

| solver | e-neutral | angular | ionization source | e-e | magnetic | tail | bulk transport |
|---|---:|---:|---:|---:|---:|---:|---:|
| `two_term` | approximate | approximate | exact | approximate | unsupported | approximate | unsupported |
| `multi_term` | approximate | approximate | unsupported | approximate | unsupported | approximate | unsupported |
| `monte_carlo` | approximate | approximate | unsupported | unsupported | approximate | approximate | unsupported |

`feature_policy.unsupported` は次の挙動を持つ。

| policy | 挙動 | 用途 |
|---|---|---|
| `fail` | 要求された未対応物理があれば実行前に失敗 | 通常の厳密検証 |
| `skip_solver` | 未対応ソルバーだけを skip し、plan に理由を残す | 比較表を作りたい場合 |
| `approximate` | 明示 fallback が許可された場合のみ degraded として実行 | 手動の妥協を記録したい場合 |

重要なのは、未対応物理を「黙って無視しない」点である。実行可否は次のような論理で決まる。

$$
\mathrm{runnable}(s) =
\begin{cases}
\mathrm{true} & \text{all requested features are supported or explicitly degraded} \\
\mathrm{false} & \text{policy is } \mathrm{skip\_solver} \\
\mathrm{error} & \text{policy is } \mathrm{fail}
\end{cases}
$$

solve plan の 1 行は、おおむね次の情報を持つ。

| 列 | 意味 |
|---|---|
| `solver` | canonical solver id |
| `runnable` | 実行対象になったか |
| `skipped` | policy により skip されたか |
| `degraded` | 近似または fallback を含むか |
| `effective_angular_scattering` | 角度散乱をどう扱ったか |
| `effective_electron_electron` | e-e をどう扱ったか |
| `effective_magnetic_field` | 磁場をどう扱ったか |
| `capability_*` | 能力表上の対応レベル |

---

## 5. direct PN `multi_term`

`multi_term` は、旧 multi-term 公開名を外し、product id として `multi_term` に統一された。現在の product method は次の 2 種類である。

| method | 入力 | 実装範囲 | metadata |
|---|---|---|---|
| `pn_closure_direct` | ordinary integral cross sections | SG reduction と同じ kinetic block を使う direct PN | `ordinary_integral_xs_closure=true`, `exact_dcs_based=false` |
| `pn_dcs` | normalized Legendre moment table | moment table を角度 moment source として direct PN に投入 | `angular_moment_source=moment_table` |

### 5.1 PN 展開

角度分布は Legendre 展開として扱う。

$$
f(E,\mu) \simeq \sum_{\ell=0}^{L} F_\ell(E) P_\ell(\mu)
$$

ここで、`lmax = L` である。`lmax: 1` は、二項近似の SG reduction と一致すべき基準経路として保護されている。

`lmax: 1` の gate では、正規化 EEDF

$$
\int_0^\infty F_0(E)\,dE = 1
$$

と平均エネルギー

$$
\langle E \rangle = \int_0^\infty E F_0(E)\,dE
$$

が `two_term` と十分近いかを確認する。実装上は `two_term` の結果をコピーせず、共有 kinetic block から direct solve する。

### 5.2 higher `lmax`

`lmax > 1` では coefficient-space の疎ブロック系を解く。現在の範囲は限定的である。

| 条件 | 現在の扱い |
|---|---|
| DC 電場 | 対応 |
| `B=0` | 対応 |
| axisymmetric `m=0` | 対応 |
| ordinary integral XS angular closure | 対応 |
| raw DCS angle table | 未対応 |
| arbitrary crossed E-B PN | 未対応 |
| `l>0` inelastic source | sink-only |

higher moment の減衰は、角度 moment $m_\ell(E)$ を用いて次の形で組み立てる。

$$
\nu_\ell(E)
= N\,v(E)\,\sigma_{\mathrm{total}}(E)\left[1-m_\ell(E)\right]
  + \nu_{\mathrm{inelastic\ loss}}(E),
\quad \ell \ge 2
$$

`ell=1` では共有 momentum frequency を使う。

$$
\nu_1(E) = \nu_m(E)
$$

このため、ordinary integral cross sections だけで走る `multi_term` は「exact DCS multi-term」ではなく、明示的に angular-closure PN solver として metadata に記録される。

### 5.3 `pn_dcs` と moment table

`pn_dcs` は raw DCS を直接読むのではなく、正規化済み Legendre moment table を読む。

```csv
energy_eV,m0,m1,m2,m3
0,1,0.0,0.0,0.0
10,1,0.2,0.04,0.008
```

入力条件は次の通り。

| 条件 | 検証 |
|---|---|
| `energy_eV` | 非負、有限、単調増加 |
| `m0` | すべて 1 |
| `m1...mL` | `[-1, 1]` の有限値 |
| extrapolation | `error` のみ |
| provenance | `dcs_derived`, `model_derived`, `unknown` |

`provenance: dcs_derived` の場合だけ、結果 metadata の `exact_dcs_based` が真になり得る。

---

## 6. 角度散乱モデル

角度散乱は solver option ではなく、`physics.angular_scattering` で要求する。

| model | closure | moment source | MC sampler |
|---|---|---|---|
| `isotropic` | `zero` | isotropic closure | 対応 |
| `momentum_power` | `power` | ordinary integral XS closure | 未定義 |
| `maxent_p1` | `maxent` | ordinary integral XS closure | 対応 |
| `moment_table` | `table` | normalized moment table | 未定義 |

ordinary XS からは第 1 moment を

$$
m_1(E) = 1 - \frac{\sigma_m(E)}{\sigma_{\mathrm{total}}(E)}
$$

として推定する。ただし、これは DCS の完全な復元ではない。`momentum_power` は

$$
m_\ell(E) = m_1(E)^\ell
$$

を仮定し、`maxent_p1` は $m_1$ を満たす maximum-entropy 型の $P_1$ 分布から高次 moment を計算する。

PN と MC の比較では、同じ角度モデルかどうかを別列で判定する。

| 比較列 | 意味 |
|---|---|
| `same_angular_model` | reference と candidate の角度 metadata が一致したか |
| `angular_model_status` | `match` / `mismatch` / `unknown` |
| `angular_sampler_treatment` | MC 側の sampler または fallback の扱い |
| `angular_model_mismatch_reason` | 比較を弱める理由 |

---

## 7. 内部 Monte Carlo backend

`monte_carlo` は外部 adapter に加えて、限定的な内部粒子 backend を持つようになった。

| backend | 実行方式 | 主な用途 |
|---|---|---|
| `external` | command または Python API を呼ぶ | 既存 MC コードとの接続 |
| `internal` | リポジトリ内の簡易粒子追跡 | 磁場 metadata と smoke 比較 |

内部 backend は DC 電場と DC 磁場を扱い、速度更新には Boris pusher を使う。電子の運動は

$$
m_e \frac{d\mathbf{v}}{dt}
= -e\left(\mathbf{E} + \mathbf{v}\times\mathbf{B}\right)
$$

であり、Boris 法では半ステップ加速、磁場回転、半ステップ加速で速度を更新する。純磁場では速度ノルムが保存されるため、磁場積分器の基本的な健全性をテストできる。

磁場要求時の扱いは次の通り。

| solver | 磁場要求の扱い |
|---|---|
| `two_term` | 未対応。policy により fail / skip / explicit fallback |
| `multi_term` | 現在の axisymmetric PN では未対応 |
| `monte_carlo` external | 出力 metadata の整合性を検査 |
| `monte_carlo` internal | `boris_lorentz_push` として実行 |

---

## 8. e-e と ionization source

e-e は schema v2 の物理 feature として整理された。

| model | 対象 | 内容 | transport |
|---|---|---|---|
| `none` | 全 solver | e-e なし | 通常 |
| `relaxation_postprocess` | Boltzmann 系 | EEDF と rates を後処理で更新 | stale と明示 |
| `fp_energy` | Boltzmann 系 | f0 energy-space の簡易 FP 緩和 | transport 再計算扱い |

`relaxation_postprocess` では輸送係数を再計算しないため、metadata に

```text
electron_electron_transport_stale = true
```

を残す。`fp_energy` は f0 のみを対象にしたエネルギー空間緩和で、完全な Landau operator ではない。

ionization source は `two_term` で明示モデルを選べる。

| energy sharing | 意味 | 現在の対応 |
|---|---|---|
| `equal` | 等分配 | 既定 |
| `primary_secondary` | secondary electron energy を明示 | `two_term` 対応 |
| `loss_only` | source なしの loss-only | `two_term` 対応 |

`multi_term` の non-default ionization source は、higher-l source term が定義されていないため、現在は policy-handled unsupported として扱う。

---

## 9. tail metrics

高エネルギー tail は、高しきい値反応の rate に強く効く。平均エネルギーだけでは tail の収束を判定できないため、結果 metadata と rates 出力に tail 指標を追加した。

tail 確率は

$$
P_{\mathrm{tail}}(E_t)
= \int_{E_t}^{\infty} F(E)\,dE
$$

反応 $r$ に対する tail rate fraction は

$$
R_{\mathrm{tail},r}
=
\frac{
  \int_{E_t}^{\infty} \sigma_r(E)\,v(E)\,F(E)\,dE
}{
  \int_0^{\infty} \sigma_r(E)\,v(E)\,F(E)\,dE
}
$$

として計算する。

| metadata | 意味 |
|---|---|
| `tail_probability` | tail threshold 以上の EEDF 確率 |
| `tail_rate_fraction_max` | 反応 rate への tail 寄与の最大値 |
| `dominant_tail_process` | tail 寄与が最大の反応 |
| `energy_grid_tail_status` | `ok` / `warning` / `insufficient` |
| `tail_fraction` | rates CSV の各反応別 tail 寄与 |

---

## 10. canonical outputs

出力は製品向けに固定され、互換目的の別名出力は削除された。

| 出力ファイル | 内容 |
|---|---|
| `<base>_summary.csv` | solver ごとの主要輸送量と compact metadata |
| `<base>_rates.csv` | 反応 rate、frequency、power loss、tail fraction |
| `<base>_eedf.csv` | normalized EEDF と EEPF |
| `<base>_solver_plan.csv` | runnable / skipped / degraded と capability |
| `<base>_comparison_summary.csv` | 比較が有効な場合の scalar / EEDF / angular 比較 |

比較で使う相対差は

$$
\Delta_{\mathrm{rel}}(x)
=
\frac{x_{\mathrm{candidate}}-x_{\mathrm{reference}}}
     {\max(|x_{\mathrm{reference}}|,10^{-300})}
$$

EEDF の L1 差は reference grid 上で

$$
\|F_c-F_r\|_1
=
\sum_i |F_c(E_i)-F_r(E_i)|\,\Delta E_i
$$

として出力する。

summary に残す metadata は意図的に絞っている。内部反復回数や過剰な診断列ではなく、第三者が比較条件を判断するための列に限定した。

| metadata group | 代表列 |
|---|---|
| solver method | `meta_solver_method`, `meta_physics_level` |
| angular | `meta_angular_model`, `meta_angular_moment_source`, `meta_exact_dcs_based` |
| direct PN | `meta_direct_pn_operator`, `meta_lmax` |
| e-e | `meta_electron_electron_treatment`, `meta_electron_electron_transport_stale` |
| magnetic | `meta_magnetic_field_treatment`, `meta_field_integrator` |
| tail | `meta_tail_probability`, `meta_energy_grid_tail_status` |

---

## 11. 外部 BOLSIG+ / MCIG 参照ベンチ

BOLSIG+ と MCIG は product solver id ではなく、benchmark reference として扱う。

| reference | 役割 | 入力 |
|---|---|---|
| BOLSIG+ | two-term 弱電離 swarm 参照 | `electron_swarm_reference_csv` または BOLSIG text |
| MCIG | Monte Carlo swarm 参照 | `electron_swarm_reference_csv` または MCIG CSV |

参照 EEDF は、どの入力形式でも正規化済み $F(E)$ にそろえる。

$$
\int F(E)\,dE = 1
$$

EEPF 入力の場合は

$$
F(E) = \mathrm{EEPF}(E)\sqrt{E}
$$

としてから正規化する。

追加された主なベンチは次の通り。

| benchmark | 目的 | 主な出力 |
|---|---|---|
| `ar_eedf_consistency` | `two_term` と `multi_term lmax=1` の direct gate | EEDF metrics, failure analysis, report |
| `ar_bolsig_plus_equivalence` | Ar 条件で BOLSIG+ と製品結果を比較 | summary, EEDF metrics, failure analysis |
| `ar_mcig_reference` | MCIG 参照との比較 | summary, confidence, failure analysis |
| `ar_bolsig_mcig_triage` | BOLSIG+ / MCIG / 製品 solver の差分原因を分類 | triage matrix, metrics, report, optional plot |

triage は単純な pass/fail ではなく、次のような原因分類を出す。

| likely cause | 解釈 |
|---|---|
| `code_regression_multi_term_lmax1` | direct PN lmax=1 gate の実装不具合 |
| `code_or_bolsig_input_mismatch` | BOLSIG+ 入力、EEDF convention、projection、rate 畳み込みの不一致 |
| `possible_two_term_approximation_limit_or_angular_model_difference` | two-term と MC の物理モデル差の可能性 |
| `multi_term_higher_l_model_issue` | higher-l damping、source/sink、tail の問題 |
| `mc_uncertainty_limited` | MC 統計が不足し、差分を断定できない |
| `angular_model_mismatch` | 角度モデル metadata が一致していない |

---

## 12. テスト拡張

テストは互換性維持ではなく、schema v2 の製品挙動を守る方向に書き直された。

| テスト領域 | 代表内容 |
|---|---|
| schema | v2 accepted、v1 / old id / unknown fields rejected |
| solve plan | feature policy、capability、skip / degraded metadata |
| outputs | canonical CSV、互換目的の別名出力が出ないこと |
| direct PN | `lmax=1` gate、higher-l smoke、moment table |
| physics | angular model、e-e、ionization source、magnetic field |
| MC | internal backend、same-as-physics metadata、Boris pusher |
| benchmark | Ar EEDF consistency、外部参照 ingest、triage |
| release shape | `electron-swarm` entry point、不要 runtime project の削除 |

`pyproject.toml` では package name を `electron-swarm` にし、CLI entry point を

```toml
[project.scripts]
electron-swarm = "electron_swarm.runner:main"
```

として公開している。`matplotlib` は本体依存から外し、plot extra に移した。

---

## 13. 削除・整理されたもの

このブランチは後方互換より製品設計を優先している。主な削除・整理は次の通り。

| 削除対象 | 理由 |
|---|---|
| 旧 standalone MC runtime project | product `monte_carlo` adapter / internal backend に統合するため |
| `swarm_comsol_exporter` | solver comparison product の中心機能ではないため |
| legacy output directories | 生成物を repository に残さないため |
| 旧 unified configs | schema v2 config に置き換えるため |
| 旧 two-term public path | `two_term` へ canonical 化 |
| 旧 multi-term public path | `multi_term` へ canonical 化 |
| 過剰な diagnostic output | 製品利用者に必要な metadata に絞るため |

---

## 14. 現在の制限

現時点で意図的に未対応として残している領域は次の通り。

| 領域 | 現在の状態 |
|---|---|
| raw DCS angle table | 未実装。moment table のみ |
| full spherical harmonics `Y_lm` | 未実装。PN は axisymmetric `m=0` |
| arbitrary crossed E-B PN | 未対応 |
| RF / time-dependent solver | 未対応 |
| finite-k hydrodynamic transport | schema validation と policy handling のみ |
| full Landau e-e operator | 未対応 |
| Coulomb particle-particle MC | 未対応 |
| generated state-resolved / superelastic chemistry | YAML 生成 helper は削除 |

この制限は単なる未完成リストではなく、比較結果を誤って読まないための境界条件である。特に ordinary integral XS だけで得た `multi_term` 結果を、DCS ベースの厳密 multi-term 解と呼ばない点が重要である。

---

## 15. レビュー時の確認ポイント

第三者レビューでは、次の順に見ると差分の意味を追いやすい。

| 順序 | 確認対象 | 見るべき点 |
|---:|---|---|
| 1 | `README.md` | 製品インターフェースが schema v2 / canonical solver になっているか |
| 2 | `docs/product_schema_v2.md` | public schema が旧 id を残していないか |
| 3 | `electron_swarm/core/config.py` | 旧 YAML の拒否と typed config の整合性 |
| 4 | `electron_swarm/orchestration/plan.py` | 未対応物理が policy 通りに処理されるか |
| 5 | `electron_swarm/solvers/multi_term/direct.py` | direct PN が `two_term` 結果コピーではないか |
| 6 | `electron_swarm/io/writers.py` | canonical output だけが出るか |
| 7 | `tests/test_product_*` | product behavior がテストで固定されているか |
| 8 | `tools/benchmark_ar_*` | 外部参照が fake data なしで扱われるか |

最小確認コマンドは次である。

```powershell
py -3 -m pytest -q
```

Monte Carlo の重い経路や外部参照を分ける場合は、pytest marker を使って段階的に実行する。

```powershell
py -3 -m pytest -q -m "not slow and not mc"
py -3 -m pytest -q -m regression
```
