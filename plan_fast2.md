# 高速化拡張計画（JIT/ベクトル化強化）  
本メモは現状コード（swarm_mc/配下）を一通り読んだうえでの高速化案。既存の安定化・統計改善（plan_fast.md）を保ちつつ、主に JIT とデータ配置最適化で wall-clock を短縮するロードマップ。デフォルト挙動は維持し、JIT は opt-in のフェーズ導入を想定。

## 0. ホットスポット予測
- 主ループ構造: `Simulation.advance_one_step` で `determine_timestep` → `determine_collisions` → `perform_collisions` → `TimeSeries.append_data`/EEDF 更新。粒子数 n_e、cross_section 数 m に対し `determine_collisions` は O(n_e * m)、`perform_collisions` は O(n_e) でメモリ再配置を伴う。
- `MonteCarlo.determine_collisions`: `searchsorted` + 線形補間 + `rand_vector` を m 行に複製→`cumsum`→比較。`collision_matrix` (m×n_e) の生成が最大の割当・演算コスト。n_e=1e5, m~50 の場合メモリ帯域支配。
- `MonteCarlo.perform_collisions`:  
  - 散乱角 (`scattering_angles`) で三角関数と乱数を 4*n_e 回呼ぶ。`unit_scattered_velocity` で near-axis 補正の分岐あり。  
  - ionization/attachment のブランチと新電子生成（`_apply_ionization_limits`）でマスク→`hstack`。Python 側で複数回の連結が発生しキャッシュ非効率。  
  - `utils.velocity_from_energy` を多数呼ぶため、ここも LUT/JIT での恩恵が大きい。
- `TimeAveragedEnergyDistribution.collect_histogram`: SST 到達後、毎ステップ `np.histogram` が走り、ステップ数 T に対し O(T*n_e_bin)。粒子数が大きい場合よりも、呼び出し回数がボトルネック。JIT or バッチ間引きで削減余地。
- `determine_timestep`: per-step 1回だがシミュレーション全体で数十万回呼ばれる。`np.interp` と log/rand が主。三角関数に比べ軽量だが、JIT で 10–20% 圧縮可能と見込む。
- I/O/ログ: `TimeSeries.append_data` は配列 append とキャッシュ無効化のみで軽量。`Output.check_sst` の polyfit は SST 直前にのみ効くので通常は非支配。

## 1. ベースライン計測を入れる（JIT前の確認）
- `Simulation.advance_one_step` 内で簡易タイマー（perf_counter）を入れ、`determine_timestep / determine_collisions / perform_collisions / collect_output_data` の4区間を集計。結果は `Output` 側でオプション保存（例: `output/perf_breakdown.csv`）。デフォルト off。
- 大規模 run 用に `run_sweep.py` から `--profile-per-step` フラグを通し、計測を opt-in。
- 実装状況: `simulation_settings.profile_per_step` を追加（デフォルト False）。有効時、`advance_one_step` で4区間の計測を取り、`<base>_perf_breakdown.csv` に保存。`run_sweep.py --profile-per-step` で上書き可。

## 2. データ配置・ベクトル化の事前整備
- cross-section テーブルを `float64` C-contiguous に強制し、`is_attachment/is_ionization/mass_ratios/thresholds` も `np.ndarray` 化した構造体にまとめる（JIT で渡しやすくする）。
- `determine_collisions` の `np.repeat` と `np.cumsum` を避け、`rand_vector` と累積和を 1-pass ベクトル比較に差し替え（配列割当を減らす）。  
  - 例: `np.searchsorted(cumprob, rand_vector, side="right")` で collision index を直接得る。
- `perform_collisions` の配列連結を two-pass で事前サイズ決定→1回の `np.empty` への埋め込みのみとし、`np.hstack` を排除（既に近いが new_e 分も含め確定サイズで埋める）。
- 実装状況: cross-section テーブルとしきい値/質量比/フラグを `float64` C-contiguous に変更（`GasMixture`）。`determine_collisions` は `np.repeat` を廃し、累積確率と乱数のブロードキャスト比較のみで index を算出。`perform_collisions` の two-pass 連結は未着手（次フェーズ）。

## 3. JIT 化ロードマップ（numba 想定、fall back 可）
### 3-1. ユーティリティ層
- `utils.velocity_from_energy / energy_from_velocity / acceleration_from_electric_field` を `@njit(cache=True, fastmath=True)` で包むヘルパを追加し、JIT が無い場合は従来関数を使用。
- `rng`：NumPy Generator は numba 非対応のため、JIT パスでは `numba.random` ベースの状態を別に持ち、`seed`/`run_label` から `SeedSequence`→`uint64` 配列を初期化。復元用の seed をログに記録。

### 3-2. 散乱角カーネル
- `scattering_angles` と `unit_scattered_velocity` を結合した `jit_scatter_directions(energies, velocities, iso, params, rng_state)` を `njit(fastmath=True)` 化。  
  - low/high 遷移窓（1–5 eV）は config 値を渡すだけの struct にし、`np.sin/np.cos` を `math` ベースに切替。  
  - near-axis 補正（sin_theta < 1e-12）も JIT 内で処理。

### 3-3. 衝突処理カーネル
- 入力: `positions (3,n)`, `velocities (3,n)`, `energies (n)`, cross-section struct（配列群）, config パラメータ（sharing, throttle 係数, max_e_max）。  
  - アウトプット: 新 pos/vel（null + survive + new_e）、衝突数/ion/att カウント。  
  - ionization share/energy loss/attachment マスクをすべて JIT 内で行い、Python 側は配列を受け取るだけにする。  
  - `max_ionization_growth_factor` と `throttle` は JIT 内で downsample 用の index を引く（`rng.choice` 相当を `numba` の random perm/selectionで置換）。

### 3-4. タイムステップ決定
- `determine_timestep` を JIT 化し、`np.interp` 部分を `searchsorted` + 線形補間に展開。`fastmath` 有効で `-log(u)` も cheap に。  
  - `timestep_min/max` の clamp は Python 側で行い、JIT 部分は純計算に専念。

### 3-5. ヒストグラム更新
- SST 後の `collect_histogram` を `njit(parallel=True)` 化した 1D bin 更新に差し替え。`np.histogram` 互換ロジックを実装し、バッファを in-place 加算。  
  - 既存 `maybe_rescale_bins` は Python のまま、カウントのみ JIT。

### 3-6. 切替制御
- Config に `use_jit: bool` を追加し、無効時は現行 numpy パスを使用。  
- JIT 依存の import 失敗時は警告を出し自動フォールバック。CI は両パスで smoke テストを行う。

### 実装状況（3.x）
- `Config.use_jit` を追加し、`to_dict` に反映。`SimulationRunner`/YAML 経由で opt-in 可能。
- `jit_kernels.py` を新設。numba が無い場合は自動フォールバック。  
  - `determine_timestep_jit`：pre-sampled rand を受け取り dt と trial_coll_freq を返却。  
  - `unit_scattered_velocity_jit`：phi/chi 用の乱数を Python 側で生成し、散乱角計算を JIT 化。  
  - `histogram_increment_jit`：EEDF ヒストグラムを in-place 加算。  
  - `jit_supported()` で numba 有無を判定。
- MonteCarlo: `determine_timestep` と `unit_scattered_velocity` で JIT パスを用意（rand は Generator で生成し再利用）。numba 無し/flag off なら従来パス。
- EEDF 収集: `collect_histogram` が `use_jit` を受け取り、numba があれば JIT で加算。
- MonteCarlo 内で velocity-from-energy を JIT ヘルパ経由に切替（エネルギー→速度の大量呼び出しを軽量化）。
- TODO: utils の JIT ラッパ、乱数 state の numba 一元化は未着手（numpy Generator から pre-sample する方針で暫定対応）。

### 4–6 ステータス
- 4 (LUT/近似): `use_velocity_lut` を追加し、max_coll_freq 計算時に energy_grid 上で速度 LUT を構築。`use_scatter_lut` を追加し、等方/異方散乱の cos_chi を逆 CDF LUT からサンプル可能に（異方は energy grid 上に 2D LUT）。`use_max_coll_lut` と `max_coll_lut_size` で max_coll_freq/period を圧縮 LUT 化し、timestep 補間を軽量化（デフォルト ON, size=256）。
- 5 (prange/並列・ログ間引き拡張): EEDF ヒストグラム加算を numba.prange で並列化。TimeSeries に `log_every_k` を追加（デフォルト1、任意に間引き可能）。sweep 並列実行を `experiment.parallel_workers` で選択可能に追加。
- 6 (検証/フォールバック強化): `use_jit` 有効かつ numba 不在時に警告し自動フォールバックするガードを追加。

### 未実装/残課題
- `perform_collisions` の計算部分: ion/attach 判定を JIT 化し、散乱角生成は JIT/LUT で並列化済み。生成・スロットリングの JIT 化と prange 並列化は未着手。
- 散乱 LUT の異方散乱対応: 実装済み（energy grid 上の 2D LUT）。低エネ窓ブレンドのパラメータ化は未着手。
- `max_coll_period` / `trial_coll_freq` のミニ LUT 化: 実装済み (`use_max_coll_lut`, `max_coll_lut_size` で圧縮)。
- RNG 再現性: numba RNG の初期 state/seed を execution.log に明示し、NumPy パスとの統計回帰テスト（平均エネルギー・w/DN を±1%で比較）。
- prange 化の拡大: 新規電子生成・スロットリングを含む `perform_collisions` 全体の並列 JIT 化と、乱数のスレッド分配。
- sweep 並列実行オプション（プロセスプール化）と GasMixtureCache の共有/再構築ポリシー。
- `use_jit`/LUT 周りの README 追記と example YAML 追加（非 JIT サンプルは現状のまま）。

## 4. 追加の LUT / 近似高速化
- 散乱角 LUT：エネルギー帯ごとに `cos_chi` の逆 CDF を 1D table にし、乱数から `np.interp` で取得するオプション。三角関数呼び出し数を半減させる（JIT パス・非 JIT パス両対応）。  
- `velocity_from_energy` LUT：`energy_grid`（cross-section gridと同じ）で `v` をプリコンピュートし、`np.interp` 再利用。`determine_collisions` と `perform_collisions` で同一 LUT を使う。  
- `max_coll_period` 近似：`determine_timestep` で使う `max_coll_period` を滑らかな区分線形に置き換え、補間を 1 乗算+1 加算で済ませるミニ LUT を導入。

## 5. 並列化とバッチ処理
- `numba.prange` を使い、電子ごとの独立計算（collision 選択、scatter、ionization energy 計算）を並列化。乱数は per-thread state を分配。  
- `TimeSeries.append_data` のログ間引き (`dn_step_skip`) をデフォルト1のまま維持しつつ、大規模 run 用に `log_every_k` を追加し、JIT と組み合わせて write コストを抑制。
- sweep 実行時のプロセス並列（run_sweep.py）をオプション化し、GasMixtureCache をプロセス間共有不可な場合は lazy 再構築。JIT とは独立。

## 6. 安全性・デバッグ
- JIT パスの再現性確認: `seed` と `run_label` をもとに numba RNG の初期 state をログ（`execution.log`）に残す。  
- `numpy` パスと JIT パスで同一 seed に対し統計量（平均エネルギー、w/DN）を比較する regression テストを追加。許容誤差を 0.5–1% に設定。  
- `use_jit` が True でも numba 未導入の場合は警告して自動で False に戻す。

## 7. 実装順（推奨）
1. 計測フックと cross-section 配列の C-contiguous 化・`determine_collisions` の再割当削減（非 JIT）。  
2. JIT ユーティリティ（velocity/energy、rng state 管理）と散乱角カーネルを導入、`use_jit` フラグを追加。  
3. 衝突処理カーネル全体の JIT 化（ionization/attachment/新電子生成まで）。  
4. `determine_timestep` JIT 化と LUT 追加。  
5. ヒストグラム JIT 化とロギング間引き。  
6. 並列化（prange）と sweep 並列実行オプション。  
7. 回帰テスト（JIT on/off）、性能ベンチ（S/M/L 粒子数 × 低/中/高 EN）を更新。

## 8. リスクと回避策
- numba 依存: インストール負担と import 失敗 → フラグで opt-in + フォールバック。  
- RNG 差異: numba の乱数と numpy Generator はビット一致しない → 統計一致をテストし、seed/entropy をログ。  
- メモリ使用: JIT カーネルで大きな一時配列を避ける実装を徹底（1-pass, prealloc）。  
- 互換性: 既存 config の挙動を変えないよう、デフォルトは非 JIT。`energy_sharing_*` やスロットリングのロジックは共通化。

## 9. 成果物イメージ
- `swarm_mc/jit_kernels.py`（numba 依存モジュール、なければ動的にスキップ）。  
- `Config` に `use_jit` 追加、`Simulation` がパスを切替。  
- プロファイル結果 (`perf_breakdown.csv`) と性能ベンチスクリプト（S/M/L テンプレ付き）。  
- ドキュメント: README に「JIT オプションの使い方」「numba が無い場合のフォールバック」を追記。
