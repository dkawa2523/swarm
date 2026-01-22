# EEDF 最適化スクリプトの使い方

本リポジトリの `EEDF_optimze.py` は、`outputs/*/eedf_table.csv` を読み込み、前処理＋混合分布フィッティングを行い、係数やグラフを出力します。前処理・モデル設定は `EEDF_optimize.yaml` で調整できます。

## 前提（Windows）

`uv` が必要です。未インストールの場合は次のいずれかで導入してください。

```powershell
winget install AstralSoftware.uv
```

または:

```powershell
scoop install uv
```

## セットアップ（Windows / PowerShell, uv）

```powershell
uv venv
uv sync
```

## 実行例
- すべての `outputs/*/eedf_table.csv` を処理:  
  `uv run python EEDF_optimze.py --root outputs --config EEDF_optimize.yaml`
- 特定ファイルのみ:  
  `uv run python EEDF_optimze.py --table outputs\\ar_n2_en_scan\\eedf_table.csv --config EEDF_optimize.yaml`
- 上限分布数や試行回数を上書き:  
  `uv run python EEDF_optimze.py --table ... --max-components 2 --num-starts 2`

出力先: 各実験フォルダの `plots/eedf_fits/` に PNG と `fit_results.csv` を保存（`--save-dir` で変更可）。

## YAML パラメータ（`EEDF_optimize.yaml`）

### preprocess
| キー | デフォルト | 意味 |
| --- | --- | --- |
| `energy_min`, `energy_max` | null | エネルギーの下限/上限でクリップ |
| `log_mad_sigma` | 4.0 | log10(EEDF) の外れ値除去（MAD スケール） |
| `log_floor` | 1e-30 | ログ取得前の固定フロア |
| `log_floor_factor` | 10.0 | 最小正 EEDF に掛ける動的フロア（最低 1 桁上からフィット） |
| `smooth.enabled` | true | log 空間の Savitzky-Golay 平滑化 |
| `smooth.window` | 11 | 平滑化窓（奇数） |
| `smooth.polyorder` | 2 | 平滑化多項式次数 |
| `resample_points` | 800 | 一様グリッドに再サンプル（0 で無効） |
| `tail_area_rel_threshold` | 1e-4 | 高エネ側累積面積がこの比率以下ならゼロ化 |
| `tail_min_energy_quantile` | 0.9 | ゼロ化を検討するエネルギー分位点 |
| `tail_rel_keep` | 1e-4 | 最大値のこの比率を下回った以降の尾部をトリム |
| `tail_keep_min_quantile` | 0.75 | トリム開始を遅らせるための最低分位点 |
| `tail_keep_min_points` | 100 | トリム後にも確保する最小データ点数 |

### model
| キー | デフォルト | 意味 |
| --- | --- | --- |
| `max_components` | 3 | 混合コンポーネントの上限（1〜3） |
| `num_starts` | 6 | コンポーネント数ごとの初期値試行回数 |
| `random_seed` | 12345 | 乱数シード |
| `alpha_bounds` | [-2.0, 4.0] | べき指数 α の下限/上限 |
| `beta_bounds` | [0.2, 4.0] | 伸張指数 β の下限/上限 |
| `amplitude_min` | 1e-30 | 振幅の下限 |
| `e_c_max_factor` | 10.0 | 特性エネルギー上限 = max(E)*factor |
| `e_c_mean_energy_factor` | 5.0 | 特性エネルギー上限 = meanE*factor |
| `shift_bounds_relative` | [-0.3, 0.3] | シフトの相対下限/上限（max(E) 比） |
| `shift_bounds_absolute` | null | 絶対値でシフト範囲を固定したい場合の [min, max] |

### selection
| キー | デフォルト | 意味 |
| --- | --- | --- |
| `criterion` | AICc | モデル選択指標（固定で AICc） |

### plot
| キー | デフォルト | 意味 |
| --- | --- | --- |
| `y_min`, `y_max` | null | y 軸の固定下限/上限（null で自動） |
| `y_margin_factor` | 0.1 | パーセンタイル範囲をこの割合だけ拡張 |
| `y_clip_percentiles` | [0.1, 99.9] | ylim 算定前に除外するパーセンタイル |

### weights
| キー | デフォルト | 意味 |
| --- | --- | --- |
| `low_energy_boost.amplitude` | 1.0 | 低エネルギー側の重み増幅（0 で無効） |
| `low_energy_boost.scale_rel_to_max` | 0.2 | 増幅の減衰スケール（max(E) 比） |

## ヒント
- 高エネルギー側の単粒子ノイズをより強く落とすときは `tail_rel_keep` を小さく、`tail_keep_min_quantile` を上げてください。
- 低エネルギー精度を上げたい場合は `low_energy_boost.amplitude` を大きくするか、`resample_points` を増やしてください。***
