# Export Spec (Swarm -> COMSOL)

## 目的
Swarm計算結果を COMSOL Plasma Module に **Interpolation Function** としてインポート可能な「Spreadsheet形式（列形式）」に変換する。

## 出力ファイル仕様（COMSOL用）

### 1) comsol_transport.csv（1D: ε̄ → μ, D）
- Data format: Spreadsheet
- Number of arguments: 1
- 列（順番固定）
  1. epsbar_eV : 平均電子エネルギー ε̄ [eV]
  2. mu_m2_V_s : 電子移動度 μ [m^2/(V·s)]
  3. De_m2_s   : 電子拡散係数 D [m^2/s]
- 先頭に `%` で始まるコメント行（任意）を許可（COMSOLが無視する）

### 2) comsol_rates.csv（1D: ε̄ → k）
- Data format: Spreadsheet
- Number of arguments: 1
- 列（順番固定）
  1. epsbar_eV
  2..N: 反応ごとの速度係数 k
    - 2体反応（e + A → ...）: [m^3/s]
    - 3体反応（e + A + M → ...）: [m^6/s]
- 反応名→列名の対応は設定ファイルで定義（prefix など）

### 3) comsol_eedf.csv（2D: ε, ε̄ → f）
- Data format: Spreadsheet
- Number of arguments: 2
- 列（順番固定）
  1. eps_eV    : 電子エネルギー ε [eV]
  2. epsbar_eV : 平均電子エネルギー ε̄ [eV]
  3. f_eV_m32  : EEDF f(ε;ε̄) [eV^(-3/2)] を推奨
- 入力EEDFが別の規格（例: eV^-1 の確率密度）で出力される場合は `processing.eedf_normalization` で変換方法を指定する。

### 4) comsol_export_report.json
- 変換に用いた入力ファイルパス
- 出力点数
- ε̄, ε の範囲
- 欠損・重複・単位変換の有無
- EEDFの正規化チェック結果（任意）

## 推奨運用
- μ, D, k は **ε̄（平均電子エネルギー）** を独立変数にして補間関数化
- EEDF は **(ε, ε̄)** の2変数補間関数として1ファイルにまとめ、COMSOL側で2D補間させる（平均エネルギーごとにファイル分割する必要はない）
