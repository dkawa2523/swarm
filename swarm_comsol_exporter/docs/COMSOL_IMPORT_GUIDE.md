# COMSOL Import Guide (Plasma Module)

このガイドは、`swarm-comsol-export` が生成するCSVを COMSOL Multiphysics (6.3/6.4) の **Plasma Module** に読み込ませるための手順をまとめたものです。

## 重要な前提
- Plasma Module の EEDF “Function” では **2D interpolation function** を使い、
  - x-data: 電子エネルギー ε [eV]
  - y-data: 平均電子エネルギー ε̄ [eV]
  を与えます。
- Plasma Module UIではエネルギーが V 表示になる場合がありますが、内部的には eV として扱われます（COMSOL公式ブログの注記参照）。

## 1) Interpolation Function を作る
Model Builder で以下を追加します。

### A. μ, D（輸送係数）
1. Global Definitions > Functions を右クリック > Interpolation
2. Settings:
   - Data source: File
   - Filename: `comsol_transport.csv`
   - Data format: Spreadsheet
   - Number of arguments: 1
3. Data Column Settings（自動認識されるはず）で以下を確認:
   - epsbar_eV 列が Argument
   - mu_m2_V_s / De_m2_s が Function values
4. Functions テーブルに2行追加して、同じファイルから2つの関数を定義:
   - 例: 関数名 `mu_e`, Position in file = 1（epsbarの次の列が第1関数列）
   - 例: 関数名 `De`,   Position in file = 2（第2関数列）

### B. 反応速度係数 k
同様に `comsol_rates.csv` を読み込み、反応ごとの関数を追加します。
- Number of arguments: 1
- Position in file を列位置に合わせて指定

### C. EEDF（2D関数）
1. Global Definitions > Functions > Interpolation
2. Settings:
   - Data source: File
   - Filename: `comsol_eedf.csv`
   - Data format: Spreadsheet
   - Number of arguments: 2
3. Data Column Settings を確認:
   - eps_eV が Argument
   - epsbar_eV が Argument
   - f_eV_m32 が Function values
4. Functions テーブル:
   - 関数名: `eedf_fun`
   - Position in file = 1（2つのArgument列の次の列）

## 2) Plasma interface へ割り当て
### EEDF の指定
Plasma interface の main node（例: “Plasma Model”）で
- Electron Energy Distribution Function: Function
- Function: `eedf_fun`

### μ, D の指定
Plasma interface の電子輸送係数設定で
- μ: `mu_e(plas.ebar)`
- D: `De(plas.ebar)`
のように設定（`plas` は Plasma interface の tag）。

## 3) 反応速度係数 k の指定
Electron impact reaction を「Lookup table / user-defined」の形式で指定し、
- k: `k_ion(plas.ebar)` のように設定します。

※Cross section を入力しない場合でも、rate coefficient を lookup table で与える運用は可能です。
