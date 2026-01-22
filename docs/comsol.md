**Swarm（電子スウォーム）計算結果を COMSOL（Plasma Module／流体モデル）へ入力する際に必要になる物理量（単位含む）と、現実的に扱いやすいデータ形式**を、COMSOL側の具体操作と紐づけてテーブルで整理
（EEDF は「平均電子エネルギーごとに別ファイル」より、**1つの2D補間関数として渡す**のがCOMSOLの想定に合います。([COMSOL Documentation][1])）

---

## 0) Swarm → COMSOL 出力ファイルの使い方（本リポジトリの exporter 前提）

Swarm の後処理で生成される COMSOL 用 CSV は以下です。

- `comsol_transport.csv`  
  1 列目 = ε̄[eV], 2 列目 = μ, 3 列目 = D  
  `%` コメント行に **列マップ** と **元列名** が入ります。
- `comsol_rates.csv`  
  1 列目 = ε̄[eV], 2 列目以降 = 反応ごとの k(ε̄)  
  列名は **断面積データの PROCESS 記述**を反映します（例: `E + Cl2 -> E + E + Cl+ + Cl`）。  
  `%` コメント行に **列マップ** が入ります。
- `comsol_eedf.csv`  
  1 列目 = ε[eV], 2 列目 = ε̄[eV], 3 列目 = f(ε,ε̄)  
  `%` コメント行に **列マップ** が入ります。
- `comsol_export_report.json`  
  列マップ、元列名、反応列一覧など **識別用メタ情報** をまとめています。

### COMSOL 側での読み込みの具体例（Column map を使う）

1. **Global Definitions → Functions → Interpolation** を作成  
2. **Data source = File / Data format = Spreadsheet** を指定  
3. **Number of arguments** を指定  
   - transport/rates は **1**  
   - eedf は **2**  
4. **Position in file** を column map と合わせて指定  
   - 例: `comsol_rates.csv` の 2 列目が “E + Cl2 -> E + E + Cl+ + Cl” に対応
5. Electron Impact Reaction では **Use lookup table** を選び、作成した関数を参照

> どの列がどの反応かは、`comsol_rates.csv` の `% Column map:` と  
> `comsol_export_report.json` の `rates.column_map` / `rates.rate_columns` で確認できます。

## 1) Swarm→COMSOLに渡す物理量・単位・推奨CSV形式（まとめ表）

> **基本方針（推奨）**
>
> * COMSOLでそのまま使えるように、**すべて「平均電子エネルギー ϕ（eV）」を軸**に揃える
> * EEDF は **f(ε, ϕ)** の **2D補間（ε: 電子エネルギー, ϕ: 平均電子エネルギー）**として1ファイル化
> * μ, D, k は 1D（ϕのみ）として、**1ファイル複数列**にまとめてもOK（COMSOL補間関数は1ファイル複数関数列を扱える）([COMSOL Documentation][2])

| 項目                                   | COMSOLでの用途（どこに効くか）                                        | 単位（COMSOL想定）                               | 依存軸          | 推奨CSV（列構成例）                             | COMSOL側の入力先（代表）                                                                                            | 根拠（該当箇所）                                                                            |
| ------------------------------------ | --------------------------------------------------------- | ------------------------------------------ | ------------ | --------------------------------------- | ---------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| 平均電子エネルギー ϕ                          | LFA（Local field approximation）などで **E/N→ϕ** 関係が必要になる場合がある | ϕ: eV（Plasma Interfaceの説明）                 | E/N または ϕ    | `EN_Vm2, meanE_eV`（必要な場合）               | **Plasma Model** または **Drift Diffusion Model** の *Mean electron energy specification* で *Use lookup table* | LFAでは **E/N と ϕ の関係提供が必要**([COMSOL Documentation][1])                               |
| 電子移動度 μe                             | 電子フラックス／導電率／輸送                                            | m²/(V·s)                                   | ϕ（eV）        | `meanE_eV, mu_m2_Vs`                    | **Drift Diffusion Model** で *Electron mobility μe*（Specify/Lookup）                                         | μe の単位、Lookupが「ϕ(eV) 対 係数」([COMSOL Documentation][3])                               |
| 電子拡散係数 De                            | 電子密度分布の拡散（数値安定にも影響）                                       | m²/s                                       | ϕ（eV）        | `meanE_eV, De_m2_s`                     | **Drift Diffusion Model** で *Electron diffusivity De*（Specify/Lookup）                                      | De の単位、Lookupが「ϕ(eV) 対 係数」([COMSOL Documentation][3])                               |
| （任意だが推奨）電子エネルギー移動度 μen / エネルギー拡散 Den | 「電子エネルギー密度」方程式側の輸送係数（モデル設定によって必要）                         | μen: m²/(V·s), Den: m²/s                   | ϕ（eV）        | `meanE_eV, mu_en_m2_Vs, Den_m2_s`       | **Drift Diffusion Model** の *Specify All* または *Use Lookup Tables*                                          | μen/Den の単位と入力項目([COMSOL Documentation][3])                                         |
| 電子衝突反応の反応速度係数 k(ϕ)（Rate coefficient） | Electron Impact Reaction の生成項・損失項                         | m³/s（COMSOLの Electron Impact Reaction で明記） | ϕ（eV）        | `meanE_eV, k_m3_s`（反応ごとにファイル or 列を分ける）  | **Electron Impact Reaction** → *Specify reaction using: Use lookup table* → *Rate coefficients*            | 「ϕ(eV) vs k（m³/s）」が明記([COMSOL Documentation][4])                                    |
| （任意）Townsend 係数 α(ϕ)                 | DC系などで Townsend 形式を使う場合                                   | m²                                         | ϕ（eV）        | `meanE_eV, alpha_m2`                    | Electron Impact Reaction → *Townsend coefficients*                                                         | 単位と形式が明記([COMSOL Documentation][4])                                                 |
| EEDF f(ε,ϕ)                          | **EEDFを外部から与える**（Function）ときに使用（※断面積からkを自動計算する場合などに特に重要）  | 典型表記：eV^(-3/2)（COMSOL blog に例）             | ε（eV）, ϕ（eV） | `eps_eV, meanE_eV, f_eV_m32`（= eV^-3/2） | Plasma Interface → *Electron Energy Distribution Function* → *Function*（2D補間）                              | x=ε(eV), y=ϕ(eV) の指定([COMSOL Documentation][1])／f(ε)の表式と単位感（eV^-3/2の例）([COMSOL][5]) |
| 電子衝突のエネルギー損失 Δe（しきい値等）               | エネルギー損失項（電子エネルギー方程式に効く）                                   | V（= eV相当として扱われる記述）                         | 反応ごと         | YAML/反応定義側に定数として持つ（CSVではなくパラメタ）         | Electron Impact Reaction → *Collision type*（Excitation/Ionization）                                         | Δe の入力が必要（単位V）([COMSOL Documentation][4])                                           |
| （中性反応）Arrhenius反応（参考）                | 中性化学・壁反応など（Swarmの範囲外）                                     | 反応次数依存                                     | T（K）, 濃度     | COMSOLの Reaction ノードのArrheniusで入力       | Heavy Species Transport / Reaction                                                                         | Arrhenius式と単位の説明([COMSOL Documentation][6])                                         |

### 重要：EEDFは「平均エネルギーごとに分ける」より 2D でまとめるのが自然

COMSOLの Plasma Interface は、EEDF を **“2次元の補間関数”として使う**場合に **x-data=電子エネルギー(eV)、y-data=平均電子エネルギー(eV)** と明記しています。([COMSOL Documentation][1])
つまり、EEDFを「条件ごとに別ファイル」にするより、**1つのCSVで (ε,ϕ,f) の点群**として持たせて補間させるのが、COMSOL側の設計に合います（Interpolation関数は最大3変数まで対応）。([COMSOL Documentation][2])

---

## 2) COMSOL側の具体操作（どこをどう設定するか）

| やりたいこと                                          | COMSOLツリー位置（代表）                                                     | 具体操作                                                                                                                                | 入力するCSVの形             | 根拠（該当箇所）                                               |
| ----------------------------------------------- | ------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- | --------------------- | ------------------------------------------------------ |
| μe(ϕ), De(ϕ) をモデルに入れる（Lookup方式）                 | **Drift Diffusion Model** → *Electron Density and Energy*           | *Electron transport properties* を **Use Lookup Tables** にして、μ/De等を **ϕ(eV)に対するテーブルとして読み込み**                                         | 2列（ϕ, 値）または 1列ϕ + 複数列 | Lookupが「ϕ(eV) vs 係数」([COMSOL Documentation][3])        |
| μe(ϕ), De(ϕ) をモデルに入れる（Interpolation Function方式） | **Global Definitions** → *Functions* → **Interpolation**            | Data source=File, Data format=Spreadsheet, Number of arguments=1（ϕ）で **μ関数・D関数**を作成 → Drift Diffusion側の入力欄に `mu_fun(ϕ)` のように参照      | 1列目=ϕ、2列目以降=関数列（複数可）  | Spreadsheet列の規約・複数関数列、最大3引数([COMSOL Documentation][2]) |
| Electron Impact Reaction に k(ϕ) を入れる            | **Heavy Species Transport / Plasma** → **Electron Impact Reaction** | *Specify reaction using* = **Use lookup table** → *Source coefficient data* で **Rate coefficients** を選び、**ϕ(eV) vs k(m³/s)** を Load | 2列（ϕ, k）              | ϕ(eV) vs k(m³/s) が明記([COMSOL Documentation][4])        |
| Townsend係数で入れる（必要な場合）                           | Electron Impact Reaction                                            | *Rate constant form* を Townsend にして、**ϕ(eV) vs α(m²)** をLoad                                                                        | 2列（ϕ, α）              | Townsend係数の単位と形式([COMSOL Documentation][4])            |
| EEDFを外部から与える（2D補間）                              | **Plasma Interface** → *Electron Energy Distribution Function*      | EEDFの選択で **Function** を選び、事前に作った **2D Interpolation function** を指定（x=ε, y=ϕ）                                                        | 3列（ε,ϕ,f）“ロング形式”推奨    | x=ε(eV), y=ϕ(eV) の指示([COMSOL Documentation][1])        |
| EEDFが本当に必要か判断する                                 | Plasma Interface（EEDF項）                                             | **断面積からkを自動計算する場合**はEEDFが必要。**kをlookupで与えるだけ**ならEEDFは必須ではない場合が多い（設定に依存）                                                             | —                     | 「断面積でsource係数定義ならEEDF選択が必要」([COMSOL Documentation][1]) |
| Local field approximation を使う場合の追加入力            | Plasma Interface / Drift Diffusion Model（Mean electron energy spec） | LFAでは **E/N と ϕ の関係**を Expression か Lookup table で与える（Boltzmann解かないなら必須）                                                            | 2列（E/N, ϕ）            | LFAの説明と必要条件([COMSOL Documentation][1])                 |

---

## 3) COMSOLの「Interpolation Function」に合わせたCSV設計（EEDF含む）

COMSOLの補間関数（Interpolation）は、ファイルが **Spreadsheet形式**なら次のルールが基本です：

* **先頭の列：引数（1〜3列まで）**、その後ろが **関数値列（複数列可）**([COMSOL Documentation][2])
* **% 行はコメント扱い**([COMSOL Documentation][2])
* 1ファイルに複数関数列を置ける（例：μとDを同じファイルに）([COMSOL Documentation][7])

### 推奨CSV例

**(A) 輸送係数（1D）：mu/D/… を1ファイルにまとめる**（引数=ϕ だけ）

```csv
% mean electron energy [eV], mu [m^2/(V*s)], De [m^2/s]
meanE_eV,mu_m2_Vs,De_m2_s
1.0,0.12,0.015
2.0,0.10,0.012
...
```

* COMSOL側では Interpolation Function を複数作り、**Position in file** で列を切り替える運用ができます（またはLookup方式）。([COMSOL Documentation][2])

**(B) 反応係数（1D）：反応ごとに1ファイル（または1ファイル複数列）**
Electron Impact Reaction の *Use lookup table* は、**ϕ(eV) vs k(m³/s)** が想定です。([COMSOL Documentation][4])

```csv
% mean electron energy [eV], k [m^3/s]
meanE_eV,k_m3_s
1.0,1.2e-16
2.0,3.4e-16
...
```

**(C) EEDF（2D）：(ε,ϕ,f) のロング形式（1ファイル）**
Plasma Interface の “Function” は **x=ε(eV), y=ϕ(eV)** を要求します。([COMSOL Documentation][1])
EEDF の表式は一般に ϕ^(-3/2) を含むため、**fの次元は eV^(-3/2) 表記がよく使われます**（COMSOL blog例）。([COMSOL][5])

```csv
% eps [eV], meanE [eV], EEDF f(eps,meanE) [eV^(-3/2)]
eps_eV,meanE_eV,f_eV_m32
0.0,2.0,0.0000
0.5,2.0,0.1234
1.0,2.0,0.1150
...
0.0,3.0,0.0000
0.5,3.0,0.1010
1.0,3.0,0.0980
...
```

---

## 4) 参照した一次情報（URLコピー用）

※ ここは **そのままコピペできるようにコードブロック**で列挙します（要望対応）。

```text
[Plasma Interface: EEDF Function (x=electron energy, y=mean electron energy), Local field approximationの説明]
https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_plasma.09.02.html

[Electron Impact Reaction: Use lookup table, Rate coefficients: mean electron energy(eV) vs k (m^3/s), Townsend係数など]
https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_heavy_species.08.06.html

[Drift Diffusion Model: μe, De, μen, Den の単位、Use Lookup Tablesが mean electron energy(eV) 軸であること]
https://doc.comsol.com/6.0/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.04.html

[COMSOL Interpolation Function: Spreadsheet形式の列ルール、1〜3引数制限、%コメントなど]
https://doc.comsol.com/5.4/doc/com.comsol.help.comsol/comsol_ref_definitions.12.047.html

[Spreadsheet Data Format: “各行が(引数..., 値...)”、複数関数列OK、%コメント]
https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/comsol_api_fileformats.50.04.html

[COMSOL Blog: EEDFの式、ε(eV), ϕ(eV)、EEDFの eV^-3/2 表記例]
https://www.comsol.com/blogs/electron-energy-distribution-function

[Application例（GUI操作例）: Electron Impact Reaction → Use lookup table → Load from File]
https://doc.comsol.com/5.6/doc/com.comsol.help.models.plasma.positive_column_2d/positive_column_2d.html
```

---

必要なら次に、あなたの **Swarm出力（μ, D, k, EEDF）を「上記CSV形式（特にEEDFの2D long形式）」へ正規化・再サンプリングする具体仕様**（ε軸の共通化、ϕ軸の選び方、外挿方針、EEDFの正規化条件など）を、前回作成した exporter（zip）仕様に合わせて追記できます。

[1]: https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_plasma.09.02.html "https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_plasma.09.02.html"
[2]: https://doc.comsol.com/5.4/doc/com.comsol.help.comsol/comsol_ref_definitions.12.047.html "https://doc.comsol.com/5.4/doc/com.comsol.help.comsol/comsol_ref_definitions.12.047.html"
[3]: https://doc.comsol.com/6.0/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.04.html "https://doc.comsol.com/6.0/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.04.html"
[4]: https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_heavy_species.08.06.html "https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_heavy_species.08.06.html"
[5]: https://www.comsol.com/blogs/electron-energy-distribution-function "https://www.comsol.com/blogs/electron-energy-distribution-function"
[6]: https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_heavy_species.08.07.html "https://doc.comsol.com/5.5/doc/com.comsol.help.plasma/plasma_ug_heavy_species.08.07.html"
[7]: https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/comsol_api_fileformats.50.04.html "https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/comsol_api_fileformats.50.04.html"

以下は **COMSOL 6.3/6.4（Plasma Module、流体モデル）に Swarm 計算結果（μ, D, k, EEDF）を入力して使う**ための、**COMSOL側の具体的な設定手順・入力形式・注意点**を、できるだけ「迷わない粒度」でまとめたものです。
（出典は COMSOL 6.3 ドキュメントを中心に示します。6.4でも基本構造は同等です。）

---

## 1. まず押さえるべき前提：COMSOL流体プラズマで Swarm データを使う “2つの流儀”

COMSOLの流体モデル（Plasma / Drift Diffusion）では、電子輸送係数や電子衝突（電子起源反応）係数を、

* **平均電子エネルギー**（mean electron energy）でパラメータ化する
* または **還元電界 E/N**（Local field approximation）でパラメータ化する

という2ルートがあります。

### A) Local energy approximation（平均電子エネルギー方程式を解く）

* COMSOLが **平均電子エネルギーを未知変数として解き**、その値で **輸送係数・反応係数を参照**します。
* 係数は **平均電子エネルギー（eV）依存のテーブル**を入れるのが自然です。
  （= Swarmで得た μ(ε̄), D(ε̄), k(ε̄) がそのまま使いやすい）

この選択は Drift Diffusion の “Mean electron energy” で **Local energy approximation** を選びます。([COMSOL Documentation][1])

### B) Local field approximation（平均電子エネルギー方程式を解かない）

* 係数は **E/N（還元電界）**で与える想定。
* その場合、COMSOLは「E/N ↔ 平均電子エネルギー」の関係を別途必要とし、**Mean Electron Energy Specification** で **式 or ルックアップテーブル**として与える必要があります。([COMSOL Documentation][1])

> 今回あなたが用意したいのは **平均電子エネルギー依存の μ, D, k と、EEDF(ε, ε̄)** なので、基本は **Local energy approximation** が最も素直です。

---

## 2. COMSOLに渡す物理量と単位・意味（あなたの想定と整合）

COMSOLの Plasma / Drift Diffusion で電子輸送・電子衝突を扱う上で、Swarm由来として最重要なのは次です。

### (1) 電子移動度 μe

* **単位**：m²/(V·s)（COMSOL表示の標準）([COMSOL Documentation][1])
* **役割**：電子ドリフト（電場による流束）を決めます
* **入力軸**：平均電子エネルギー ε̄（eV）で与えるのが自然（Local energy approximation）

### (2) 電子拡散係数 De

* **単位**：m²/s ([COMSOL Documentation][1])
* **役割**：電子拡散流束を決めます
* **入力軸**：平均電子エネルギー ε̄（eV）

> COMSOLでは **μだけ指定して De を Einstein で自動計算**するモードもあります（Maxwellian仮定）。
> ただし SwarmでDeまで出しているなら **Specify all / Use lookup tables** で **μとDe両方**を入れた方が、Maxwellian仮定に縛られにくいです。([COMSOL Documentation][1])

### (3) 電子衝突（電子起源）反応速度係数 kf（例：電離、付着、励起…）

* COMSOLの「Electron Impact Reaction」で、断面積から計算する代わりに **平均電子エネルギー vs 速度係数**として **lookup table**で与えられます。([COMSOL Documentation][1])
* **単位**：典型は **m³/s**（e + 中性 の 2体反応の rate coefficient）([COMSOL Documentation][1])
  （※反応次数で変わるので、GUIでの反応式に従い整合させます。Rate constant 設定では “次数に依存”と明記されています。）([COMSOL Documentation][1])
* **入力軸**：平均電子エネルギー ε̄（eV）

### (4) EEDF f(ε, ε̄)

* COMSOLは「EEDF を選ぶ」機能があり、**Function（2D補間関数）**として与える場合、
  **x 軸 = 電子エネルギー ε（eV）**、**y 軸 = 平均電子エネルギー ε̄（eV）**で与える、と明記されています。([COMSOL Documentation][1])
* ただし重要：
  **EEDFは“断面積データからソース係数を計算する”場合に必要**とされます。([COMSOL Documentation][1])
  つまり **あなたが k(ε̄) を直接与える運用なら、計算上はEEDFが必須ではない**ことが多いです（ただし後述の用途あり）。

---

## 3. COMSOLに入力するデータ形式：基本は “Interpolation Function” と “Lookup Table” の2通り

COMSOL側の取り込みは大きく2方式あります。

### 方式1：Global Definitions の **Interpolation Function**（推奨：汎用・再利用性が高い）

* ファイル（.csv/.txt/.dat）から補間関数を作り、あらゆる設定欄で **式として呼び出す**方式です。
* COMSOLは **.csv も読み込み可**で、区切りはスペース/カンマ/セミコロン/タブ等に対応します。([COMSOL Documentation][2])
* Spreadsheet形式では **「各行が1点」**で、
  1列目=第1引数、2列目=第2引数、…、最後列=関数値、というのが基本です。([COMSOL Documentation][3])
  複数関数を1ファイルにまとめることも可能（引数列の後ろに値列を複数並べる）。([COMSOL Documentation][3])

### 方式2：Plasma/DriftDiffusion/ElectronImpactReaction ノードの **Use lookup tables（内蔵テーブル）**

* そのノードの設定欄で **“Use lookup tables”** を選び、
  **平均電子エネルギー（eV）に対して各係数のテーブル**を直接ロード/入力する方式です。([COMSOL Documentation][1])
* 反応のほうも **Use lookup table** を選ぶと、
  **mean electron energy (eV) vs rate coefficient (m³/s)** の表をファイルロードできる、と明記されています。([COMSOL Documentation][1])

> 実務的には、
>
> * **μ, De は “Use lookup tables”**（Plasma Model/Drift Diffusion Model）で入れる
> * **k は Electron Impact Reaction の “Use lookup table”** で入れる
> * **EEDFは Interpolation Function（2D）** で入れる
>   が、COMSOLの機能構造と整合して扱いやすいです。

---

## 4. COMSOL側での具体的な入力手順（GUI想定）

以下は「Plasma interface を使う（例：ICPリアクタ例題のような構成）」の想定で書きます。
（Drift Diffusion interface単体でも同様で、ノード名が dd になるだけです。）

---

### 4.1 事前：平均電子エネルギーモデルを決める

#### ✅ Local energy approximation（推奨）

* Drift Diffusion interface では Electron Properties の Mean electron energy で
  **Local energy approximation** を選ぶ。([COMSOL Documentation][1])

#### ⚠️ Local field approximation を使う場合

* “Mean Electron Energy Specification” で **E/N と ε̄ の関係**を **式またはlookup table**で与える必要があります。([COMSOL Documentation][1])
  （Swarm側がE/Nベースならこのルートも可能ですが、今回の設計は ε̄ベースが中心。）

---

### 4.2 電子移動度 μe と電子拡散係数 De の入力

#### 入力先（Plasma interface）

* Model Builder → **Plasma（例：plas）→ Plasma Model → Electron Density and Energy**
  ここに “Electron transport properties list” があり、
  **Specify mobility only / Specify all / Use lookup tables** 等を選べます。([COMSOL Documentation][1])

#### A) 最も直球：Use lookup tables で入れる

1. Electron transport properties で **Use lookup tables** を選択
2. **平均電子エネルギー（eV）に対して**、μe, De, Den, μen… のテーブルをロード/入力
   （「transport properties as listed above versus mean electron energy (eV)」と明記）([COMSOL Documentation][1])

**メリット**：Plasma Moduleの想定する形に直結、分かりやすい
**注意**：Den（電子エネルギー拡散係数）や μen（電子エネルギー移動度）も必要になる場合があります（モデル設定で）。無ければ Einstein 関係で代替される場合もありますが、正確性は設定に依存します。([COMSOL Documentation][1])

#### B) Interpolation Function で入れる（柔軟・再利用性が高い）

1. Global Definitions → Functions → **Interpolation**
2. Data source = File、Data format = Spreadsheet
   .csv を読める（.txt/.csv/.dat）([COMSOL Documentation][2])
3. μe 用の 1D 関数 `mu_e(epsbar)` を作る
4. Electron transport properties を **Specify all** にして、
   μe欄に `mu_e(<平均電子エネルギー変数>)` を指定
5. De 欄にも `De_e(<平均電子エネルギー変数>)` のように指定

> 平均電子エネルギー変数は、Plasma interface名（例：plas）に依存して `<name>.<variable_name>` 形式で参照します。([COMSOL Documentation][1])
> COMSOLブログ例では `int2(plas.ebar)` のように **`plas.ebar` を引数に使う**例が示されています。([COMSOL][4])

---

### 4.3 電子起源反応（電離/付着/励起など）の k(ε̄) の入力

#### 入力先

* Model Builder → Plasma → Plasma Model → **Electron Impact Reaction**（反応ごとに追加）

#### 反応係数を “lookup table” で与える

1. Electron Impact Reaction ノードで反応式（例：`e+Ar=>e+e+Ar+` など）を設定
2. “Specify reaction using” を **Use lookup table** にする
3. “Source coefficient data” が出るので、**Rate coefficients** を選択
4. **平均電子エネルギー（eV） vs 反応速度係数（m³/s）** の表を入力/ファイルロード
   （この形式が明記されています）([COMSOL Documentation][1])

#### 重要：Collision type と Energy loss（閾値）

Electron Impact Reaction には Collision type（Elastic/Excitation/Attachment/Ionization）があり、
**Excitation を選ぶと energy loss Δε（しきい値）を入れる**必要があります。([COMSOL Documentation][1])
Ionization も同様に反応に応じたエネルギー損失を設定するのが一般的です。

> これは「kを入れたら終わり」ではなく、**電子エネルギー方程式の収支（電子のエネルギー損失）**に影響します。
> Swarmが “反応係数” を出しても、反応のしきい値エネルギーは別途（反応定義として）入れる必要が出やすい点です。

#### “Rate constant（式）”で入れる場合

* “Specify reaction using = Rate constant” にすると、Forward rate constant に式を入れます。
* このとき **単位は反応次数に依存**します。([COMSOL Documentation][1])
  （lookup table の方が単位事故が少ないです。）

---

### 4.4 EEDF（ε, ε̄）の入力と、COMSOL内での扱い

#### COMSOLが求めるEEDF “Function” の形式

Plasma Moduleでは EEDF を **Function（2D補間関数）**として与える場合、

* **x-data：電子エネルギー ε（eV）**
* **y-data：平均電子エネルギー ε̄（eV）**

であることが明記されています。([COMSOL Documentation][1])

また、**断面積データからソース係数を計算する場合はEEDF選択が必要**と説明されています。([COMSOL Documentation][1])

#### 実際の入力手順（推奨）

1. Global Definitions → Functions → **Interpolation** を追加
2. Data source = File、Data format = Spreadsheet
3. **2引数関数**（Interpolation は 1〜3 引数サポート）([COMSOL Documentation][2])
4. CSVの列を「第1引数=ε, 第2引数=ε̄, 値=f」に割り当てる
   （Spreadsheet形式は “各行が1点、列が引数と関数値”）([COMSOL Documentation][3])
5. Plasma Model → Electron Energy Distribution Function で **Function** を選択し、この2D関数を指定する

#### 重要：EEDFを“平均エネルギーごとにファイル分割”する意味は薄い

COMSOLが求めるのは **(ε, ε̄) → f** の2D関数なので、
**E/Nや平均エネルギーごとに別ファイル**に分けるより、
**1つの2Dテーブル（散布点形式でもOK）として持たせた方が、COMSOL側で自然に補間できます。**（この前提が docs に書かれている）([COMSOL Documentation][1])

---

## 5. 推奨する “CSVの具体形” （COMSOLが読みやすい形）

### 5.1 μ と D（1D：平均電子エネルギー依存）

例（CSV、ヘッダは `%` コメントにするとCOMSOLで安全）：([COMSOL Documentation][2])

```text
% ebar_eV, mu_m2_Vs, De_m2_s
0.1, 0.35, 0.012
0.2, 0.40, 0.015
...
```

* 1ファイルにまとめてもOK（引数列の後ろに複数の値列を置ける）([COMSOL Documentation][3])
* あるいは `ebar, mu` の2列ファイルを別々でもOK

### 5.2 反応係数 k（1D：平均電子エネルギー依存）

Electron Impact Reaction の Use lookup table が想定する形は：

* **mean electron energy (eV) vs rate coefficient (m³/s)**([COMSOL Documentation][1])

例：

```text
% ebar_eV, k_m3_s
0.1, 1.2e-16
0.2, 3.4e-16
...
```

### 5.3 EEDF（2D：ε と ε̄）

Spreadsheet形式は「各行=1点、列=引数…値」でOKです。([COMSOL Documentation][3])
例：

```text
% eps_eV, ebar_eV, f
0.0,  1.0,  ...
0.1,  1.0,  ...
...
0.0,  2.0,  ...
0.1,  2.0,  ...
...
```

---

## 6. COMSOL側で “補間が正しく入ったか” を確認する方法（必須）

### 6.1 Interpolation Function のテスト

* Global Definitions → Functions → Interpolation を右クリック → **Plot** で
  1D/2Dプロットして形がおかしくないか確認
* 外挿（extrapolation）は “Constant” にするのが事故りにくい（範囲外で暴れない）。
  Interpolationヘルプでは外挿に Constant / Specific value 等があることが説明されています。([COMSOL Documentation][2])

### 6.2 反応係数の確認

* Results → Derived Values → Global Evaluation などで
  `k_fun(plas.ebar)` を評価し、期待するオーダーかチェック

---

## 7. よくある落とし穴（Swarm→COMSOL連携で高頻度）

### 7.1 “平均電子エネルギーの単位”混乱

* COMSOL内部では mean electron energy を **V として扱う箇所**があり（初期値 e0 が V 等）、([COMSOL Documentation][1])
  ブログでも「electron energy の単位は V だが eV として扱われる」と注意があります。([COMSOL][4])
* 実務では **テーブルは eV で作り、関数引数も eV として設定**し、COMSOLの式欄でも eV で整合させるのが安全です（COMSOLは eV 単位も扱えます）。

### 7.2 μだけ入れてDeがEinstein計算になってしまう

* “Specify mobility only” を選ぶと、Deなどが Maxwellian仮定で Einstein 関係から計算されます。([COMSOL Documentation][1])
* SwarmのDeを使いたいなら **Specify all** か **Use lookup tables** に寄せる。

### 7.3 k(ε̄) を入れても “energy loss（しきい値）” を入れ忘れる

* Excitationなどは Δε（しきい値）を入れる設計です。([COMSOL Documentation][1])
* ここが抜けると **平均電子エネルギーの収支が壊れる**ことがあります。

### 7.4 EEDFを入れたのに反応係数が変わらない

* そもそも **k を lookup table で与えているなら EEDFは反応係数計算に使われません**。
* EEDFが必要なのは主に「断面積からkを計算する」ルート（cross section data）です。([COMSOL Documentation][1])
  → ただし **EEDFを後処理・可視化**に使う（COMSOL上で参照・比較する）目的は十分あります。

---

## 8. まとめ：あなたの用途に対する “COMSOL側の最適構成” 提案

あなたの条件（断面積入力なし、Swarmで μ, D, k, EEDF を作る）だと、COMSOL側は下記が最も安定です。

1. **Mean electron energy = Local energy approximation**（ε̄で係数参照）([COMSOL Documentation][1])
2. 電子輸送（μ,De）は

   * **Use lookup tables**（Plasma Model / Drift Diffusion Model）([COMSOL Documentation][1])
     もしくは
   * Interpolation Function + Specify all（式で呼ぶ）
3. 電子起源反応は

   * Electron Impact Reaction → **Use lookup table** → **Rate coefficients**（ε̄[eV] vs k[m³/s]）([COMSOL Documentation][1])
   * 反応ごとに Collision type と Δε を適切に設定([COMSOL Documentation][1])
4. EEDFは

   * Global Interpolation Function（2D）として `f(ε, ε̄)` を作り、
   * 必要なら Plasma Model の EEDF を Function にする（cross sectionを使う時に必須）([COMSOL Documentation][1])
   * ただし **kをテーブル入力する運用では必須ではない**（用途は主に可視化・整合チェック）

---

必要なら、次のステップとして「あなたの exporter が出す CSV（μ/D/k/EEDF）を前提に、COMSOL側の Model Builder のノード構造（Plasma interface / dd interface）別に、**“このCSVをこのノードにこう割り当てる”**を **チェックリスト形式**（GUI操作順）で落とし込みます。
また、**mph（Python）で同じ設定を自動で埋め込む場合**に、どのノード・どのプロパティを触るべきかも併記できます。

[1]: https://doc.comsol.com/6.3/doc/com.comsol.help.plasma/PlasmaModuleUsersGuide.pdf "https://doc.comsol.com/6.3/doc/com.comsol.help.plasma/PlasmaModuleUsersGuide.pdf"
[2]: https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_ref_definitions.21.063.html "https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_ref_definitions.21.063.html"
[3]: https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_api_fileformats.54.04.html "https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_api_fileformats.54.04.html"
[4]: https://www.comsol.com/blogs/the-boltzmann-equation-two-term-approximation-interface "The Boltzmann Equation, Two-Term Approximation Interface | COMSOL Blog"
