# 固定ガス条件での電子衝突経路

Swarmは、指定した組成・温度・密度・DC E/Nで電子衝突の速度係数を計算します。
反応経路出力では、出典付きカタログの反応式とその係数を結び付け、直接の生成・消費経路を
順位付けします。重粒子反応、壁、組成の時間発展、電子密度は利用先のプラズマモデルが扱います。
[LoKIの電子・化学計算の役割分担](https://nprime.tecnico.ulisboa.pt/loki/tools.html)

## 実行

リポジトリのルートから、実断面積を使うAr 90%・N₂ 10%の例を実行できます。

```powershell
py -3 -m swarm_workflow.cli run examples/reaction_pathways_ar_n2.yaml --plan --json
py -3 -m swarm_workflow.cli run examples/reaction_pathways_ar_n2.yaml
```

使う断面積ファイルとガス条件は既存の `cross_sections` と `conditions` で指定します。
`output.export` でsolverとカタログを指定します。生成種を絞る場合は `focus_species` を、
微小な経路を順位から外す場合は `minimum_activity_s_inv` を追加します。
カタログは断面積データに付随する定義で、同じデータを使う別の組成・E/Nでも使えます。

```yaml
output:
  export:
    solver: two_term
    quantities: [reaction_pathways]
    pathways:
      catalog: cross_sections/ar_n2.reactions.yaml
      focus_species: Ar+
      minimum_activity_s_inv: 10.0
```

保存済みのproduct manifestまたはworkflow SQLiteから、solverを再実行せずに
着目種と活動下限を変えられます。

```powershell
py -3 -m swarm_workflow.cli export-reaction-rates outputs/examples/ar_n2_pathways_manifest.json --solver two_term --pathway-catalog examples/cross_sections/ar_n2.reactions.yaml --focus-species Ar+ --output outputs/argon_ion_paths
```

## 出力の見方

| ファイル | 読む内容 |
| --- | --- |
| `*_reaction_rates.csv` | 元の反応別係数 `k`、対象種のモル分率、数値状態、利用できるMCの区間推定 |
| `*_reaction_pathways.csv` | 反応式、化学量論のJSON、カテゴリ、分解能、出典、反応頻度、順位、判定理由 |
| `*_manifest.json` | 組成・密度・E/N、元の計算の由来、カタログと断面積のSHA-256、規則、全判定件数 |

経路表の1行は、係数表の同じ `series_id, solver, E_over_N_Td, species, process,
process_type, threshold_eV` に対応します。反応式を読むだけでなく、`reactants` と
`products` のJSONから正味化学量論を取り出せます。周波数は
`event_frequency_s_inv = n_gas × x_target × k` です。着目種がある場合は、その種の
正味生成数を掛けた符号付きの `focus_source_per_electron_s_inv` も出力します。
体積あたりの反応源には利用先でその状態の電子密度を掛けます。

順位の範囲は、同じ計算条件で式が同定され、数値根拠が確認できた行です。通常は
カテゴリ内で頻度を比較し、着目種がある場合は生成と消費を分けて正味反応源の大きさを
比較します。大小は順位とスコアで示し、同点は同順位です。総括断面積と内訳は
重なる可能性があるため、各チャネルの頻度を個別に保持します。

| 判定 | 規則 |
| --- | --- |
| `ranked` | 収束した決定論計算の係数、または独立MC反復の95%区間の正の下端を採用。スコアが正で指定下限以上 |
| `below_minimum_activity` / `zero_activity` | 根拠付きスコアが下限未満、またはゼロ。係数と根拠は保持 |
| `does_not_change_focus` | 着目種の正味化学量論がゼロ |
| `unmapped_channel` | 対応するカタログ反応が未登録。係数と電子1個あたりの頻度を保持 |
| `scattering_only` / `unsupported_collision_family` | 散乱、または化学反応頻度を定義できない衝突型。係数表で照合可能 |
| `unranked_numerical_status` | 収束が確認できない決定論点、後処理で変わった点など。係数表に元の状態を記録 |
| `mc_censored` / `mc_requires_independent_replicas` | 全反復ゼロ観測、または独立反復不足 |
| `mc_interval_unavailable` / `mc_interval_includes_zero` | 区間の根拠が不足、または区間がゼロを含む |

MCの下端を使った順位は、反応ごとの保守的な選別指標です。反応間の差の有意性には、
反復を対応させた差の検定が別途必要です。順位が付いた行がゼロの場合は
`no_ranked_channels` と記録し、判定表を保存します。

## カタログの追加

[Ar/N₂カタログ](../examples/cross_sections/ar_n2.reactions.yaml)を例に、次を記述します。
初めて使う断面積では、先に `quantities: [reaction_rates]` で対象種・衝突型・
`process` 名・閾値の一覧を取得できます。

- `datasets` に、使用する断面積ファイルすべてのSHA-256と出典を登録します。
  照合対象はファイルの内容です。ファイル集合の不一致はエラーになります。
- `species` に、種ID・元素数・電荷を明示します。電子のID `e` と電荷−1は共通です。
  励起状態や集約状態は、組成が同じでも独立の種IDにします。
- `reactions` に、対象種、断面積の正確な `process` 名、閾値、化学カテゴリ、
  生成物の整数個数、分解能 `resolved` または `lumped`、出典内の位置を記載します。
  閾値未宣言の付着チャネルには `null` を使います。式の左辺は電子1個と対象種1個です。

共通処理は元素・電荷、種ID、重複、データ由来を確認します。カテゴリは計算側の衝突型と
生成電子・重粒子の個数を規定します。励起と解離はともにSwarmでは `excitation` 型を
使えますが、解離は複数の重粒子を必要とします。電離、解離性電離、付着、解離性付着、
超弾性もそれぞれ確認します。個々の生成種・状態・分岐は、出典を調べたカタログの責務です。
正の割合を持つ対象種の登録反応が係数に見つからない場合は、ラベルや閾値のずれをエラーにします。

## 実データでの確認

300 K、13.3 Pa、Ar 90%・N₂ 10%で、20、50、100 Tdのtwo-term計算を行いました。
各点49チャネルのうち18経路に順位を付け、29行を未登録、2行を散乱として保持しました。
励起の首位は3点とも `e + N2 → e + N2(v=1)` でした。50 Tdでは
`1.7332e6 s⁻¹`、Arの基底状態からの直接電離は `6.1999e3 s⁻¹` でした。
純Arとの組成スイープと、保存結果からのAr⁺生成経路の再選別も確認しました。

ArはCOMSOL例題の基底状態の励起・電離、N₂は既存のBiagi Magboltz 8.97断面積の
振動励起15件と電子励起の集約状態1件を使いました。N₂の電離は生成枝を確定できる情報が
得られるまで未登録としています。20 Tdでは閾値格子を細分化したtwo-termが収束しないため、
この例の格子設定には細分化を無効化する指定を明記しています。
50、100 Tdでは両方の設定が収束し、N₂の振動励起首位とArの直接電離の係数差は最大0.0082%でした。
[COMSOLのAr反応表](https://doc.comsol.com/6.4/doc/com.comsol.help.models.plasma.boltzmann_argon/boltzmann_argon.html)
[LXCat Biagi](https://lxcat.net/Biagi)
