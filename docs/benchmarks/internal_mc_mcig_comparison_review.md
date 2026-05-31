# Internal MC と MCIG EEDF 比較レビュー

このメモは、Ar 単ガスの MCIG GUI actual と本コード Internal MC の
EEDF 比較について、最近の改善内容、観測された効果、物理精度への影響を
整理したものです。MCIG は有用な Monte Carlo reference ですが、ここでは
「絶対的な正解」や fitting target とは扱いません。

## 対象

- gas: Ar
- field: B=0, DC
- E/N: 100, 200, 400, 600 Td
- 本コード MC: `fixed_particle_single_daughter`
- warmup: 600 collisions
- production: 1200 collisions
- particles: 768
- EEDF estimator: time-sampled null-clock histogram

最新比較出力:

- curves:
  `outputs/benchmarks/internal_mc_fixed_warmup_hp_100_200_400_600/comparison/ar_mcig_two_term_internal_mc_fixed_warmup_hp_curves.csv`
- graph:
  `outputs/benchmarks/latest_mcig_internal_mc_eedf/ar_mcig_vs_internal_mc_latest_eedf.png`
- metrics:
  `outputs/benchmarks/latest_mcig_internal_mc_eedf/ar_mcig_vs_internal_mc_latest_metrics.csv`

## 再現手順

この比較は、外部 MCIG を本コード内で solver として実行するものではなく、
MCIG GUI などで得た EEDF 出力を reference として読み込み、本コードの
Internal MC 出力と同じ `F(E) [1/eV]` に正規化して比較します。

### 1. 依存関係

通常の開発環境を用意します。

```powershell
py -3 -m pip install -e ".[dev,plot]"
```

plot extra がない場合でも CSV は生成できますが、EEDF 図を再生成するには
`matplotlib` が必要です。

### 2. 本コード Internal MC の再実行

今回使った high-precision warmup MC config は以下です。

```powershell
outputs\benchmarks\internal_mc_fixed_warmup_hp_100_200_400_600\mc_internal_fixed_warmup_hp.yaml
```

重要な設定:

```yaml
run:
  solvers:
    - id: monte_carlo
  e_over_n_Td: [100, 200, 400, 600]
physics:
  angular_scattering:
    model: isotropic
  ionization:
    energy_sharing: equal
solvers:
  monte_carlo:
    backend: internal
    angular_scattering: same_as_physics
    population_model: fixed_particle_single_daughter
    particles: 768
    warmup_collisions: 600
    max_collisions: 1200
    seed: 20260601
```

Internal MC audit を再実行するには:

```powershell
py -3 tools\benchmark_internal_mc_audit.py `
  --config outputs\benchmarks\internal_mc_fixed_warmup_hp_100_200_400_600\mc_internal_fixed_warmup_hp.yaml `
  --output-dir outputs\benchmarks\internal_mc_fixed_warmup_hp_check
```

出力:

- `outputs\benchmarks\internal_mc_fixed_warmup_hp_check\internal_mc_audit_summary.csv`
- `outputs\benchmarks\internal_mc_fixed_warmup_hp_check\internal_mc_tail_uncertainty.csv`
- `outputs\benchmarks\internal_mc_fixed_warmup_hp_check\internal_mc_audit_report.md`

この audit は product runner を増やさず、開発時に MC energy balance と
tail uncertainty を確認するための薄い wrapper です。

### 3. MCIG reference の準備

MCIG は benchmark target であり、本コードが MCIG を正解として
fitting するわけではありません。MCIG GUI などで Ar 100/200/400/600 Td
の EEDF を保存し、以下のいずれかの形で比較用 CSV に変換します。

推奨 canonical format:

```csv
case_id,E_over_N_Td,energy_eV,eedf_eV_inv,mean_energy_eV
ar_100,100,0.01,0.00123,6.42
...
```

または EEPF 入力:

```csv
case_id,E_over_N_Td,energy_eV,eepf_eV_m32,mean_energy_eV
ar_100,100,0.01,0.0123,6.42
...
```

EEPF の場合は比較 ingest 側で `F(E)=EEPF*sqrt(E)` に変換し、
`integral F(E)dE = 1` に正規化します。MCIG 側の angular scattering model、
growth model、uncertainty が分かる場合は metadata として記録してください。
分からない場合、その比較は `unknown reference setting` を含む評価として
扱います。

今回の最新比較 CSV では、MCIG reference は `MCIG GUI actual` という
source 名で入っています。

### 4. 比較曲線 CSV の作成

今回使った curves CSV は以下です。

```powershell
outputs\benchmarks\internal_mc_fixed_warmup_hp_100_200_400_600\comparison\ar_mcig_two_term_internal_mc_fixed_warmup_hp_curves.csv
```

列:

- `E_over_N_Td`
- `source`
- `energy_eV`
- `eedf_eV_inv`
- `counts`
- `bin_width_eV`
- `effective_sample_count`
- `relative_standard_error`

`source` は少なくとも以下を含みます。

- `MCIG GUI actual`
- `internal MC fixed warmup HP (768, warmup 600 + production 1200)`
- `two_term native SG hold`

### 5. EEDF グラフの再生成

以下の Python snippet で、最新比較図と metrics 抽出を再生成できます。

```powershell
@'
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

src = Path("outputs/benchmarks/internal_mc_fixed_warmup_hp_100_200_400_600/comparison/ar_mcig_two_term_internal_mc_fixed_warmup_hp_curves.csv")
metrics_src = Path("outputs/benchmarks/internal_mc_fixed_warmup_hp_100_200_400_600/comparison/ar_mcig_two_term_internal_mc_fixed_warmup_hp_metrics.csv")
out_dir = Path("outputs/benchmarks/latest_mcig_internal_mc_eedf")
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(src)
metrics = pd.read_csv(metrics_src)

mcig = "MCIG GUI actual"
internal = "internal MC fixed warmup HP (768, warmup 600 + production 1200)"
two_term = "two_term native SG hold"
styles = {
    mcig: ("MCIG GUI actual", "#111111", "-", 2.2),
    internal: ("Internal MC fixed warmup HP", "#d62728", "-", 1.7),
    two_term: ("two-term native SG", "#6f6f6f", "--", 1.25),
}

fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6))
for ax, en in zip(axes.ravel(), sorted(df["E_over_N_Td"].unique())):
    sub_en = df[df["E_over_N_Td"] == en]
    for source in [mcig, internal, two_term]:
        sub = sub_en[sub_en["source"] == source].sort_values("energy_eV")
        x = sub["energy_eV"].to_numpy(float)
        y = sub["eedf_eV_inv"].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        label, color, ls, lw = styles[source]
        ax.plot(x[mask], y[mask], label=label, color=color, linestyle=ls, linewidth=lw)
        if source == internal and "relative_standard_error" in sub:
            rse = sub.loc[mask, "relative_standard_error"].to_numpy(float)
            finite = np.isfinite(rse) & (rse > 0) & (rse < 1)
            ax.fill_between(
                x[mask][finite],
                np.maximum(y[mask][finite] * (1 - rse[finite]), 1e-40),
                y[mask][finite] * (1 + rse[finite]),
                color=color,
                alpha=0.16,
                linewidth=0,
            )
    row = metrics[
        (metrics["E_over_N_Td"] == en)
        & (metrics["reference"] == mcig)
        & (metrics["candidate"] == internal)
    ].iloc[0]
    ax.set_title(
        f"Ar {en:g} Td, L1={row['eedf_relative_l1']:.3f}, "
        f"mean diff={100*row['mean_energy_relative_difference']:.1f}%"
    )
    ax.set_yscale("log")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("EEDF F(E) (1/eV)")
    ax.grid(True, which="both", alpha=0.22)

handles, labels = axes.ravel()[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(out_dir / "ar_mcig_vs_internal_mc_latest_eedf.png", dpi=180)
fig.savefig(out_dir / "ar_mcig_vs_internal_mc_latest_eedf.svg")

metrics[
    (metrics["reference"] == mcig)
    & (metrics["candidate"] == internal)
].to_csv(out_dir / "ar_mcig_vs_internal_mc_latest_metrics.csv", index=False)
'@ | py -3 -
```

再生成されるファイル:

- `outputs/benchmarks/latest_mcig_internal_mc_eedf/ar_mcig_vs_internal_mc_latest_eedf.png`
- `outputs/benchmarks/latest_mcig_internal_mc_eedf/ar_mcig_vs_internal_mc_latest_eedf.svg`
- `outputs/benchmarks/latest_mcig_internal_mc_eedf/ar_mcig_vs_internal_mc_latest_metrics.csv`

### 6. 評価に使うファイル

第三者が結果を確認する場合は、最低限この 4 ファイルを見ます。

1. MC config:
   `outputs/benchmarks/internal_mc_fixed_warmup_hp_100_200_400_600/mc_internal_fixed_warmup_hp.yaml`
2. curves CSV:
   `outputs/benchmarks/internal_mc_fixed_warmup_hp_100_200_400_600/comparison/ar_mcig_two_term_internal_mc_fixed_warmup_hp_curves.csv`
3. metrics CSV:
   `outputs/benchmarks/latest_mcig_internal_mc_eedf/ar_mcig_vs_internal_mc_latest_metrics.csv`
4. audit summary:
   `outputs/benchmarks/internal_mc_fixed_warmup_hp_check/internal_mc_audit_summary.csv`

## 評価方法

この benchmark は「MCIG と完全一致するか」ではなく、以下を順番に見る
triage です。

### 1. EEDF convention を確認

全 source の EEDF は以下に統一します。

```text
F(E), unit 1/eV
integral F(E)dE = 1
```

MCIG が EEPF を出している場合は、必ず EEDF に変換してから比較します。
EEDF/EEPF の取り違えは tail と低エネルギー側の両方を大きく歪めます。

### 2. Internal MC audit status を確認

まず `internal_mc_audit_summary.csv` で以下を見る。

- `mc_energy_balance_status`
- `mc_tail_comparison_status`
- `mc_tail_uncertainty_status`
- `mc_min_tail_bin_count`
- `mc_max_resolved_energy_eV`
- `mc_null_collision_acceptance_fraction`
- `mc_max_collision_to_trial_ratio`

`mc_energy_balance_status != ok` の場合、MCIG との差を物理差として議論する
前に、Internal MC の数値 audit を確認します。

`mc_tail_comparison_status != ok` の場合、高エネルギー tail の比較は
定量評価ではなく診断扱いにします。

### 3. EEDF relative L1 を見る

`ar_mcig_vs_internal_mc_latest_metrics.csv` の `eedf_relative_l1` を見る。

目安:

- `< 0.03`: かなり近い
- `0.03 - 0.06`: 実用的には近いが設定差の確認が必要
- `0.06 - 0.12`: 高 E/N では model difference または tail estimator 差を疑う
- `> 0.12`: convention、angular model、growth model、または MC 統計設計を再確認

この閾値は acceptance test ではなく、benchmark triage の目安です。

### 4. Mean energy difference を見る

`mean_energy_relative_difference` は EEDF shape よりも頑健な指標です。

目安:

- `< 2%`: よく一致
- `2 - 5%`: 設定差を確認
- `5 - 10%`: angular / ionization / growth model 差を疑う
- `> 10%`: reference setting または MC estimator を要再確認

### 5. Tail を見る

高エネルギー tail は、MCIG 側と本コード側で sampling uncertainty が異なる
ため、EEDF L1 より慎重に解釈します。

特に本コードの EEDF は 80 eV 以降で bin が粗くなるため、600 Td の
80-100 eV 付近の滑らかさは物理的なノイズ低下ではなく bin averaging の
効果です。tail の細かい構造を評価する場合は、development-only fine-bin
estimator で再集計してください。

### 6. MCIG との差を分類

差が大きい場合は、以下の順に疑います。

1. EEDF/EEPF convention mismatch
2. MCIG angular scattering model が不明または不一致
3. MCIG growth model / ionization secondary treatment の違い
4. Internal MC tail count / ESS 不足
5. null-collision majorant or high-energy extrapolation issue
6. cross-section input mismatch
7. 本コードの実装 bug

MCIG は benchmark target であり、差があること自体を即座に本コード bug と
見なさないでください。

### 7. 再現時の期待値

同じ config、同じ seed、同じ MCIG reference CSV を使えば、最新比較では
おおむね以下になります。

| E/N (Td) | EEDF relative L1 | mean difference |
|---:|---:|---:|
| 100 | 0.031 | 1.9% |
| 200 | 0.045 | 3.9% |
| 400 | 0.085 | 7.4% |
| 600 | 0.108 | 10.6% |

Internal MC audit は全点で:

```text
mc_energy_balance_status = ok
mc_tail_comparison_status = ok
```

になることを期待します。

## 改善内容

### 1. Ionization energy sharing を schema 設定に接続

以前の Internal MC は、電離後の追跡粒子 energy を実質的に
`E - I` とする loss-only 的な挙動に寄っていました。これは単一粒子を
追跡する fixed-particle MC では高エネルギー tail を過大評価しやすい
モデルです。

現在は `physics.ionization.energy_sharing` に従います。

- `equal`: 追跡電子 energy は `(E - I) / 2`
- `primary_secondary`: primary または secondary を 1/2 確率で追跡
- `loss_only`: 診断用として明示された場合のみ

製品 metadata には以下を出します。

- `mc_population_model=fixed_particle_single_daughter`
- `ionization_branching_model=single_daughter_sampling`
- `secondary_electron_tracking=false`

効果:

- 電離後に過大な primary energy を持ち続ける bias が減る。
- 高 E/N で Internal MC の tail が下がり、BOLSIG+/MCIG reference に
  近づく。
- ただし full branching ではないため、実粒子群の増殖そのものはまだ
  表現していない。

### 2. Warmup collisions を導入

MC の初期速度分布からの過渡を EEDF estimator に混ぜると、特に高 E/N で
tail と mean energy が初期条件に影響されます。

現在は `warmup_collisions` を product knob として残し、warmup 中の
flight は EEDF/transport sampling に入れません。production sampling は
warmup 後に開始します。

効果:

- 100 Td では two-term との EEDF L1 が大きく改善した。
- MCIG との比較でも、低エネルギー側と mean energy のずれが減った。
- warmup は物理モデルを変えるものではなく、定常 swarm 分布を測るための
  estimator 改善です。

### 3. Null-collision majorant を fail-fast 化

過去の MC では、もし実行中に `nu(E) > nu_trial` が起きても確率を
`min(prob, 1)` 的に丸めると、tail 側で見かけ上の衝突不足や数値加熱を
隠す危険がありました。

現在は trial collision frequency を energy limit まで評価し、実行中に
majorant violation があれば fail-fast します。

効果:

- 高エネルギー tail の破綻を黙殺しない。
- tail の過大評価が数値 majorant 不足由来かどうか切り分けやすい。

### 4. EEDF binwise uncertainty を出力

Internal MC の EEDF CSV には以下の列を出します。

- `sample_count`
- `effective_sample_count`
- `relative_standard_error`

また summary metadata には compact に以下を出します。

- `mc_tail_uncertainty_status`
- `mc_min_tail_bin_count`
- `mc_max_resolved_energy_eV`
- `mc_tail_comparison_status`

効果:

- tail の差を solver mismatch と即断せず、統計的に弱い bin を識別できる。
- MCIG に uncertainty がない場合でも、本コード MC 側の tail 信頼度を
  可視化できる。

### 5. Full branching 実験を product API から外した

一時的に検討した weighted full branching は、低エネルギー側の EEDF を
悪化させ、product としては解釈が複雑でした。現在は product schema から
外し、roadmap/diagnostic 扱いに戻しています。

効果:

- 製品 MC の意味が明確になった。
- MCIG matching correction のような ad hoc 経路を避けた。
- 本体コードは fixed-particle single-daughter に集中できる。

## 最新比較結果

MCIG GUI actual と Internal MC fixed warmup HP の比較:

| E/N (Td) | EEDF relative L1 | MCIG mean (eV) | Internal MC mean (eV) | mean difference |
|---:|---:|---:|---:|---:|
| 100 | 0.0314 | 6.4276 | 6.5468 | 1.86% |
| 200 | 0.0453 | 7.7577 | 8.0622 | 3.92% |
| 400 | 0.0854 | 10.1091 | 10.8569 | 7.40% |
| 600 | 0.1080 | 12.4614 | 13.7822 | 10.60% |

Internal MC audit status:

| E/N (Td) | energy balance | tail comparison | warmup | production |
|---:|---|---|---:|---:|
| 100 | ok | ok | 600 | 1200 |
| 200 | ok | ok | 600 | 1200 |
| 400 | ok | ok | 600 | 1200 |
| 600 | ok | ok | 600 | 1200 |

## 物理的な精度への影響

### 改善された点

- 電離 energy sharing の過大 tail bias が減った。
- 初期過渡を捨てることで、定常 EEDF に近い推定になった。
- majorant violation を fail-fast にしたため、数値加熱を隠しにくくなった。
- tail uncertainty を出すことで、統計不足と物理差を分けやすくなった。

### まだ同一視できない点

MCIG と本コード Internal MC は、完全に同じ物理モデルとは限りません。
以下が一致していない、または reference 側から読めない場合、差は実装 bug
ではなく model difference の可能性があります。

- angular scattering model
- growth / nonconservative swarm model
- ionization secondary treatment
- EEDF sampling clock
- rate convolution convention
- MCIG 側の binwise uncertainty

特に 400, 600 Td では Internal MC の mean energy が MCIG より高く、
EEDF L1 も増えます。この差は、単純な統計ノイズだけではなく、angular
model、growth model、ionization treatment、または reference setting の
違いが効いている可能性があります。

## Bin 幅の注意

現在の Internal MC EEDF histogram は、80 eV 以降で energy bin が粗く
なります。

- 0-80 eV: 0.5 eV bin
- 80 eV 以降: 10 eV bin

そのため 600 Td の 80-100 eV 付近では、ノイズが急に小さくなったように
見えます。これは衝突計算の精度が急に上がったためではなく、広い bin に
よる平均化です。

影響:

- 粒子軌道計算、衝突判定、mean energy には直接影響しない。
- EEDF 曲線、EEDF L1、tail probability、tail shape 評価には影響する。

高エネルギー tail の細かい構造を評価する場合は、development benchmark
側で fine-bin estimator を使って再集計する必要があります。

## 評価

現状の Internal MC は、100-200 Td では MCIG とかなり近く、製品の
fixed-particle MC として妥当な範囲に入っています。400-600 Td では差が
大きくなりますが、energy balance と tail comparison status は ok であり、
少なくとも明白な数値破綻は検出されていません。

したがって現時点の結論は以下です。

- 低から中 E/N: product MC として実用的な妥当性あり。
- 高 E/N: tail と mean energy の差は残る。MCIG と合わせ込むべきではなく、
  angular model、growth model、reference uncertainty を切り分けるべき。
- 本体 product API は fixed-particle single-daughter + warmup に絞った
  ままが安全。

## 次に検討すべきこと

1. 80 eV 以降も fine-bin で出す development-only EEDF estimator
2. seed ensemble による MC confidence interval
3. MCIG 出力から angular model と uncertainty をより明示的に ingest
4. 同じ angular model での MCIG / Internal MC 比較
5. full branching は product default ではなく、別モデルとして方程式、
   observable、normalization を固定してから再検討
