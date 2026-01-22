# Integration with Swarm repository

本ツールは「Swarm計算の出力ファイル」を入力として動作するため、Swarm本体を壊さずに後処理として組み込めます。

## 推奨フロー
1. Swarm計算を通常通り実行
2. Swarmの出力ディレクトリ（run dir）に、transport / rates / eedf のCSV（またはEEDFの分割CSV群）が残る
3. 本ツールを実行して COMSOL用CSVを生成

## Swarm側への最小変更（推奨）
- Swarmの実行YAMLに `comsol_export.enabled: true` を追加
- Swarmの最後に1行追加してフックを呼ぶ

例（擬似コード）:

```python
from swarm_comsol_export.hook import maybe_export_comsol

# swarm_main(...)
# ... swarm calc ...
maybe_export_comsol(cfg, run_dir=output_dir)
```

cfg は dict でも YAML パスでもOK（hook参照）。
