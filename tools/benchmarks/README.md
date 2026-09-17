# Benchmark tools

This package owns benchmark orchestration, EEDF comparison, and adapters for
external reference solvers. Each script can be run directly from the repository
root, for example:

```powershell
py -3 tools\benchmarks\run_product_benchmarks.py
```

Reference parsers and command runners are isolated in `references/`; they do
not participate in product solver execution.
