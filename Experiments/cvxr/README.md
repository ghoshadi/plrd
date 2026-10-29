# Reproducing the numbers in the `optimizer = "cvxr"` PR

Run every script from the repository root, for example

```sh
Rscript Experiments/cvxr/table2.R
```

The scripts source the package's own `R/` files, so they test this checkout as it stands,
without installing it. `common.R` holds the shared code: `plrd.optim()`'s grid and bias
operator, rebuilt exactly as `plrd.optim()` builds them, and the scoring used throughout.
Every solution is judged only by the weights it returns, against the problem
`plrd.optim()` itself sets up:

- the constraints, checked at the bin centers;
- the true worst-case bias $B\lVert Av\rVert_1$, with $v$ the weight per bin;
- the fold objective $\hat\sigma^2\sum\gamma_i^2 + \text{bias}^2$.

No solver's self-report is trusted.

| script | reproduces | time |
|---|---|---|
| `table2.R` | Table 2 through `plrd()` under the three optimizer paths; Issue 1's table; Issue 2's numbers on Matsudaira and Jacob & Lefgren; the fold-objective comparison; the timings | ~1 min |
| `audit_paths.R` | Issue 2's table: the three `plrd.optim()` paths on the `test_plrd.R` design, seeds 1–3, plus a centered control | ~30 s |
| `spline_df_sweep.R` | Issue 1: the error in the reported worst-case bias as `spline.df` varies | ~15 s |
| `fold_noise.R` | the optimizer's effect on senate against the seed-to-seed variation of $\hat\tau$ | ~10 s |
| `plot_settings.R` | the seven `test_plot.R` designs: every fold re-solved with PIQP and ECOS | ~15 s |
| `tolerance.R` | why `R/cvxr.R` sets `tol_gap_rel = 1e-7`: the folds that report `optimal_inaccurate` under CLARABEL's defaults, with solver logs | ~20 s |
| `solver_comparison.R` | CLARABEL vs PIQP on 120 fold problems, scaled and unscaled, with timings | ~5 min |

Packages beyond the package's own Imports: `piqp` and `ECOSolveR` for `plot_settings.R` and
`tolerance.R`, and `piqp` and `bench` for `solver_comparison.R`.

**Data.** `table2.R` and `fold_noise.R` read the CSV files already in
`Experiments/Datasets/`, prepared as in `Experiments/Datasets/Results_for_real_datasets.R`.
The other scripts use simulated data. This directory adds no data.

**Output.** Only `solver_comparison.R` writes files, to `$CVXR_EXPT_OUT` if set and
otherwise to a temporary directory. Nothing is written inside the repository.

Timings depend on the machine, so compare ratios rather than seconds. The numbers in the
PR were produced with R 4.6.1, CVXR 1.9.3 and clarabel 0.11.3 (piqp 0.6.4.1, ECOSolveR
0.6.2, bench 1.1.4) on macOS.
