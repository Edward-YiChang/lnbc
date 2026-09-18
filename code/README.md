# Reproducing the manuscript results

This directory contains the R and shell code needed to reproduce the simulation
and malaria analyses. Run commands from this directory unless a command changes
into `real data/`. The code assumes that R and the required packages are
already available; setup and package-installation scripts are not part of this
bundle.

The main simulation uses R 4.5.0 and the packages used by the methods and
reporting scripts (`glmnet`, `mgcv`, `spatstat.univar`, `rootSolve`, and
`data.table`).

## Source files

| Path | Purpose |
|---|---|
| `bootstrap_core.R` | Shared data generation, fitting, bootstrap, ROC, scalar CI, and joint-region code for the unified simulation. |
| `run_all.sh` | Runs or submits one unified simulation job per outer replication. |
| `methods/EMBC/method.R` | Proposed EMBC estimator. |
| `methods/NP/method_np.r` | NP-inc estimator. |
| `methods/LBNP/method_lbnp.r` | LBNP and LBNP-log estimators. |
| `methods/LBNP/original/FUN.R` | Supporting functions used by the LBNP implementation. |
| `simulation/make_parameters.R` | Writes the true values and distribution parameters used for main Table 1. |
| `simulation/make_tables.R` | Summarizes retained scalar results and writes the six simulation LaTeX tables plus summary CSV files. |
| `simulation/make_plots.R` | Reads retained ROC results and writes the five manuscript ROC-panel PNGs. |
| `gof/gof.r` | Goodness-of-fit fitting and bootstrap functions. |
| `gof/test.r` | GOF simulation helpers and scenario definitions. |
| `gof/gof_symm.R` | One outer replication of the lognormal GOF simulation. |
| `gof/gof_symm_weibull.R` | One outer replication of the Weibull GOF simulation. |
| `gof/run_gof.slurm` | SLURM wrapper for the lognormal or Weibull GOF runners. |
| `gof/summarize_gof_symm.R` | Combines GOF replication files into CSV summaries. |
| `gof/make_gof_symm_table.R` | Converts a GOF summary into a LaTeX rejection-rate table. |
| `gof/real_data_gof.R` | Runs the malaria-data model-based GOF bootstrap. |
| `real data/malaria.txt` | The two-column malaria data used by the real-data analysis. |
| `real data/malaria.R` | Fits the four methods and saves malaria point, bootstrap, and ROC results. |
| `real data/summarize_scaled_results.R` | Writes the malaria point-estimate table, ROC plot, and transformation-parameter summary. |
| `real data/summarize_gof.R` | Writes the weighted Sup3 malaria GOF result. |

## Unified simulation

The default grid is `n = 100, 300, 500`, `J = 0.3, 0.5`,
`pi0 = pi1 = 0.75, 0.90`, and distributions `lognorm`, `weibull`, and
`weibull2` (reported as Lognormal, Weibull, and Gamma with `kappa = 0, 1/2,
1`). The default methods are `EMBC`, `EMBC_fixk`, `NP-inc`, `LBNP`, and
`LBNP-log`. The final run uses 1,000 outer replications and 500 bootstrap
samples per replication.

Submit the full study on SLURM with:

```bash
bash run_all.sh --final --reps 1000 --account XXX
```

Run the same workflow directly on a machine with enough resources with:

```bash
bash run_all.sh --local --final --reps 1000 --cpus 10
```

For a small smoke run, use one setting and one replication:

```bash
bash run_all.sh --local --pilot --reps 1 --rep-id 1 \
  --dist lognorm --n 20 --J 3 --pi 9
```

Raw simulation files are written below `result/point_est/`,
`result/coverage/`, `result/coverage_joint/`, and `result/ROC_plot/`. Keep
these files because the reporting scripts read them directly. Run each
replication ID only once in a given output directory: the unified runner
appends rows, so rerunning an ID creates duplicate rows. Use a fresh output
directory for an independent rerun.

After the unified run, regenerate the simulation outputs with:

```bash
Rscript simulation/make_parameters.R
Rscript simulation/make_tables.R
Rscript simulation/make_plots.R
```

`make_tables.R` can take a different raw-result directory as its first
argument; `make_plots.R` can do the same. The scripts require a complete set
of retained raw files and stop before replacing any table when an input is
missing.

## Main-text outputs

| Manuscript item | Generated file | How it is produced |
|---|---|---|
| Table 1: parameters and true accuracy measures | `simulation/latex/table1_parameters.csv` | `Rscript simulation/make_parameters.R` |
| Table 2: point estimation, `pi = 0.90` | `simulation/latex/point_main_pi90.tex` | `Rscript simulation/make_tables.R` |
| Table 3: coverage and average length, `pi = 0.90` | `simulation/latex/coverage_main_pi90.tex` | `Rscript simulation/make_tables.R`; remove the generated `Ours_fixk` columns before including the four-method paper table |
| Figure 1: ROC panels for `n = 300`, `J = 0.5`, `pi = 0.90` | `simulation/plots/n300J5pi9_roc_band_panel.png` | `Rscript simulation/make_plots.R` |
| Table 4: malaria estimates and 95% bootstrap CIs | `real data/biomarker_div_100000/malaria_div100000_point_estimates.tex` | `Rscript summarize_scaled_results.R biomarker_div_100000` from `real data/` |
| Figure 2: malaria ROC curves | `real data/biomarker_div_100000/malaria_div100000_roc.pdf` and `.png` | `Rscript summarize_scaled_results.R biomarker_div_100000` from `real data/` |
| Section 6: estimated transformation parameter and CI | `real data/biomarker_div_100000/kappa_summary.csv` | `Rscript summarize_scaled_results.R biomarker_div_100000` from `real data/` |
| Section 6: malaria GOF result | `real data/biomarker_div_100000/gof_sup3.csv` | GOF commands in the real-data section below |

Table 1 is intentionally a CSV export of the unrounded values. Format or
include it in the manuscript as needed; `make_parameters.R` does not create a
LaTeX wrapper for it. The simulation table script also writes
`point_appendix_kappa.tex`, a four-method point-estimate export for both pi
values that is not a numbered manuscript table.

## Supplementary simulation outputs

The supplement's additional simulation section is Section S5.3.

| Supplement item | Generated file | How it is produced |
|---|---|---|
| Table S1: estimated-k versus fixed-k coverage | `simulation/latex/coverage_appendix_kappa.tex` | `Rscript simulation/make_tables.R`; combine the generated pi=0.90 and pi=0.75 blocks side by side, retain only `Ours` and `Ours_fixk`, and relabel them `Our` and `Our-F` |
| Table S2: lognormal GOF size and power | `gof/result_symm/summary/gof_rejection_table.tex` | Lognormal GOF run, then the two GOF summary commands below |
| Table S3: point estimation, `pi = 0.75` | `simulation/latex/point_appendix_pi75.tex` | `Rscript simulation/make_tables.R` |
| Table S4: coverage and average length, `pi = 0.75` | `simulation/latex/coverage_appendix_pi75.tex` | `Rscript simulation/make_tables.R`; remove the generated `Ours_fixk` columns before including the four-method paper table |
| Figure S1: `J = 0.3`, `pi = 0.75` | `simulation/plots/J3pi75_roc_band_panel.png` | `Rscript simulation/make_plots.R` |
| Figure S2: `J = 0.3`, `pi = 0.90` | `simulation/plots/J3pi9_roc_band_panel.png` | `Rscript simulation/make_plots.R` |
| Figure S3: `J = 0.5`, `pi = 0.75` | `simulation/plots/J5pi75_roc_band_panel.png` | `Rscript simulation/make_plots.R` |
| Figure S4: `J = 0.5`, `pi = 0.90` | `simulation/plots/J5pi9_roc_band_panel.png` | `Rscript simulation/make_plots.R` |

The four ROC figures contain panels for all three sample sizes and all three
reported distributions. The GOF table contains several tests; for Table S2,
use only rows with `test = sup3_quantile` in `gof_summary.csv`. This is the
weighted statistic displayed as Sup3-Q in the generated table. The supplement
mentions Weibull and Gamma GOF patterns narratively, but only the lognormal GOF
table is a reported supplementary table. This bundle has a Weibull GOF runner
and no Gamma GOF runner.

### GOF simulation for Table S2

Submit the 1,000 lognormal outer replications with:

```bash
sbatch --array=1-1000 --account=XXX --cpus-per-task=10 \
  --mem=8G --time=1-12:00:00 gof/run_gof.slurm lognormal
```

On a non-SLURM machine, run one replication with:

```bash
Rscript gof/gof_symm.R --rep-id 1 --B 500 --cores 10 --scenario both
```

Repeat the command for replication IDs 1 through 1,000, using each ID once.
Then summarize and write the table:

```bash
Rscript gof/summarize_gof_symm.R --out-dir gof/result_symm
Rscript gof/make_gof_symm_table.R gof/result_symm
```

The summary directory contains `gof_results.csv`, `gof_summary.csv`,
`gof_failures.csv`, and `gof_rejection_table.tex`.

The Weibull diagnostic run, which is not a numbered table in the supplement,
uses the same commands with `gof_symm_weibull.R`,
`gof/result_symm_weibull`, and the `weibull` argument to `run_gof.slurm`.

## Malaria real-data analysis

The included `real data/malaria.txt` has two headerless columns: nominal group
`R` and biomarker `T`. The manuscript analysis removes zero biomarker values,
divides `T` by 100,000, and uses `pi0 = 1` and `pi1 = 0.677`.

Run the four-method analysis and generate Table 4, Figure 2, and the kappa
summary with:

```bash
cd "real data"
Rscript malaria.R 500 10 "EMBC,NP-inc,LBNP,LBNP-log" \
  overwrite scale:100000 biomarker_div_100000
Rscript summarize_scaled_results.R biomarker_div_100000
cd ..
```

The analysis writes `result.rds`, `estimates.csv`, and `roc_bands.csv` before
the summary script writes the manuscript table, ROC PDF/PNG, and
`kappa_summary.csv` under `real data/biomarker_div_100000/`.

For the reported weighted Sup3 GOF result, run from `code/`:

```bash
Rscript gof/real_data_gof.R --variant div100000 --B 500 \
  --cores 10 --seed 20260818
cd "real data"
Rscript summarize_gof.R
cd ..
```

This writes the bootstrap files under `gof/real_data_result/div100000/` and
writes `real data/biomarker_div_100000/gof_sup3.csv`. The reported statistic is
the weighted `sup3 = (n0/n)D0 + (n1/n)D1` result.

## Output conventions

All paths above are relative to `code/`. The generated LaTeX files are complete
table environments and can be included in the manuscript. The generated PNG
and PDF files are the figure artifacts. Preserve generated raw results so the
reporting scripts can regenerate the manuscript outputs.
