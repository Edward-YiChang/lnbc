# Reproducing the manuscript results

This folder contains the R code and malaria data for the main text and
Supplementary Material. It covers the point-estimation and coverage studies,
ROC figures, lognormal and Weibull goodness-of-fit (GOF) studies, and malaria
analysis. Gamma (k=1) GOF will be added separately. The theoretical derivations
in supplementary Sections S1–S4 do not require separate numerical runs.

Run commands from this folder unless a change of directory is shown. Full
simulation and bootstrap runs should be executed on a computing server.
Generated results are written locally; no precomputed results are required.

## 1. Software and installation

Use R 4.5.0 and a Unix-like system with bash. Required R packages are `glmnet`,
`mgcv`, `spatstat.univar`, `rootSolve`, and `data.table`.

```bash
bash setup_server.sh
```

The installer installs missing packages from CRAN. On Alliance, the scripts use
`StdEnv/2023` and `r/4.5.0`. For another server, load R and use the direct-run
commands below. LaTeX tables require `booktabs`, `multirow`, and `graphicx`;
define `\newcommand{\kk}{\kappa}` for the simulation tables.

## 2. Results-to-code guide

Paths in this table are relative to this folder. Sections below give the commands.

| Manuscript result | Generation step | Output |
|---|---|---|
| Main Table 1: simulation parameters and true accuracy measures | Section 3, parameter command | `simulation/latex/table1_parameters.csv` |
| Main Table 2: point estimation, pi=0.90 | Section 3, table command | `simulation/latex/point_main_pi90.tex` |
| Main Table 3: coverage and interval length, pi=0.90 | Section 3, table command | `simulation/latex/coverage_main_pi90.tex` |
| Main Figure 1: ROC curves, n=300, J=0.5, pi=0.90 | Section 3, plot command | `simulation/plots/n300J5pi9_roc_band_panel.png` |
| Supplement Table S1 / Section S5.1: fixed-k comparison | Section 3, table command | `simulation/latex/coverage_appendix_kappa.tex`: Ours and Ours_fixk columns for both pi values |
| Supplement Table S2: point estimation, pi=0.75 | Section 3, table command | `simulation/latex/point_appendix_pi75.tex` |
| Supplement Table S3: coverage, pi=0.75 | Section 3, table command | `simulation/latex/coverage_appendix_pi75.tex` |
| Supplement Section S5.2: GOF size and power | Section 4 | `gof/result_symm/summary/` and `gof/result_symm_weibull/summary/` |
| Supplement Section S5.4 / Figures S1–S4: additional ROC panels | Section 3, plot command | Four files `simulation/plots/J{3,5}pi{75,9}_roc_band_panel.png` |
| Main Table 4: malaria point estimates and CIs | Section 5 | `real data/biomarker_div_100000/malaria_div100000_point_estimates.tex` |
| Main Figure 2: malaria ROC curves | Section 5 | `real data/biomarker_div_100000/malaria_div100000_roc.pdf` (also PNG) |
| Main Section 6: estimated k and its CI | Section 5 | `real data/biomarker_div_100000/kappa_summary.csv` |
| Main Section 6: malaria GOF result | Section 6 | `real data/biomarker_div_100000/gof_sup3.csv` |

Coverage tables include all five methods. For the fixed-k comparison in
Table S1, use the Ours and Ours_fixk columns; this output contains one table
per pi value. The four supplementary ROC files are identified by J and pi. `J3`/`J5` mean 0.3/0.5; `pi75`/`pi9` mean 0.75/0.90.

## 3. Point estimation, coverage, and ROC simulations

### Methods and distributions

| Manuscript label | Internal method key | Reports |
|---|---|---|
| Ours | `EMBC` | Point, coverage, ROC |
| Ours_fixk | `EMBC_fixk` | Coverage |
| NP | `NP-inc` | Point, coverage, ROC |
| LBNP | `LBNP` | Point, coverage, ROC |
| LBNP-log | `LBNP-log` | Point, coverage, ROC |

NP always denotes `NP-inc` in the reports, including the malaria analysis.
Ours_fixk estimates k from the original data and holds that estimate fixed in
the bootstrap. LBNP-log applies the method's log transformation.

| Manuscript distribution | Internal key |
|---|---|
| Lognormal (k=0) | `lognorm` |
| Weibull (k=1/2) | `weibull` |
| Gamma (k=1) | `weibull2` |

`weibull2` generates a Weibull distribution of shape 1, equivalently Gamma of
shape 1. Do not replace it with the legacy `gamma` key. In Table 1, lognormal
parameters are log-mean and log-variance; Weibull parameters are scale and shape;
Gamma parameters are shape and scale.

The study uses n0=n1 in {100,300,500}, J in {0.3,0.5}, and pi0=pi1 in
{0.75,0.90}, giving 36 settings per method. Each setting has 1,000 outer
replications and 500 nonparametric bootstrap samples. CIs are 95% percentile
intervals. The ROC grid is 0,0.01,...,1.

### Run the study

On Alliance/SLURM:

```bash
bash run_all.sh --final --reps 1000 --account YOUR_ACCOUNT
```

The default grid includes all five methods and all three distributions. One job
is submitted per outer replication. Defaults are 8 nodes, 10 cores per node,
16 GB per node, and two days; adjust `--nodes`, `--cpus`, `--mem`, and `--time`
for your allocation. On a non-SLURM server, run directly:

```bash
bash run_all.sh --local --final --reps 1000 --cpus 10
```

Raw outputs go to `result/point_est/`, `result/coverage/`, `result/ROC_plot/`,
and `result/coverage_joint/`, separated by internal method key. The joint-region
files are produced by the shared pipeline but are not used in these tables.
Run each replication ID only once into a given output directory: this runner
appends rows. Use a fresh directory for an independent rerun. The bootstrap
base seed is 123; data/scenario seeds are derived in `bootstrap_core.R`.

### Generate manuscript outputs

After the simulations finish:

```bash
Rscript simulation/make_parameters.R
Rscript simulation/make_tables.R
Rscript simulation/make_plots.R
```

The parameter command does not require simulation outputs. It writes unrounded
values for Table 1; display them to two decimal places in the manuscript.
The table command writes six LaTeX files plus `point_summary.csv` and
`coverage_summary.csv` in `simulation/latex/`. The plot command writes five PNGs.
Both reporting scripts can alternatively take a raw-result directory as their
first argument. They accept fresh `_roc.csv` and retained `_roc_band.csv` names,
preferring `_roc.csv` when both exist.

Point tables report 100 x relative bias and 100 x MSE. Coverage tables report
CP in percent and average interval length. In the summary CSVs, `rb`, `mse`,
and `cp` are unscaled; `al` is interval length. The reported metrics are ROC(0.2),
AUC, J, sensitivity (`se`), and specificity (`sp`).

For the equal-weight averages across the 36 settings discussed in Section S5.1,
after generating the tables run:

```bash
Rscript -e 'x <- read.csv("simulation/latex/coverage_summary.csv"); x <- subset(x, method %in% c("EMBC", "EMBC_fixk")); x$cp <- 100*x$cp; print(aggregate(cp ~ method + metric, x, mean))'
```

Nonfinite estimates/lengths and missing coverage flags are excluded from their
respective means. Available replication and missing-value counts are recorded
in the CSVs. LBNP failure columns count missing estimates among recorded
replications. The retained study used fewer than 1,000 recorded replications in
some files (984 for Ours_fixk); completing all 1,000 can change the reported
Monte Carlo values. R/package versions can also affect numerical results.

## 4. GOF simulation results: supplementary Section S5.2

Both runners include homogeneous and heterogeneous scenarios:

| Runner | Homogeneous | Heterogeneous | Output directory |
|---|---|---|---|
| `gof/gof_symm.R` | Lognormal log-variances (1,1) | Log-variances (1,4) | `gof/result_symm/` |
| `gof/gof_symm_weibull.R` | Weibull shapes (0.5,0.5) | Shapes (0.5,0.25) | `gof/result_symm_weibull/` |

Each replication covers n in {100,300,500}, pi in {0.9,0.75}, J labels in
{0.3,0.5}, and both scenarios: 24 settings per distribution. Lognormal means
are (0,0.77) or (0,1.35); Weibull scales are (0.5,2.6793) or (0.5,9.7345).
Heterogeneous cases retain these means/scales, so their J labels refer to the
homogeneous baseline. Gamma (k=1) GOF is not included yet.

Submit from this folder:

```bash
sbatch --array=1-1000 --account=YOUR_ACCOUNT --cpus-per-task=10 --mem=8G --time=1-12:00:00 gof/run_gof.slurm lognormal
sbatch --array=1-1000 --account=YOUR_ACCOUNT --cpus-per-task=10 --mem=8G --time=1-12:00:00 gof/run_gof.slurm weibull
```

Resource requests are examples. Each task uses B=500 successful model-based
bootstrap fits, alpha=0.05, and base seed 20260818. At most 2B bootstrap attempts
are allowed; unsuccessful replications are recorded as failures. On a non-SLURM
server, run the following for each replication ID 1 through 1000:

```bash
Rscript gof/gof_symm.R --rep-id 1 --B 500 --cores 10 --scenario both
Rscript gof/gof_symm_weibull.R --rep-id 1 --B 500 --cores 10 --scenario both
```

Summarize after completion:

```bash
Rscript gof/summarize_gof_symm.R --out-dir gof/result_symm
Rscript gof/make_gof_symm_table.R gof/result_symm
Rscript gof/summarize_gof_symm.R --out-dir gof/result_symm_weibull
Rscript gof/make_gof_symm_table.R gof/result_symm_weibull
```

Each output directory receives a `summary/` folder with `gof_results.csv`,
`gof_summary.csv`, `gof_failures.csv`, and `gof_rejection_table.tex`. The five
implemented GOF variants are retained. The manuscript's weighted statistic
is `sup3_quantile`, displayed as Sup3-Q in these comparison tables:
Dn = (n0/n)D0 + (n1/n)D1. The `_quantile` tests use bootstrap critical values;
the normal-approximation residual test uses its p-value. Summaries average saved
rejection flags over nonmissing flags and report failure counts.

GOF runners replace a scenario file when its replication ID is rerun. Keep the
two distributions and different B/seed settings in separate directories.

## 5. Malaria Table 4, Figure 2, and transformation estimate

The included `real data/malaria.txt` has 408 rows and two columns without a
header: nominal group R (0=dry season, 1=wet season) and parasite density T.
Removing zero densities leaves 81 dry-season and 211 wet-season observations.
The analysis divides T by 100,000 and uses pi0=1 and pi1=0.677, as in Section 6.
The data are copied unchanged from the repository's malaria analysis. The study
is described by Kitua et al. (1996), Tropical Medicine & International Health
1:475–484, and Qin and Leung (2005), Biometrics 61:456–464.

On a computing server, run:

```bash
cd "real data"
Rscript malaria.R 500 10 "EMBC,NP-inc,LBNP,LBNP-log" overwrite scale:100000 biomarker_div_100000
Rscript summarize_scaled_results.R biomarker_div_100000
cd ..
```

The first command fits the four methods and uses 500 nonparametric bootstrap
samples with seed 123. It writes `result.rds`, `estimates.csv`, and
`roc_bands.csv` under `real data/biomarker_div_100000/`. The second produces
Table 4, Figure 2 (PDF and PNG), and `kappa_summary.csv`. Table and k intervals
use successful bootstrap fits with valid ROC curves. NP-inc is displayed as NP.

## 6. Malaria GOF: weighted sup3 statistic

Run from this folder on a computing server:

```bash
Rscript gof/real_data_gof.R --variant div100000 --B 500 --cores 10 --seed 20260818
cd "real data"
Rscript summarize_gof.R
cd ..
```

This performs a separate model-based bootstrap on the same rescaled malaria
data and saves the fit/bootstrap under `gof/real_data_result/div100000/`.
The summary selects sup3 = (n0/n)D0 + (n1/n)D1 and writes
`real data/biomarker_div_100000/gof_sup3.csv`. The p-value is the proportion of
bootstrap sup3 values at least as large as observed sup3, with no plus-one
correction; reject when p <= 0.05, as in Section 3.3.

For the retained bootstrap draws, this convention gives 299/500 = 0.598.
Use the weighted sup3 result for the manuscript GOF analysis.

## 7. File organization

- `bootstrap_core.R` and `methods/`: common bootstrap and estimator code.
- `simulation/`: parameter, table, and figure generation.
- `gof/`: GOF implementation, simulation runners, and reporting.
  `gof/test.r` defines the shared simulation functions required by both runners.
- `real data/`: malaria input data, analysis, and reporting.
- `run_all.sh`, `setup_server.sh`, `install_packages.R`: execution/setup.

Generated simulation results are not bundled. All commands regenerate outputs
from the supplied code and, for the real-data analysis, the supplied malaria data.
