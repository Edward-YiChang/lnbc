# Reproducibility code

This directory contains the code needed to run the simulation study, the
malaria analysis, and the reported ROC plots. It assumes that R and the
packages used by the analysis are already available. Generated raw results
under `result/` are the inputs to the reporting scripts.

All paths below are relative to `code/`; run the real-data scripts from the
`real data/` directory.

## Files

| Path | Purpose |
|---|---|
| `bootstrap_core.R` | Shared data generation, model fitting, bootstrap, confidence intervals, joint regions, and ROC calculations. |
| `run_all.sh` | Runs the unified simulation grid. |
| `methods/EMBC/method.R` | EMBC estimator. |
| `methods/NP/method_np.r` | NP-inc estimator. |
| `methods/LBNP/method_lbnp.r` | LBNP and LBNP-log estimators. |
| `methods/LBNP/original/FUN.R` | Supporting LBNP functions. |
| `simulation/make_parameters.R` | Writes true parameter and accuracy values for Table 1 to `simulation/summary/table1_parameters.csv`. |
| `simulation/summarize_results.R` | Summarizes retained simulation results into numeric point-estimate and coverage CSV files. |
| `simulation/make_plots.R` | Produces the five simulation ROC-panel figures. |
| `real data/malaria.txt` | Headerless malaria data with group and biomarker columns. |
| `real data/malaria.R` | Fits the four methods and saves malaria estimates, bootstrap results, and ROC bands. |
| `real data/summarize_scaled_results.R` | Produces malaria ROC plots, numeric point estimates, and the kappa summary. |
| `gof/gof.r`, `gof/test.r` | Goodness-of-fit functions and scenario definitions. |
| `gof/gof_symm.R`, `gof/gof_symm_weibull.R` | Lognormal and Weibull GOF replication runners. |
| `gof/run_gof.slurm` | Batch wrapper for GOF runners. |
| `gof/summarize_gof_symm.R` | Combines GOF replications into numeric CSV summaries. |
| `gof/real_data_gof.R`, `real data/summarize_gof.R` | Computes and summarizes the malaria weighted Sup3 GOF result. |

## Simulation outputs

Run the unified simulation with `run_all.sh`. It stores raw point estimates,
confidence intervals, joint regions, and ROC results under `result/`. Use
`simulation/summarize_results.R` to obtain numeric summaries in
`simulation/summary/`, and use `simulation/make_plots.R` to obtain ROC figures
in `simulation/plots/`.

The default grid uses `n = 100, 300, 500`, `J = 0.3, 0.5`,
`pi0 = pi1 = 0.75, 0.90`, and Lognormal, Weibull, and Gamma distributions
(`kappa = 0, 1/2, 1`). Retain the raw result files because all summary and
plot files are regenerated from them.

## Manuscript tables

| Item | Numeric source | Contents |
|---|---|---|
| Table 1 | `simulation/summary/table1_parameters.csv` | Simulation parameters and true accuracy measures. |
| Table 2 | `simulation/summary/point_summary.csv` | Point-estimation summaries for `pi = 0.90`. |
| Table 3 | `simulation/summary/coverage_summary.csv` | Coverage probabilities and average interval lengths for `pi = 0.90`. |
| Table 4 | `real data/biomarker_div_100000/malaria_div100000_point_estimates.csv` | Malaria point estimates and percentile bootstrap intervals. |
| Table S1 | `simulation/summary/coverage_summary.csv` | Estimated-k versus fixed-k coverage; select `EMBC` and `EMBC_fixk`. |
| Table S2 | `gof/result_symm/summary/gof_summary.csv` | Lognormal GOF rejection rates; use the `sup3_quantile` rows. |
| Table S3 | `simulation/summary/point_summary.csv` | Point-estimation summaries for `pi = 0.75`. |
| Table S4 | `simulation/summary/coverage_summary.csv` | Coverage probabilities and average interval lengths for `pi = 0.75`. |

The CSV files contain the numeric values used in the manuscript. Formatting
them as tables is left to the manuscript source. For Tables 2 and 3 and their
supplementary counterparts, report `100 * rb`, `100 * mse`, and `100 * cp`
to match the manuscript's scales; Tables 3 and S4 should exclude `EMBC_fixk`
unless the fixed-k comparison is being shown.

## Manuscript figures

`simulation/make_plots.R` produces these files in `simulation/plots/`:

| Item | Figure file |
|---|---|
| Figure 1 | `n300J5pi9_roc_band_panel.png` |
| Figure S1 | `J3pi75_roc_band_panel.png` |
| Figure S2 | `J3pi9_roc_band_panel.png` |
| Figure S3 | `J5pi75_roc_band_panel.png` |
| Figure S4 | `J5pi9_roc_band_panel.png` |

The four supplementary panels contain all three sample sizes. Figure 1 is the
`n = 300`, `J = 0.5`, `pi = 0.90` panel.

For the malaria analysis, `real data/malaria.R` creates the saved result and
`real data/summarize_scaled_results.R` produces:

- Figure 2: `real data/biomarker_div_100000/malaria_div100000_roc.pdf` and
  `malaria_div100000_roc.png`.
- Table 4 source: `malaria_div100000_point_estimates.csv`.
- Transformation summary: `kappa_summary.csv`.

The GOF scripts produce Table S2's numeric inputs in
`gof/result_symm/summary/`. The real-data GOF scripts produce
`real data/biomarker_div_100000/gof_sup3.csv`.
