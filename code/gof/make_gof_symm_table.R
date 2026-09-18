#!/usr/bin/env Rscript

# Convert the symmetric-sample GOF summary to a manuscript-ready LaTeX table.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript gof/make_gof_symm_table.R [RESULT_DIR]")
result_dir <- if (length(args)) args[1] else file.path("gof", "result_symm")
input_file <- file.path(result_dir, "summary", "gof_summary.csv")
output_file <- file.path(result_dir, "summary", "gof_rejection_table.tex")

if (!file.exists(input_file)) {
  stop("Run this script from the repository root; cannot find ", input_file)
}

x <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "config_id", "distribution", "n", "pi", "J", "scenario", "test",
  "b0", "b1", "reps", "successful_reps", "failed_reps", "rejection_rate",
  "mean_bootstrap_successes"
)
missing <- setdiff(required, names(x))
if (length(missing) > 0L) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

distribution <- unique(x$distribution)
if (length(distribution) != 1L || !distribution %in% c("lognormal", "weibull")) {
  stop("Provide one lognormal or Weibull summary directory at a time.")
}
shape_label <- if (distribution == "lognormal") "log-variances" else "shapes"
dist_label <- if (distribution == "lognormal") "lognormal" else "Weibull"

test_order <- c(
  "residual_pvalue", "residual_quantile", "sup1_quantile",
  "sup2_quantile", "sup3_quantile"
)
unknown_tests <- setdiff(unique(x$test), test_order)
if (length(unknown_tests) > 0L) {
  stop("Unknown test(s): ", paste(unknown_tests, collapse = ", "))
}
if (anyDuplicated(paste(x$config_id, x$scenario, x$test, sep = "|"))) {
  stop("Duplicate configuration/scenario/test rows found in ", input_file)
}

format_rate <- function(value, bold = FALSE) {
  if (length(value) != 1L || !is.finite(value)) return("--")
  out <- formatC(value, format = "f", digits = 3)
  if (bold) paste0("\\textbf{", out, "}") else out
}

scenario_titles <- c(homo = "empirical size", hete = "empirical power")
present_scenarios <- names(scenario_titles)[names(scenario_titles) %in% x$scenario]
if (length(present_scenarios) == 0L) stop("No recognized scenarios found.")

body <- character()
for (scenario in present_scenarios) {
  z <- x[x$scenario == scenario, , drop = FALSE]
  variances <- unique(z[c("b0", "b1")])
  if (nrow(variances) != 1L) {
    stop("Multiple parameter pairs found for scenario ", scenario)
  }
  scenario_label <- paste0(
    if (scenario == "homo") "Panel A: Homogeneous" else "Panel B: Heterogeneous",
    " ", shape_label, " $(", format(variances$b0, trim = TRUE), ",",
    format(variances$b1, trim = TRUE), ")$ (", scenario_titles[[scenario]], ")"
  )
  configs <- unique(z[c("config_id", "n", "pi", "J")])
  configs <- configs[order(configs$n, -configs$pi, configs$J), , drop = FALSE]

  body <- c(
    body,
    paste0("\\multicolumn{9}{l}{\\textit{", scenario_label, "}} \\\\"),
    "\\addlinespace[2pt]"
  )
  for (i in seq_len(nrow(configs))) {
    one <- z[z$config_id == configs$config_id[i], , drop = FALSE]
    if (!setequal(one$test, test_order)) {
      stop("Incomplete test set for ", configs$config_id[i], "/", scenario)
    }
    one <- one[match(test_order, one$test), , drop = FALSE]
    successful <- unique(one$successful_reps)
    if (length(successful) != 1L) {
      stop("Test-specific successful replication counts for ",
           configs$config_id[i], "/", scenario)
    }

    rates <- one$rejection_rate
    bold <- rep(FALSE, length(rates))
    if (scenario == "hete" && any(is.finite(rates))) {
      bold <- is.finite(rates) & abs(rates - max(rates, na.rm = TRUE)) < 1e-12
    }
    cells <- mapply(format_rate, rates, bold, USE.NAMES = FALSE)
    body <- c(
      body,
      paste(
        format(configs$n[i], trim = TRUE),
        formatC(configs$pi[i], format = "f", digits = 2),
        formatC(configs$J[i], format = "f", digits = 1),
        format(successful, trim = TRUE),
        paste(cells, collapse = " & "),
        sep = " & "
      ) |>
        paste0(" \\\\")
    )
  }
  if (scenario != tail(present_scenarios, 1L)) body <- c(body, "\\midrule")
}

bootstrap_counts <- unique(x$mean_bootstrap_successes)
bootstrap_note <- if (length(bootstrap_counts) == 1L &&
                      is.finite(bootstrap_counts)) {
  paste0(" Each test uses ", format(bootstrap_counts, trim = TRUE),
         " successful bootstrap fits.")
} else {
  ""
}

lines <- c(
  paste0("% Auto-generated from ", input_file),
  "% Requires \\usepackage{booktabs}",
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4.5pt}",
  paste0("\\caption{Empirical rejection rates for the goodness-of-fit tests with symmetric ", dist_label, " samples.}"),
  paste0("\\label{tab:gof-symm-", distribution, "-rejection}"),
  "\\begin{tabular}{rccrrrrrr}",
  "\\toprule",
  paste0(
    "$n_0=n_1$ & $\\pi_0=\\pi_1$ & $J$ & Successful & Residual-$p$ & ",
    "Residual-Q & Sup1-Q & Sup2-Q & Sup3-Q \\\\"
  ),
  "& & & replications & & & & & \\\\ ",
  "\\midrule",
  body,
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{2pt}",
  "\\begin{minipage}{0.96\\textwidth}",
  "\\footnotesize",
  paste0(
    "\\textit{Notes:} The nominal level is 0.05. Residual-$p$ denotes the ",
    "normal-approximation residual test, Residual-Q denotes the bootstrap-quantile ",
    "residual test, and Sup1-Q--Sup3-Q denote the three bootstrap-quantile ",
    "supremum tests. Boldface identifies the largest empirical power within each ",
    "heterogeneous row. Rejection rates exclude failed replications. The J labels refer to the homogeneous baseline settings.", bootstrap_note
  ),
  "\\end{minipage}",
  "\\end{table}"
)

writeLines(lines, output_file)
cat("Wrote", output_file, "\n")
