#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd())
if (!file.exists(file.path(repo_root, "bootstrap_core.R"))) {
  stop("Run this script from the repository root: Rscript simulation/make_tables.R")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.")
}

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else file.path(repo_root, "result")
output_dir <- file.path(repo_root, "simulation", "latex")

methods <- c("EMBC", "NP-inc", "LBNP", "LBNP-log")
coverage_methods <- c("EMBC", "EMBC_fixk", "NP-inc", "LBNP", "LBNP-log")
dists <- c("lognorm", "weibull", "weibull2")
n_values <- c(100, 300, 500)
J_values <- c("3", "5")
pi_values <- c("75", "9")
metric_order <- c("ROC_0.20", "AUC", "J", "se", "sp")
roc_suffix <- c(
  EMBC = "roc_band.csv",
  `NP-inc` = "roc.csv",
  LBNP = "roc_band.csv",
  `LBNP-log` = "roc_band.csv"
)

result_path <- function(type, method, dist, n, J, pi) {
  stem <- paste0("n", n, "J", J, "pi", pi, dist, "_final_B500_reps1000_")
  suffix <- if (type == "coverage") "ci.csv" else roc_suffix[[method]]
  path <- file.path(result_dir, type, method, paste0(stem, suffix))
  unified <- file.path(result_dir, type, method, paste0(stem, "roc.csv"))
  if (type == "ROC_plot" && file.exists(unified)) unified else path
}

required <- expand.grid(
  method = methods,
  dist = dists,
  n = n_values,
  J = J_values,
  pi = pi_values,
  stringsAsFactors = FALSE
)
required$coverage_path <- mapply(
  result_path, "coverage", required$method, required$dist,
  required$n, required$J, required$pi, USE.NAMES = FALSE
)
required$roc_path <- mapply(
  result_path, "ROC_plot", required$method, required$dist,
  required$n, required$J, required$pi, USE.NAMES = FALSE
)

stopifnot(nrow(required) == 144L)
all_paths <- c(required$coverage_path, required$roc_path)
missing_paths <- all_paths[!file.exists(all_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing ", length(missing_paths), " required raw inputs:\n",
    paste(missing_paths, collapse = "\n")
  )
}

coverage_columns <- c(
  "method", "dist", "n", "J_key", "pi_key", "pi", "B", "total_reps",
  "rep", "metric", "estimate", "ci_length", "true", "covered"
)
roc_columns <- c(
  "method", "dist", "n", "J_key", "pi_key", "pi", "B", "total_reps",
  "rep", "fpr", "tpr", "ci_low", "ci_high", "true_tpr", "covered"
)

bad_schema <- character()
for (path in required$coverage_path) {
  columns <- names(read.csv(path, nrows = 0L, check.names = FALSE))
  if (!all(coverage_columns %in% columns)) bad_schema <- c(bad_schema, path)
}
for (path in required$roc_path) {
  columns <- names(read.csv(path, nrows = 0L, check.names = FALSE))
  if (!all(roc_columns %in% columns)) bad_schema <- c(bad_schema, path)
}
if (length(bad_schema) > 0L) {
  stop("Required columns are missing from:\n", paste(bad_schema, collapse = "\n"))
}
message("Validated 144 scalar and 144 ROC raw-result files.")

safe_mean <- function(x) {
  if (length(x) == 0L || !any(is.finite(x))) return(NA_real_)
  mean(x[is.finite(x)])
}

summarize_rows <- function(x, metric_name = NULL) {
  if (!is.null(metric_name)) x[, metric := metric_name]
  point <- x[, {
    truth_value <- true[is.finite(true)][1L]
    if (length(truth_value) == 0L) truth_value <- NA_real_
    list(
      n_reps_available = data.table::uniqueN(rep),
      true = truth_value,
      mean_est = safe_mean(estimate),
      rb = if (is.finite(truth_value) && truth_value != 0) {
        safe_mean((estimate - truth_value) / truth_value)
      } else {
        safe_mean(estimate - truth_value)
      },
      mse = safe_mean((estimate - truth_value)^2),
      n_missing = sum(!is.finite(estimate))
    )
  }, by = .(method, dist, n, J_key, pi_key, pi, metric, B, total_reps)]
  coverage <- x[, list(
    n_reps_available = data.table::uniqueN(rep),
    n_valid = sum(!is.na(covered)),
    n_missing = sum(is.na(covered)),
    cp = safe_mean(as.numeric(covered)),
    al = safe_mean(ci_length)
  ), by = .(method, dist, n, J_key, pi_key, pi, metric, B, total_reps)]
  list(point = point, coverage = coverage)
}

point_parts <- vector("list", 2L * nrow(required))
coverage_parts <- vector("list", 2L * nrow(required))
part <- 0L

for (i in seq_len(nrow(required))) {
  scalar <- data.table::fread(
    required$coverage_path[[i]],
    select = coverage_columns,
    showProgress = FALSE
  )
  scalar <- scalar[metric %in% c("AUC", "J", "se", "sp")]
  if (scalar[, anyDuplicated(paste(rep, metric, sep = ":"))] != 0L) {
    stop("Duplicated replication/metric rows in ", required$coverage_path[[i]])
  }
  scalar_summary <- summarize_rows(scalar)
  part <- part + 1L
  point_parts[[part]] <- scalar_summary$point
  coverage_parts[[part]] <- scalar_summary$coverage

  roc <- data.table::fread(
    required$roc_path[[i]],
    select = roc_columns,
    showProgress = FALSE
  )
  roc <- roc[abs(fpr - 0.2) < 1e-12]
  if (nrow(roc) == 0L) stop("ROC(0.2) is missing from ", required$roc_path[[i]])
  if (roc[, anyDuplicated(rep)] != 0L) {
    stop("Duplicated ROC(0.2) replication rows in ", required$roc_path[[i]])
  }
  roc[, `:=`(
    estimate = tpr,
    true = true_tpr,
    ci_length = ci_high - ci_low
  )]
  roc_summary <- summarize_rows(roc, "ROC_0.20")
  part <- part + 1L
  point_parts[[part]] <- roc_summary$point
  coverage_parts[[part]] <- roc_summary$coverage
}

point <- data.table::rbindlist(point_parts)
coverage <- data.table::rbindlist(coverage_parts)
expected_rows <- length(methods) * length(dists) * length(n_values) *
  length(J_values) * length(pi_values) * length(metric_order)
if (nrow(point) != expected_rows || nrow(coverage) != expected_rows) {
  stop("Expected ", expected_rows, " summary rows per result type; got ",
       nrow(point), " point and ", nrow(coverage), " coverage rows.")
}


# Fixed-k coverage is stored for all five reported metrics in scalar CI files.
fixed_parts <- list()
for (dist in dists) for (n in n_values) for (J in J_values) for (pi in pi_values) {
  path <- result_path("coverage", "EMBC_fixk", dist, n, J, pi)
  if (!file.exists(path)) stop("Missing fixed-k input: ", path)
  x <- data.table::fread(path, select = coverage_columns, showProgress = FALSE)
  x <- x[metric %in% metric_order]
  if (anyDuplicated(paste(x$rep, x$metric))) stop("Duplicate rows in ", path)
  if (!setequal(x$metric, metric_order)) stop("Missing metrics in ", path)
  fixed_parts[[length(fixed_parts) + 1L]] <- summarize_rows(x)$coverage
}
coverage <- data.table::rbindlist(c(list(coverage), fixed_parts))
stopifnot(nrow(coverage) == expected_rows / length(methods) * length(coverage_methods))

metric_tex <- c(
  AUC = "AUC",
  ROC_0.20 = "$\\mathrm{ROC}(0.2)$",
  se = "$\\eta$",
  sp = "$\\tau$",
  J = "$J$"
)
dist_tex <- c(
  lognorm = "\\shortstack{Lognormal\\\\$(\\kk=0)$}",
  weibull = "\\shortstack{Weibull\\\\$(\\kk=1/2)$}",
  weibull2 = "\\shortstack{Gamma\\\\$(\\kk=1)$}"
)
profiles <- list(
  main_pi90 = list(dists = dists, pis = "9"),
  appendix_pi75 = list(dists = dists, pis = "75"),
  appendix_kappa = list(dists = dists, pis = c("75", "9"))
)

pick <- function(x, method, dist, n, J, pi, metric) {
  method_value <- method
  dist_value <- dist
  n_value <- n
  J_value <- J
  pi_value <- pi
  metric_value <- metric
  y <- x[
    method == method_value & dist == dist_value & n == n_value &
      as.character(J_key) == J_value & as.character(pi_key) == pi_value &
      metric == metric_value
  ]
  if (nrow(y) == 0L) return(NULL)
  y[1L]
}
num_fixed <- function(x, digits) {
  if (length(x) == 0L || !is.finite(x)) "--" else formatC(x, format = "f", digits = digits)
}
num_point <- function(x) {
  if (length(x) == 0L || !is.finite(x)) return("--")
  num_fixed(x, if (abs(x) >= 10) 1 else 2)
}
pi_caption <- function(pi) if (pi == "9") "0.90" else "0.75"

point_cells <- function(dist, n, J, pi, metric) {
  unlist(lapply(methods, function(method) {
    y <- pick(point, method, dist, n, J, pi, metric)
    if (is.null(y)) {
      return(if (method %in% c("LBNP", "LBNP-log")) c("--", "--", "--") else c("--", "--"))
    }
    cells <- c(num_point(100 * y$rb), num_point(100 * y$mse))
    if (method %in% c("LBNP", "LBNP-log")) {
      fail <- as.integer(y$n_missing)
      cells <- c(cells, if (fail > 0L) paste0("\\textbf{", fail, "}") else "0")
    }
    cells
  }))
}
coverage_cells <- function(dist, n, J, pi, metric) {
  unlist(lapply(coverage_methods, function(method) {
    y <- pick(coverage, method, dist, n, J, pi, metric)
    if (is.null(y)) c("--", "--") else c(num_fixed(100 * y$cp, 1), num_fixed(y$al, 2))
  }))
}

make_table <- function(profile_name, profile, pi, type = c("point", "coverage")) {
  type <- match.arg(type)
  is_point <- type == "point"
  caption <- if (is_point) {
    paste0(
      "Point-estimate \\%RB and MSE for $\\pi_0=\\pi_1=", pi_caption(pi),
      "$, with MSE reported as $100\\times\\mathrm{MSE}$."
    )
  } else {
    paste0(
      "Coverage probability and average confidence-interval length for ",
      "$\\pi_0=\\pi_1=", pi_caption(pi), "$."
    )
  }
  tabular <- if (is_point) "llc|*{10}{c}|*{10}{c}" else "llc|*{10}{c}|*{10}{c}"
  row_end <- "\\\\"
  headers <- if (is_point) c(
    paste0("& & & \\multicolumn{10}{c}{$J_0=0.3$} & \\multicolumn{10}{c}{$J_0=0.5$} ", row_end),
    "\\cline{4-13}\\cline{14-23}",
    paste0("Distribution & $n_0=n_1$ & & \\multicolumn{2}{c}{Ours} & \\multicolumn{2}{c}{NP} & \\multicolumn{3}{c}{LBNP} & \\multicolumn{3}{c}{LBNP-log} & \\multicolumn{2}{c}{Ours} & \\multicolumn{2}{c}{NP} & \\multicolumn{3}{c}{LBNP} & \\multicolumn{3}{c}{LBNP-log} ", row_end),
    "\\cline{4-5}\\cline{6-7}\\cline{8-10}\\cline{11-13}\\cline{14-15}\\cline{16-17}\\cline{18-20}\\cline{21-23}",
    paste0("& & & \\%RB & MSE & \\%RB & MSE & \\%RB & MSE & Fail & \\%RB & MSE & Fail & \\%RB & MSE & \\%RB & MSE & \\%RB & MSE & Fail & \\%RB & MSE & Fail ", row_end))
  else c(
    paste0("& & & \\multicolumn{10}{c}{$J_0=0.3$} & \\multicolumn{10}{c}{$J_0=0.5$} ", row_end),
    "\\cline{4-13}\\cline{14-23}",
    paste0("Distribution & $n_0=n_1$ & & ", paste(rep(c(
      "\\multicolumn{2}{c}{Ours}", "\\multicolumn{2}{c}{Ours\\_fixk}",
      "\\multicolumn{2}{c}{NP}", "\\multicolumn{2}{c}{LBNP}",
      "\\multicolumn{2}{c}{LBNP-log}"), 2), collapse = " & "), " ", row_end),
    paste0(paste0("\\cline{", seq(4, 22, 2), "-", seq(5, 23, 2), "}", collapse = "")),
    paste0("& & & ", paste(rep("CP (\\%) & AL", 10), collapse = " & "), " ", row_end))
  lines <- c(
    "% Auto-generated from result/ by simulation/make_tables.R",
    "\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{tab:qe_", type, "_", profile_name, "_pi", pi, "}"),
    "\\small", "\\setlength{\\tabcolsep}{2pt}", "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{", tabular, "}"), "\\toprule", headers, "\\midrule"
  )
  for (d_i in seq_along(profile$dists)) {
    dist <- profile$dists[[d_i]]
    for (n_i in seq_along(n_values)) {
      n <- n_values[[n_i]]
      for (m_i in seq_along(metric_order)) {
        metric <- metric_order[[m_i]]
        cells <- if (is_point) {
          c(point_cells(dist, n, "3", pi, metric), point_cells(dist, n, "5", pi, metric))
        } else {
          c(coverage_cells(dist, n, "3", pi, metric), coverage_cells(dist, n, "5", pi, metric))
        }
        dist_cell <- if (n_i == 1L && m_i == 1L) {
          paste0("\\multirow{", length(n_values) * length(metric_order), "}{*}{", dist_tex[[dist]], "}")
        } else ""
        n_cell <- if (m_i == 1L) paste0("\\multirow{", length(metric_order), "}{*}{", n, "}") else ""
        lines <- c(lines, paste0(paste(c(dist_cell, n_cell, metric_tex[[metric]], cells), collapse = " & "), " ", row_end))
      }
      if (n_i < length(n_values)) lines <- c(lines, "\\cline{2-23}")
    }
    if (d_i < length(profile$dists)) lines <- c(lines, "\\midrule")
  }
  c(lines, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table}")
}

tables <- list()
for (profile_name in names(profiles)) {
  profile <- profiles[[profile_name]]
  tables[[paste0("point_", profile_name, ".tex")]] <- unlist(
    lapply(profile$pis, function(pi) make_table(profile_name, profile, pi, "point")),
    use.names = FALSE
  )
  tables[[paste0("coverage_", profile_name, ".tex")]] <- unlist(
    lapply(profile$pis, function(pi) make_table(profile_name, profile, pi, "coverage")),
    use.names = FALSE
  )
}

missing_lines <- lapply(tables, function(lines) lines[grepl("--", lines, fixed = TRUE)])
missing_lines <- missing_lines[lengths(missing_lines) > 0L]
if (length(missing_lines) > 0L) {
  details <- unlist(Map(
    function(name, lines) paste0(name, ": ", lines),
    names(missing_lines), missing_lines
  ), use.names = FALSE)
  stop(
    "At least one table contains a missing-value placeholder; no files were replaced:\n",
    paste(details, collapse = "\n")
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(point, file.path(output_dir, "point_summary.csv"))
data.table::fwrite(coverage, file.path(output_dir, "coverage_summary.csv"))
for (name in names(tables)) writeLines(tables[[name]], file.path(output_dir, name))
message("Regenerated six LaTeX tables in ", normalizePath(output_dir))
