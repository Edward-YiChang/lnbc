#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd())
if (!file.exists(file.path(repo_root, "bootstrap_core.R"))) {
  stop("Run this script from the code directory.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.")
}

args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else {
  file.path(repo_root, "result")
}
output_dir <- file.path(repo_root, "simulation", "summary")

methods <- c("EMBC", "NP-inc", "LBNP", "LBNP-log")
coverage_methods <- c("EMBC", "EMBC_fixk", "NP-inc", "LBNP", "LBNP-log")
dists <- c("lognorm", "weibull", "weibull2")
n_values <- c(100, 300, 500)
J_values <- c("3", "5")
pi_values <- c("75", "9")
metric_order <- c("ROC_0.20", "AUC", "J", "se", "sp")
roc_suffix <- c(
  EMBC = "roc_band.csv", `NP-inc` = "roc.csv",
  LBNP = "roc_band.csv", `LBNP-log` = "roc_band.csv"
)

result_path <- function(type, method, dist, n, J, pi) {
  stem <- paste0("n", n, "J", J, "pi", pi, dist,
                 "_final_B500_reps1000_")
  suffix <- if (type == "coverage") "ci.csv" else roc_suffix[[method]]
  path <- file.path(result_dir, type, method, paste0(stem, suffix))
  unified <- file.path(result_dir, type, method, paste0(stem, "roc.csv"))
  if (type == "ROC_plot" && file.exists(unified)) unified else path
}

required <- expand.grid(
  method = methods, dist = dists, n = n_values,
  J = J_values, pi = pi_values, stringsAsFactors = FALSE
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
if (length(missing_paths)) {
  stop("Missing required raw inputs:\n", paste(missing_paths, collapse = "\n"))
}

coverage_columns <- c(
  "method", "dist", "n", "J_key", "pi_key", "pi", "B", "total_reps",
  "rep", "metric", "estimate", "ci_length", "true", "covered"
)
roc_columns <- c(
  "method", "dist", "n", "J_key", "pi_key", "pi", "B", "total_reps",
  "rep", "fpr", "tpr", "ci_low", "ci_high", "true_tpr", "covered"
)

bad_schema <- c(
  vapply(required$coverage_path, function(path) {
    !all(coverage_columns %in% names(read.csv(path, nrows = 0L,
                                                check.names = FALSE)))
  }, logical(1)),
  vapply(required$roc_path, function(path) {
    !all(roc_columns %in% names(read.csv(path, nrows = 0L,
                                        check.names = FALSE)))
  }, logical(1))
)
if (any(bad_schema)) stop("Required raw-result columns are missing.")

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

summarize_rows <- function(x, metric_name = NULL) {
  if (!is.null(metric_name)) x[, metric := metric_name]
  point <- x[, {
    truth_value <- true[is.finite(true)][1L]
    if (!length(truth_value)) truth_value <- NA_real_
    list(
      n_reps_available = data.table::uniqueN(rep), true = truth_value,
      mean_est = safe_mean(estimate),
      rb = if (is.finite(truth_value) && truth_value != 0) {
        safe_mean((estimate - truth_value) / truth_value)
      } else safe_mean(estimate - truth_value),
      mse = safe_mean((estimate - truth_value)^2),
      n_missing = sum(!is.finite(estimate))
    )
  }, by = .(method, dist, n, J_key, pi_key, pi, metric, B, total_reps)]
  coverage <- x[, list(
    n_reps_available = data.table::uniqueN(rep),
    n_valid = sum(!is.na(covered)), n_missing = sum(is.na(covered)),
    cp = safe_mean(as.numeric(covered)), al = safe_mean(ci_length)
  ), by = .(method, dist, n, J_key, pi_key, pi, metric, B, total_reps)]
  list(point = point, coverage = coverage)
}

point_parts <- list()
coverage_parts <- list()
for (i in seq_len(nrow(required))) {
  scalar <- data.table::fread(required$coverage_path[[i]],
                              select = coverage_columns, showProgress = FALSE)
  scalar <- scalar[metric %in% c("AUC", "J", "se", "sp")]
  if (scalar[, anyDuplicated(paste(rep, metric, sep = ":"))] != 0L) {
    stop("Duplicated replication/metric rows in ", required$coverage_path[[i]])
  }
  summary <- summarize_rows(scalar)
  point_parts[[length(point_parts) + 1L]] <- summary$point
  coverage_parts[[length(coverage_parts) + 1L]] <- summary$coverage

  roc <- data.table::fread(required$roc_path[[i]],
                           select = roc_columns, showProgress = FALSE)
  roc <- roc[abs(fpr - 0.2) < 1e-12]
  if (!nrow(roc)) stop("ROC(0.2) is missing from ", required$roc_path[[i]])
  if (roc[, anyDuplicated(rep)] != 0L) {
    stop("Duplicated ROC(0.2) replication rows in ", required$roc_path[[i]])
  }
  roc[, `:=`(estimate = tpr, true = true_tpr, ci_length = ci_high - ci_low)]
  summary <- summarize_rows(roc, "ROC_0.20")
  point_parts[[length(point_parts) + 1L]] <- summary$point
  coverage_parts[[length(coverage_parts) + 1L]] <- summary$coverage
}

point <- data.table::rbindlist(point_parts)
coverage <- data.table::rbindlist(coverage_parts)
expected_rows <- length(methods) * length(dists) * length(n_values) *
  length(J_values) * length(pi_values) * length(metric_order)
if (nrow(point) != expected_rows || nrow(coverage) != expected_rows) {
  stop("Unexpected number of summary rows.")
}

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
if (nrow(coverage) != expected_rows / length(methods) * length(coverage_methods)) {
  stop("Unexpected number of coverage summary rows.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(point, file.path(output_dir, "point_summary.csv"))
data.table::fwrite(coverage, file.path(output_dir, "coverage_summary.csv"))
message("Wrote numeric simulation summaries to ", normalizePath(output_dir))
