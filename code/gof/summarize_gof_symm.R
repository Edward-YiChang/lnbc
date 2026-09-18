#!/usr/bin/env Rscript

# Combine symmetric-sample GOF replication files into result_symm/summary.

usage <- function() {
  cat(
    "Usage: Rscript gof/summarize_gof_symm.R [--out-dir DIR]\n\n",
    "Options:\n",
    "  --out-dir DIR  Directory containing scenario/rep_####.csv files.\n",
    "                 Default: gof/result_symm relative to this script.\n",
    "  --help         Show this help.\n",
    sep = ""
  )
}

script_arg <- commandArgs()[grep("^--file=", commandArgs())][1L]
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
out_dir <- file.path(script_dir, "result_symm")

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  arg <- args[i]
  if (arg %in% c("--help", "-h")) {
    usage()
    quit(status = 0)
  }
  if (arg == "--out-dir") {
    if (i == length(args)) stop("--out-dir requires a directory.")
    out_dir <- args[i + 1L]
    i <- i + 2L
    next
  }
  stop("Unknown option: ", arg, "\n\n", capture.output(usage()))
}

out_dir <- normalizePath(out_dir, mustWork = FALSE)
if (!dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir)

result_files <- list.files(
  out_dir,
  pattern = "^rep_[0-9]+[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(result_files) == 0L) {
  stop("No replication files named rep_####.csv found under ", out_dir)
}

required <- c(
  "config_id", "scenario", "rep", "distribution", "n", "n0", "n1",
  "pi", "pi0", "pi1", "J", "a0", "a1", "b0", "b1", "sdlog0",
  "sdlog1", "test", "statistic", "p_value", "reject", "n_boot_ok",
  "error"
)
read_result <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop("Missing column(s) in ", path, ": ", paste(missing, collapse = ", "))
  }
  x$source_file <- path
  x
}

results <- do.call(rbind, lapply(result_files, read_result))
row.names(results) <- NULL
results$scenario <- as.character(results$scenario)
results$test <- as.character(results$test)
results$rep <- as.integer(results$rep)
results$reject <- as.logical(results$reject)
results$error <- as.character(results$error)

row_key <- paste(
  results$config_id, results$scenario, results$rep, results$test, sep = "|"
)
if (anyDuplicated(row_key)) {
  stop("Duplicate configuration/scenario/replication/test rows were found.")
}

group_key <- interaction(
  results$config_id, results$scenario, results$test,
  drop = TRUE, lex.order = TRUE
)
groups <- split(results, group_key)
gof_summary <- do.call(rbind, lapply(groups, function(x) {
  successful <- !is.na(x$reject)
  data.frame(
    config_id = x$config_id[1L],
    distribution = x$distribution[1L],
    n = x$n[1L],
    n0 = x$n0[1L],
    n1 = x$n1[1L],
    pi = x$pi[1L],
    pi0 = x$pi0[1L],
    pi1 = x$pi1[1L],
    J = x$J[1L],
    a0 = x$a0[1L],
    a1 = x$a1[1L],
    b0 = x$b0[1L],
    b1 = x$b1[1L],
    sdlog0 = x$sdlog0[1L],
    sdlog1 = x$sdlog1[1L],
    scenario = x$scenario[1L],
    test = x$test[1L],
    reps = length(unique(x$rep)),
    successful_reps = sum(successful),
    failed_reps = sum(!successful),
    rejection_rate = if (any(successful)) mean(x$reject[successful]) else NA_real_,
    mean_p_value = if (any(successful)) mean(x$p_value[successful]) else NA_real_,
    mean_statistic = if (any(successful)) mean(x$statistic[successful]) else NA_real_,
    mean_bootstrap_successes = if (any(successful)) {
      mean(x$n_boot_ok[successful])
    } else {
      NA_real_
    },
    row.names = NULL
  )
}))
row.names(gof_summary) <- NULL
gof_summary <- gof_summary[
  order(gof_summary$n, -gof_summary$pi, gof_summary$J,
        gof_summary$scenario, gof_summary$test),
  , drop = FALSE
]

failures <- results[
  !is.na(results$error) & nzchar(results$error),
  , drop = FALSE
]

summary_dir <- file.path(out_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
results_file <- file.path(summary_dir, "gof_results.csv")
summary_file <- file.path(summary_dir, "gof_summary.csv")
failures_file <- file.path(summary_dir, "gof_failures.csv")
write.csv(results, results_file, row.names = FALSE, na = "")
write.csv(gof_summary, summary_file, row.names = FALSE, na = "")
write.csv(failures, failures_file, row.names = FALSE, na = "")

cat("Read", length(result_files), "replication file(s) from", out_dir, "\n\n")
print(gof_summary, row.names = FALSE)
cat(
  "\nWrote:\n", results_file, "\n", summary_file, "\n", failures_file,
  "\n", sep = ""
)

