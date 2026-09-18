#!/usr/bin/env Rscript

# Apply the EMBC goodness-of-fit tests to one malaria-data transformation.

usage <- function() {
  cat(
    "Usage: Rscript gof/real_data_gof.R --variant NAME [options]\n\n",
    "Variants: original, div100000, div200000\n\n",
    "Options:\n",
    "  --variant NAME   Required data transformation.\n",
    "  --B N            Successful bootstrap fits. Default: 500.\n",
    "  --cores N        Parallel bootstrap workers. Default: 32.\n",
    "  --alpha X        Test level. Default: 0.05.\n",
    "  --seed N         Base seed. Default: 20260818.\n",
    "  --pi0 X          Group-0 purity. Default: 1.\n",
    "  --pi1 X          Group-1 purity. Default: 0.677.\n",
    "  --repo-root DIR  Repository root.\n",
    "  --out-dir DIR    Output root. Default: gof/real_data_result.\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs()[grep("^--file=", commandArgs())][1]
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

settings <- list(
  variant = NULL,
  B = 500L,
  cores = 32L,
  alpha = 0.05,
  seed = 20260818L,
  pi0 = 1,
  pi1 = 0.677,
  repo_root = normalizePath(file.path(script_dir, ".."), mustWork = FALSE),
  out_dir = file.path(script_dir, "real_data_result")
)

i <- 1L
while (i <= length(args)) {
  arg <- args[i]
  if (arg %in% c("--help", "-h")) {
    usage()
    quit(status = 0)
  }
  if (i == length(args)) stop(arg, " requires a value.")
  value <- args[i + 1L]
  if (arg == "--variant") settings$variant <- value
  else if (arg == "--B") settings$B <- as.integer(value)
  else if (arg == "--cores") settings$cores <- as.integer(value)
  else if (arg == "--alpha") settings$alpha <- as.numeric(value)
  else if (arg == "--seed") settings$seed <- as.integer(value)
  else if (arg == "--pi0") settings$pi0 <- as.numeric(value)
  else if (arg == "--pi1") settings$pi1 <- as.numeric(value)
  else if (arg == "--repo-root") settings$repo_root <- value
  else if (arg == "--out-dir") settings$out_dir <- value
  else stop("Unknown option: ", arg)
  i <- i + 2L
}

variants <- c("original", "div100000", "div200000")
if (is.null(settings$variant) || !settings$variant %in% variants) {
  stop("--variant must be one of: ", paste(variants, collapse = ", "), ".")
}
if (is.na(settings$B) || settings$B < 2L) stop("--B must be at least 2.")
if (is.na(settings$cores) || settings$cores < 1L) {
  stop("--cores must be positive.")
}

repo_root <- normalizePath(settings$repo_root)
source(file.path(repo_root, "gof", "gof.r"))

raw <- read.table(
  file.path(repo_root, "real data", "malaria.txt"),
  header = FALSE,
  col.names = c("R", "T")
)
data <- raw[raw$T != 0, , drop = FALSE]

divisor <- switch(
  settings$variant,
  original = 1,
  div100000 = 100000,
  div200000 = 200000
)
data$T <- data$T / divisor
if (any(!is.finite(data$T)) || any(data$T <= 0)) {
  stop("The transformed nonzero biomarker values must all be finite and positive.")
}

gof <- lnbc_gof_test(
  data = list(T = data$T, R = data$R),
  pi0 = settings$pi0,
  pi1 = settings$pi1,
  B = settings$B,
  alpha = settings$alpha,
  seed = settings$seed,
  cores = settings$cores,
  repo_root = repo_root,
  verbose = TRUE
)

tests <- data.frame(
  variant = settings$variant,
  divisor = divisor,
  n0 = sum(data$R == 0),
  n1 = sum(data$R == 1),
  pi0 = settings$pi0,
  pi1 = settings$pi1,
  B = settings$B,
  gof$tests,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

config <- settings
config$repo_root <- repo_root
config$divisor <- divisor
config$n0 <- sum(data$R == 0)
config$n1 <- sum(data$R == 1)

variant_dir <- file.path(settings$out_dir, settings$variant)
dir.create(variant_dir, recursive = TRUE, showWarnings = FALSE)
rds_file <- file.path(variant_dir, "gof_result.rds")
csv_file <- file.path(variant_dir, "gof_tests.csv")
rds_tmp <- paste0(rds_file, ".tmp-", Sys.getpid())
csv_tmp <- paste0(csv_file, ".tmp-", Sys.getpid())

saveRDS(list(config = config, tests = tests, gof = gof), rds_tmp)
write.csv(tests, csv_tmp, row.names = FALSE, na = "")
if (!file.rename(rds_tmp, rds_file)) stop("Could not write ", rds_file)
if (!file.rename(csv_tmp, csv_file)) stop("Could not write ", csv_file)

cat("Wrote real-data GOF results for", settings$variant, "to", variant_dir, "\n")
