#!/usr/bin/env Rscript

# Run one outer replication of the symmetric-sample GOF study.
#
# Each replication covers
#   n0 = n1 = 100, 300, 500;
#   pi0 = pi1 = 0.90, 0.75;
#   J = 0.3, 0.5; and
#   homo/heterogeneous lognormal distributions with log-variances (1, 1)
#   and (1, 4), respectively.

usage <- function() {
  cat(
    "Usage: Rscript gof/gof_symm.R --rep-id N [options]\n\n",
    "Runs one outer replication over the complete symmetric GOF grid.\n",
    "When --rep-id is omitted, SLURM_ARRAY_TASK_ID is used.\n\n",
    "Options:\n",
    "  --rep-id N       Outer replication number.\n",
    "  --B N            Successful bootstrap fits per test. Default: 500.\n",
    "  --cores N        Parallel bootstrap workers. Default: SLURM allocation,\n",
    "                   or 1 outside SLURM.\n",
    "  --alpha X        Test level. Default: 0.05.\n",
    "  --seed N         Base simulation seed. Default: 20260818.\n",
    "  --n X,Y,...     Symmetric sample sizes from 100,300,500.\n",
    "                   Default: 100,300,500.\n",
    "  --pi X,Y,...    Equal pi0=pi1 values from 0.9,0.75.\n",
    "                   Default: 0.9,0.75.\n",
    "  --J X,Y,...     J values from 0.3,0.5. Default: 0.3,0.5.\n",
    "  --scenario S    both, homo, or hete. Default: both.\n",
    "  --repo-root DIR  Repository root. Default: parent of this script.\n",
    "  --out-dir DIR    Output root. Default: gof/result_symm.\n",
    "  --help           Show this help.\n",
    sep = ""
  )
}

script_arg <- commandArgs()[grep("^--file=", commandArgs())][1L]
script_path <- sub("^--file=", "", script_arg)
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
default_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (is.na(default_cores) || default_cores < 1L) default_cores <- 1L

settings <- list(
  rep_id = suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))),
  B = 500L,
  cores = default_cores,
  alpha = 0.05,
  seed = 20260818,
  n = c(100L, 300L, 500L),
  pi = c(0.90, 0.75),
  J = c(0.3, 0.5),
  scenario = "both",
  repo_root = file.path(script_dir, ".."),
  out_dir = file.path(script_dir, "result_symm")
)

parse_values <- function(value, type, option) {
  values <- strsplit(value, ",", fixed = TRUE)[[1L]]
  numeric_values <- suppressWarnings(as.numeric(values))
  if (length(numeric_values) == 0L || any(!is.finite(numeric_values))) {
    stop(option, " requires a comma-separated list of numbers.")
  }
  if (type == "integer" && any(numeric_values != floor(numeric_values))) {
    stop(option, " requires whole numbers.")
  }
  parsed <- if (type == "integer") as.integer(numeric_values) else numeric_values
  unique(parsed)
}

args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args)) {
  arg <- args[i]
  if (arg %in% c("--help", "-h")) {
    usage()
    quit(status = 0)
  }
  if (i == length(args)) stop(arg, " requires a value.")
  value <- args[i + 1L]
  if (arg == "--rep-id") settings$rep_id <- as.integer(value)
  else if (arg == "--B") settings$B <- as.integer(value)
  else if (arg == "--cores") settings$cores <- as.integer(value)
  else if (arg == "--alpha") settings$alpha <- as.numeric(value)
  else if (arg == "--seed") settings$seed <- as.numeric(value)
  else if (arg %in% c("--n", "--n-values")) {
    settings$n <- parse_values(value, "integer", arg)
  }
  else if (arg %in% c("--pi", "--pi-values")) {
    settings$pi <- parse_values(value, "numeric", arg)
  }
  else if (arg %in% c("--J", "--J-values")) {
    settings$J <- parse_values(value, "numeric", arg)
  }
  else if (arg == "--scenario") settings$scenario <- value
  else if (arg == "--repo-root") settings$repo_root <- value
  else if (arg == "--out-dir") settings$out_dir <- value
  else stop("Unknown option: ", arg, "\n\n", capture.output(usage()))
  i <- i + 2L
}

if (length(settings$rep_id) != 1L || is.na(settings$rep_id) ||
    settings$rep_id < 1L) {
  stop("Provide a positive --rep-id or set SLURM_ARRAY_TASK_ID.")
}
if (length(settings$B) != 1L || is.na(settings$B) || settings$B < 2L) {
  stop("--B must be an integer of at least 2.")
}
if (length(settings$cores) != 1L || is.na(settings$cores) ||
    settings$cores < 1L) {
  stop("--cores must be a positive integer.")
}
if (length(settings$alpha) != 1L || !is.finite(settings$alpha) ||
    settings$alpha <= 0 || settings$alpha >= 1) {
  stop("--alpha must be in (0, 1).")
}
if (length(settings$seed) != 1L || !is.finite(settings$seed)) {
  stop("--seed must be finite.")
}
if (any(!settings$n %in% c(100L, 300L, 500L))) {
  stop("--n supports only 100, 300, and 500.")
}
if (any(!settings$pi %in% c(0.90, 0.75))) {
  stop("--pi supports only 0.9 and 0.75.")
}
if (any(!settings$J %in% c(0.3, 0.5))) {
  stop("--J supports only 0.3 and 0.5.")
}
if (!settings$scenario %in% c("both", "homo", "hete")) {
  stop("--scenario must be both, homo, or hete.")
}

repo_root <- normalizePath(settings$repo_root, mustWork = TRUE)
out_dir <- normalizePath(settings$out_dir, mustWork = FALSE)
test_file <- file.path(repo_root, "gof", "test.r")
if (!file.exists(test_file)) stop("Cannot find ", test_file, ".")
source(test_file)

n_grid <- settings$n
pi_grid <- settings$pi
J_grid <- settings$J
scenarios <- switch(
  settings$scenario,
  both = c("homo", "hete"),
  homo = "homo",
  hete = "hete"
)
distribution <- "lognormal"
setting_count <- length(n_grid) * length(pi_grid) * length(J_grid) *
  length(scenarios)
cat(
  "Running", setting_count, "setting(s) for replication", settings$rep_id,
  "\n",
  "n:", paste(n_grid, collapse = ","),
  "| pi:", paste(pi_grid, collapse = ","),
  "| J:", paste(J_grid, collapse = ","),
  "| scenario:", paste(scenarios, collapse = ","), "\n"
)

for (n in n_grid) {
  for (pi in pi_grid) {
    for (J in J_grid) {
      for (scenario in scenarios) {
        long_scenario <- if (scenario == "homo") {
          "homogeneous"
        } else {
          "heterogeneous"
        }
        log_variance0 <- 1
        log_variance1 <- if (scenario == "homo") 1 else 4
        parameters <- lnbc_gof_lognormal_parameters(J, long_scenario)

        scenario_offset <- if (scenario == "hete") 500000000 else 0
        config_offset <- 1000000 * n + 10 * round(100 * pi) + round(10 * J)
        rep_seed <- as.integer(
          (settings$seed + scenario_offset + config_offset +
             1009 * settings$rep_id) %% .Machine$integer.max
        )
        if (rep_seed < 1L) rep_seed <- 1L

        result <- lnbc_run_gof_normal_scenario(
          scenario = long_scenario,
          distribution = distribution,
          reps = 1L,
          B = settings$B,
          n0 = n,
          n1 = n,
          pi0 = pi,
          pi1 = pi,
          mu0 = parameters$a0,
          mu1 = parameters$a1,
          homogeneous_sd = sqrt(log_variance0),
          heterogeneous_sd = sqrt(c(log_variance0, log_variance1)),
          alpha = settings$alpha,
          seed = rep_seed,
          cores = settings$cores,
          repo_root = repo_root,
          verbose = TRUE
        )

        config_id <- sprintf("n%d_pi%s_J%s", n, pi, J)
        result$scenario <- scenario
        result$rep <- settings$rep_id
        result$distribution <- distribution
        result$n <- n
        result$n0 <- n
        result$n1 <- n
        result$pi <- pi
        result$pi0 <- pi
        result$pi1 <- pi
        result$J <- J
        result$a0 <- parameters$a0
        result$a1 <- parameters$a1
        result$b0 <- log_variance0
        result$b1 <- log_variance1
        result$sdlog0 <- sqrt(log_variance0)
        result$sdlog1 <- sqrt(log_variance1)
        result$config_id <- config_id

        scenario_dir <- file.path(out_dir, config_id, scenario)
        dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
        out_file <- file.path(
          scenario_dir, sprintf("rep_%04d.csv", settings$rep_id)
        )
        tmp_file <- paste0(out_file, ".tmp-", Sys.getpid())
        write.csv(result, tmp_file, row.names = FALSE, na = "")
        if (!file.rename(tmp_file, out_file)) {
          stop("Could not move completed result to ", out_file)
        }
        cat("Wrote", out_file, "\n")
      }
    }
  }
}
