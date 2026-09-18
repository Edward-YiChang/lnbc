# Simulation study for the GOF tests. This file only defines functions; it does
# not run the study when sourced.

.gof_test_source_dir <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile)))
    }
  }
  getwd()
}

.gof_test_dir <- .gof_test_source_dir()
if (!exists("lnbc_gof_test", mode = "function")) {
  source(file.path(.gof_test_dir, "gof.r"))
}

# Parameters from the kappa = 0 lognormal settings used in the GOF study.
# The heterogeneous alternative changes only b1 from 1 to 2.
lnbc_gof_lognormal_parameters <- function(
    J = c(0.3, 0.5), scenario = c("homogeneous", "heterogeneous")) {
  scenario <- match.arg(scenario)
  if (length(J) != 1 || !is.finite(J)) {
    stop("J must be one finite value.")
  }
  J_key <- format(round(J, 10), trim = TRUE, scientific = FALSE)
  a1 <- switch(
    J_key,
    "0.3" = 0.77,
    "0.5" = 1.35,
    stop("Unsupported J value: ", J, ". Choose 0.3 or 0.5.")
  )
  data.frame(
    J = J,
    AUC = if (J_key == "0.3") 0.71 else 0.83,
    eta_star0 = if (J_key == "0.3") 0.65 else 0.75,
    tau_star0 = if (J_key == "0.3") 0.65 else 0.75,
    a0 = 0,
    b0 = 1,
    a1 = a1,
    b1 = if (scenario == "heterogeneous") 2 else 1,
    row.names = NULL
  )
}

lnbc_gof_normal_data <- function(n0, n1, pi0, pi1, mu0, mu1, sd0, sd1,
                                 seed = NULL) {
  n0 <- as.integer(n0)
  n1 <- as.integer(n1)
  if (length(n0) != 1 || length(n1) != 1 || any(is.na(c(n0, n1))) ||
      n0 < 1 || n1 < 1) {
    stop("n0 and n1 must be positive integers.")
  }
  .gof_check_probability(pi0, "pi0")
  .gof_check_probability(pi1, "pi1")
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }
  if (length(mu0) != 1 || length(mu1) != 1 ||
      length(sd0) != 1 || length(sd1) != 1 ||
      !all(is.finite(c(mu0, mu1, sd0, sd1))) || sd0 <= 0 || sd1 <= 0) {
    stop("Normal means must be finite and standard deviations must be positive.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  g0 <- rbinom(n0, size = 1, prob = 1 - pi0)
  g1 <- rbinom(n1, size = 1, prob = pi1)
  G <- c(g0, g1)
  means <- ifelse(G == 0, mu0, mu1)
  sds <- ifelse(G == 0, sd0, sd1)
  log_survival_at_zero <- pnorm(
    0, mean = means, sd = sds, lower.tail = FALSE, log.p = TRUE
  )
  u <- pmax(runif(n0 + n1), .Machine$double.xmin)
  T <- qnorm(
    log_survival_at_zero + log(u),
    mean = means, sd = sds, lower.tail = FALSE, log.p = TRUE
  )

  list(
    T = T,
    R = c(rep.int(0L, n0), rep.int(1L, n1)),
    G = G
  )
}

lnbc_gof_lognormal_data <- function(n0, n1, pi0, pi1,
                                    meanlog0, meanlog1, sdlog0, sdlog1,
                                    seed = NULL) {
  n0 <- as.integer(n0)
  n1 <- as.integer(n1)
  if (length(n0) != 1 || length(n1) != 1 || any(is.na(c(n0, n1))) ||
      n0 < 1 || n1 < 1) {
    stop("n0 and n1 must be positive integers.")
  }
  .gof_check_probability(pi0, "pi0")
  .gof_check_probability(pi1, "pi1")
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }
  if (length(meanlog0) != 1 || length(meanlog1) != 1 ||
      length(sdlog0) != 1 || length(sdlog1) != 1 ||
      !all(is.finite(c(meanlog0, meanlog1, sdlog0, sdlog1))) ||
      sdlog0 <= 0 || sdlog1 <= 0) {
    stop("Lognormal meanlog values must be finite and sdlog values must be positive.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  g0 <- rbinom(n0, size = 1, prob = 1 - pi0)
  g1 <- rbinom(n1, size = 1, prob = pi1)
  G <- c(g0, g1)
  meanlog <- ifelse(G == 0, meanlog0, meanlog1)
  sdlog <- ifelse(G == 0, sdlog0, sdlog1)
  T <- rlnorm(n0 + n1, meanlog = meanlog, sdlog = sdlog)

  list(
    T = T,
    R = c(rep.int(0L, n0), rep.int(1L, n1)),
    G = G
  )
}

lnbc_gof_weibull_data <- function(n0, n1, pi0, pi1, a0, a1, b0, b1,
                                  seed = NULL) {
  if (length(n0) != 1L || length(n1) != 1L ||
      any(!is.finite(c(n0, n1))) || any(c(n0, n1) < 1) ||
      any(c(n0, n1) != floor(c(n0, n1)))) {
    stop("n0 and n1 must be positive integers.")
  }
  .gof_check_probability(pi0, "pi0")
  .gof_check_probability(pi1, "pi1")
  if (pi0 + pi1 <= 1) stop("pi0 + pi1 must be greater than 1.")
  parameters <- list(a0, a1, b0, b1)
  if (!all(vapply(parameters, function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }, logical(1)))) {
    stop("Weibull a0, a1, b0, and b1 must be positive finite scalars.")
  }
  if (!is.null(seed)) set.seed(seed)
  G <- c(rbinom(n0, 1, 1 - pi0), rbinom(n1, 1, pi1))
  T <- rweibull(n0 + n1, shape = ifelse(G == 0, b0, b1),
               scale = ifelse(G == 0, a0, a1))
  list(T = T, R = c(rep.int(0L, n0), rep.int(1L, n1)), G = G)
}

lnbc_run_gof_normal_scenario <- function(
    scenario = c("homogeneous", "heterogeneous"),
    reps = 1000, B = 499, n0 = 100, n1 = 300,
    pi0 = 0.9, pi1 = 0.9, mu0 = 10, mu1 = 12,
    homogeneous_sd = 1, heterogeneous_sd = c(1, 2),
    alpha = 0.05, seed = 20260818, fit_control = list(),
    cores = 1,
    repo_root = normalizePath(file.path(.gof_test_dir, "..")),
    verbose = FALSE, distribution = c("normal", "lognormal", "weibull"),
    weibull_parameters = NULL) {
  scenario <- match.arg(scenario)
  distribution <- match.arg(distribution)
  if (distribution == "weibull") {
    if (!is.list(weibull_parameters) || length(weibull_parameters) != 4L ||
        !setequal(names(weibull_parameters), c("a0", "a1", "b0", "b1")) ||
        !all(vapply(weibull_parameters, function(x) {
          is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
        }, logical(1)))) {
      stop("weibull_parameters must be a named list of positive finite scalars a0, a1, b0, b1.")
    }
  } else if (!is.null(weibull_parameters)) {
    stop("weibull_parameters is only supported for distribution = 'weibull'.")
  }
  reps <- as.integer(reps)
  if (length(reps) != 1 || is.na(reps) || reps < 1) {
    stop("reps must be a positive integer.")
  }
  if (length(homogeneous_sd) != 1 || !is.finite(homogeneous_sd) ||
      homogeneous_sd <= 0) {
    stop("homogeneous_sd must be one positive value.")
  }
  if (length(heterogeneous_sd) != 2 || any(!is.finite(heterogeneous_sd)) ||
      any(heterogeneous_sd <= 0)) {
    stop("heterogeneous_sd must contain two positive values.")
  }

  if (scenario == "homogeneous") {
    sd0 <- homogeneous_sd
    sd1 <- homogeneous_sd
  } else {
    sd0 <- heterogeneous_sd[1]
    sd1 <- heterogeneous_sd[2]
  }

  set.seed(seed)
  data_seeds <- sample.int(.Machine$integer.max, reps)
  gof_seeds <- sample.int(.Machine$integer.max, reps)
  test_names <- c("sup1_quantile", "sup2_quantile", "sup3_quantile",
                  "residual_quantile", "residual_pvalue")
  results <- vector("list", reps)

  for (rep_id in seq_len(reps)) {
    one <- tryCatch({
      data <- if (distribution == "normal") {
        lnbc_gof_normal_data(
          n0 = n0, n1 = n1, pi0 = pi0, pi1 = pi1,
          mu0 = mu0, mu1 = mu1, sd0 = sd0, sd1 = sd1,
          seed = data_seeds[rep_id]
        )
      } else if (distribution == "weibull") {
        do.call(lnbc_gof_weibull_data, c(
          list(n0 = n0, n1 = n1, pi0 = pi0, pi1 = pi1,
               seed = data_seeds[rep_id]), weibull_parameters
        ))
      } else {
        lnbc_gof_lognormal_data(
          n0 = n0, n1 = n1, pi0 = pi0, pi1 = pi1,
          meanlog0 = mu0, meanlog1 = mu1,
          sdlog0 = sd0, sdlog1 = sd1,
          seed = data_seeds[rep_id]
        )
      }
      gof <- lnbc_gof_test(
        data = data, pi0 = pi0, pi1 = pi1, B = B, alpha = alpha,
        seed = gof_seeds[rep_id], fit_control = fit_control,
        cores = cores, repo_root = repo_root, verbose = verbose
      )
      out <- gof$tests
      out$n_boot_ok <- gof$n_boot_ok
      out$error <- NA_character_
      out
    }, error = function(e) {
      data.frame(
        test = test_names,
        statistic = NA_real_,
        critical_low = NA_real_,
        critical_high = NA_real_,
        p_value = NA_real_,
        reject = NA,
        n_boot_ok = 0L,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
    one$scenario <- scenario
    one$rep <- rep_id
    one$data_seed <- data_seeds[rep_id]
    one$gof_seed <- gof_seeds[rep_id]
    results[[rep_id]] <- one

    if (isTRUE(verbose) && (rep_id == 1 || rep_id %% 10 == 0 || rep_id == reps)) {
      cat("GOF", scenario, "replication", rep_id, "of", reps, "\n")
    }
  }

  out <- do.call(rbind, results)
  row.names(out) <- NULL
  out[, c("scenario", "rep", "test", "statistic", "critical_low",
          "critical_high", "p_value", "reject", "n_boot_ok", "error",
          "data_seed", "gof_seed")]
}

.gof_simulation_summary <- function(results) {
  groups <- split(results, interaction(results$scenario, results$test, drop = TRUE))
  out <- lapply(groups, function(x) {
    successful <- !is.na(x$reject)
    rejection_rate <- if (any(successful)) mean(x$reject[successful]) else NA_real_
    data.frame(
      scenario = x$scenario[1],
      test = x$test[1],
      reps = nrow(x),
      successful_reps = sum(successful),
      rejection_rate = rejection_rate,
      mean_bootstrap_successes = if (any(successful)) {
        mean(x$n_boot_ok[successful])
      } else {
        NA_real_
      },
      row.names = NULL
    )
  })
  summary <- do.call(rbind, out)
  row.names(summary) <- NULL
  summary[order(summary$scenario, summary$test), ]
}

lnbc_run_gof_study <- function(
    reps = 1000, B = 499, n0 = 100, n1 = 300,
    pi0 = 0.9, pi1 = 0.9, mu0 = 10, mu1 = 12,
    homogeneous_sd = 1, heterogeneous_sd = c(1, 2),
    alpha = 0.05, seed = 20260818, fit_control = list(),
    cores = 1,
    repo_root = normalizePath(file.path(.gof_test_dir, "..")),
    verbose = FALSE, distribution = c("normal", "lognormal")) {
  distribution <- match.arg(distribution)
  common <- list(
    reps = reps, B = B, n0 = n0, n1 = n1,
    pi0 = pi0, pi1 = pi1, mu0 = mu0, mu1 = mu1,
    homogeneous_sd = homogeneous_sd, heterogeneous_sd = heterogeneous_sd,
    alpha = alpha, fit_control = fit_control, repo_root = repo_root,
    cores = cores, verbose = verbose, distribution = distribution
  )
  homogeneous <- do.call(
    lnbc_run_gof_normal_scenario,
    c(list(scenario = "homogeneous", seed = seed), common)
  )
  heterogeneous <- do.call(
    lnbc_run_gof_normal_scenario,
    c(list(scenario = "heterogeneous", seed = seed + 1), common)
  )
  results <- rbind(homogeneous, heterogeneous)
  row.names(results) <- NULL

  list(
    results = results,
    summary = .gof_simulation_summary(results),
    settings = list(
      reps = reps, B = B, n0 = n0, n1 = n1,
      pi0 = pi0, pi1 = pi1, mu0 = mu0, mu1 = mu1,
      homogeneous_sd = homogeneous_sd,
      heterogeneous_sd = heterogeneous_sd,
      alpha = alpha, seed = seed, cores = cores, distribution = distribution
    )
  )
}

# Example (intentionally not run):
# study <- lnbc_run_gof_study(reps = 1000, B = 499)
