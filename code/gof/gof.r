# Goodness-of-fit tests for the EMBC Box-Cox density-ratio model.

.gof_source_dir <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile)))
    }
  }
  getwd()
}

.gof_dir <- .gof_source_dir()
.gof_embc_env <- new.env(parent = globalenv())

.gof_check_probability <- function(x, name) {
  if (length(x) != 1 || !is.finite(x) || x <= 0 || x > 1) {
    stop(name, " must be a number in (0, 1].")
  }
}

.gof_check_data <- function(data) {
  if (!is.list(data) || is.null(data$T) || is.null(data$R)) {
    stop("data must be a list with T and R.")
  }
  if (length(data$T) != length(data$R)) {
    stop("T and R must have the same length.")
  }
  if (!all(data$R %in% c(0, 1))) {
    stop("R must contain only 0/1 labels.")
  }
  if (sum(data$R == 0) == 0 || sum(data$R == 1) == 0) {
    stop("R must contain at least one 0 and one 1.")
  }
  if (!all(is.finite(data$T)) || any(data$T <= 0)) {
    stop("T must contain only finite positive values for the Box-Cox model.")
  }
  list(T = as.numeric(data$T), R = as.integer(data$R))
}

.gof_get_embc_fit_fun <- function(repo_root) {
  if (exists("EMBC_fixpi", mode = "function", inherits = TRUE)) {
    return(get("EMBC_fixpi", mode = "function", inherits = TRUE))
  }

  method_file <- file.path(repo_root, "methods", "EMBC", "method.R")
  if (!file.exists(method_file)) {
    stop("Cannot find the EMBC implementation at ", method_file, ".")
  }
  if (!exists("EMBC_fixpi", envir = .gof_embc_env, mode = "function", inherits = FALSE)) {
    sys.source(method_file, envir = .gof_embc_env)
  }
  get("EMBC_fixpi", envir = .gof_embc_env, inherits = FALSE)
}

.gof_fit_model <- function(data, pi0, pi1, fit_fun, fit_control) {
  if (!is.list(fit_control)) {
    stop("fit_control must be a list.")
  }
  reserved <- intersect(names(fit_control), c("data", "pi0", "pi1"))
  if (length(reserved) > 0) {
    stop("fit_control cannot contain: ", paste(reserved, collapse = ", "), ".")
  }
  do.call(fit_fun, c(list(data = data, pi0 = pi0, pi1 = pi1), fit_control))
}

.gof_extract_fit <- function(fit, n) {
  core <- if (is.list(fit) && !is.null(fit$best)) fit$best else fit
  required <- c("alpha", "beta", "k", "pi")
  if (!is.list(core) || !all(required %in% names(core))) {
    stop("fit must contain alpha, beta, k, and empirical-likelihood weights pi.")
  }
  if (length(core$alpha) != 1 || length(core$beta) != 1 || length(core$k) != 1 ||
      !all(is.finite(c(core$alpha, core$beta, core$k)))) {
    stop("The fitted alpha, beta, and k must be finite scalars.")
  }
  if (length(core$pi) != n || any(!is.finite(core$pi)) || any(core$pi < 0) ||
      sum(core$pi) <= 0) {
    stop("The fitted empirical-likelihood weights pi are invalid.")
  }
  core
}

.gof_boxcox <- function(x, k) {
  if (abs(k) < 1e-8) {
    return(log(x))
  }
  expm1(k * log(x)) / k
}

.gof_normalize_weights <- function(w) {
  w / sum(w)
}

.gof_density_weights <- function(data, fit) {
  core <- .gof_extract_fit(fit, length(data$T))
  eta <- core$alpha + core$beta * .gof_boxcox(data$T, core$k)
  if (any(!is.finite(eta))) {
    stop("The fitted log-density ratio is not finite at every observation.")
  }

  p <- .gof_normalize_weights(core$pi)
  log_q <- log(p) + eta
  log_q_normalizer <- max(log_q) + log(sum(exp(log_q - max(log_q))))
  q <- exp(log_q - log_q_normalizer)
  list(p = p, q = q, eta = eta - log_q_normalizer)
}

.gof_logspace_add <- function(x, y) {
  largest <- pmax(x, y)
  largest + log(exp(x - largest) + exp(y - largest))
}

.gof_reference_probability <- function(data, eta, pi0, pi1) {
  n0 <- sum(data$R == 0)
  n1 <- sum(data$R == 1)
  log_f0_star <- .gof_logspace_add(log(pi0), log1p(-pi0) + eta)
  log_f1_star <- .gof_logspace_add(log1p(-pi1), log(pi1) + eta)
  plogis(log(n1) - log(n0) + log_f1_star - log_f0_star)
}

# Compute the observed GOF statistics for a supplied EMBC fit.
lnbc_gof_statistics <- function(data, fit, pi0, pi1) {
  data <- .gof_check_data(data)
  .gof_check_probability(pi0, "pi0")
  .gof_check_probability(pi1, "pi1")
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }

  weights <- .gof_density_weights(data, fit)
  r <- .gof_reference_probability(data, weights$eta, pi0, pi1)
  residual <- sum((data$R - r)^2 - r * (1 - r))

  f0_star_weights <- pi0 * weights$p + (1 - pi0) * weights$q
  f1_star_weights <- (1 - pi1) * weights$p + pi1 * weights$q
  ord <- order(data$T)
  last_at_value <- !duplicated(data$T[ord], fromLast = TRUE)
  empirical_f0_star <- cumsum(data$R[ord] == 0) / sum(data$R == 0)
  empirical_f1_star <- cumsum(data$R[ord] == 1) / sum(data$R == 1)
  fitted_f0_star <- cumsum(f0_star_weights[ord])
  fitted_f1_star <- cumsum(f1_star_weights[ord])
  D0 <- max(abs(
    empirical_f0_star[last_at_value] - fitted_f0_star[last_at_value]
  ))
  D1 <- max(abs(
    empirical_f1_star[last_at_value] - fitted_f1_star[last_at_value]
  ))
  n0 <- sum(data$R == 0)
  n1 <- sum(data$R == 1)
  n <- n0 + n1

  # Three combined supremum statistics. D0 and D1 are retained in the
  # returned vector for diagnostics, but only these combinations are tested.
  sup1 <- max(D0, D1)
  sup2 <- max(n0 * D0, n1 * D1)
  sup3 <- n0 / n * D0 + n1 / n * D1

  c(
    residual = residual,
    D0 = D0,
    D1 = D1,
    sup1 = sup1,
    sup2 = sup2,
    sup3 = sup3
  )
}

.gof_bootstrap_sample <- function(data, fit, pi0, pi1) {
  weights <- .gof_density_weights(data, fit)
  f0_star_weights <- pi0 * weights$p + (1 - pi0) * weights$q
  f1_star_weights <- (1 - pi1) * weights$p + pi1 * weights$q
  n0 <- sum(data$R == 0)
  n1 <- sum(data$R == 1)

  index0 <- sample.int(length(data$T), n0, replace = TRUE, prob = f0_star_weights)
  index1 <- sample.int(length(data$T), n1, replace = TRUE, prob = f1_star_weights)
  list(
    T = c(data$T[index0], data$T[index1]),
    R = c(rep.int(0L, n0), rep.int(1L, n1))
  )
}

# Parametric/model-based bootstrap calibration of the proposed GOF tests.
lnbc_gof_test <- function(data, pi0, pi1, B = 499, alpha = 0.05,
                          seed = NULL, fit = NULL, fit_control = list(),
                          max_bootstrap_attempts = 2 * B,
                          cores = 1,
                          repo_root = normalizePath(file.path(.gof_dir, "..")),
                          verbose = FALSE) {
  data <- .gof_check_data(data)
  .gof_check_probability(pi0, "pi0")
  .gof_check_probability(pi1, "pi1")
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }
  B <- as.integer(B)
  if (length(B) != 1 || is.na(B) || B < 2) {
    stop("B must be an integer of at least 2.")
  }
  max_bootstrap_attempts <- as.integer(max_bootstrap_attempts)
  if (length(max_bootstrap_attempts) != 1 || is.na(max_bootstrap_attempts) ||
      max_bootstrap_attempts < B) {
    stop("max_bootstrap_attempts must be an integer at least as large as B.")
  }
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number in (0, 1).")
  }
  if (!is.null(seed) && (length(seed) != 1 || !is.finite(seed))) {
    stop("seed must be NULL or a finite scalar.")
  }
  cores <- as.integer(cores)
  if (length(cores) != 1 || is.na(cores) || cores < 1) {
    stop("cores must be a positive integer.")
  }
  if (.Platform$OS.type != "unix") {
    cores <- 1L
  }
  cores <- min(cores, B)

  fit_fun <- .gof_get_embc_fit_fun(repo_root)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(fit)) {
    fit <- .gof_fit_model(data, pi0, pi1, fit_fun, fit_control)
  }
  observed <- lnbc_gof_statistics(data, fit, pi0, pi1)
  attempt_seeds <- sample.int(.Machine$integer.max, max_bootstrap_attempts)
  bootstrap_statistics <- matrix(
    NA_real_, nrow = B, ncol = length(observed),
    dimnames = list(paste0("boot_", seq_len(B)), names(observed))
  )
  attempt_errors <- rep.int(NA_character_, max_bootstrap_attempts)
  attempt_ok <- rep.int(FALSE, max_bootstrap_attempts)
  successful_seeds <- integer(B)
  attempt <- 0L
  n_boot_ok <- 0L

  run_bootstrap_attempt <- function(attempt_index) {
    set.seed(attempt_seeds[attempt_index])
    tryCatch({
      bootstrap_data <- .gof_bootstrap_sample(data, fit, pi0, pi1)
      bootstrap_fit <- .gof_fit_model(
        bootstrap_data, pi0, pi1, fit_fun, fit_control
      )
      list(
        statistics = lnbc_gof_statistics(
          bootstrap_data, bootstrap_fit, pi0, pi1
        ),
        error = NA_character_
      )
    }, error = function(e) {
      list(
        statistics = c(
          residual = NA_real_, D0 = NA_real_, D1 = NA_real_,
          sup1 = NA_real_, sup2 = NA_real_, sup3 = NA_real_
        ),
        error = conditionMessage(e)
      )
    })
  }

  while (n_boot_ok < B && attempt < max_bootstrap_attempts) {
    # Put every currently required attempt into one dynamically scheduled
    # worker pool. This keeps fast workers busy instead of waiting at a
    # cores-sized batch barrier for the slowest EMBC fit. If any attempts
    # fail, the next iteration schedules only the required replacements.
    batch_size <- min(
      B - n_boot_ok,
      max_bootstrap_attempts - attempt
    )
    batch_indices <- attempt + seq_len(batch_size)
    if (cores > 1L && batch_size > 1L) {
      results <- parallel::mclapply(
        batch_indices,
        run_bootstrap_attempt,
        mc.cores = min(cores, batch_size),
        mc.preschedule = FALSE,
        mc.set.seed = FALSE
      )
    } else {
      results <- lapply(batch_indices, run_bootstrap_attempt)
    }

    for (batch_position in seq_along(batch_indices)) {
      attempt <- batch_indices[batch_position]
      result <- results[[batch_position]]
      finite <- all(is.finite(result$statistics))
      attempt_ok[attempt] <- finite
      if (finite) {
        n_boot_ok <- n_boot_ok + 1L
        bootstrap_statistics[n_boot_ok, ] <- result$statistics
        successful_seeds[n_boot_ok] <- attempt_seeds[attempt]
      } else {
        attempt_errors[attempt] <- if (is.na(result$error)) {
          "Bootstrap fit returned non-finite GOF statistics."
        } else {
          result$error
        }
      }
    }

    if (isTRUE(verbose)) {
      cat("GOF bootstrap", n_boot_ok, "successful fits of", B,
          "after", attempt, "attempts\n")
    }
  }

  if (n_boot_ok != B) {
    stop("Only ", n_boot_ok, " of ", B,
         " required bootstrap fits succeeded in ", attempt, " attempts.")
  }

  residual_boot <- bootstrap_statistics[, "residual"]
  supremum_names <- c("sup1", "sup2", "sup3")
  supremum_critical <- vapply(supremum_names, function(name) {
    unname(quantile(
      bootstrap_statistics[, name], probs = 1 - alpha,
      type = 1, names = FALSE
    ))
  }, numeric(1))
  supremum_p <- vapply(supremum_names, function(name) {
    mean(bootstrap_statistics[, name] >= observed[name])
  }, numeric(1))

  residual_critical <- unname(quantile(
    residual_boot, probs = c(alpha / 2, 1 - alpha / 2),
    type = 1, names = FALSE
  ))
  residual_lower_p <- (1 + sum(residual_boot <= observed["residual"])) /
    (n_boot_ok + 1)
  residual_upper_p <- (1 + sum(residual_boot >= observed["residual"])) /
    (n_boot_ok + 1)
  residual_quantile_p <- min(1, 2 * min(residual_lower_p, residual_upper_p))

  residual_mean <- mean(residual_boot)
  residual_sd <- sd(residual_boot)
  if (is.finite(residual_sd) && residual_sd > 0) {
    residual_z <- observed["residual"] / residual_sd
    residual_pvalue <- 2 * pnorm(-abs(residual_z))
  } else {
    residual_z <- NA_real_
    residual_pvalue <- NA_real_
  }
  normal_critical <- qnorm(1 - alpha / 2)

  tests <- data.frame(
    test = c("sup1_quantile", "sup2_quantile", "sup3_quantile",
             "residual_quantile", "residual_pvalue"),
    statistic = c(observed[supremum_names], observed["residual"], residual_z),
    critical_low = c(rep(NA_real_, 3), residual_critical[1], -normal_critical),
    critical_high = c(supremum_critical, residual_critical[2], normal_critical),
    p_value = c(supremum_p, residual_quantile_p, residual_pvalue),
    reject = c(
      observed[supremum_names] > supremum_critical,
      observed["residual"] < residual_critical[1] ||
        observed["residual"] > residual_critical[2],
      if (is.na(residual_pvalue)) NA else residual_pvalue < alpha
    ),
    row.names = NULL
  )

  list(
    tests = tests,
    observed = observed,
    residual_bootstrap_mean = residual_mean,
    residual_bootstrap_sd = residual_sd,
    bootstrap_statistics = bootstrap_statistics,
    bootstrap_ok = attempt_ok[seq_len(attempt)],
    bootstrap_errors = attempt_errors[seq_len(attempt)],
    n_boot_ok = n_boot_ok,
    n_bootstrap_attempts = attempt,
    cores = cores,
    B = B,
    alpha = alpha,
    fit = fit,
    seeds = list(
      seed = seed,
      bootstrap_attempts = attempt_seeds[seq_len(attempt)],
      bootstrap_successes = successful_seeds
    )
  )
}
