cat("LNBC bootstrap_core\n")

lnbc_metrics <- c("AUC", "opt_cut", "se", "sp", "J")
lnbc_methods <- c("EMBC", "NP", "NP-inc", "LBNP", "LBNP-log", "EMBC_fixk")

.lnbc_bootstrap_state <- new.env(parent = emptyenv())

lnbc_repo_root <- function() {
  repo_root <- getwd()
  if (file.exists(file.path(repo_root, "bootstrap_core.R"))) {
    return(repo_root)
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      repo_root <- dirname(normalizePath(frames[[i]]$ofile))
      break
    }
  }
  if (!file.exists(file.path(repo_root, "bootstrap_core.R"))) {
    stop("Run this script from the LNBC repository root.")
  }
  repo_root
}

lnbc_source_method <- function(method, repo_root = lnbc_repo_root()) {
  method <- match.arg(method, lnbc_methods)

  if (method %in% c("EMBC", "EMBC_fixk") && !isTRUE(.lnbc_bootstrap_state$EMBC)) {
    source(file.path(repo_root, "EMBC", "method.R"), local = .GlobalEnv)
    .lnbc_bootstrap_state$EMBC <- TRUE
  }
  if (method %in% c("NP", "NP-inc") && !isTRUE(.lnbc_bootstrap_state$NP)) {
    source(file.path(repo_root, "NP", "method_np.r"), local = .GlobalEnv)
    .lnbc_bootstrap_state$NP <- TRUE
  }
  if (method %in% c("LBNP", "LBNP-log") && !isTRUE(.lnbc_bootstrap_state$LBNP)) {
    source(file.path(repo_root, "LBNP", "method_lbnp.r"), local = .GlobalEnv)
    .lnbc_bootstrap_state$LBNP <- TRUE
  }
  invisible(TRUE)
}

lnbc_pi_value <- function(pi_key) {
  pi_key <- as.character(pi_key)
  if (pi_key == "9") return(0.9)
  if (pi_key == "75") return(0.75)
  value <- suppressWarnings(as.numeric(pi_key))
  if (is.na(value)) stop("Unknown pi key: ", pi_key)
  if (value > 1) value <- value / 100
  if (!is.finite(value) || value <= 0 || value > 1) {
    stop("pi must be in (0, 1].")
  }
  value
}

lnbc_replication_seeds <- function(data_reps, seed_pool_size = 1000) {
  data_reps <- as.integer(data_reps)
  if (is.na(data_reps) || data_reps < 1) {
    stop("data_reps must be a positive integer.")
  }
  seed_pool_size <- max(as.integer(seed_pool_size), data_reps)
  set.seed(123)
  sample(1:1e8, seed_pool_size)[seq_len(data_reps)]
}

lnbc_scenario <- function(dist, J_key) {
  J_key <- as.character(J_key)
  switch(
    dist,
    lognorm = list(kind = "lognorm", label = "Lognormal", a0 = 0, b0 = 1,
                   a1 = if (J_key == "3") 0.77 else 1.35, b1 = 1),
    gamma = list(kind = "gamma", label = "Gamma", a0 = 1.5, b0 = 1,
                 a1 = if (J_key == "3") 2.47 else 3.39, b1 = 1),
    weibull = list(kind = "weibull", label = "Weibull", a0 = 0.5, b0 = 0.5,
                   a1 = if (J_key == "3") 2.6793 else 9.7345, b1 = 0.5),
    gamma2 = list(kind = "gamma2", label = "Gamma unequal-rate", a0 = 1, b0 = 0.5,
                  a1 = 1, b1 = if (J_key == "3") 1.157433 else 2.201749),
    weibull2 = list(kind = "weibull2", label = "Weibull b0=b1=1", a0 = 1, b0 = 1,
                    a1 = if (J_key == "3") 2.314865 else 4.403498, b1 = 1),
    stop("Unknown distribution: ", dist)
  )
}

lnbc_true_lognorm <- function(a0, a1, b0, b1) {
  opt <- optimize(
    function(z) -(plnorm(exp(z), a0, sqrt(b0)) - plnorm(exp(z), a1, sqrt(b1))),
    c(min(a0, a1) - 10 * sqrt(max(b0, b1)), max(a0, a1) + 10 * sqrt(max(b0, b1)))
  )
  cstar <- exp(opt$minimum)
  se <- 1 - plnorm(cstar, a1, sqrt(b1))
  sp <- plnorm(cstar, a0, sqrt(b0))
  c(AUC = pnorm((a1 - a0) / sqrt(b0 + b1)), opt_cut = cstar,
    se = se, sp = sp, J = se + sp - 1)
}

lnbc_true_gamma_equal <- function(a0, a1, rate) {
  cstar <- exp((lgamma(a1) - lgamma(a0) - (a1 - a0) * log(rate)) / (a1 - a0))
  se <- 1 - pgamma(cstar, shape = a1, rate = rate)
  sp <- pgamma(cstar, shape = a0, rate = rate)
  c(AUC = pbeta(0.5, shape1 = a0, shape2 = a1), opt_cut = cstar,
    se = se, sp = sp, J = se + sp - 1)
}

lnbc_true_gamma_unequal <- function(a, b0, b1) {
  cstar <- a * log(b1 / b0) / (b1 - b0)
  if (b1 > b0) {
    auc <- pbeta(b1 / (b0 + b1), shape1 = a, shape2 = a)
    se <- pgamma(cstar, shape = a, rate = b1)
    sp <- 1 - pgamma(cstar, shape = a, rate = b0)
  } else {
    auc <- 1 - pbeta(b1 / (b0 + b1), shape1 = a, shape2 = a)
    se <- 1 - pgamma(cstar, shape = a, rate = b1)
    sp <- pgamma(cstar, shape = a, rate = b0)
  }
  c(AUC = auc, opt_cut = cstar, se = se, sp = sp, J = se + sp - 1)
}

lnbc_true_weibull <- function(a0, a1, shape) {
  d <- a0^(-shape) - a1^(-shape)
  cstar <- (-shape * log(a0 / a1) / d)^(1 / shape)
  se <- 1 - pweibull(cstar, shape = shape, scale = a1)
  sp <- pweibull(cstar, shape = shape, scale = a0)
  c(AUC = a1^shape / (a0^shape + a1^shape), opt_cut = cstar,
    se = se, sp = sp, J = se + sp - 1)
}

lnbc_true_values <- function(dist, J_key) {
  s <- lnbc_scenario(dist, J_key)
  if (dist == "lognorm") return(lnbc_true_lognorm(s$a0, s$a1, s$b0, s$b1))
  if (dist == "gamma") return(lnbc_true_gamma_equal(s$a0, s$a1, s$b0))
  if (dist == "gamma2") return(lnbc_true_gamma_unequal(s$a0, s$b0, s$b1))
  lnbc_true_weibull(s$a0, s$a1, s$b0)
}

lnbc_make_data <- function(n, dist, J_key, pi_key, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  s <- lnbc_scenario(dist, J_key)
  pi_val <- lnbc_pi_value(pi_key)

  n0 <- as.integer(n)
  n1 <- as.integer(n)
  if (is.na(n0) || n0 < 1) stop("n must be a positive integer.")
  n00 <- rbinom(1, n0, pi_val)
  n11 <- rbinom(1, n1, pi_val)

  if (s$kind == "lognorm") {
    x <- c(rlnorm(n00, meanlog = s$a0, sdlog = sqrt(s$b0)),
           rlnorm(n0 - n00, meanlog = s$a1, sdlog = sqrt(s$b1)))
    y <- c(rlnorm(n11, meanlog = s$a1, sdlog = sqrt(s$b1)),
           rlnorm(n1 - n11, meanlog = s$a0, sdlog = sqrt(s$b0)))
  } else if (s$kind == "gamma") {
    x <- c(rgamma(n00, shape = s$a0, rate = s$b0),
           rgamma(n0 - n00, shape = s$a1, rate = s$b1))
    y <- c(rgamma(n11, shape = s$a1, rate = s$b1),
           rgamma(n1 - n11, shape = s$a0, rate = s$b0))
  } else if (s$kind == "gamma2") {
    x <- c(rgamma(n00, shape = s$a0, rate = s$b0),
           rgamma(n0 - n00, shape = s$a1, rate = s$b1))
    y <- c(rgamma(n11, shape = s$a1, rate = s$b1),
           rgamma(n1 - n11, shape = s$a0, rate = s$b0))
  } else {
    x <- c(rweibull(n00, shape = s$b0, scale = s$a0),
           rweibull(n0 - n00, shape = s$b1, scale = s$a1))
    y <- c(rweibull(n11, shape = s$b1, scale = s$a1),
           rweibull(n1 - n11, shape = s$b0, scale = s$a0))
  }

  list(
    T = c(x, y),
    R = c(rep(0, n0), rep(1, n1)),
    G = c(rep(0, n00), rep(1, n0 - n00), rep(1, n11), rep(0, n1 - n11)),
    true_param = lnbc_true_values(dist, J_key),
    scenario = s
  )
}

lnbc_check_data <- function(data) {
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
  if (!all(is.finite(data$T))) {
    stop("T must be finite.")
  }
  if (!is.null(data$G) && length(data$G) != length(data$T)) {
    stop("G must have the same length as T when provided.")
  }
  data
}

lnbc_empty_estimate <- function() {
  setNames(rep(NA_real_, length(lnbc_metrics)), lnbc_metrics)
}

lnbc_np_transform <- function(x, dist) {
  if (dist == "real_logged") {
    x_range <- range(x)
    if (!all(is.finite(x_range)) || diff(x_range) <= 0) {
      stop("Logged real-data biomarker must have a nonzero finite range.")
    }
    return(list(
      value = (x - x_range[1]) / diff(x_range),
      inverse = function(z) z * diff(x_range) + x_range[1]
    ))
  }
  if (dist == "real") {
    if (any(x <= 0)) stop("Real-data NP analysis requires positive biomarker values.")
    x_log <- log(x)
    x_range <- range(x_log)
    if (!all(is.finite(x_range)) || diff(x_range) <= 0) {
      stop("Real-data biomarker must have a nonzero finite range.")
    }
    return(list(
      value = (x_log - x_range[1]) / diff(x_range),
      inverse = function(z) exp(z * diff(x_range) + x_range[1])
    ))
  }
  if (dist == "gamma2") {
    return(list(value = 1 / (1 + x), inverse = function(z) (1 - z) / z))
  }
  list(value = x / (1 + x), inverse = function(z) z / (1 - z))
}

lnbc_lbnp_transform <- function(x, dist, transform = "identity") {
  transform <- match.arg(transform, c("identity", "log"))
  if (transform == "identity") {
    if (dist == "gamma2") return(list(value = -x, inverse = function(z) -z))
    return(list(value = x, inverse = function(z) z))
  }
  if (dist == "gamma2") return(list(value = -log(x), inverse = function(z) exp(-z)))
  list(value = log(x), inverse = exp)
}

lnbc_default_roc_grid <- function(n_points = 101) {
  n_points <- as.integer(n_points)
  if (is.na(n_points) || n_points < 2) {
    stop("n_points must be an integer of at least 2.")
  }
  seq(0, 1, length.out = n_points)
}

lnbc_check_roc_grid <- function(fpr_grid) {
  fpr_grid <- as.numeric(fpr_grid)
  if (length(fpr_grid) < 2 || any(!is.finite(fpr_grid)) ||
      any(fpr_grid < 0) || any(fpr_grid > 1)) {
    stop("fpr_grid must contain at least two finite values in [0, 1].")
  }
  sort(unique(fpr_grid))
}

lnbc_empty_roc_curve <- function(fpr_grid) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  data.frame(fpr = fpr_grid, tpr = rep(NA_real_, length(fpr_grid)))
}

lnbc_sanitize_roc_curve <- function(roc, fpr_grid) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  if (is.null(roc)) return(lnbc_empty_roc_curve(fpr_grid))
  if (is.numeric(roc)) {
    roc <- data.frame(fpr = fpr_grid, tpr = as.numeric(roc))
  }
  if (!is.data.frame(roc) || !"tpr" %in% names(roc)) {
    return(lnbc_empty_roc_curve(fpr_grid))
  }
  if (!"fpr" %in% names(roc)) roc$fpr <- fpr_grid
  out <- data.frame(fpr = fpr_grid, tpr = rep(NA_real_, length(fpr_grid)))
  keep <- is.finite(roc$fpr) & is.finite(roc$tpr)
  if (any(keep)) {
    out$tpr <- approx(
      x = pmin(pmax(roc$fpr[keep], 0), 1),
      y = pmin(pmax(roc$tpr[keep], 0), 1),
      xout = fpr_grid,
      ties = max,
      rule = 2
    )$y
  }
  out$tpr[out$fpr <= 0] <- 0
  out$tpr[out$fpr >= 1] <- 1
  out$tpr <- pmin(pmax(out$tpr, 0), 1)
  out
}

lnbc_weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) {
    if (!is.finite(p)) return(NA_real_)
    if (p <= 0) return(x[1])
    if (p >= 1) return(x[length(x)])
    x[which(cw >= p)[1]]
  }, numeric(1))
}

lnbc_weighted_ecdf_eval <- function(x, w, q) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(q)))
  x <- x[ok]
  w <- w[ok] / sum(w[ok])
  vapply(q, function(z) {
    if (!is.finite(z)) return(NA_real_)
    sum(w[x <= z])
  }, numeric(1))
}

lnbc_weighted_roc_curve <- function(score, w0, w1, fpr_grid) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  cutoff <- lnbc_weighted_quantile(score, w0, 1 - fpr_grid)
  tpr <- 1 - lnbc_weighted_ecdf_eval(score, w1, cutoff)
  lnbc_sanitize_roc_curve(data.frame(fpr = fpr_grid, tpr = tpr), fpr_grid)
}

lnbc_embc_roc_curve <- function(fit, data, fpr_grid) {
  T <- data$T
  score <- fit$alpha + fit$beta * boxcox(T, fit$k)
  w0 <- fit$pi
  w1 <- fit$pi * exp(score)
  lnbc_weighted_roc_curve(score, w0, w1, fpr_grid)
}

lnbc_np_roc_curve <- function(dat, pi0, pi1, fpr_grid) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  if (!is.data.frame(dat) || !all(c("R", "biomarker") %in% names(dat))) {
    stop("dat must be a data.frame with columns R and biomarker.")
  }
  R <- dat$R
  Y <- dat$biomarker
  if (!all(R %in% c(0, 1)) || any(!is.finite(Y))) {
    stop("Invalid R/biomarker values for NP ROC.")
  }
  denom <- pi1 + pi0 - 1
  if (!is.finite(denom) || denom <= 0) {
    stop("pi0 + pi1 must be greater than 1.")
  }

  EP <- sort(unique(Y))
  F1 <- rep(0, length(EP) + 1)
  F0 <- rep(0, length(EP) + 1)
  for (i in seq_along(EP)) {
    G1_t <- sum(Y <= EP[i] & R == 1) / sum(R == 1)
    G0_t <- sum(Y <= EP[i] & R == 0) / sum(R == 0)
    F1[i + 1] <- pi0 / denom * G1_t - (1 - pi1) / denom * G0_t
    F0[i + 1] <- pi1 / denom * G0_t - (1 - pi0) / denom * G1_t
  }

  EP_ext <- c(0, EP, 1)
  delta_EP <- diff(EP_ext)
  F_points <- c(F1, F0)
  F_points <- F_points * (F_points >= 0) * (F_points <= 1)
  F_points <- sort(unique(F_points))
  if (length(F_points) == 0) F_points <- c(0, 1)

  Q1 <- Q0 <- numeric(length(F_points))
  for (i in seq_along(F_points)) {
    Q1[i] <- sum((F1 <= F_points[i]) * delta_EP)
    Q0[i] <- sum((F0 <= F_points[i]) * delta_EP)
  }

  roc_fun <- function(s) {
    tem <- sum((F0 <= (1 - s)) * delta_EP)
    idx <- max(c(which(Q1 <= tem), 1))
    1 - F_points[idx]
  }
  tpr <- vapply(fpr_grid, roc_fun, numeric(1))
  lnbc_sanitize_roc_curve(data.frame(fpr = fpr_grid, tpr = tpr), fpr_grid)
}

lnbc_lbnp_roc_curve <- function(fit, fpr_grid) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  roc <- tryCatch(
    lbnp.roc(fpr_grid, fit),
    error = function(e) rep(NA_real_, length(fpr_grid))
  )
  lnbc_sanitize_roc_curve(roc, fpr_grid)
}

lnbc_true_class_cdf <- function(x, dist, J_key, class_id) {
  s <- lnbc_scenario(dist, J_key)
  if (class_id == 0) {
    a <- s$a0
    b <- s$b0
  } else {
    a <- s$a1
    b <- s$b1
  }
  if (dist == "lognorm") return(plnorm(x, meanlog = a, sdlog = sqrt(b)))
  if (dist %in% c("gamma", "gamma2")) return(pgamma(x, shape = a, rate = b))
  pweibull(x, shape = b, scale = a)
}

lnbc_true_class_quantile <- function(p, dist, J_key, class_id) {
  s <- lnbc_scenario(dist, J_key)
  if (class_id == 0) {
    a <- s$a0
    b <- s$b0
  } else {
    a <- s$a1
    b <- s$b1
  }
  p <- pmin(pmax(p, 0), 1)
  if (dist == "lognorm") return(qlnorm(p, meanlog = a, sdlog = sqrt(b)))
  if (dist %in% c("gamma", "gamma2")) return(qgamma(p, shape = a, rate = b))
  qweibull(p, shape = b, scale = a)
}

lnbc_true_roc_curve <- function(dist, J_key, fpr_grid = lnbc_default_roc_grid()) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  if (dist == "gamma2") {
    cutoff <- lnbc_true_class_quantile(fpr_grid, dist, J_key, class_id = 0)
    tpr <- lnbc_true_class_cdf(cutoff, dist, J_key, class_id = 1)
  } else {
    cutoff <- lnbc_true_class_quantile(1 - fpr_grid, dist, J_key, class_id = 0)
    tpr <- 1 - lnbc_true_class_cdf(cutoff, dist, J_key, class_id = 1)
  }
  lnbc_sanitize_roc_curve(data.frame(fpr = fpr_grid, tpr = tpr), fpr_grid)
}

lnbc_fit_method_once <- function(method, data, pi0, pi1, dist, seed = NULL,
                                 repo_root = lnbc_repo_root(),
                                 lbnp_k = 10,
                                 lbnp_nu_seq = 3 * 10^seq(from = -2, to = 2.5, length.out = 10),
                                 lbnp_maxit = 500,
                                 lbnp_thres = 1e-5,
                                 lbnp_transform = NULL,
                                 embc_fixed_k = NULL,
                                 roc_grid = NULL) {
  method <- match.arg(method, lnbc_methods)
  data <- lnbc_check_data(data)
  if (!is.null(roc_grid)) roc_grid <- lnbc_check_roc_grid(roc_grid)
  lnbc_source_method(method, repo_root = repo_root)
  if (!is.null(seed)) set.seed(seed)

  out <- tryCatch({
    if (method %in% c("EMBC", "EMBC_fixk")) {
      embc_data <- data
      inverse_cutoff <- identity
      if (dist == "real") {
        if (any(embc_data$T <= 0)) {
          stop("Real-data EMBC analysis requires positive biomarker values.")
        }
        embc_data$T <- log(embc_data$T)
        inverse_cutoff <- exp
      }
      if (method == "EMBC_fixk" && !is.null(embc_fixed_k)) {
        fit <- EMBC_fixpi_fixk(
          embc_data, pi0 = pi0, pi1 = pi1, k = embc_fixed_k
        )
      } else {
        fit <- EMBC_fixpi(embc_data, pi0 = pi0, pi1 = pi1)
      }
      estimate <- get_est(fit, embc_data, EM = TRUE)[lnbc_metrics]
      estimate["opt_cut"] <- inverse_cutoff(estimate["opt_cut"])
      roc <- if (!is.null(roc_grid)) {
        lnbc_embc_roc_curve(fit$best, embc_data, roc_grid)
      } else NULL
      list(estimate = estimate, ok = all(is.finite(estimate)),
           conv = TRUE, error = NA_character_, nu = NA_real_, roc = roc,
           embc_k = unname(fit$best$k))
    } else if (method %in% c("NP", "NP-inc")) {
      tr <- lnbc_np_transform(data$T, dist)
      dat <- data.frame(R = data$R, biomarker = tr$value)
      fit <- if (method == "NP-inc") {
        np.inc.roc(dat, pi0 = pi0, pi1 = pi1)
      } else {
        np.roc(dat, pi0 = pi0, pi1 = pi1,
               paper_youden = dist %in% c("real", "real_logged"))
      }
      estimate <- c(AUC = fit$AUC[1], opt_cut = tr$inverse(fit$opt_cut[1]),
                    se = fit$se[1], sp = fit$sp[1], J = fit$J[1])
      roc <- if (is.null(roc_grid)) {
        NULL
      } else if (method == "NP-inc") {
        np.inc.roc.curve(dat, pi0 = pi0, pi1 = pi1, fpr = roc_grid)
      } else {
        lnbc_np_roc_curve(dat, pi0 = pi0, pi1 = pi1, fpr_grid = roc_grid)
      }
      list(estimate = estimate[lnbc_metrics], ok = all(is.finite(estimate)),
           conv = TRUE, error = NA_character_, nu = NA_real_, roc = roc)
    } else {
      transform <- lbnp_transform
      if (is.null(transform)) {
        transform <- if (method == "LBNP-log") "log" else "identity"
      }
      tr <- lnbc_lbnp_transform(data$T, dist, transform = transform)
      fit <- lbnp.fit.cv5(
        data.frame(R = data$R, biomarker = tr$value),
        pi0 = pi0,
        pi1 = pi1,
        nu.seq = lbnp_nu_seq,
        k = lbnp_k,
        maxit = lbnp_maxit,
        thres = lbnp_thres,
        seed = if (is.null(seed)) 123 else seed
      )
      if (!isTRUE(fit$conv)) {
        list(estimate = lnbc_empty_estimate(), ok = FALSE, conv = FALSE,
             error = "LBNP did not converge", nu = fit$nu,
             roc = if (!is.null(roc_grid)) lnbc_empty_roc_curve(roc_grid) else NULL)
      } else {
        estimate <- c(AUC = fit$stats$AUC[1],
                      opt_cut = tr$inverse(fit$stats$opt_cut[1]),
                      se = fit$stats$se[1],
                      sp = fit$stats$sp[1],
                      J = fit$stats$J[1])
        roc <- if (!is.null(roc_grid)) lnbc_lbnp_roc_curve(fit, roc_grid) else NULL
        list(estimate = estimate[lnbc_metrics], ok = all(is.finite(estimate)),
             conv = TRUE, error = NA_character_, nu = fit$nu, roc = roc)
      }
    }
  }, error = function(e) {
    list(estimate = lnbc_empty_estimate(), ok = FALSE, conv = FALSE,
         error = conditionMessage(e), nu = NA_real_,
         roc = if (!is.null(roc_grid)) lnbc_empty_roc_curve(roc_grid) else NULL)
  })

  out$estimate <- as.numeric(out$estimate)
  names(out$estimate) <- lnbc_metrics
  if (!is.null(roc_grid)) {
    out$roc <- lnbc_sanitize_roc_curve(out$roc, roc_grid)
  }
  out
}

lnbc_stratified_boot_sample <- function(data, seed = NULL) {
  data <- lnbc_check_data(data)
  if (!is.null(seed)) set.seed(seed)

  id0 <- which(data$R == 0)
  id1 <- which(data$R == 1)
  boot_id <- c(sample(id0, length(id0), replace = TRUE),
               sample(id1, length(id1), replace = TRUE))

  out <- list(T = data$T[boot_id], R = data$R[boot_id])
  if (!is.null(data$G)) out$G <- data$G[boot_id]
  if (!is.null(data$true_param)) out$true_param <- data$true_param
  if (!is.null(data$scenario)) out$scenario <- data$scenario
  out
}

lnbc_percentile_ci <- function(point_est, boot_mat, conf_level = 0.95,
                               log_cutoff = TRUE) {
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.")
  }
  alpha <- (1 - conf_level) / 2
  probs <- c(alpha, 1 - alpha)
  boot_mat <- as.matrix(boot_mat)

  rows <- lapply(lnbc_metrics, function(metric) {
    estimate <- unname(point_est[metric])
    vals <- as.numeric(boot_mat[, metric])
    usable <- is.finite(vals)
    transformation <- "identity"

    if (metric == "opt_cut" && isTRUE(log_cutoff)) {
      usable <- usable & vals > 0
      vals <- log(vals[usable])
      transformation <- "log"
    } else {
      vals <- vals[usable]
    }

    if (!is.finite(estimate) || length(vals) == 0) {
      ci <- c(NA_real_, NA_real_)
    } else {
      ci <- as.numeric(quantile(vals, probs = probs, names = FALSE, na.rm = FALSE))
      if (transformation == "log") ci <- exp(ci)
    }

    data.frame(
      metric = metric,
      estimate = estimate,
      ci_low = ci[1],
      ci_high = ci[2],
      ci_length = ci[2] - ci[1],
      n_boot = nrow(boot_mat),
      n_valid = length(vals),
      n_missing = nrow(boot_mat) - length(vals),
      transformation = transformation,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lnbc_bootstrap_method <- function(data, method, B = 5, pi0, pi1, dist,
                                  seed = 123, conf_level = 0.95,
                                  repo_root = lnbc_repo_root(),
                                  lbnp_k = 10,
                                  lbnp_nu_seq = 3 * 10^seq(from = -2, to = 2.5, length.out = 10),
                                  lbnp_maxit = 500,
                                  lbnp_thres = 1e-5,
                                  lbnp_transform = NULL,
                                  cores = 1,
                                  verbose = TRUE) {
  B <- as.integer(B)
  if (is.na(B) || B < 1) stop("B must be a positive integer.")
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1) cores <- 1
  if (.Platform$OS.type != "unix") cores <- 1
  cores <- min(cores, B)

  data <- lnbc_check_data(data)
  set.seed(seed)
  point_seed <- sample.int(.Machine$integer.max, 1)
  boot_seeds <- sample.int(.Machine$integer.max, B)
  fit_seeds <- sample.int(.Machine$integer.max, B)

  if (isTRUE(verbose)) {
    cat("Fitting point estimate:", method, "\n")
  }
  point <- lnbc_fit_method_once(
    method, data, pi0 = pi0, pi1 = pi1, dist = dist,
    seed = point_seed, repo_root = repo_root, lbnp_k = lbnp_k,
    lbnp_nu_seq = lbnp_nu_seq, lbnp_transform = lbnp_transform
  )

  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(lnbc_metrics),
                     dimnames = list(paste0("boot_", seq_len(B)), lnbc_metrics))
  boot_status <- data.frame(
    boot = seq_len(B),
    ok = rep(FALSE, B),
    conv = rep(FALSE, B),
    nu = rep(NA_real_, B),
    error = rep(NA_character_, B),
    stringsAsFactors = FALSE
  )

  fit_boot <- function(b) {
    if (isTRUE(verbose)) {
      cat("Bootstrap", b, "of", B, "\n")
    }
    boot_data <- lnbc_stratified_boot_sample(data, seed = boot_seeds[b])
    fit <- lnbc_fit_method_once(
      method, boot_data, pi0 = pi0, pi1 = pi1, dist = dist,
      seed = fit_seeds[b], repo_root = repo_root, lbnp_k = lbnp_k,
      lbnp_nu_seq = lbnp_nu_seq, lbnp_transform = lbnp_transform
    )
    list(b = b, fit = fit)
  }

  if (cores > 1) {
    if (isTRUE(verbose)) {
      cat("Running", B, "bootstrap fits with", cores, "cores\n")
    }
    boot_fits <- parallel::mclapply(seq_len(B), fit_boot, mc.cores = cores,
                                    mc.preschedule = FALSE)
  } else {
    boot_fits <- lapply(seq_len(B), fit_boot)
  }

  for (boot_result in boot_fits) {
    b <- boot_result$b
    fit <- boot_result$fit
    boot_mat[b, ] <- fit$estimate
    boot_status$ok[b] <- isTRUE(fit$ok)
    boot_status$conv[b] <- isTRUE(fit$conv)
    boot_status$nu[b] <- fit$nu
    boot_status$error[b] <- fit$error
  }

  ci <- lnbc_percentile_ci(point$estimate, boot_mat, conf_level = conf_level)
  list(
    method = method,
    dist = dist,
    pi0 = pi0,
    pi1 = pi1,
    B = B,
    conf_level = conf_level,
    n0 = sum(data$R == 0),
    n1 = sum(data$R == 1),
    point = point$estimate,
    point_status = data.frame(ok = point$ok, conv = point$conv, nu = point$nu,
                              error = point$error, stringsAsFactors = FALSE),
    bootstrap = boot_mat,
    bootstrap_status = boot_status,
    ci = ci,
    seeds = list(point_seed = point_seed, boot_seeds = boot_seeds, fit_seeds = fit_seeds)
  )
}

lnbc_roc_percentile_band <- function(point_roc, boot_roc, conf_level = 0.95) {
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.")
  }
  if (!is.data.frame(point_roc) || !all(c("fpr", "tpr") %in% names(point_roc))) {
    stop("point_roc must be a data.frame with fpr and tpr.")
  }
  boot_roc <- as.matrix(boot_roc)
  if (ncol(boot_roc) != nrow(point_roc)) {
    stop("boot_roc must have one column for each point_roc row.")
  }

  alpha <- (1 - conf_level) / 2
  probs <- c(alpha, 1 - alpha)
  ci <- t(vapply(seq_len(ncol(boot_roc)), function(j) {
    vals <- as.numeric(boot_roc[, j])
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) return(c(NA_real_, NA_real_))
    as.numeric(quantile(vals, probs = probs, names = FALSE, na.rm = FALSE))
  }, numeric(2)))
  n_valid <- colSums(is.finite(boot_roc))

  data.frame(
    fpr = point_roc$fpr,
    tpr = point_roc$tpr,
    ci_low = ci[, 1],
    ci_high = ci[, 2],
    n_boot = nrow(boot_roc),
    n_valid = n_valid,
    n_missing = nrow(boot_roc) - n_valid,
    stringsAsFactors = FALSE
  )
}

lnbc_roc_band_method <- function(data, method, B = 100, pi0, pi1, dist,
                                 seed = 123, conf_level = 0.95,
                                 fpr_grid = lnbc_default_roc_grid(),
                                 repo_root = lnbc_repo_root(),
                                 lbnp_k = 10,
                                 lbnp_nu_seq = 3 * 10^seq(from = -2, to = 2.5, length.out = 10),
                                 lbnp_transform = NULL,
                                 cores = 1,
                                 verbose = TRUE) {
  method <- match.arg(method, lnbc_methods)
  B <- as.integer(B)
  if (is.na(B) || B < 1) stop("B must be a positive integer.")
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1) cores <- 1
  if (.Platform$OS.type != "unix") cores <- 1
  cores <- min(cores, B)

  data <- lnbc_check_data(data)
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  set.seed(seed)
  point_seed <- sample.int(.Machine$integer.max, 1)
  boot_seeds <- sample.int(.Machine$integer.max, B)
  fit_seeds <- sample.int(.Machine$integer.max, B)

  if (isTRUE(verbose)) {
    cat("Fitting ROC point curve:", method, "\n")
  }
  point <- lnbc_fit_method_once(
    method, data, pi0 = pi0, pi1 = pi1, dist = dist,
    seed = point_seed, repo_root = repo_root, lbnp_k = lbnp_k,
    lbnp_nu_seq = lbnp_nu_seq, lbnp_maxit = lbnp_maxit,
    lbnp_thres = lbnp_thres, lbnp_transform = lbnp_transform,
    roc_grid = fpr_grid
  )

  boot_roc <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(fpr_grid),
    dimnames = list(paste0("boot_", seq_len(B)), paste0("fpr_", fpr_grid))
  )
  boot_status <- data.frame(
    boot = seq_len(B),
    ok = rep(FALSE, B),
    conv = rep(FALSE, B),
    roc_ok = rep(FALSE, B),
    nu = rep(NA_real_, B),
    error = rep(NA_character_, B),
    stringsAsFactors = FALSE
  )

  fit_boot <- function(b) {
    if (isTRUE(verbose)) {
      cat("ROC bootstrap", b, "of", B, "\n")
    }
    boot_data <- lnbc_stratified_boot_sample(data, seed = boot_seeds[b])
    fit <- lnbc_fit_method_once(
      method, boot_data, pi0 = pi0, pi1 = pi1, dist = dist,
      seed = fit_seeds[b], repo_root = repo_root, lbnp_k = lbnp_k,
      lbnp_nu_seq = lbnp_nu_seq, lbnp_maxit = lbnp_maxit,
      lbnp_thres = lbnp_thres, lbnp_transform = lbnp_transform,
      roc_grid = fpr_grid
    )
    list(b = b, fit = fit)
  }

  if (cores > 1) {
    if (isTRUE(verbose)) {
      cat("Running", B, "ROC bootstrap fits with", cores, "cores\n")
    }
    boot_fits <- parallel::mclapply(seq_len(B), fit_boot, mc.cores = cores,
                                    mc.preschedule = FALSE)
  } else {
    boot_fits <- lapply(seq_len(B), fit_boot)
  }

  for (boot_result in boot_fits) {
    b <- boot_result$b
    fit <- boot_result$fit
    boot_roc[b, ] <- fit$roc$tpr
    boot_status$ok[b] <- isTRUE(fit$ok)
    boot_status$conv[b] <- isTRUE(fit$conv)
    boot_status$roc_ok[b] <- all(is.finite(fit$roc$tpr))
    boot_status$nu[b] <- fit$nu
    boot_status$error[b] <- fit$error
  }

  band <- lnbc_roc_percentile_band(point$roc, boot_roc, conf_level = conf_level)
  list(
    method = method,
    dist = dist,
    pi0 = pi0,
    pi1 = pi1,
    B = B,
    conf_level = conf_level,
    n0 = sum(data$R == 0),
    n1 = sum(data$R == 1),
    fpr_grid = fpr_grid,
    point = point$estimate,
    point_roc = point$roc,
    point_status = data.frame(
      ok = point$ok,
      conv = point$conv,
      roc_ok = all(is.finite(point$roc$tpr)),
      nu = point$nu,
      error = point$error,
      stringsAsFactors = FALSE
    ),
    bootstrap_roc = boot_roc,
    bootstrap_status = boot_status,
    band = band,
    seeds = list(point_seed = point_seed, boot_seeds = boot_seeds, fit_seeds = fit_seeds)
  )
}

lnbc_plot_roc_bands <- function(band_df, true_roc = NULL, png_file = NULL,
                                pdf_file = NULL, title = "ROC confidence bands",
                                width = 7, height = 6, dpi = 150) {
  if (!is.data.frame(band_df) || !all(c("method", "fpr", "tpr", "ci_low", "ci_high") %in% names(band_df))) {
    stop("band_df must contain method, fpr, tpr, ci_low, and ci_high.")
  }
  methods <- unique(band_df$method)
  palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
               "#66a61e", "#e6ab02", "#a6761d", "#666666")
  cols <- setNames(rep(palette, length.out = length(methods)), methods)

  draw_plot <- function() {
    old_mar <- par("mar")
    on.exit(par(mar = old_mar), add = TRUE)
    par(mar = c(4.2, 4.5, 3.2, 1.2))
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
         xlab = "False positive rate", ylab = "Sensitivity", main = title)
    grid(col = "gray88", lty = 1)
    abline(0, 1, col = "gray70", lty = 3)

    for (method in methods) {
      x <- band_df[band_df$method == method, ]
      x <- x[order(x$fpr), ]
      ok <- is.finite(x$fpr) & is.finite(x$ci_low) & is.finite(x$ci_high)
      if (sum(ok) >= 2) {
        polygon(
          c(x$fpr[ok], rev(x$fpr[ok])),
          c(x$ci_low[ok], rev(x$ci_high[ok])),
          col = grDevices::adjustcolor(cols[[method]], alpha.f = 0.18),
          border = NA
        )
      }
    }

    if (!is.null(true_roc)) {
      lines(true_roc$fpr, true_roc$tpr, col = "black", lty = 2, lwd = 2)
    }

    for (method in methods) {
      x <- band_df[band_df$method == method, ]
      x <- x[order(x$fpr), ]
      lines(x$fpr, x$tpr, col = cols[[method]], lwd = 2)
    }

    legend_labels <- methods
    legend_cols <- unname(cols[methods])
    legend_lty <- rep(1, length(methods))
    legend_lwd <- rep(2, length(methods))
    if (!is.null(true_roc)) {
      legend_labels <- c(legend_labels, "Truth")
      legend_cols <- c(legend_cols, "black")
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 2)
    }
    legend("bottomright", legend = legend_labels, col = legend_cols,
           lty = legend_lty, lwd = legend_lwd, bg = "white", cex = 0.85)
  }

  out <- c()
  if (!is.null(png_file)) {
    dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(png_file, width = width, height = height, units = "in", res = dpi)
    draw_plot()
    grDevices::dev.off()
    out <- c(out, png = png_file)
  }
  if (!is.null(pdf_file)) {
    dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(pdf_file, width = width, height = height)
    draw_plot()
    grDevices::dev.off()
    out <- c(out, pdf = pdf_file)
  }
  invisible(out)
}

lnbc_run_roc_band_study <- function(methods = lnbc_methods,
                                    dist,
                                    n,
                                    J_key,
                                    pi_key,
                                    B = 100,
                                    rep_id = 1,
                                    seed = 123,
                                    data_seed = NULL,
                                    conf_level = 0.95,
                                    fpr_grid = lnbc_default_roc_grid(),
                                    repo_root = lnbc_repo_root(),
                                    out_dir = file.path(repo_root, "ROC_plot"),
                                    mode = "custom",
                                    cores = 1,
                                    save_plot = TRUE,
                                    verbose = TRUE) {
  methods <- match.arg(methods, lnbc_methods, several.ok = TRUE)
  n <- as.integer(n)
  B <- as.integer(B)
  rep_id <- as.integer(rep_id)
  if (is.na(n) || n < 1) stop("n must be a positive integer.")
  if (is.na(B) || B < 1) stop("B must be a positive integer.")
  if (is.na(rep_id) || rep_id < 1) stop("rep_id must be a positive integer.")

  pi_val <- lnbc_pi_value(pi_key)
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  if (is.null(data_seed)) {
    data_seed <- lnbc_replication_seeds(rep_id)[rep_id]
  }
  data <- lnbc_make_data(n = n, dist = dist, J_key = J_key,
                         pi_key = pi_key, seed = data_seed)
  true_roc <- lnbc_true_roc_curve(dist, J_key, fpr_grid)
  truth <- lnbc_true_values(dist, J_key)

  results <- list()
  for (method in methods) {
    if (isTRUE(verbose)) {
      cat("ROC band study method", method, "dist", dist, "n", n,
          "J", J_key, "pi", pi_key, "\n")
    }
    method_seed <- lnbc_grid_boot_seed(
      rep_id = rep_id,
      method = method,
      dist = dist,
      n = n,
      J_key = J_key,
      pi_key = pi_key,
      B = B,
      seed = seed
    )
    results[[method]] <- lnbc_roc_band_method(
      data = data,
      method = method,
      B = B,
      pi0 = pi_val,
      pi1 = pi_val,
      dist = dist,
      seed = method_seed,
      conf_level = conf_level,
      fpr_grid = fpr_grid,
      repo_root = repo_root,
      cores = cores,
      verbose = verbose
    )
  }

  band_rows <- lapply(results, function(res) {
    x <- res$band
    x$method <- res$method
    x
  })
  band_df <- do.call(rbind, band_rows)
  rownames(band_df) <- NULL
  band_df$dist <- dist
  band_df$n <- n
  band_df$J_key <- as.character(J_key)
  band_df$pi_key <- as.character(pi_key)
  band_df$pi <- pi_val
  band_df$B <- B
  band_df$rep <- rep_id
  band_df$data_seed <- data_seed
  band_df$conf_level <- conf_level
  band_df <- band_df[, c("method", "dist", "n", "J_key", "pi_key", "pi",
                         "B", "rep", "data_seed", "conf_level",
                         "fpr", "tpr", "ci_low", "ci_high",
                         "n_boot", "n_valid", "n_missing")]

  point_status <- do.call(rbind, lapply(results, function(res) {
    x <- res$point_status
    x$method <- res$method
    x
  }))
  rownames(point_status) <- NULL
  boot_status <- do.call(rbind, lapply(results, function(res) {
    x <- res$bootstrap_status
    x$method <- res$method
    x
  }))
  rownames(boot_status) <- NULL

  run_dir <- file.path(out_dir, paste0(mode, "_B", B))
  prefix <- file.path(
    run_dir,
    paste0("n", n, "J", J_key, "pi", pi_key, dist,
           "_rep", rep_id, "_", mode, "_B", B)
  )
  out <- list(
    call = list(methods = methods, dist = dist, n = n,
                J_key = as.character(J_key), pi_key = as.character(pi_key),
                B = B, rep_id = rep_id, seed = seed, data_seed = data_seed,
                conf_level = conf_level, cores = as.integer(cores)),
    truth = truth,
    true_roc = true_roc,
    scenario = lnbc_scenario(dist, J_key),
    data = data,
    results = results,
    band = band_df,
    point_status = point_status,
    bootstrap_status = boot_status,
    files = list(
      rds = paste0(prefix, "_roc.rds"),
      band_csv = paste0(prefix, "_roc_band.csv"),
      point_status_csv = paste0(prefix, "_roc_point_status.csv"),
      boot_status_csv = paste0(prefix, "_roc_boot_status.csv"),
      png = paste0(prefix, "_roc.png"),
      pdf = paste0(prefix, "_roc.pdf")
    )
  )

  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, out$files$rds)
  write.csv(band_df, out$files$band_csv, row.names = FALSE)
  write.csv(point_status, out$files$point_status_csv, row.names = FALSE)
  write.csv(boot_status, out$files$boot_status_csv, row.names = FALSE)

  if (isTRUE(save_plot)) {
    plot_title <- paste0("ROC bands: ", dist, ", n=", n, ", J=", J_key,
                         ", pi=", pi_key, ", B=", B, ", rep=", rep_id)
    lnbc_plot_roc_bands(
      band_df = band_df,
      true_roc = true_roc,
      png_file = out$files$png,
      pdf_file = out$files$pdf,
      title = plot_title
    )
  }

  cat("Saved ", out$files$rds, "\n", sep = "")
  cat("Saved ", out$files$band_csv, "\n", sep = "")
  if (isTRUE(save_plot)) {
    cat("Saved ", out$files$png, "\n", sep = "")
    cat("Saved ", out$files$pdf, "\n", sep = "")
  }
  out
}

lnbc_add_roc_truth <- function(band_df, true_roc) {
  if (!is.data.frame(band_df) || !"fpr" %in% names(band_df)) {
    stop("band_df must contain fpr.")
  }
  if (!is.data.frame(true_roc) || !all(c("fpr", "tpr") %in% names(true_roc))) {
    stop("true_roc must contain fpr and tpr.")
  }
  band_df$true_tpr <- approx(
    x = true_roc$fpr,
    y = true_roc$tpr,
    xout = band_df$fpr,
    ties = max,
    rule = 2
  )$y
  band_df$covered <- with(band_df, ifelse(
    is.finite(ci_low) & is.finite(ci_high) & is.finite(true_tpr),
    ci_low <= true_tpr & true_tpr <= ci_high,
    NA
  ))
  band_df
}

lnbc_by_rep_roc_file <- function(out_dir, method, dist, n, J_key, pi_key,
                                 B, total_reps, mode = "final") {
  file.path(
    out_dir,
    method,
    paste0("n", n, "J", J_key, "pi", pi_key, dist,
           "_", mode, "_B", B, "_reps", total_reps, "_roc_band.csv")
  )
}

lnbc_by_rep_roc_plot_file <- function(plot_dir, dist, n, J_key, pi_key,
                                      rep_id, B, total_reps, mode = "final",
                                      ext = "png") {
  file.path(
    plot_dir,
    paste0("n", n, "J", J_key, "pi", pi_key, dist,
           "_rep", rep_id, "_", mode, "_B", B,
           "_reps", total_reps, "_roc.", ext)
  )
}

lnbc_run_roc_replication_grid <- function(rep_id,
                                          methods = lnbc_methods,
                                          dists = c("lognorm", "gamma", "weibull", "gamma2", "weibull2"),
                                          n_vals = c(100, 300, 500),
                                          J_keys = c("3", "5"),
                                          pi_keys = c("75", "9"),
                                          B = 500,
                                          total_reps = 1000,
                                          seed = 123,
                                          conf_level = 0.95,
                                          fpr_grid = lnbc_default_roc_grid(),
                                          repo_root = lnbc_repo_root(),
                                          out_dir = NULL,
                                          plot_dir = NULL,
                                          mode = "final",
                                          cores = 1,
                                          append = TRUE,
                                          save_plot = FALSE,
                                          verbose = TRUE) {
  methods <- match.arg(methods, lnbc_methods, several.ok = TRUE)
  rep_id <- as.integer(rep_id)
  if (is.na(rep_id) || rep_id < 1) stop("rep_id must be a positive integer.")
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  data_seed <- lnbc_replication_seeds(rep_id)[rep_id]

  rows <- list()
  row_i <- 1

  for (dist in dists) {
    for (J_key in J_keys) {
      for (n in n_vals) {
        for (pi_key in pi_keys) {
          pi_val <- lnbc_pi_value(pi_key)
          true_roc <- lnbc_true_roc_curve(dist, J_key, fpr_grid)
          data <- lnbc_make_data(n = n, dist = dist, J_key = J_key,
                                 pi_key = pi_key, seed = data_seed)
          scenario_rows <- list()

          for (method in methods) {
            if (isTRUE(verbose)) {
              cat("rep", rep_id, "method", method, "dist", dist,
                  "n", n, "J", J_key, "pi", pi_key, "\n")
            }

            boot_seed <- lnbc_grid_boot_seed(
              rep_id = rep_id,
              method = method,
              dist = dist,
              n = n,
              J_key = J_key,
              pi_key = pi_key,
              B = B,
              seed = seed
            )
            res <- lnbc_roc_band_method(
              data = data,
              method = method,
              B = B,
              pi0 = pi_val,
              pi1 = pi_val,
              dist = dist,
              seed = boot_seed,
              conf_level = conf_level,
              fpr_grid = fpr_grid,
              repo_root = repo_root,
              cores = cores,
              verbose = FALSE
            )

            band <- lnbc_add_roc_truth(res$band, true_roc)
            boot_conv <- res$bootstrap_status$ok &
              res$bootstrap_status$conv &
              res$bootstrap_status$roc_ok
            band$rep <- rep_id
            band$data_seed <- data_seed
            band$boot_seed <- boot_seed
            band$method <- method
            band$dist <- dist
            band$n <- as.integer(n)
            band$J_key <- as.character(J_key)
            band$pi_key <- as.character(pi_key)
            band$pi <- pi_val
            band$B <- as.integer(B)
            band$total_reps <- as.integer(total_reps)
            band$conf_level <- conf_level
            band$n_boot_fail <- band$n_boot - band$n_valid
            band$point_ok <- res$point_status$ok[1]
            band$point_conv <- res$point_status$conv[1]
            band$point_roc_ok <- res$point_status$roc_ok[1]
            band$point_error <- res$point_status$error[1]
            band$n_boot_ok <- sum(boot_conv, na.rm = TRUE)
            band$boot_converge_prob <- mean(boot_conv, na.rm = TRUE)

            band <- band[, c("method", "dist", "n", "J_key", "pi_key", "pi",
                             "B", "total_reps", "rep", "data_seed", "boot_seed",
                             "conf_level", "fpr", "tpr", "ci_low", "ci_high",
                             "true_tpr", "covered", "n_boot", "n_valid",
                             "n_missing", "n_boot_fail", "point_ok",
                             "point_conv", "point_roc_ok", "point_error",
                             "n_boot_ok", "boot_converge_prob")]

            if (isTRUE(append) && !is.null(out_dir)) {
              path <- lnbc_by_rep_roc_file(
                out_dir = out_dir,
                method = method,
                dist = dist,
                n = n,
                J_key = J_key,
                pi_key = pi_key,
                B = B,
                total_reps = total_reps,
                mode = mode
              )
              lnbc_append_csv_locked(band, path)
            }

            rows[[row_i]] <- band
            scenario_rows[[method]] <- band
            row_i <- row_i + 1
          }

          if (isTRUE(save_plot) && length(scenario_rows) > 0) {
            if (is.null(plot_dir)) {
              plot_dir <- if (is.null(out_dir)) {
                file.path(repo_root, "ROC_plot", paste0(mode, "_B", B, "_reps", total_reps), "plots")
              } else {
                file.path(out_dir, "plots")
              }
            }
            scenario_df <- do.call(rbind, scenario_rows)
            plot_title <- paste0("ROC bands: ", dist, ", n=", n, ", J=", J_key,
                                 ", pi=", pi_key, ", B=", B, ", rep=", rep_id)
            png_file <- lnbc_by_rep_roc_plot_file(
              plot_dir = plot_dir,
              dist = dist,
              n = n,
              J_key = J_key,
              pi_key = pi_key,
              rep_id = rep_id,
              B = B,
              total_reps = total_reps,
              mode = mode,
              ext = "png"
            )
            pdf_file <- lnbc_by_rep_roc_plot_file(
              plot_dir = plot_dir,
              dist = dist,
              n = n,
              J_key = J_key,
              pi_key = pi_key,
              rep_id = rep_id,
              B = B,
              total_reps = total_reps,
              mode = mode,
              ext = "pdf"
            )
            lnbc_plot_roc_bands(
              band_df = scenario_df,
              true_roc = true_roc,
              png_file = png_file,
              pdf_file = pdf_file,
              title = plot_title
            )
          }
        }
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

lnbc_grid_boot_seed <- function(rep_id, method, dist, n, J_key, pi_key,
                                B, seed = 123) {
  # Pair these methods on the same point-fit seed and bootstrap samples. Their
  # only intended difference is whether k is re-estimated in a resample.
  seed_method <- if (identical(method, "EMBC_fixk")) "EMBC" else method
  method_idx <- match(seed_method, lnbc_methods)
  dist_idx <- match(dist, c("lognorm", "gamma", "weibull", "gamma2", "weibull2"))
  if (is.na(method_idx)) stop("Unknown method: ", method)
  if (is.na(dist_idx)) stop("Unknown distribution: ", dist)

  pi_num <- as.integer(round(100 * lnbc_pi_value(pi_key)))
  raw <- as.numeric(seed) +
    as.integer(rep_id) * 1000003 +
    method_idx * 100003 +
    dist_idx * 10007 +
    as.integer(n) * 101 +
    as.integer(J_key) * 1009 +
    pi_num * 13 +
    as.integer(B)

  as.integer((raw %% (.Machine$integer.max - 1)) + 1)
}

lnbc_by_rep_ci_file <- function(out_dir, method, dist, n, J_key, pi_key,
                                B, total_reps, mode = "final") {
  file.path(
    out_dir,
    method,
    paste0("n", n, "J", J_key, "pi", pi_key, dist,
           "_", mode, "_B", B, "_reps", total_reps, "_ci.csv")
  )
}

lnbc_append_csv_locked <- function(x, path, lock_timeout = 3600,
                                   lock_poll = 0.2) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(path, ".lock")
  start <- Sys.time()
  while (!dir.create(lock_dir, showWarnings = FALSE)) {
    if (as.numeric(difftime(Sys.time(), start, units = "secs")) > lock_timeout) {
      stop("Timed out waiting for lock: ", lock_dir)
    }
    Sys.sleep(lock_poll)
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)

  has_file <- file.exists(path) && file.info(path)$size > 0
  utils::write.table(
    x,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !has_file,
    append = has_file,
    qmethod = "double"
  )
  invisible(path)
}

lnbc_run_replication_grid <- function(rep_id,
                                      methods = lnbc_methods,
                                      dists = c("lognorm", "gamma", "weibull", "gamma2", "weibull2"),
                                      n_vals = c(100, 300, 500),
                                      J_keys = c("3", "5"),
                                      pi_keys = c("75", "9"),
                                      B = 500,
                                      total_reps = 1000,
                                      seed = 123,
                                      conf_level = 0.95,
                                      repo_root = lnbc_repo_root(),
                                      out_dir = NULL,
                                      mode = "final",
                                      cores = 1,
                                      append = TRUE,
                                      verbose = TRUE) {
  rep_id <- as.integer(rep_id)
  if (is.na(rep_id) || rep_id < 1) stop("rep_id must be a positive integer.")
  data_seed <- lnbc_replication_seeds(rep_id)[rep_id]

  rows <- list()
  row_i <- 1

  for (dist in dists) {
    for (J_key in J_keys) {
      for (n in n_vals) {
        for (pi_key in pi_keys) {
          pi_val <- lnbc_pi_value(pi_key)
          truth <- lnbc_true_values(dist, J_key)
          data <- lnbc_make_data(n = n, dist = dist, J_key = J_key,
                                 pi_key = pi_key, seed = data_seed)

          for (method in methods) {
            if (isTRUE(verbose)) {
              cat("rep", rep_id, "method", method, "dist", dist,
                  "n", n, "J", J_key, "pi", pi_key, "\n")
            }

            boot_seed <- lnbc_grid_boot_seed(
              rep_id = rep_id,
              method = method,
              dist = dist,
              n = n,
              J_key = J_key,
              pi_key = pi_key,
              B = B,
              seed = seed
            )
            res <- lnbc_bootstrap_method(
              data = data,
              method = method,
              B = B,
              pi0 = pi_val,
              pi1 = pi_val,
              dist = dist,
              seed = boot_seed,
              conf_level = conf_level,
              repo_root = repo_root,
              cores = cores,
              verbose = FALSE
            )

            ci <- lnbc_add_coverage(res$ci, truth)
            boot_conv <- res$bootstrap_status$ok & res$bootstrap_status$conv
            ci$rep <- rep_id
            ci$data_seed <- data_seed
            ci$boot_seed <- boot_seed
            ci$method <- method
            ci$dist <- dist
            ci$n <- as.integer(n)
            ci$J_key <- as.character(J_key)
            ci$pi_key <- as.character(pi_key)
            ci$pi <- pi_val
            ci$B <- as.integer(B)
            ci$total_reps <- as.integer(total_reps)
            ci$conf_level <- conf_level
            ci$n_boot_fail <- ci$n_boot - ci$n_valid
            ci$point_ok <- res$point_status$ok[1]
            ci$point_conv <- res$point_status$conv[1]
            ci$point_error <- res$point_status$error[1]
            ci$n_boot_ok <- sum(boot_conv, na.rm = TRUE)
            ci$boot_converge_prob <- mean(boot_conv, na.rm = TRUE)

            if (isTRUE(append) && !is.null(out_dir)) {
              path <- lnbc_by_rep_ci_file(
                out_dir = out_dir,
                method = method,
                dist = dist,
                n = n,
                J_key = J_key,
                pi_key = pi_key,
                B = B,
                total_reps = total_reps,
                mode = mode
              )
              lnbc_append_csv_locked(ci, path)
            }

            rows[[row_i]] <- ci
            row_i <- row_i + 1
          }
        }
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

lnbc_add_coverage <- function(ci_df, truth) {
  ci_df$true <- as.numeric(truth[ci_df$metric])
  ci_df$covered <- with(ci_df, ifelse(
    is.finite(estimate) & is.finite(ci_low) & is.finite(ci_high) & is.finite(true),
    ci_low <= true & true <= ci_high,
    NA
  ))
  ci_df
}

lnbc_coverage_summary <- function(ci_df) {
  if (nrow(ci_df) == 0) return(data.frame())
  mean_or_na <- function(x) {
    if (!any(!is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }
  groups <- split(ci_df, interaction(ci_df$method, ci_df$dist, ci_df$n, ci_df$J_key,
                                     ci_df$pi_key, ci_df$metric, drop = TRUE))
  rows <- lapply(groups, function(x) {
    x1 <- x[1, ]
    data.frame(
      method = x1$method,
      dist = x1$dist,
      n = x1$n,
      J_key = x1$J_key,
      pi_key = x1$pi_key,
      metric = x1$metric,
      n_reps = nrow(x),
      n_covered_available = sum(!is.na(x$covered)),
      coverage = mean_or_na(x$covered),
      avg_length = mean_or_na(x$ci_length),
      avg_n_valid_boot = mean_or_na(x$n_valid),
      total_boot_fail = sum(x$n_boot - x$n_valid, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

run_bootstrap_study <- function(method, dist, n, J_key, pi_key,
                                B = 5, data_reps = 1, seed = 123,
                                conf_level = 0.95,
                                repo_root = lnbc_repo_root(),
                                out_file = NULL,
                                cores = 1,
                                verbose = TRUE) {
  method <- match.arg(method, lnbc_methods)
  data_reps <- as.integer(data_reps)
  if (is.na(data_reps) || data_reps < 1) {
    stop("data_reps must be a positive integer.")
  }
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1) cores <- 1
  if (.Platform$OS.type != "unix") cores <- 1
  cores <- min(cores, data_reps)

  pi_val <- lnbc_pi_value(pi_key)
  truth <- lnbc_true_values(dist, J_key)

  data_seeds <- lnbc_replication_seeds(data_reps)
  set.seed(seed)
  boot_study_seeds <- sample.int(.Machine$integer.max, data_reps)

  run_replication <- function(r) {
    if (isTRUE(verbose)) {
      cat("Data replication", r, "of", data_reps, "\n")
    }
    data <- lnbc_make_data(n = n, dist = dist, J_key = J_key,
                           pi_key = pi_key, seed = data_seeds[r])
    res <- lnbc_bootstrap_method(
      data = data,
      method = method,
      B = B,
      pi0 = pi_val,
      pi1 = pi_val,
      dist = dist,
      seed = boot_study_seeds[r],
      conf_level = conf_level,
      repo_root = repo_root,
      cores = 1,
      verbose = isTRUE(verbose) && cores == 1
    )

    ci <- lnbc_add_coverage(res$ci, truth)
    ci$rep <- r
    ci$method <- method
    ci$dist <- dist
    ci$n <- as.integer(n)
    ci$J_key <- as.character(J_key)
    ci$pi_key <- as.character(pi_key)
    ci$pi <- pi_val
    ci$n_boot_fail <- ci$n_boot - ci$n_valid
    ci$point_ok <- res$point_status$ok[1]
    ci$point_error <- res$point_status$error[1]

    list(replication = res, ci = ci)
  }

  if (cores > 1) {
    if (isTRUE(verbose)) {
      cat("Running", data_reps, "data replications with", cores, "cores\n")
    }
    rep_results <- parallel::mclapply(seq_len(data_reps), run_replication,
                                      mc.cores = cores, mc.preschedule = FALSE)
  } else {
    rep_results <- lapply(seq_len(data_reps), run_replication)
  }

  reps <- lapply(rep_results, `[[`, "replication")
  ci_rows <- lapply(rep_results, `[[`, "ci")

  ci_summary <- do.call(rbind, ci_rows)
  rownames(ci_summary) <- NULL
  coverage_summary <- lnbc_coverage_summary(ci_summary)

  out <- list(
    call = list(method = method, dist = dist, n = as.integer(n),
                J_key = as.character(J_key), pi_key = as.character(pi_key),
                B = as.integer(B), data_reps = data_reps, seed = seed,
                conf_level = conf_level, cores = as.integer(cores)),
    truth = truth,
    scenario = lnbc_scenario(dist, J_key),
    data_seeds = data_seeds,
    boot_study_seeds = boot_study_seeds,
    replications = reps,
    ci_summary = ci_summary,
    coverage_summary = coverage_summary
  )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, out_file)
    write.csv(ci_summary, sub("\\.rds$", "_ci.csv", out_file), row.names = FALSE)
    write.csv(coverage_summary, sub("\\.rds$", "_coverage.csv", out_file), row.names = FALSE)
    cat("Saved", out_file, "\n")
  }

  out
}

## Unified bootstrap workflow -------------------------------------------------

lnbc_source_method <- function(method, repo_root = lnbc_repo_root()) {
  method <- match.arg(method, lnbc_methods)

  method_file <- function(...) {
    candidates <- c(
      file.path(repo_root, "methods", ...),
      file.path(repo_root, ...)
    )
    candidates[file.exists(candidates)][1]
  }

  if (method %in% c("EMBC", "EMBC_fixk") && !isTRUE(.lnbc_bootstrap_state$EMBC)) {
    path <- method_file("EMBC", "method.R")
    if (is.na(path)) stop("Cannot find EMBC method.R.")
    source(path, local = .GlobalEnv)
    .lnbc_bootstrap_state$EMBC <- TRUE
  }
  if (method %in% c("NP", "NP-inc") && !isTRUE(.lnbc_bootstrap_state$NP)) {
    path <- method_file("NP", "method_np.r")
    if (is.na(path)) stop("Cannot find NP method_np.r.")
    source(path, local = .GlobalEnv)
    .lnbc_bootstrap_state$NP <- TRUE
  }
  if (method %in% c("LBNP", "LBNP-log") && !isTRUE(.lnbc_bootstrap_state$LBNP)) {
    path <- method_file("LBNP", "method_lbnp.r")
    if (is.na(path)) stop("Cannot find LBNP method_lbnp.r.")
    source(path, local = .GlobalEnv)
    .lnbc_bootstrap_state$LBNP <- TRUE
  }
  invisible(TRUE)
}

lnbc_roc_metric_names <- function(fpr_grid) {
  paste0("ROC_", sprintf("%.2f", lnbc_check_roc_grid(fpr_grid)))
}

lnbc_all_metric_names <- function(fpr_grid = NULL) {
  if (is.null(fpr_grid)) return(lnbc_metrics)
  c(lnbc_metrics, lnbc_roc_metric_names(fpr_grid))
}

lnbc_estimate_with_roc <- function(fit, fpr_grid = NULL) {
  out <- fit$estimate[lnbc_metrics]
  if (!is.null(fpr_grid)) {
    roc_names <- lnbc_roc_metric_names(fpr_grid)
    # A failed fit must remain missing. Re-sanitizing an empty curve can use
    # its forced endpoints (0, 0) and (1, 1) to create an artificial diagonal.
    if (!isTRUE(fit$ok) || is.null(fit$roc) || !all(is.finite(fit$roc$tpr))) {
      out <- c(out, setNames(rep(NA_real_, length(roc_names)), roc_names))
    } else {
      roc <- lnbc_sanitize_roc_curve(fit$roc, fpr_grid)
      out <- c(out, setNames(roc$tpr, roc_names))
    }
  }
  out
}

lnbc_truth_vector <- function(dist, J_key, fpr_grid = NULL) {
  out <- lnbc_true_values(dist, J_key)[lnbc_metrics]
  if (!is.null(fpr_grid)) {
    true_roc <- lnbc_true_roc_curve(dist, J_key, fpr_grid)
    out <- c(out, setNames(true_roc$tpr, lnbc_roc_metric_names(fpr_grid)))
  }
  out
}

lnbc_percentile_ci <- function(point_est, boot_mat, conf_level = 0.95,
                               log_cutoff = TRUE,
                               metrics = colnames(as.matrix(boot_mat))) {
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.")
  }
  alpha <- (1 - conf_level) / 2
  probs <- c(alpha, 1 - alpha)
  boot_mat <- as.matrix(boot_mat)
  metrics <- intersect(as.character(metrics), colnames(boot_mat))

  rows <- lapply(metrics, function(metric) {
    estimate <- unname(point_est[metric])
    vals <- as.numeric(boot_mat[, metric])
    usable <- is.finite(vals)
    transformation <- "identity"

    if (metric == "opt_cut" && isTRUE(log_cutoff)) {
      usable <- usable & vals > 0
      vals <- log(vals[usable])
      transformation <- "log"
    } else {
      vals <- vals[usable]
    }

    if (!is.finite(estimate) || length(vals) == 0) {
      ci <- c(NA_real_, NA_real_)
    } else {
      ci <- as.numeric(quantile(vals, probs = probs, names = FALSE, na.rm = FALSE))
      if (transformation == "log") ci <- exp(ci)
    }

    data.frame(
      metric = metric,
      estimate = estimate,
      ci_low = ci[1],
      ci_high = ci[2],
      ci_length = ci[2] - ci[1],
      n_boot = nrow(boot_mat),
      n_valid = length(vals),
      n_missing = nrow(boot_mat) - length(vals),
      transformation = transformation,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lnbc_default_ci_metrics <- function(fpr_grid = lnbc_default_roc_grid()) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  c("AUC", "opt_cut", "J", "se", "sp", "ROC_0.10", "ROC_0.20")[
    c("AUC", "opt_cut", "J", "se", "sp", "ROC_0.10", "ROC_0.20") %in%
      lnbc_all_metric_names(fpr_grid)
  ]
}

lnbc_bootstrap_method <- function(data, method, B = 5, pi0, pi1, dist,
                                  seed = 123, conf_level = 0.95,
                                  repo_root = lnbc_repo_root(),
                                  lbnp_k = 10,
                                  lbnp_nu_seq = 3 * 10^seq(from = -2, to = 2.5, length.out = 10),
                                  lbnp_maxit = 500,
                                  lbnp_thres = 1e-5,
                                  lbnp_transform = NULL,
                                  roc_grid = NULL,
                                  ci_metrics = NULL,
                                  cores = 1,
                                  verbose = TRUE) {
  method <- match.arg(method, lnbc_methods)
  B <- as.integer(B)
  if (is.na(B) || B < 1) stop("B must be a positive integer.")
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1) cores <- 1
  if (.Platform$OS.type != "unix") cores <- 1
  cores <- min(cores, B)
  if (!is.null(roc_grid)) roc_grid <- lnbc_check_roc_grid(roc_grid)
  metric_names <- lnbc_all_metric_names(roc_grid)
  if (is.null(ci_metrics)) ci_metrics <- metric_names

  data <- lnbc_check_data(data)
  set.seed(seed)
  point_seed <- sample.int(.Machine$integer.max, 1)
  boot_seeds <- sample.int(.Machine$integer.max, B)
  fit_seeds <- sample.int(.Machine$integer.max, B)

  if (isTRUE(verbose)) cat("Fitting point estimate:", method, "\n")
  point_fit <- lnbc_fit_method_once(
    method, data, pi0 = pi0, pi1 = pi1, dist = dist,
    seed = point_seed, repo_root = repo_root, lbnp_k = lbnp_k,
    lbnp_nu_seq = lbnp_nu_seq, lbnp_maxit = lbnp_maxit,
    lbnp_thres = lbnp_thres, lbnp_transform = lbnp_transform,
    roc_grid = roc_grid
  )
  point <- lnbc_estimate_with_roc(point_fit, roc_grid)
  embc_fixed_k <- if (method == "EMBC_fixk") {
    if (isTRUE(point_fit$ok) && length(point_fit$embc_k) == 1 &&
        is.finite(point_fit$embc_k)) point_fit$embc_k else NA_real_
  } else NULL

  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(metric_names),
                     dimnames = list(paste0("boot_", seq_len(B)), metric_names))
  boot_status <- data.frame(
    boot = seq_len(B),
    ok = rep(FALSE, B),
    conv = rep(FALSE, B),
    roc_ok = if (is.null(roc_grid)) NA else rep(FALSE, B),
    nu = rep(NA_real_, B),
    embc_k = rep(NA_real_, B),
    error = rep(NA_character_, B),
    stringsAsFactors = FALSE
  )

  fit_boot <- function(b) {
    if (isTRUE(verbose)) cat("Bootstrap", b, "of", B, "\n")
    boot_data <- lnbc_stratified_boot_sample(data, seed = boot_seeds[b])
    fit <- lnbc_fit_method_once(
      method, boot_data, pi0 = pi0, pi1 = pi1, dist = dist,
      seed = fit_seeds[b], repo_root = repo_root, lbnp_k = lbnp_k,
      lbnp_nu_seq = lbnp_nu_seq, lbnp_maxit = lbnp_maxit,
      lbnp_thres = lbnp_thres, lbnp_transform = lbnp_transform,
      embc_fixed_k = embc_fixed_k,
      roc_grid = roc_grid
    )
    list(b = b, fit = fit, estimate = lnbc_estimate_with_roc(fit, roc_grid))
  }

  if (cores > 1) {
    if (isTRUE(verbose)) cat("Running", B, "bootstrap fits with", cores, "cores\n")
    boot_fits <- parallel::mclapply(seq_len(B), fit_boot, mc.cores = cores,
                                    mc.preschedule = FALSE)
  } else {
    boot_fits <- lapply(seq_len(B), fit_boot)
  }

  for (boot_result in boot_fits) {
    b <- boot_result$b
    fit <- boot_result$fit
    boot_mat[b, ] <- boot_result$estimate
    boot_status$ok[b] <- isTRUE(fit$ok)
    boot_status$conv[b] <- isTRUE(fit$conv)
    if (!is.null(roc_grid)) {
      boot_status$roc_ok[b] <- !is.null(fit$roc) && all(is.finite(fit$roc$tpr))
    }
    boot_status$nu[b] <- fit$nu
    boot_status$embc_k[b] <- if (is.null(fit$embc_k)) NA_real_ else fit$embc_k
    boot_status$error[b] <- fit$error
  }

  ci <- lnbc_percentile_ci(point, boot_mat, conf_level = conf_level,
                           metrics = ci_metrics)
  list(
    method = method,
    dist = dist,
    pi0 = pi0,
    pi1 = pi1,
    B = B,
    conf_level = conf_level,
    n0 = sum(data$R == 0),
    n1 = sum(data$R == 1),
    fpr_grid = roc_grid,
    point = point,
    point_status = data.frame(
      ok = point_fit$ok,
      conv = point_fit$conv,
      roc_ok = if (is.null(roc_grid)) NA else !is.null(point_fit$roc) && all(is.finite(point_fit$roc$tpr)),
      nu = point_fit$nu,
      embc_k = if (is.null(point_fit$embc_k)) NA_real_ else point_fit$embc_k,
      error = point_fit$error,
      stringsAsFactors = FALSE
    ),
    point_roc = point_fit$roc,
    bootstrap = boot_mat,
    bootstrap_status = boot_status,
    ci = ci,
    seeds = list(point_seed = point_seed, boot_seeds = boot_seeds, fit_seeds = fit_seeds)
  )
}

lnbc_logit <- function(x, eps = 1e-6) {
  x <- pmin(pmax(x, eps), 1 - eps)
  log(x / (1 - x))
}

lnbc_try_solve <- function(mat, ridge = 1e-8) {
  inv <- tryCatch(solve(mat), error = function(e) NULL)
  if (!is.null(inv)) return(list(inv = inv, ridge = 0, ok = TRUE))
  inv <- tryCatch(solve(mat + diag(ridge, nrow(mat))), error = function(e) NULL)
  if (!is.null(inv)) return(list(inv = inv, ridge = ridge, ok = TRUE))
  list(inv = matrix(NA_real_, nrow(mat), ncol(mat)), ridge = ridge, ok = FALSE)
}

lnbc_joint_logit_region <- function(point_est, boot_mat, truth,
                                    conf_level = 0.95, mc_draws = 10000,
                                    seed = NULL, eps = 1e-6) {
  required <- c("se", "sp")
  point_pair <- as.numeric(point_est[required])
  true_pair <- as.numeric(truth[required])
  boot_pair <- as.matrix(boot_mat[, required, drop = FALSE])
  ok <- is.finite(boot_pair[, 1]) & is.finite(boot_pair[, 2]) &
    boot_pair[, 1] > 0 & boot_pair[, 1] < 1 &
    boot_pair[, 2] > 0 & boot_pair[, 2] < 1
  n_valid <- sum(ok)
  crit <- stats::qchisq(conf_level, df = 2)

  fail <- function(status) {
    data.frame(
      coverage = NA,
      area = NA_real_,
      n_boot = nrow(boot_mat),
      n_valid = n_valid,
      mc_draws = as.integer(mc_draws),
      chi_sq_cutoff = crit,
      ridge = NA_real_,
      status = status,
      stringsAsFactors = FALSE
    )
  }

  if (!all(is.finite(point_pair)) || !all(is.finite(true_pair)) || n_valid < 3) {
    return(fail("insufficient finite se/sp bootstrap pairs"))
  }

  g_boot <- cbind(lnbc_logit(boot_pair[ok, 1], eps), lnbc_logit(boot_pair[ok, 2], eps))
  v_boot <- stats::cov(g_boot)
  if (!all(is.finite(v_boot))) return(fail("non-finite covariance"))
  inv <- lnbc_try_solve(v_boot)
  if (!isTRUE(inv$ok)) return(fail("singular covariance"))

  center <- lnbc_logit(point_pair, eps)
  qfun <- function(mu) {
    diff <- sweep(cbind(lnbc_logit(mu[, 1], eps), lnbc_logit(mu[, 2], eps)),
                  2, center, "-")
    rowSums((diff %*% inv$inv) * diff)
  }

  true_q <- qfun(matrix(true_pair, nrow = 1))
  if (!is.null(seed)) set.seed(seed)
  draws <- cbind(stats::runif(mc_draws), stats::runif(mc_draws))
  area <- mean(qfun(draws) <= crit, na.rm = TRUE)

  data.frame(
    coverage = ifelse(is.finite(true_q), true_q <= crit, NA),
    area = area,
    n_boot = nrow(boot_mat),
    n_valid = n_valid,
    mc_draws = as.integer(mc_draws),
    chi_sq_cutoff = crit,
    ridge = inv$ridge,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

lnbc_add_coverage <- function(ci_df, truth) {
  ci_df$true <- as.numeric(truth[ci_df$metric])
  ci_df$covered <- with(ci_df, ifelse(
    is.finite(estimate) & is.finite(ci_low) & is.finite(ci_high) & is.finite(true),
    ci_low <= true & true <= ci_high,
    NA
  ))
  ci_df
}

lnbc_point_rows <- function(res, truth) {
  data.frame(
    metric = names(res$point),
    estimate = as.numeric(res$point),
    true = as.numeric(truth[names(res$point)]),
    stringsAsFactors = FALSE
  )
}

lnbc_roc_rows <- function(res, truth, fpr_grid) {
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  metric <- lnbc_roc_metric_names(fpr_grid)
  band <- lnbc_roc_percentile_band(res$point_roc, res$bootstrap[, metric, drop = FALSE],
                                   conf_level = res$conf_level)
  band$metric <- metric
  band$true_tpr <- as.numeric(truth[metric])
  band$covered <- with(band, ifelse(
    is.finite(ci_low) & is.finite(ci_high) & is.finite(true_tpr),
    ci_low <= true_tpr & true_tpr <= ci_high,
    NA
  ))
  band[, c("metric", "fpr", "tpr", "ci_low", "ci_high", "true_tpr",
           "covered", "n_boot", "n_valid", "n_missing")]
}

lnbc_result_file <- function(out_dir, family, method, dist, n, J_key, pi_key,
                             B, total_reps, mode = "final", suffix) {
  file.path(
    out_dir,
    family,
    method,
    paste0("n", n, "J", J_key, "pi", pi_key, dist,
           "_", mode, "_B", B, "_reps", total_reps, suffix)
  )
}

lnbc_run_unified_replication_grid <- function(rep_id,
                                              methods = lnbc_methods,
                                              dists = c("lognorm", "gamma", "weibull", "gamma2", "weibull2"),
                                              n_vals = c(100, 300, 500),
                                              J_keys = c("3", "5"),
                                              pi_keys = c("75", "9"),
                                              B = 500,
                                              total_reps = 1000,
                                              seed = 123,
                                              conf_level = 0.95,
                                              fpr_grid = lnbc_default_roc_grid(),
                                              mc_draws = 10000,
                                              repo_root = lnbc_repo_root(),
                                              out_dir = file.path(repo_root, "result"),
                                              mode = "final",
                                              cores = 1,
                                              append = TRUE,
                                              save_rds = FALSE,
                                              verbose = TRUE) {
  methods <- match.arg(methods, lnbc_methods, several.ok = TRUE)
  rep_id <- as.integer(rep_id)
  if (is.na(rep_id) || rep_id < 1) stop("rep_id must be a positive integer.")
  fpr_grid <- lnbc_check_roc_grid(fpr_grid)
  data_seed <- lnbc_replication_seeds(rep_id)[rep_id]

  rows <- list(point = list(), coverage = list(), joint = list(), roc = list())
  row_i <- c(point = 1L, coverage = 1L, joint = 1L, roc = 1L)

  for (dist in dists) {
    for (J_key in J_keys) {
      for (n in n_vals) {
        for (pi_key in pi_keys) {
          pi_val <- lnbc_pi_value(pi_key)
          truth <- lnbc_truth_vector(dist, J_key, fpr_grid)
          data <- lnbc_make_data(n = n, dist = dist, J_key = J_key,
                                 pi_key = pi_key, seed = data_seed)

          for (method in methods) {
            if (isTRUE(verbose)) {
              cat("rep", rep_id, "method", method, "dist", dist,
                  "n", n, "J", J_key, "pi", pi_key, "\n")
            }

            boot_seed <- lnbc_grid_boot_seed(
              rep_id = rep_id, method = method, dist = dist, n = n,
              J_key = J_key, pi_key = pi_key, B = B, seed = seed
            )
            res <- lnbc_bootstrap_method(
              data = data,
              method = method,
              B = B,
              pi0 = pi_val,
              pi1 = pi_val,
              dist = dist,
              seed = boot_seed,
              conf_level = conf_level,
              repo_root = repo_root,
              roc_grid = fpr_grid,
              ci_metrics = lnbc_default_ci_metrics(fpr_grid),
              cores = cores,
              verbose = FALSE
            )
            boot_conv <- res$bootstrap_status$ok & res$bootstrap_status$conv
            if ("roc_ok" %in% names(res$bootstrap_status)) {
              boot_conv <- boot_conv & (is.na(res$bootstrap_status$roc_ok) | res$bootstrap_status$roc_ok)
            }
            common <- list(
              method = method, dist = dist, n = as.integer(n),
              J_key = as.character(J_key), pi_key = as.character(pi_key),
              pi = pi_val, B = as.integer(B), total_reps = as.integer(total_reps),
              rep = rep_id, data_seed = data_seed, boot_seed = boot_seed,
              conf_level = conf_level, point_ok = res$point_status$ok[1],
              point_conv = res$point_status$conv[1],
              point_roc_ok = res$point_status$roc_ok[1],
              point_error = res$point_status$error[1],
              n_boot_ok = sum(boot_conv, na.rm = TRUE),
              boot_converge_prob = mean(boot_conv, na.rm = TRUE)
            )

            point <- lnbc_point_rows(res, truth)
            for (nm in names(common)) point[[nm]] <- common[[nm]]
            point <- point[, c(names(common), "metric", "estimate", "true")]

            ci <- lnbc_add_coverage(res$ci, truth)
            ci$n_boot_fail <- ci$n_boot - ci$n_valid
            for (nm in names(common)) ci[[nm]] <- common[[nm]]
            ci <- ci[, c(names(common), "metric", "estimate", "ci_low", "ci_high",
                         "ci_length", "true", "covered", "n_boot", "n_valid",
                         "n_missing", "n_boot_fail", "transformation")]

            joint <- lnbc_joint_logit_region(
              point_est = res$point,
              boot_mat = res$bootstrap,
              truth = truth,
              conf_level = conf_level,
              mc_draws = mc_draws,
              seed = boot_seed + 991L
            )
            for (nm in names(common)) joint[[nm]] <- common[[nm]]
            joint <- joint[, c(names(common), "coverage", "area", "n_boot",
                               "n_valid", "mc_draws", "chi_sq_cutoff",
                               "ridge", "status")]

            roc <- lnbc_roc_rows(res, truth, fpr_grid)
            roc$n_boot_fail <- roc$n_boot - roc$n_valid
            for (nm in names(common)) roc[[nm]] <- common[[nm]]
            roc <- roc[, c(names(common), "metric", "fpr", "tpr", "ci_low",
                           "ci_high", "true_tpr", "covered", "n_boot",
                           "n_valid", "n_missing", "n_boot_fail")]

            if (isTRUE(append) && !is.null(out_dir)) {
              lnbc_append_csv_locked(
                point,
                lnbc_result_file(out_dir, "point_est", method, dist, n, J_key,
                                 pi_key, B, total_reps, mode, "_point.csv")
              )
              lnbc_append_csv_locked(
                ci,
                lnbc_result_file(out_dir, "coverage", method, dist, n, J_key,
                                 pi_key, B, total_reps, mode, "_ci.csv")
              )
              lnbc_append_csv_locked(
                joint,
                lnbc_result_file(out_dir, "coverage_joint", method, dist, n, J_key,
                                 pi_key, B, total_reps, mode, "_joint.csv")
              )
              lnbc_append_csv_locked(
                roc,
                lnbc_result_file(out_dir, "ROC_plot", method, dist, n, J_key,
                                 pi_key, B, total_reps, mode, "_roc.csv")
              )
              if (isTRUE(save_rds)) {
                rds_file <- lnbc_result_file(out_dir, "raw", method, dist, n, J_key,
                                             pi_key, B, total_reps, mode,
                                             paste0("_rep", rep_id, ".rds"))
                dir.create(dirname(rds_file), recursive = TRUE, showWarnings = FALSE)
                saveRDS(res, rds_file)
              }
            }

            rows$point[[as.integer(row_i["point"])]] <- point
            rows$coverage[[as.integer(row_i["coverage"])]] <- ci
            rows$joint[[as.integer(row_i["joint"])]] <- joint
            rows$roc[[as.integer(row_i["roc"])]] <- roc
            row_i <- row_i + c(point = 1L, coverage = 1L, joint = 1L, roc = 1L)
          }
        }
      }
    }
  }

  lapply(rows, function(x) {
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
  })
}
