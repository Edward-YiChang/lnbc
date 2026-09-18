# Main proposed NP ROC method from the NP paper code.
# The original NP folder contains simulation scripts; this file exposes the
# single-biomarker proposed estimator without running those simulations.

.np_check_inputs <- function(dat, pi0, pi1, scale01 = FALSE) {
  if (!is.data.frame(dat) || !all(c("R", "biomarker") %in% names(dat))) {
    stop("dat must be a data.frame with columns R and biomarker.")
  }
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }

  R <- dat$R
  Y <- dat$biomarker
  if (!all(R %in% c(0, 1))) {
    stop("R must contain only 0/1 labels.")
  }
  if (sum(R == 0) == 0 || sum(R == 1) == 0) {
    stop("R must contain at least one 0 and one 1.")
  }
  if (!all(is.finite(Y))) {
    stop("biomarker must be finite.")
  }

  if (isTRUE(scale01)) {
    y_range <- range(Y)
    if (diff(y_range) == 0) {
      stop("biomarker cannot be scaled because all values are equal.")
    }
    Y <- (Y - y_range[1]) / diff(y_range)
  } else if (min(Y) < 0 || max(Y) > 1) {
    stop("np.roc assumes biomarker is scaled to [0, 1]; set scale01 = TRUE or scale before calling.")
  }

  data.frame(R = R, biomarker = Y)
}

.np_trapz <- function(x, y) {
  idx <- 2:length(x)
  as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1])) / 2
}

np.auc <- function(dat, pi0 = 1, pi1 = 1, scale01 = FALSE) {
  dat <- .np_check_inputs(dat, pi0, pi1, scale01)
  R <- dat$R
  Y <- dat$biomarker

  apparent_auc <- mean(outer(
    Y[R == 1],
    Y[R == 0],
    FUN = function(y1, y0) (y1 > y0) + 0.5 * (y1 == y0)
  ))

  auc <- (apparent_auc + (pi1 + pi0) / 2 - 1) / (pi1 + pi0 - 1)
  c(AUC = auc, apparent_AUC = apparent_auc)
}

np.roc <- function(dat, pi0 = 1, pi1 = 1, t = 0.2,
                   s.range = c(0.1, 0.3), scale01 = FALSE,
                   paper_youden = FALSE) {
  dat <- .np_check_inputs(dat, pi0, pi1, scale01)
  R <- dat$R
  Y <- dat$biomarker
  n <- length(R)
  denom <- pi1 + pi0 - 1

  auc <- np.auc(dat, pi0 = pi0, pi1 = pi1)

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

  Q1 <- Q0 <- numeric(length(F_points))
  for (i in seq_along(F_points)) {
    Q1[i] <- sum((F1 <= F_points[i]) * delta_EP)
    Q0[i] <- sum((F0 <= F_points[i]) * delta_EP)
  }

  q0 <- sum((F0 <= (1 - t)) * delta_EP)
  q1 <- 1 - F_points[max(c(which(Q1 <= q0), 1))]

  t_eval <- sort(unique(c(Q0, Q1)))
  delta_F <- numeric(length(t_eval))
  for (i in seq_along(t_eval)) {
    q0_i <- F_points[max(c(which(Q0 <= t_eval[i]), 1))]
    q1_i <- F_points[max(c(which(Q1 <= t_eval[i]), 1))]
    delta_F[i] <- q0_i - q1_i
  }

  ROC.fun <- function(s) {
    tem <- sum((F0 <= (1 - s)) * delta_EP)
    1 - F_points[max(c(which(Q1 <= tem), 1))]
  }
  pAUC <- try(
    integrate(Vectorize(ROC.fun), s.range[1], s.range[2])$value /
      (s.range[2] - s.range[1]),
    silent = TRUE
  )
  if (inherits(pAUC, "try-error")) {
    x_eval <- seq(from = s.range[1], to = s.range[2], length.out = 1000)
    pAUC <- .np_trapz(x_eval, Vectorize(ROC.fun)(x_eval)) /
      (s.range[2] - s.range[1])
  }

  cutoff <- np.cutoff.di(dat, pi0 = pi0, pi1 = pi1)
  G1_cut <- sum(Y <= cutoff$cutoff & R == 1) / sum(R == 1)
  G0_cut <- sum(Y <= cutoff$cutoff & R == 0) / sum(R == 0)
  F1_cut <- pi0 / denom * G1_cut - (1 - pi1) / denom * G0_cut
  F0_cut <- pi1 / denom * G0_cut - (1 - pi0) / denom * G1_cut
  se <- 1 - F1_cut
  sp <- F0_cut

  J_est <- if (isTRUE(paper_youden)) max(delta_F) else se + sp - 1

  data.frame(
    ROC = q1,
    AUC = unname(auc["AUC"]),
    apparent_AUC = unname(auc["apparent_AUC"]),
    cutoff = cutoff$cutoff,
    opt_cut = cutoff$cutoff,
    Youden = cutoff$Youden,
    J = J_est,
    se = se,
    sp = sp,
    pAUC = pAUC,
    n = n
  )
}

np.cutoff.di <- function(dat, pi0 = 1, pi1 = 1, scale01 = FALSE) {
  dat <- .np_check_inputs(dat, pi0, pi1, scale01)
  R <- dat$R
  Y <- dat$biomarker

  F0 <- ecdf(Y[R == 0])
  F1 <- ecdf(Y[R == 1])
  eval <- sort(unique(Y))
  delta_F <- F0(eval) - F1(eval)

  data.frame(
    cutoff = eval[which.max(delta_F)],
    Youden = max(delta_F) / (pi0 + pi1 - 1)
  )
}

# Increasing-rearrangement NP estimator from Section 4.2 of Sun et al.
#
# The plug-in estimates of F1 and F0 are differences of empirical CDFs and
# need not be increasing in a finite sample.  On the unit biomarker domain,
# Section 4.2 rearranges each estimate through
#   Q(y) = integral I{F(u) <= y} du,  F+(t) = inf{y: Q(y) >= t}.
# The implementation below evaluates these integrals exactly over the
# empirical-CDF intervals.  Clipping commutes with increasing rearrangement
# and makes the finite-sample result a proper CDF on [0, 1].
.np_inc_rearrangement <- function(values, interval_weights) {
  keep <- is.finite(values) & is.finite(interval_weights) & interval_weights > 0
  values <- pmin(pmax(values[keep], 0), 1)
  interval_weights <- interval_weights[keep]
  if (length(values) == 0 || sum(interval_weights) <= 0) {
    stop("Cannot rearrange an empty CDF estimate.")
  }
  interval_weights <- interval_weights / sum(interval_weights)
  ord <- order(values)
  sorted_values <- values[ord]
  sorted_weights <- interval_weights[ord]
  cumulative_weights <- cumsum(sorted_weights)

  Q <- function(y) {
    y <- as.numeric(y)
    vapply(y, function(z) {
      if (!is.finite(z)) return(NA_real_)
      if (z < 0) return(0)
      if (z >= 1) return(1)
      sum(interval_weights[values <= z])
    }, numeric(1))
  }
  F_plus <- function(t) {
    t <- as.numeric(t)
    vapply(t, function(p) {
      if (!is.finite(p)) return(NA_real_)
      if (p <= 0) return(0)
      if (p >= 1) return(1)
      sorted_values[which(cumulative_weights >= p)[1]]
    }, numeric(1))
  }

  list(
    Q = Q,
    F_plus = F_plus,
    values = values,
    weights = interval_weights,
    knots = sort(unique(c(0, cumulative_weights, 1)))
  )
}

.np_inc_components <- function(dat, pi0 = 1, pi1 = 1, scale01 = FALSE) {
  dat <- .np_check_inputs(dat, pi0, pi1, scale01)
  R <- dat$R
  Y <- dat$biomarker
  denom <- pi0 + pi1 - 1

  # Empirical CDFs are constant between these knots. Values at individual
  # knots have Lebesgue measure zero in the rearrangement integral.
  domain_knots <- sort(unique(c(0, Y[Y > 0 & Y < 1], 1)))
  interval_weights <- diff(domain_knots)
  interval_midpoints <- (domain_knots[-length(domain_knots)] +
                           domain_knots[-1]) / 2
  G1 <- vapply(interval_midpoints, function(u) mean(Y[R == 1] <= u), numeric(1))
  G0 <- vapply(interval_midpoints, function(u) mean(Y[R == 0] <= u), numeric(1))
  raw_F1 <- (pi0 * G1 - (1 - pi1) * G0) / denom
  raw_F0 <- (pi1 * G0 - (1 - pi0) * G1) / denom

  rearranged1 <- .np_inc_rearrangement(raw_F1, interval_weights)
  rearranged0 <- .np_inc_rearrangement(raw_F0, interval_weights)
  roc <- function(s) {
    s <- pmin(pmax(as.numeric(s), 0), 1)
    out <- 1 - rearranged1$F_plus(rearranged0$Q(1 - s))
    out[s <= 0] <- 0
    out[s >= 1] <- 1
    pmin(pmax(out, 0), 1)
  }

  list(
    dat = dat,
    raw_F1 = raw_F1,
    raw_F0 = raw_F0,
    interval_weights = interval_weights,
    F1_plus = rearranged1$F_plus,
    F0_plus = rearranged0$F_plus,
    Q1 = rearranged1$Q,
    Q0 = rearranged0$Q,
    roc = roc,
    cutoff_knots = sort(unique(c(rearranged1$knots, rearranged0$knots)))
  )
}

np.inc.roc.curve <- function(dat, pi0 = 1, pi1 = 1,
                             fpr = seq(0, 1, by = 0.01),
                             scale01 = FALSE) {
  components <- .np_inc_components(dat, pi0, pi1, scale01)
  data.frame(fpr = as.numeric(fpr), tpr = components$roc(fpr))
}

np.inc.roc <- function(dat, pi0 = 1, pi1 = 1, t = 0.2,
                       s.range = c(0.1, 0.3), scale01 = FALSE) {
  components <- .np_inc_components(dat, pi0, pi1, scale01)
  dat <- components$dat

  # Both rearranged CDFs are step functions. Their joint set of probability
  # knots and the intervening midpoints cover every distinct Youden value.
  knots <- components$cutoff_knots
  cutoff_grid <- sort(unique(c(
    knots,
    if (length(knots) > 1) (knots[-1] + knots[-length(knots)]) / 2 else numeric()
  )))
  F1_cut <- components$F1_plus(cutoff_grid)
  F0_cut <- components$F0_plus(cutoff_grid)
  J_grid <- F0_cut - F1_cut
  best <- which.max(J_grid)

  auc <- np.auc(dat, pi0 = pi0, pi1 = pi1)
  x_auc <- seq(0, 1, length.out = 2001)
  auc_integrated <- .np_trapz(x_auc, components$roc(x_auc))
  x_pauc <- seq(s.range[1], s.range[2], length.out = 1001)
  pAUC <- .np_trapz(x_pauc, components$roc(x_pauc)) /
    (s.range[2] - s.range[1])

  data.frame(
    ROC = components$roc(t),
    AUC = unname(auc["AUC"]),
    apparent_AUC = unname(auc["apparent_AUC"]),
    AUC_integrated = auc_integrated,
    cutoff = cutoff_grid[best],
    opt_cut = cutoff_grid[best],
    Youden = J_grid[best],
    J = J_grid[best],
    se = 1 - F1_cut[best],
    sp = F0_cut[best],
    pAUC = pAUC,
    n = nrow(dat)
  )
}

# Short alias for comparison scripts.
NP <- np.roc
NP.inc <- np.inc.roc
