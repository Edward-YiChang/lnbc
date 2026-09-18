# Main proposed LBNP method from the downloaded GitHub code.
# Source this file to load DRest0.EM and use lbnp.fit for a single biomarker.

.lbnp_source_dir <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile)))
    }
  }
  getwd()
}

.lbnp_dir <- .lbnp_source_dir()
source(file.path(.lbnp_dir, "original", "FUN.R"))

.lbnp_check_packages <- function() {
  pkgs <- c("mgcv", "spatstat.univar", "rootSolve")
  missing <- pkgs[!vapply(pkgs, require, logical(1), quietly = TRUE, character.only = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "))
  }
}

.lbnp_prepare_data <- function(dat) {
  if (!is.data.frame(dat) || !all(c("R", "biomarker") %in% names(dat))) {
    stop("dat must be a data.frame with columns R and biomarker.")
  }
  if (!all(dat$R %in% c(0, 1))) {
    stop("R must contain only 0/1 labels.")
  }
  if (sum(dat$R == 0) == 0 || sum(dat$R == 1) == 0) {
    stop("R must contain at least one 0 and one 1.")
  }
  if (!all(is.finite(dat$biomarker))) {
    stop("biomarker must be finite.")
  }
  data.frame(R = dat$R, biomarker = dat$biomarker)
}

.lbnp_cut_stats <- function(dat, fit, cutoff, ELtuneType = 0) {
  if (is.na(cutoff)) {
    return(c(opt_cut = NA, J = NA, se = NA, sp = NA))
  }
  if (ELtuneType == 0) {
    Fg0 <- ewcdf(dat$biomarker, weights = fit$weight.Dis$G0)
    Fg1 <- ewcdf(dat$biomarker, weights = fit$weight.Dis$G1)
  } else {
    Fg0 <- ewcdf(dat$biomarker, weights = fit$weight.Dis.est$G0)
    Fg1 <- ewcdf(dat$biomarker, weights = fit$weight.Dis.est$G1)
  }
  sp <- Fg0(cutoff)
  se <- 1 - Fg1(cutoff)
  c(opt_cut = cutoff, J = se + sp - 1, se = se, sp = sp)
}

.lbnp_add_cut_stats <- function(dat, fit, stats, ELtuneType = 0) {
  cut_stats <- .lbnp_cut_stats(dat, fit, stats$cutoff, ELtuneType)
  stats$opt_cut <- cut_stats["opt_cut"]
  stats$J <- cut_stats["J"]
  stats$se <- cut_stats["se"]
  stats$sp <- cut_stats["sp"]
  stats
}

lbnp.fit <- function(dat, pi0, pi1,
                     nu = 1, k = 10, ord = 4, ord.pen = 2,
                     maxit = 500, thres = 1e-4,
                     knots = NULL, ini.G = NULL,
                     s = 0.2, s.range = c(0.1, 0.3),
                     ELtuneType = 0, eval.pts = NULL,
                     alt.auc = FALSE) {
  .lbnp_check_packages()
  dat <- .lbnp_prepare_data(dat)
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }

  fit <- DRest0.EM(
    dat = dat,
    pi0 = pi0,
    pi1 = pi1,
    nu = nu,
    k = k,
    ord = ord,
    ord.pen = ord.pen,
    maxit = maxit,
    thres = thres,
    knots = knots,
    ini.G = ini.G
  )

  summary <- sum.DRest0.EM(
    dat = dat,
    EMout = fit,
    s = s,
    ELtuneType = ELtuneType,
    eval.pts = eval.pts,
    s.range = s.range,
    alt.auc = alt.auc
  )
  summary$stats.df <- .lbnp_add_cut_stats(dat, fit, summary$stats.df, ELtuneType)

  list(
    fit = fit,
    data = dat,
    stats = summary$stats.df,
    IC = summary$IC.df,
    h = summary$h.df,
    conv = fit$conv,
    pi0 = pi0,
    pi1 = pi1,
    nu = nu
  )
}

lbnp.fit.cv5 <- function(dat, pi0, pi1,
                         nu.seq = 3 * 10^seq(from = -2, to = 2.5, length.out = 10),
                         k = 10, ord = 4, ord.pen = 2,
                         maxit = 500, thres = 1e-5,
                         knots = "q", ini.G = NULL,
                         nfolds = 5, seed = 123,
                         s = 0.2, s.range = c(0.1, 0.3),
                         ELtuneType = 0, eval.pts = NULL,
                         alt.auc = FALSE) {
  .lbnp_check_packages()
  dat <- .lbnp_prepare_data(dat)
  if (pi0 + pi1 <= 1) {
    stop("pi0 + pi1 must be greater than 1.")
  }

  cvout <- cv.DRest.EM(
    dat = dat,
    pi0 = pi0,
    pi1 = pi1,
    nfolds = nfolds,
    seed = seed,
    nu = nu.seq,
    k = k,
    ord = ord,
    ord.pen = ord.pen,
    maxit = maxit,
    thres = thres,
    knots = knots,
    ini.G = ini.G
  )

  tune.id <- cvout$tune.id
  if (length(tune.id) == 0 || is.na(tune.id) || !isTRUE(cvout$res$conv[tune.id])) {
    return(list(
      fit = NULL,
      data = dat,
      stats = data.frame(AUC = NA, cutoff = NA, opt_cut = NA, J = NA, se = NA, sp = NA),
      cv = cvout,
      conv = FALSE,
      pi0 = pi0,
      pi1 = pi1,
      nu = NA
    ))
  }

  fit <- DRest0.EM(
    dat = dat,
    pi0 = pi0,
    pi1 = pi1,
    nu = nu.seq[tune.id],
    k = k,
    ord = ord,
    ord.pen = ord.pen,
    maxit = maxit,
    thres = thres,
    knots = knots,
    ini.G = ini.G
  )

  if (!isTRUE(fit$conv)) {
    return(list(
      fit = fit,
      data = dat,
      stats = data.frame(AUC = NA, cutoff = NA, opt_cut = NA, J = NA, se = NA, sp = NA),
      cv = cvout,
      conv = FALSE,
      pi0 = pi0,
      pi1 = pi1,
      nu = nu.seq[tune.id]
    ))
  }

  summary <- sum.DRest0.EM(
    dat = dat,
    EMout = fit,
    s = s,
    ELtuneType = ELtuneType,
    eval.pts = eval.pts,
    s.range = s.range,
    alt.auc = alt.auc
  )
  summary$stats.df <- .lbnp_add_cut_stats(dat, fit, summary$stats.df, ELtuneType)

  list(
    fit = fit,
    data = dat,
    stats = summary$stats.df,
    IC = summary$IC.df,
    h = summary$h.df,
    cv = cvout,
    conv = TRUE,
    pi0 = pi0,
    pi1 = pi1,
    nu = nu.seq[tune.id]
  )
}

lbnp.roc <- function(s.seq, object, dat = NULL, ELtuneType = 0) {
  .lbnp_check_packages()
  if (is.list(object) && !is.null(object$fit)) {
    fit <- object$fit
    if (is.null(dat)) {
      dat <- object$data
    }
  } else {
    fit <- object
  }
  dat <- .lbnp_prepare_data(dat)
  est.ROC(s.seq = s.seq, dat = dat, EMout = fit, ELtuneType = ELtuneType)
}

# Short aliases for comparison scripts.
LBNP <- lbnp.fit
LBNP.cv5 <- lbnp.fit.cv5
LBNP.roc <- lbnp.roc
