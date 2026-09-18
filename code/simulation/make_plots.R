#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd())
if (!file.exists(file.path(repo_root, "bootstrap_core.R"))) {
  stop("Run this script from the repository root: Rscript simulation/make_plots.R")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.")
}

args <- commandArgs(trailingOnly = TRUE)
input_root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else file.path(repo_root, "result")
result_dir <- file.path(input_root, "ROC_plot")
output_dir <- file.path(repo_root, "simulation", "plots")

roc_methods <- c("EMBC", "NP-inc", "LBNP", "LBNP-log")
roc_labels <- c(
  EMBC = "Ours",
  `NP-inc` = "NP",
  LBNP = "LBNP",
  `LBNP-log` = "LBNP-log"
)
roc_cols <- c(
  EMBC = "#0072B2",
  `NP-inc` = "#D55E00",
  LBNP = "#009E73",
  `LBNP-log` = "#CC79A7"
)
roc_suffix <- c(
  EMBC = "roc_band.csv",
  `NP-inc` = "roc.csv",
  LBNP = "roc_band.csv",
  `LBNP-log` = "roc_band.csv"
)
roc_lwd <- 1.25
axis_title_cex <- 1

dist_values <- c("lognorm", "weibull", "weibull2")
dist_labels <- list(
  lognorm = expression(paste("Lognormal (", kappa, " = 0)")),
  weibull = expression(paste("Weibull (", kappa, " = 1/2)")),
  weibull2 = expression(paste("Gamma (", kappa, " = 1)"))
)
n_values <- c(100, 300, 500)
J_values <- c("3", "5")
pi_values <- c("75", "9")

roc_path <- function(method, dist, n_val, J_key, pi_key) {
  path <- file.path(
    result_dir,
    method,
    paste0(
      "n", n_val, "J", J_key, "pi", pi_key, dist,
      "_final_B500_reps1000_", roc_suffix[[method]]
    )
  )
  unified <- sub("roc_band\\.csv$", "roc.csv", path)
  if (file.exists(unified)) unified else path
}

required <- expand.grid(
  method = roc_methods,
  dist = dist_values,
  n = n_values,
  J = J_values,
  pi = pi_values,
  stringsAsFactors = FALSE
)
required$path <- mapply(
  roc_path,
  required$method,
  required$dist,
  required$n,
  required$J,
  required$pi,
  USE.NAMES = FALSE
)

stopifnot(nrow(required) == 144L, length(unique(required$path)) == 144L)
missing_paths <- required$path[!file.exists(required$path)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing ", length(missing_paths), " of 144 required ROC inputs:\n",
    paste(missing_paths, collapse = "\n")
  )
}

required_columns <- c("fpr", "tpr", "true_tpr")
bad_schema <- vapply(required$path, function(path) {
  columns <- names(read.csv(path, nrows = 0L, check.names = FALSE))
  !all(required_columns %in% columns)
}, logical(1))
if (any(bad_schema)) {
  stop(
    "Required columns are missing from:\n",
    paste(required$path[bad_schema], collapse = "\n")
  )
}
message("Validated all 144 required raw ROC inputs.")

read_curve <- function(path) {
  x <- data.table::fread(
    path,
    select = required_columns,
    showProgress = FALSE
  )
  x <- x[
    is.finite(fpr) & is.finite(tpr) & is.finite(true_tpr),
    lapply(.SD, mean, na.rm = TRUE),
    by = fpr,
    .SDcols = c("tpr", "true_tpr")
  ]
  as.data.frame(x[order(fpr)])
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plots <- expand.grid(J = J_values, pi = pi_values, stringsAsFactors = FALSE)
plots$main <- FALSE
plots <- rbind(plots, data.frame(J = "5", pi = "9", main = TRUE))
for (i in seq_len(nrow(plots))) {
    J_key <- plots$J[i]
    pi_key <- plots$pi[i]
    main <- plots$main[i]
    plot_n <- if (main) 300 else n_values
    plot_path <- file.path(output_dir, paste0(
      if (main) "n300" else "", "J", J_key, "pi", pi_key, "_roc_band_panel.png"
    ))
    grDevices::png(plot_path, width = if (main) 12 else 10,
                  height = if (main) 5 else 10.5, units = "in", res = 160)
    panels <- length(plot_n) * 3L
    graphics::layout(
      rbind(matrix(seq_len(panels), nrow = length(plot_n), byrow = TRUE),
            rep(panels + 1L, 3L)),
      heights = c(rep(1, length(plot_n)), 0.22)
    )
    graphics::par(oma = c(0.2, 5.2, 0.5, 0.5))

    for (n_val in plot_n) {
      for (dist_index in seq_along(dist_values)) {
        dist <- dist_values[[dist_index]]
        graphics::par(mar = c(2.5, 2.5, 2.0, 0.7))
        curves <- lapply(roc_methods, function(method) {
          read_curve(roc_path(method, dist, n_val, J_key, pi_key))
        })
        names(curves) <- roc_methods
        graphics::plot(
          NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
          xaxs = "i", yaxs = "i", xlab = "", ylab = "", cex.axis = 0.9
        )
        graphics::grid(col = "gray88")
        graphics::abline(0, 1, col = "gray70", lty = 3)
        truth <- curves[[1L]]
        graphics::lines(
          truth$fpr, truth$true_tpr, col = "black", lty = 2, lwd = roc_lwd
        )
        for (method in roc_methods) {
          z <- curves[[method]]
          graphics::lines(z$fpr, z$tpr, col = roc_cols[[method]], lwd = roc_lwd)
        }
        if (n_val == plot_n[[1L]]) {
          graphics::mtext(dist_labels[[dist]], side = 3, line = 0.6, cex = 0.95)
        }
        if (dist_index == 1L && !main) {
          graphics::mtext(
            bquote(paste(n[0], " = ", n[1], " = ", .(n_val))), side = 2, line = 3.5, cex = 0.9
          )
        }
        rm(curves)
        gc(FALSE)
      }
    }

    graphics::par(mar = c(0, 0, 0, 0), cex = 1)
    graphics::plot.new()
    graphics::text(0.5, 0.87, "False positive rate", cex = axis_title_cex)
    legend_labels <- c(unname(roc_labels[roc_methods]), "Truth")
    legend_cols <- c(unname(roc_cols[roc_methods]), "black")
    legend_lty <- c(rep(1, length(roc_methods)), 2)
    legend_x <- seq(0.15, 0.82, length.out = length(legend_labels))
    graphics::segments(
      legend_x - 0.035, 0.28, legend_x, 0.28,
      col = legend_cols, lty = legend_lty, lwd = roc_lwd
    )
    graphics::text(
      legend_x + 0.01, 0.28, legend_labels, adj = c(0, 0.5), cex = 0.82
    )
    graphics::mtext(
      "True positive rate", outer = TRUE, side = 2, line = 3.5,
      at = 0.55, cex = axis_title_cex
    )
    grDevices::dev.off()
    message("Wrote ", plot_path)
}
message("Wrote four combined ROC figures and the main n300J5pi9 figure.")
