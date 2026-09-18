args <- commandArgs(trailingOnly = TRUE)
result_dir <- if (length(args) >= 1) args[1] else "biomarker_div_100000"

result_file <- file.path(result_dir, "result.rds")
if (!file.exists(result_file)) {
  stop("Cannot find scaled real-data result: ", result_file)
}

saved <- readRDS(result_file)
scale_divisor <- saved$config$scale_divisor
is_original <- identical(saved$config$preprocess, "original") ||
  identical(saved$config$biomarker_transform, "original unscaled T")
if (is_original) {
  analysis_tag <- "original"
  file_stem <- "malaria_original"
  caption_scale <- "the original unscaled biomarker values"
  plot_title <- "Malaria ROC curves: original biomarker scale"
} else {
  if (length(scale_divisor) != 1 || !is.finite(scale_divisor) ||
      scale_divisor <= 0) {
    stop("Result is neither an original-scale result nor a valid scaled result.")
  }
  divisor_tag <- format(scale_divisor, scientific = FALSE, trim = TRUE)
  divisor_label <- format(scale_divisor, scientific = FALSE, trim = TRUE,
                          big.mark = ",")
  analysis_tag <- paste0("div", divisor_tag)
  file_stem <- paste0("malaria_", analysis_tag)
  caption_scale <- paste0(
    "the biomarker divided by ",
    gsub(",", "{,}", divisor_label, fixed = TRUE)
  )
  plot_title <- paste0("Malaria ROC curves: biomarker / ", divisor_label)
}
pdf_file <- if (length(args) >= 2) args[2] else {
  file.path(result_dir, paste0(file_stem, "_roc.pdf"))
}

method_order <- c("EMBC", "NP-inc", "LBNP", "LBNP-log")
method_order <- method_order[method_order %in% names(saved$methods)]
display_names <- setNames(method_order, method_order)
display_names["EMBC"] <- "Ours"
display_names["NP-inc"] <- "NP"
metric_names <- c("ROC_0.20", "AUC", "J", "se", "sp")
metric_labels <- c("$\\operatorname{ROC}(0.2)$", "AUC", "$J$",
                   "Sensitivity", "Specificity")
probs <- c((1 - saved$config$conf_level) / 2,
           1 - (1 - saved$config$conf_level) / 2)

successful_fit <- function(z) {
  z$bootstrap_status$ok & z$bootstrap_status$conv &
    z$bootstrap_status$roc_ok
}

point_value <- function(z, metric) {
  if (metric == "ROC_0.20") {
    i <- which(abs(z$roc$fpr - 0.2) < 1e-12)
    if (length(i) != 1) stop("ROC grid does not contain 0.2.")
    return(z$roc$tpr[i])
  }
  value <- z$estimates$estimate[z$estimates$metric == metric]
  if (length(value) != 1) stop("Missing point estimate for ", metric, ".")
  value
}

fmt <- function(x) sprintf("%.3f", x)
table_cells <- t(vapply(method_order, function(method) {
  z <- saved$methods[[method]]
  ok <- successful_fit(z)
  vapply(metric_names, function(metric) {
    values <- z$bootstrap[ok, metric]
    interval <- quantile(values, probs = probs, names = FALSE, na.rm = TRUE)
    sprintf("%s (%s, %s)", fmt(point_value(z, metric)),
            fmt(interval[1]), fmt(interval[2]))
  }, character(1))
}, character(length(metric_names))))
colnames(table_cells) <- metric_labels
rownames(table_cells) <- display_names[method_order]

valid_counts <- vapply(method_order, function(method) {
  sum(successful_fit(saved$methods[[method]]))
}, integer(1))
total_counts <- vapply(method_order, function(method) {
  nrow(saved$methods[[method]]$bootstrap_status)
}, integer(1))

tex_lines <- c(
  "% Requires \\usepackage{booktabs,graphicx}",
  "\\begin{table}[htbp]",
  "\\centering",
  paste0(
    "\\caption{Malaria-data point estimates using ", caption_scale,
    ". Parentheses contain 95\\% percentile bootstrap confidence intervals.}"
  ),
  paste0("\\label{tab:malaria-", analysis_tag, "}"),
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  paste(c("Method", metric_labels), collapse = " & ") |> paste0(" \\\\"),
  "\\midrule",
  vapply(seq_along(method_order), function(i) {
    paste(c(display_names[method_order[i]], table_cells[i, ]), collapse = " & ") |>
      paste0(" \\\\")
  }, character(1)),
  "\\bottomrule",
  "\\end{tabular}",
  "}%",
  "\\par\\smallskip",
  "\\begin{minipage}{0.98\\linewidth}",
  paste0(
    "\\footnotesize \\textit{Note:} $\\operatorname{ROC}(0.2)$ is the ",
    "sensitivity at false-positive rate 0.2. Confidence intervals use only ",
    "successful fits with valid ROC curves. Valid bootstrap fits: ",
    paste0(display_names[method_order], " ", valid_counts, "/", total_counts,
           collapse = "; "), "."
  ),
  "\\end{minipage}",
  "\\end{table}"
)
tex_file <- file.path(
  result_dir,
  paste0(file_stem, "_point_estimates.tex")
)
writeLines(tex_lines, tex_file)

fpr_grid <- saved$config$fpr_grid
roc_curves <- lapply(method_order, function(method) {
  z <- saved$methods[[method]]
  data.frame(fpr = fpr_grid, tpr = z$roc$tpr)
})
names(roc_curves) <- method_order

colors <- c(EMBC = "#0072B2", NP = "#D55E00", `NP-inc` = "#E69F00",
            LBNP = "#009E73", `LBNP-log` = "#CC79A7")
roc_lwd <- 1.25
plot_methods <- method_order
plot_display_names <- display_names
plot_colors <- colors
if (identical(analysis_tag, "div100000")) {
  plot_methods <- setdiff(plot_methods, "NP")
  plot_display_names["NP-inc"] <- "NP"
  plot_colors["NP-inc"] <- colors["NP"]
  plot_title <- "ROC curve for Malaria data"
}
x_axis_label <- if (identical(analysis_tag, "div100000")) {
  "False positive rate"
} else {
  "False-positive rate (1 - specificity)"
}
y_axis_label <- if (identical(analysis_tag, "div100000")) {
  "True positive rate"
} else {
  "True-positive rate (sensitivity)"
}
draw_roc <- function() {
  old <- par(no.readonly = TRUE)
  on.exit({
    layout(1)
    par(old)
  })
  layout(matrix(1:2, ncol = 1), heights = c(1, 0.18))
  par(mar = c(3.2, 4.8, 3.4, 1.1), las = 1)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = "", ylab = y_axis_label, axes = FALSE)
  axis(1, at = seq(0, 1, 0.2))
  axis(2, at = seq(0, 1, 0.2))
  box()
  grid(col = "gray88", lty = 1)
  abline(0, 1, col = "black", lty = 2, lwd = roc_lwd)

  for (method in plot_methods) {
    z <- roc_curves[[method]]
    lines(z$fpr, z$tpr, col = plot_colors[[method]], lwd = roc_lwd)
  }
  title(plot_title, line = 1.4, cex.main = 1.05, font.main = 1)

  par(mar = c(0, 0, 0, 0), las = 1)
  plot.new()
  text(0.5, 0.88, x_axis_label, cex = 1)
  legend(
    "center",
    legend = unname(plot_display_names[plot_methods]),
    col = unname(plot_colors[plot_methods]),
    lwd = rep(roc_lwd, length(plot_methods)),
    lty = rep(1, length(plot_methods)),
    bty = "n", ncol = length(plot_methods), cex = 0.82, x.intersp = 0.8,
    y.intersp = 1.0
  )
}

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
pdf(pdf_file, width = 9, height = 6, useDingbats = FALSE, version = "1.4")
draw_roc()
dev.off()

png_file <- file.path(result_dir, paste0(file_stem, "_roc.png"))
png_type <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
  "quartz"
} else if (capabilities("cairo")) {
  "cairo"
} else {
  "Xlib"
}
png(png_file, width = 2160, height = 1440, res = 240,
    type = png_type)
draw_roc()
dev.off()

cat("Wrote", tex_file, pdf_file, png_file, sep = "\n")
cat("\n")

# Transformation estimate and percentile CI reported in Section 6.
z <- saved$methods$EMBC
ok <- successful_fit(z) & is.finite(z$bootstrap_status$embc_k)
k_interval <- quantile(z$bootstrap_status$embc_k[ok], probs, names = FALSE)
write.csv(data.frame(estimate = z$point_status$embc_k,
                     ci_low = k_interval[1], ci_high = k_interval[2]),
          file.path(result_dir, "kappa_summary.csv"), row.names = FALSE)
