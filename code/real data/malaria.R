source("../bootstrap_core.R")

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 500L
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
if (is.na(B) || B < 1) stop("B must be a positive integer.")
if (is.na(cores) || cores < 1) stop("cores must be a positive integer.")

default_methods <- c("EMBC", "NP-inc", "LBNP", "LBNP-log")
methods <- if (length(args) >= 3) {
  trimws(strsplit(args[3], ",", fixed = TRUE)[[1]])
} else {
  default_methods
}
unknown_methods <- setdiff(methods, default_methods)
if (length(unknown_methods) > 0 || length(methods) == 0) {
  stop("Method argument must be one or more of: ",
       paste(default_methods, collapse = ", "),
       ". Separate multiple methods with commas.")
}
save_mode <- if (length(args) >= 4) tolower(trimws(args[4])) else "overwrite"
if (!save_mode %in% c("merge", "overwrite")) {
  stop("Save mode must be either 'merge' or 'overwrite'.")
}
preprocess <- if (length(args) >= 5) tolower(trimws(args[5])) else "scale:100000"
output_dir <- if (length(args) >= 6) trimws(args[6]) else "biomarker_div_100000"
if (!nzchar(output_dir)) stop("Output directory must not be empty.")

scale_match <- regexec("^scale:([0-9]+(?:[.][0-9]+)?)$", preprocess)
scale_parts <- regmatches(preprocess, scale_match)[[1]]
if (preprocess %in% c("log", "original")) {
  scale_divisor <- NA_real_
} else if (length(scale_parts) == 2) {
  scale_divisor <- as.numeric(scale_parts[2])
  if (!is.finite(scale_divisor) || scale_divisor <= 0) {
    stop("Scale divisor must be positive and finite.")
  }
} else {
  stop("Preprocessing must be 'original', 'log', or 'scale:<positive divisor>'.")
}
metrics <- c("AUC", "J", "se", "sp")
conf_level <- 0.95
fpr_grid <- seq(0, 1, by = 0.01)

data <- read.table("malaria.txt", header = FALSE, col.names = c("R", "T"))
data <- data[data$T != 0, , drop = FALSE]
if (preprocess == "log") {
  data$T <- log(data$T)
  # LBNP-log applies log() once more, so all methods require transformed T>0.
  data <- data[data$T > 0, , drop = FALSE]
  inclusion <- "original T != 0 and transformed log(T) > 0"
  biomarker_transform <- "shared input log(T)"
  method_transforms <- c(
    EMBC = "shared log(T) input",
    NP = "shared log(T) input, min-max scaled to [0,1]",
    `NP-inc` = "shared log(T) input, min-max scaled to [0,1]",
    LBNP = "shared log(T) input",
    `LBNP-log` = "log of the shared log(T) input"
  )
} else if (preprocess == "original") {
  inclusion <- "original T != 0"
  biomarker_transform <- "original unscaled T"
  method_transforms <- c(
    EMBC = biomarker_transform,
    NP = paste0(biomarker_transform, ", min-max scaled to [0,1]"),
    `NP-inc` = paste0(biomarker_transform, ", min-max scaled to [0,1]"),
    LBNP = biomarker_transform,
    `LBNP-log` = paste0("log of ", biomarker_transform)
  )
} else {
  data$T <- data$T / scale_divisor
  inclusion <- "original T != 0"
  biomarker_transform <- paste0("shared input T/", format(scale_divisor, scientific = FALSE))
  method_transforms <- c(
    EMBC = biomarker_transform,
    NP = paste0(biomarker_transform, ", min-max scaled to [0,1]"),
    `NP-inc` = paste0(biomarker_transform, ", min-max scaled to [0,1]"),
    LBNP = biomarker_transform,
    `LBNP-log` = paste0("log of ", biomarker_transform)
  )
}

# Match real-application.R in the authors' repository.
lbnp_nu_seq <- 10^seq(from = -3, to = 2, length.out = 11)

roc_metrics <- lnbc_roc_metric_names(fpr_grid)
results <- setNames(vector("list", length(methods)), methods)

for (method in methods) {
  fit <- lnbc_bootstrap_method(
    data = data,
    method = method,
    B = B,
    pi0 = 1,
    pi1 = 0.677,
    dist = "real_logged",
    seed = 123,
    conf_level = conf_level,
    repo_root = normalizePath(".."),
    lbnp_k = 10,
    lbnp_nu_seq = lbnp_nu_seq,
    lbnp_maxit = 1000,
    lbnp_thres = 1e-5,
    lbnp_transform = NULL,
    roc_grid = fpr_grid,
    ci_metrics = metrics,
    cores = cores,
    verbose = TRUE
  )

  results[[method]] <- list(
    estimates = fit$ci,
    roc = lnbc_roc_percentile_band(
      fit$point_roc,
      fit$bootstrap[, roc_metrics, drop = FALSE],
      conf_level = conf_level
    ),
    point_status = fit$point_status,
    bootstrap_status = fit$bootstrap_status,
    bootstrap = fit$bootstrap[, c(metrics, roc_metrics), drop = FALSE],
    seeds = fit$seeds
  )
}

config <- list(
    B = B,
    cores = cores,
    conf_level = conf_level,
    pi0 = 1,
    pi1 = 0.677,
    lbnp_k = 10,
    lbnp_nu_seq = lbnp_nu_seq,
    lbnp_maxit = 1000,
    lbnp_thres = 1e-5,
    fpr_grid = fpr_grid,
    inclusion = inclusion,
    biomarker_transform = biomarker_transform,
    preprocess = preprocess,
    scale_divisor = scale_divisor,
    method_transforms = method_transforms,
    paper_proposed_mapping = "EMBC",
    np_scaling = paste0(biomarker_transform, ", then min-max to [0,1]"),
    np_youden = "grid-based estimator reported in the authors' FUN-np.R"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(output_dir, "result.rds")
estimate_file <- file.path(output_dir, "estimates.csv")
roc_file <- file.path(output_dir, "roc_bands.csv")

if (save_mode == "merge" && file.exists(result_file)) {
  existing <- readRDS(result_file)
  compatibility_keys <- c(
    "B", "conf_level", "pi0", "pi1", "fpr_grid", "inclusion",
    "biomarker_transform"
  )
  incompatible <- compatibility_keys[!vapply(compatibility_keys, function(key) {
    !is.null(existing$config[[key]]) &&
      isTRUE(all.equal(existing$config[[key]], config[[key]], check.attributes = FALSE))
  }, logical(1))]
  if (length(incompatible) > 0) {
    stop(
      "Cannot merge with result.rds because these settings differ or are missing: ",
      paste(incompatible, collapse = ", "),
      ". Use save mode 'overwrite' or preserve/rename the existing result first."
    )
  }
  merged_results <- existing$methods
  merged_results[methods] <- results
  merged_order <- default_methods[default_methods %in% names(merged_results)]
  merged_results <- merged_results[merged_order]
  results <- merged_results
  methods <- merged_order
  config$merged_methods <- methods
}

result <- list(
  config = config,
  methods = results
)

saveRDS(result, file = result_file)

# Rebuild the CSV files from the final method collection so merge mode cannot
# leave duplicate or stale rows.
estimate_table <- do.call(rbind, lapply(methods, function(method) {
  cbind(method = method, results[[method]]$estimates)
}))
rownames(estimate_table) <- NULL
roc_table <- do.call(rbind, lapply(methods, function(method) {
  cbind(method = method, results[[method]]$roc)
}))
rownames(roc_table) <- NULL
write.csv(estimate_table, file = estimate_file, row.names = FALSE)
write.csv(roc_table, file = roc_file, row.names = FALSE)
cat(
  "Saved", paste(methods, collapse = ", "), "using", save_mode,
  "mode in", normalizePath(output_dir), "with", biomarker_transform, "preprocessing.\n"
)
