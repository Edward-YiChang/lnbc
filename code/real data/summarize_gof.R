# Run from real data/ after ../gof/real_data_gof.R --variant div100000.
saved <- readRDS("../gof/real_data_result/div100000/gof_result.rds")
x <- saved$gof
config <- saved$config
stopifnot(config$divisor == 100000)
observed <- unname(x$observed["sup3"])
boot <- x$bootstrap_statistics[, "sup3"]
stopifnot(all(is.finite(boot)), length(boot) == config$B,
          isTRUE(all.equal(observed, unname(
            config$n0 / (config$n0 + config$n1) * x$observed["D0"] +
            config$n1 / (config$n0 + config$n1) * x$observed["D1"]))))
p_value <- mean(boot >= observed)
result <- data.frame(statistic = "sup3", Dn = observed,
                     n0 = config$n0, n1 = config$n1, B = length(boot),
                     exceedances = sum(boot >= observed), p_value = p_value,
                     alpha = config$alpha, reject = p_value <= config$alpha)
dir.create("biomarker_div_100000", showWarnings = FALSE)
write.csv(result, "biomarker_div_100000/gof_sup3.csv", row.names = FALSE)
print(result, row.names = FALSE)
