# Main-text Table 1: true parameters and diagnostic accuracy measures.
source("bootstrap_core.R")
labels <- c(lognorm = "Lognormal (k=0)", weibull = "Weibull (k=1/2)",
            weibull2 = "Gamma (k=1)")
rows <- list()
for (dist in names(labels)) for (J_key in c("3", "5")) {
  s <- lnbc_scenario(dist, J_key)
  truth <- lnbc_true_values(dist, J_key)
  # Express Weibull shape 1 as Gamma(shape, scale) for the manuscript.
  parameters <- if (dist == "weibull2") c(1, s$a0, 1, s$a1) else {
    c(s$a0, s$b0, s$a1, s$b1)
  }
  rows[[length(rows) + 1L]] <- data.frame(
    distribution = labels[[dist]], AUC = unname(truth["AUC"]),
    J = unname(truth["J"]), se = unname(truth["se"]), sp = unname(truth["sp"]),
    a0 = parameters[1], b0 = parameters[2], a1 = parameters[3], b1 = parameters[4]
  )
}
dir.create("simulation/latex", recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), "simulation/latex/table1_parameters.csv",
          row.names = FALSE)
