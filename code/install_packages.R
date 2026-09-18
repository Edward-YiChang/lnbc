cat("LNBC install_packages\n")

repo <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library"))
pkgs <- c("glmnet", "mgcv", "spatstat.univar", "rootSolve")
install_order <- c("RcppEigen", "glmnet", "spatstat.utils",
                   "spatstat.univar", "mgcv", "rootSolve")

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

cat("Library:", lib, "\n")
cat("CRAN repo:", repo, "\n")
cat("R version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("R executable:", file.path(R.home("bin"), "R"), "\n")
cat(".libPaths():\n")
cat(paste0("  ", .libPaths(), collapse = "\n"), "\n", sep = "")
cat("MAKEFLAGS:", Sys.getenv("MAKEFLAGS", unset = ""), "\n")
cat("R_MAKEVARS_USER:", Sys.getenv("R_MAKEVARS_USER", unset = ""), "\n")
cat("Required packages:", paste(pkgs, collapse = ", "), "\n")

is_available <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

for (pkg in install_order) {
  if (is_available(pkg)) {
    path <- tryCatch(find.package(pkg)[1], error = function(e) NA_character_)
    cat(sprintf("Already available: %-16s %s\n", pkg, path))
    next
  }
  cat("Installing:", pkg, "\n")
  tryCatch(
    install.packages(pkg, repos = repo, lib = lib, dependencies = TRUE,
                     Ncpus = 1),
    error = function(e) {
      cat("Install error for ", pkg, ": ", conditionMessage(e), "\n", sep = "")
    }
  )
}

cat("Verifying packages...\n")
failed <- character()
for (pkg in pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-16s %s\n", pkg, if (ok) "OK" else "FAILED"))
  if (!ok) failed <- c(failed, pkg)
}

if (length(failed) > 0) {
  stop("Failed to install/load package(s): ", paste(failed, collapse = ", "))
}

cat("Package setup complete.\n")
