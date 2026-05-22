# Suppress R CMD check NOTEs for variables created by RTMB::getAll()
utils::globalVariables(c(
  "log_r", "log_K", "log_m", "log_sigma_obs",
  "log_sigma_proc", "proc_dev"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "SurplusProductionModel v",
    utils::packageVersion("SurplusProductionModel"),
    " loaded."
  )
  packageStartupMessage("Pella-Tomlinson surplus production model.")
}

.onLoad <- function(libname, pkgname) {
  # Check for RTMB availability
  if (!requireNamespace("RTMB", quietly = TRUE)) {
    warning("RTMB package not available. Install with: install.packages('RTMB')")
  }
}

# Clean up on unload
.onUnload <- function(libpath) {
  # Clean up any global variables or connections if needed
}
