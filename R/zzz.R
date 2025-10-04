# Byte-compile on load
.onLoad <- function(libname, pkgname) {
  # Package startup message
  packageStartupMessage("SurplusProductionModel v", 
                        utils::packageVersion("SurplusProductionModel"), 
                        " loaded.")
  packageStartupMessage("Pella-Tomlinson surplus production model for Antarctic toothfish.")
  
  # Check for RTMB availability
  if (!requireNamespace("RTMB", quietly = TRUE)) {
    warning("RTMB package not available. Install with: install.packages('RTMB')")
  }
}

# Clean up on unload
.onUnload <- function(libpath) {
  # Clean up any global variables or connections if needed
}