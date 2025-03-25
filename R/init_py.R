#' Initialize Python Environment for GRNBoost2
#'
#' Sets the Python path and imports required Python modules (arboreto, pandas, numpy)
#' for running GRNBoost2.
#'
#' @param python_path Path to the Python executable (e.g., "/usr/bin/python3").
#' @param required Logical. If TRUE, throws an error if Python is not available.
#'
#' @return A list containing Python modules: `arboreto`, `pandas`, and `numpy`.
#'
#' @export
init_py <- function(python_path = "/usr/bin/python3", required = TRUE) {
  reticulate::use_python(python_path, required = required)

  modules <- list(
    arboreto = reticulate::import("arboreto.algo", delay_load = TRUE),
    pandas   = reticulate::import("pandas", delay_load = TRUE),
    numpy    = reticulate::import("numpy", delay_load = TRUE)
  )
  
  message("Python modules successfully loaded.")
  return(modules)
}

