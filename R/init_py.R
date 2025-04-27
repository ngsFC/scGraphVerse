#' Initialize Python Environment for GRNBoost2
#'
#' Sets up the Python environment and loads required Python modules
#' (\code{arboreto}, \code{pandas}, and \code{numpy}) needed for running GRNBoost2.
#'
#' @param python_path Character string. Path to the Python executable to use 
#'   (e.g., \code{"/usr/bin/python3"}).
#' @param required Logical. If \code{TRUE}, throws an error if Python is not available
#'   or if the specified path is invalid. Default is \code{TRUE}.
#'
#' @return
#' A list containing three Python module objects:
#' \itemize{
#'   \item \code{arboreto}: GRNBoost2 algorithm module.
#'   \item \code{pandas}: Data handling module.
#'   \item \code{numpy}: Numerical operations module.
#' }
#'
#' @details
#' This function uses the \pkg{reticulate} package to bind R to a specified Python interpreter
#' and lazily loads the modules needed for GRNBoost2 network inference. 
#' If any module is missing or incompatible, an informative error will be raised
#' (if \code{required = TRUE}).
#'
#' @importFrom reticulate use_python import
#' @export
#'
#' @examples
#' if (interactive() && requireNamespace("reticulate", quietly = TRUE)) {
#'   # Only run if reticulate is installed and session is interactive
#'   # Initialize Python environment (adjust python_path to your system if needed)
#'   grnboost_modules <- init_py()
#'   
#'   # Access numpy to generate random numbers
#'   grnboost_modules$numpy$random$rand(5)
#' }

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
