#' Initialize Python Environment for GRNBoost2
#'
#' Sets up the Python environment and lazily loads modules required for
#' running GRNBoost2: \code{arboreto}, \code{pandas}, and \code{numpy}.
#'
#' @param python_path Character string. Path to the Python executable,
#'   e.g., \code{"/usr/bin/python3"}.
#' @param required Logical. If \code{TRUE}, errors if Python is not
#'   available or path is invalid. Default: \code{TRUE}.
#'
#' @return A list with three Python module objects:
#'   \itemize{
#'     \item \code{arboreto}: GRNBoost2 algorithm module.
#'     \item \code{pandas}: Data handling module.
#'     \item \code{numpy}: Numerical operations module.
#'   }
#'
#' @details Uses \pkg{reticulate} to bind R to the specified Python
#' interpreter and lazily import modules needed for GRNBoost2. If a module
#' is missing or incompatible, an informative error is raised (when
#' \code{required = TRUE}).
#'
#' @importFrom reticulate use_python import
#' @export
#'
#' @examples
#' # Initialize Python environment (adjust python_path as needed)
#' modules <- init_py()
init_py <- function(
    python_path = "/usr/bin/python3",
    required = TRUE) {
    reticulate::use_python(python_path, required = required)

    modules <- list(
        arboreto = reticulate::import(
            "arboreto.algo",
            delay_load = TRUE
        ),
        pandas = reticulate::import(
            "pandas",
            delay_load = TRUE
        ),
        numpy = reticulate::import(
            "numpy",
            delay_load = TRUE
        )
    )

    message("Python modules successfully loaded.")
    return(modules)
}
