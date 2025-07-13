#' Initialize Python Environment for GRNBoost2
#'
#' Sets up the Python environment and lazily loads modules required for
#' running GRNBoost2: \code{arboreto}, \code{pandas}, and \code{numpy}.
#' Automatically installs missing Python packages if requested.
#'
#' @param python_path Character string. Path to the Python executable,
#'   e.g., \code{"/usr/bin/python3"}.
#' @param required Logical. If \code{TRUE}, errors if Python is not
#'   available or path is invalid. Default: \code{TRUE}.
#' @param install_missing Logical. If \code{TRUE}, automatically installs
#'   missing Python packages. Default: \code{FALSE}.
#' @param install_method Character string. Installation method when
#'   \code{install_missing = TRUE}. Options: "auto", "conda", "pip".
#'   Default: "auto".
#' @param verbose Logical. If \code{TRUE}, shows installation progress.
#'   Default: \code{TRUE}.
#'
#' @return A list with three Python module objects:
#'   \itemize{
#'     \item \code{arboreto}: GRNBoost2 algorithm module.
#'     \item \code{pandas}: Data handling module.
#'     \item \code{numpy}: Numerical operations module.
#'   }
#'
#' @details Uses \pkg{reticulate} to bind R to the specified Python
#' interpreter and lazily import modules needed for GRNBoost2. If
#' \code{install_missing = TRUE}, automatically installs the 'arboreto'
#' package using the specified method if not found.
#'
#' @importFrom reticulate use_python import py_module_available py_install conda_install
#' @export
#'
#' @examples
#' # Initialize Python environment (adjust python_path as needed)
#' modules <- init_py()
#' 
#' # Initialize with automatic installation of missing packages
#' modules <- init_py(install_missing = TRUE)
init_py <- function(
    python_path = "/usr/bin/python3",
    required = TRUE,
    install_missing = FALSE,
    install_method = "auto",
    verbose = TRUE) {
    reticulate::use_python(python_path, required = required)
    
    if (install_missing && !reticulate::py_module_available("arboreto")) {
        if (verbose) {
            message("Installing Python package 'arboreto' for GRNBoost2...")
            message("This may take a few minutes.")
        }
        
        tryCatch({
            if (install_method == "auto") {
                if (verbose) message("Trying conda installation...")
                tryCatch({
                    reticulate::conda_install("arboreto", channel = "bioconda")
                    if (verbose) message("Successfully installed via conda.")
                }, error = function(e) {
                    if (verbose) message("Conda failed, trying pip...")
                    reticulate::py_install("arboreto")
                    if (verbose) message("Successfully installed via pip.")
                })
            } else if (install_method == "conda") {
                reticulate::conda_install("arboreto", channel = "bioconda")
                if (verbose) message("Successfully installed via conda.")
            } else if (install_method == "pip") {
                reticulate::py_install("arboreto")
                if (verbose) message("Successfully installed via pip.")
            } else {
                stop("Invalid install_method. Use 'auto', 'conda', or 'pip'.")
            }
        }, error = function(e) {
            if (verbose) {
                message("Installation failed: ", e$message)
                message("Try manual installation: pip install arboreto")
            }
            if (required) stop("Failed to install required Python package 'arboreto'.")
        })
    }

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

    if (verbose) message("Python modules successfully loaded.")
    return(modules)
}
