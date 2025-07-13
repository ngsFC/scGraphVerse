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

#' Setup GRNBoost2 Python Dependencies
#'
#' Automatically installs the Python package 'arboreto' required for GRNBoost2.
#' This function provides a convenient way to set up the Python environment
#' without manual command-line installation.
#'
#' @param method Character string specifying installation method.
#'   Options: "auto" (default), "pip", "conda", "virtualenv".
#' @param conda Character string. Name of conda environment to use.
#'   If NULL, uses the default environment. Only used when method = "conda".
#' @param python_version Character string. Python version to use when
#'   creating new environments. Default: "3.8".
#' @param force Logical. If TRUE, reinstalls arboreto even if already present.
#'   Default: FALSE.
#' @param verbose Logical. If TRUE, shows installation progress messages.
#'   Default: TRUE.
#'
#' @return Logical. TRUE if installation successful, FALSE otherwise.
#'
#' @details This function attempts to install the 'arboreto' Python package
#' required for GRNBoost2 network inference. It first checks if arboreto is
#' already available and installs it only if needed (unless force = TRUE).
#'
#' Installation methods:
#' \itemize{
#'   \item \code{"auto"}: Tries conda first, falls back to pip
#'   \item \code{"pip"}: Uses pip to install from PyPI
#'   \item \code{"conda"}: Uses conda to install from conda-forge/bioconda
#'   \item \code{"virtualenv"}: Creates new virtual environment with arboreto
#' }
#'
#' @importFrom reticulate py_module_available py_install conda_install
#' @export
#'
#' @examples
#' \dontrun{
#' # Automatic installation (recommended)
#' setup_grnboost2()
#'
#' # Force reinstallation using conda
#' setup_grnboost2(method = "conda", force = TRUE)
#'
#' # Install in specific conda environment
#' setup_grnboost2(method = "conda", conda = "my_env")
#' }
setup_grnboost2 <- function(method = "auto", conda = NULL, 
                           python_version = "3.8", force = FALSE, 
                           verbose = TRUE) {
    
    # Check if already installed (unless forcing)
    if (!force && reticulate::py_module_available("arboreto")) {
        if (verbose) message("arboreto is already installed.")
        return(TRUE)
    }
    
    if (verbose) {
        message("Installing Python package 'arboreto' for GRNBoost2...")
        message("This may take a few minutes.")
    }
    
    success <- FALSE
    
    tryCatch({
        if (method == "auto") {
            # Try conda first, fallback to pip
            if (verbose) message("Trying conda installation...")
            tryCatch({
                reticulate::conda_install("arboreto", channel = "bioconda")
                success <- TRUE
                if (verbose) message("Successfully installed via conda.")
            }, error = function(e) {
                if (verbose) message("Conda failed, trying pip...")
                reticulate::py_install("arboreto")
                success <<- TRUE
                if (verbose) message("Successfully installed via pip.")
            })
        } else if (method == "conda") {
            reticulate::conda_install("arboreto", envname = conda, 
                                     channel = "bioconda")
            success <- TRUE
            if (verbose) message("Successfully installed via conda.")
        } else if (method == "pip") {
            reticulate::py_install("arboreto")
            success <- TRUE
            if (verbose) message("Successfully installed via pip.")
        } else if (method == "virtualenv") {
            env_name <- paste0("scgraphverse_py", gsub("\\.", "", python_version))
            reticulate::virtualenv_create(env_name, python = python_version)
            reticulate::virtualenv_install(env_name, "arboreto")
            reticulate::use_virtualenv(env_name)
            success <- TRUE
            if (verbose) message("Successfully installed in virtual environment: ", env_name)
        } else {
            stop("Invalid method. Use 'auto', 'conda', 'pip', or 'virtualenv'.")
        }
    }, error = function(e) {
        if (verbose) {
            message("Installation failed with error: ", e$message)
            message("\nTroubleshooting suggestions:")
            message("1. Ensure Python is properly installed")
            message("2. Try: reticulate::py_install('arboreto')")
            message("3. Manual install: pip install arboreto")
            message("4. See vignette('python-setup', package = 'scGraphVerse')")
        }
        success <<- FALSE
    })
    
    # Verify installation
    if (success) {
        Sys.sleep(1)  # Brief pause for module loading
        if (reticulate::py_module_available("arboreto")) {
            if (verbose) message("Installation verified successfully!")
            return(TRUE)
        } else {
            if (verbose) message("Installation completed but module not detected. You may need to restart R.")
            return(FALSE)
        }
    }
    
    return(FALSE)
}
