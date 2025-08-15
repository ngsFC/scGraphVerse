#' Initialize Python Environment for GRNBoost2
#'
#' Sets up the Python environment and lazily loads modules required for
#' running GRNBoost2: \code{arboreto}, \code{pandas}, and \code{numpy}.
#' Automatically installs missing Python packages if requested.
#'
#' @param python_path Character string. Path to the Python executable,
#'   e.g., \code{"/usr/bin/python3"}. For optimal GRNBoost2 compatibility,
#'   Python 3.8.10 is recommended.
#' @param required Logical. If \code{TRUE}, errors if Python is not
#'   available or path is invalid. Default: \code{TRUE}.
#' @param install_missing Logical. If \code{TRUE}, automatically installs
#'   missing Python packages. Default: \code{FALSE}.
#' @param install_method Character string. Installation method when
#'   \code{install_missing = TRUE}. Options: "auto", "conda", "pip".
#'   Default: "auto".
#' @param create_conda_env Logical. If \code{TRUE}, creates a conda 
#'   environment with Python 3.8.10 for optimal GRNBoost2 compatibility.
#'   If conda is not available, automatically installs Miniconda first.
#'   Default: \code{FALSE}.
#' @param conda_env_name Character string. Name for the conda environment.
#'   Default: "scgraphverse-py38".
#' @param install_conda Logical. If \code{TRUE}, automatically installs 
#'   Miniconda when conda is not found. Default: \code{TRUE}.
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
#' @importFrom reticulate use_python import py_module_available py_install
#'   conda_install
#' @export
#'
#' @examples
#' # Initialize Python environment (adjust python_path as needed)
#' modules <- init_py()
#'
#' # Initialize with automatic installation of missing packages
#' modules <- init_py(install_missing = TRUE)
#' 
#' # Create conda environment with Python 3.8.10 for optimal compatibility
#' # This will automatically install Miniconda if not found
#' modules <- init_py(create_conda_env = TRUE, install_missing = TRUE)
#' 
#' # Use existing conda environment
#' modules <- init_py(create_conda_env = TRUE, 
#'                    conda_env_name = "my-py38-env",
#'                    install_missing = TRUE)
#'                    
#' # Disable automatic Miniconda installation
#' modules <- init_py(create_conda_env = TRUE, 
#'                    install_conda = FALSE,
#'                    install_missing = TRUE)
init_py <- function(
    python_path = "/usr/bin/python3",
    required = TRUE,
    install_missing = FALSE,
    install_method = "auto",
    create_conda_env = FALSE,
    conda_env_name = "scgraphverse-py38",
    install_conda = TRUE,
    verbose = TRUE) {
    
    # Create conda environment with Python 3.8.10 if requested
    if (create_conda_env) {
        if (verbose) {
            message("Creating conda environment '", conda_env_name, "' with Python 3.8.10...")
        }
        
        tryCatch({
            # Check if conda is available
            conda_available <- tryCatch({
                !reticulate::conda_binary() == ""
            }, error = function(e) FALSE)
            
            if (!conda_available && install_conda) {
                if (verbose) {
                    message("Conda not found. Installing Miniconda...")
                    message("This may take several minutes...")
                }
                
                # Install Miniconda using reticulate
                reticulate::install_miniconda()
                
                if (verbose) message("✓ Miniconda installed successfully")
                conda_available <- TRUE
            }
            
            if (conda_available) {
                # Check if environment already exists
                existing_envs <- reticulate::conda_list()
                if (!conda_env_name %in% existing_envs$name) {
                    if (verbose) {
                        message("Creating conda environment with Python 3.8.10...")
                    }
                    # Create new environment with Python 3.8.10
                    reticulate::conda_create(
                        envname = conda_env_name,
                        python_version = "3.8.10"
                    )
                    if (verbose) message("✓ Conda environment created successfully")
                } else {
                    if (verbose) message("✓ Conda environment already exists")
                }
                
                # Use the conda environment
                reticulate::use_condaenv(conda_env_name, required = TRUE)
                if (verbose) message("✓ Using conda environment: ", conda_env_name)
                
            } else {
                stop("Conda installation failed or install_conda=FALSE. Please install Anaconda/Miniconda manually.")
            }
        }, error = function(e) {
            if (verbose) {
                message("Failed to create conda environment: ", e$message)
                message("Falling back to specified python_path")
            }
            reticulate::use_python(python_path, required = required)
        })
    } else {
        reticulate::use_python(python_path, required = required)
    }
    
    # Check Python version compatibility for GRNBoost2
    if (verbose) {
        tryCatch({
            py_version <- reticulate::py_config()$version
            if (verbose) message("Using Python version: ", py_version)
            
            # Extract major.minor.patch version
            version_parts <- strsplit(py_version, "\\.")[[1]]
            major <- as.numeric(version_parts[1])
            minor <- as.numeric(version_parts[2])
            
            # Check if not using recommended Python 3.8.10
            if (!(major == 3 && minor == 8)) {
                message("INFO: Python ", major, ".", minor, " detected.")
                message("For optimal GRNBoost2 compatibility, Python 3.8.10 is recommended")
                message("On macOS: pyenv install 3.8.10 && pyenv global 3.8.10")
                message("On Windows: Download Python 3.8.10 from python.org/downloads/release/python-3810/")
            } else {
                message("✓ Using Python 3.8.x - optimal for GRNBoost2")
            }
        }, error = function(e) {
            if (verbose) message("Could not check Python version")
        })
    }

    if (install_missing && !reticulate::py_module_available("arboreto")) {
        if (verbose) {
            message("Installing Python package 'arboreto' for GRNBoost2...")
            message("This may take a few minutes.")
        }

        tryCatch(
            {
                if (install_method == "auto") {
                    if (verbose) message("Trying conda installation...")
                    tryCatch(
                        {
                            reticulate::conda_install("arboreto", 
                                                        channel = "bioconda")
                            if (verbose) message("Successfully installed.")
                        },
                        error = function(e) {
                            if (verbose) message("Conda failed, trying pip...")
                            reticulate::py_install("arboreto")
                            if (verbose) message("Successfully installed.")
                        }
                    )
                } else if (install_method == "conda") {
                    reticulate::conda_install("arboreto", channel = "bioconda")
                    if (verbose) message("Successfully installed via conda.")
                } else if (install_method == "pip") {
                    reticulate::py_install("arboreto")
                    if (verbose) message("Successfully installed via pip.")
                } else {
                    stop("Invalid install_method. Use 'auto',
                        'conda', or 'pip'.")
                }
            },
            error = function(e) {
                if (verbose) {
                    message("Installation failed: ", e$message)
                    message("Try manual installation: pip install arboreto")
                }
                if (required) {
                    stop("Failed to install Python package 'arboreto'.")
                }
            }
        )
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
