# Initialize Python Environment for GRNBoost2

Sets up the Python environment and lazily loads modules required for
running GRNBoost2: `arboreto`, `pandas`, and `numpy`. Automatically
installs missing Python packages if requested.

## Usage

``` r
init_py(
  python_path = "/usr/bin/python3",
  required = TRUE,
  install_missing = FALSE,
  install_method = "auto",
  verbose = TRUE
)
```

## Arguments

- python_path:

  Character string. Path to the Python executable, e.g.,
  `"/usr/bin/python3"`. For optimal GRNBoost2 compatibility, Python
  3.8.x is strongly recommended.

- required:

  Logical. If `TRUE`, errors if Python is not available or path is
  invalid. Default: `TRUE`.

- install_missing:

  Logical. If `TRUE`, automatically installs missing Python packages.
  Default: `FALSE`.

- install_method:

  Character string. Installation method when `install_missing = TRUE`.
  Options: "auto", "conda", "pip". Default: "auto".

- verbose:

  Logical. If `TRUE`, shows installation progress. Default: `TRUE`.

## Value

A list with three Python module objects:

- `arboreto`: GRNBoost2 algorithm module.

- `pandas`: Data handling module.

- `numpy`: Numerical operations module.

## Details

Uses reticulate to bind R to the specified Python interpreter and lazily
import modules needed for GRNBoost2. If `install_missing = TRUE`,
automatically installs the 'arboreto' package using the specified method
if not found.

## Examples

``` r
# Initialize Python environment (handles missing modules gracefully)
tryCatch(
    {
        modules <- init_py(required = FALSE)
    },
    error = function(e) {
        message("Python environment not available: ", e$message)
    }
)
#> Using Python version: 3.12
#> Could not check Python version
#> Python modules not available:ModuleNotFoundError: No module named 'arboreto'
#> Run `reticulate::py_last_error()` for details.
```
