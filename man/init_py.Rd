% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_py.R
\name{init_py}
\alias{init_py}
\title{Initialize Python Environment for GRNBoost2}
\usage{
init_py(python_path = "/usr/bin/python3", required = TRUE)
}
\arguments{
\item{python_path}{Character string. Path to the Python executable,
e.g., \code{"/usr/bin/python3"}.}

\item{required}{Logical. If \code{TRUE}, errors if Python is not
available or path is invalid. Default: \code{TRUE}.}
}
\value{
A list with three Python module objects:
\itemize{
\item \code{arboreto}: GRNBoost2 algorithm module.
\item \code{pandas}: Data handling module.
\item \code{numpy}: Numerical operations module.
}
}
\description{
Sets up the Python environment and lazily loads modules required for
running GRNBoost2: \code{arboreto}, \code{pandas}, and \code{numpy}.
}
\details{
Uses \pkg{reticulate} to bind R to the specified Python
interpreter and lazily import modules needed for GRNBoost2. If a module
is missing or incompatible, an informative error is raised (when
\code{required = TRUE}).
}
\examples{
# Initialize Python environment (adjust python_path as needed)
modules <- init_py()
}
