#' @importFrom cli rule
.onAttach <- function(libname, pkgname) {
    msg <- paste(
        cli::rule(center = "\u2728 Welcome to scGraphVerse \u2728", line = 2),
        "",
        "\U0001F680 Network analysis tools for single-cell data.",
        "",
        "This work is supported by the project:",
        "National Centre for HPC, Big Data and Quantum Computing",
        "Funded by European Union \u2013 Next Generation EU \u2013 CN00000013",
        "CUP: B93C22000620006",
        "",
        cli::rule(line = 2),
        sep = "\n"
    )
    packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
}

.onUnload <- function(libpath) {
}
