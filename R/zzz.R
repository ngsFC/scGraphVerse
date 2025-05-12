.onAttach <- function(libname, pkgname) {
    msg <- paste(
        cli::rule(center = "✨ Welcome to scGraphVerse ✨", line = 2),
        "",
        "🚀 Network analysis tools for single-cell data.",
        "",
        "This work is supported by the project:",
        "National Centre for HPC, Big Data and Quantum Computing",
        "Funded by European Union – Next Generation EU – CN00000013 🇪🇺",
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
    library.dynam.unload("scGraphVerse", libpath)
}
