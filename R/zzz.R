# Package load hook to load compiled C code
.onLoad <- function(libname, pkgname) {
    library.dynam("scGraphVerse", pkgname, libname)
}

.onUnload <- function(libpath) {
    library.dynam.unload("scGraphVerse", libpath)
}
