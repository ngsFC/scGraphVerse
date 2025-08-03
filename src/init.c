#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "rf.h"

// Register native routines for .C calls
static const R_CMethodDef cMethods[] = {
    {"regRF", (DL_FUNC) &regRF, 44},
    {NULL, NULL, 0}
};

void R_init_scGraphVerse(DllInfo *info) {
    // Register C methods
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    
    // Print debug info (will appear when package loads)
    Rprintf("scGraphVerse: Registered %d C methods\n", 1);
}