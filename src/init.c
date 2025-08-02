#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "rf.h"

// Register native routines
static const R_CallMethodDef callMethods[] = {
    {"regTree", (DL_FUNC) &regTree, 19},
    {NULL, NULL, 0}
};

void R_init_scGraphVerse(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}