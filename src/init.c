#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Full declaration of regRF - matches original JRF signature */
void regRF(
    double *x, double *y, double *weight, 
    int *xdim, int *sampsize, int *totsize,
    int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
    int *cat, int *maxcat, int *jprint, int *doProx, int *oobprox,
    int *biasCorr, double *yptr, double *errimp, double *impmat,
    double *impSD, double *prox, int *treeSize, int *nodestatus,
    int *lDaughter, int *rDaughter, double *avnode, int *mbest,
    double *upper, double *mse, int *keepf, int *replace,
    int *testdat, double *xts, int *nts, double *yts, int *labelts,
    double *yTestPred, double *proxts, double *msets, double *coef,
    int *nout, int *inbag, int *nclasses
);

static const R_CMethodDef CEntries[] = {
    {"regRF", (DL_FUNC) &regRF, 44},  /* matches original JRF signature */
    {NULL, NULL, 0}
};

void R_init_scGraphVerse(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

