#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Full declaration of regRF */
void regRF(
    double *x, double *y, double *ww,
    int *xdim, int *nsample, int *nodesize, int *nrnodes,
    int *ntree, int *mtry, int *imp, int *cat, int *maxcat,
    int *do_trace, int *proximity, int *oob_prox, int *corr_bias,
    double *ypred, double *impout, double *impmat, double *impSD,
    double *prox, int *ndbigtree, int *nodestatus, int *leftDaughter,
    int *rightDaughter, double *nodepred, int *bestvar, double *xbestsplit,
    double *mse, int *keep, int *replace, int *testdat, double *xts,
    int *ntest, double *yts, int *labelts, double *ytestpred,
    double *proxts, double *msets, double *coef, int *oob_times,
    int *inbag, int *nclasses
);

static const R_CMethodDef CEntries[] = {
    {"regRF", (DL_FUNC) &regRF, 44},  /* matches your R call */
    {NULL, NULL, 0}
};

void R_init_scGraphVerse(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

