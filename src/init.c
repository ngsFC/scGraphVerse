#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void regRF(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, int *, double *, double *, int *, int *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, int *, int *);

static const R_CMethodDef CEntries[] = {
    {"regRF", (DL_FUNC) &regRF, 44},
    {NULL, NULL, 0}
};

void R_init_scGraphVerse(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}