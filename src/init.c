
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP corr(SEXP, SEXP);
extern SEXP maxstatpermdist(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP newcorr(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"corr",            (DL_FUNC) &corr,            2},
    {"maxstatpermdist", (DL_FUNC) &maxstatpermdist, 7},
    {"newcorr",         (DL_FUNC) &newcorr,         2},
    {NULL, NULL, 0}
};

void R_init_maxstat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
