#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(adapt)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(curbf)(void *, void *);
extern void F77_NAME(fdpss)(void *, void *, void *, void *, void *);
extern void F77_NAME(jkcoh1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mw2jkw)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mw2wta)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mweave)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(northf)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(quicksinef)(void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"adapt",      (DL_FUNC) &F77_NAME(adapt),       9},
    {"curbf",      (DL_FUNC) &F77_NAME(curbf),       2},
    {"fdpss",      (DL_FUNC) &F77_NAME(fdpss),       5},
    {"jkcoh1",     (DL_FUNC) &F77_NAME(jkcoh1),     19},
    {"mw2jkw",     (DL_FUNC) &F77_NAME(mw2jkw),     20},
    {"mw2wta",     (DL_FUNC) &F77_NAME(mw2wta),     15},
    {"mweave",     (DL_FUNC) &F77_NAME(mweave),     10},
    {"northf",     (DL_FUNC) &F77_NAME(northf),      6},
    {"quicksinef", (DL_FUNC) &F77_NAME(quicksinef),  7},
    {NULL, NULL, 0}
};

void R_init_multitaper(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
