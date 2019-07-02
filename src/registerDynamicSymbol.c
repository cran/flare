#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void SFGen(void *, void *, void *, void *);
extern void slim_dantzig_ladm_scr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void slim_dantzig_ladm_scr2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void slim_lad_ladm_scr_btr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void slim_lasso_ladm_scr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void slim_lq_ladm_scr_btr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void slim_sqrt_ladm_scr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sugm_clime_ladm_scr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sugm_tiger_ladm_scr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"SFGen",                  (DL_FUNC) &SFGen,                   4},
    {"slim_dantzig_ladm_scr",  (DL_FUNC) &slim_dantzig_ladm_scr,  21},
    {"slim_dantzig_ladm_scr2", (DL_FUNC) &slim_dantzig_ladm_scr2, 18},
    {"slim_lad_ladm_scr_btr",  (DL_FUNC) &slim_lad_ladm_scr_btr,  21},
    {"slim_lasso_ladm_scr",    (DL_FUNC) &slim_lasso_ladm_scr,    20},
    {"slim_lq_ladm_scr_btr",   (DL_FUNC) &slim_lq_ladm_scr_btr,   22},
    {"slim_sqrt_ladm_scr",     (DL_FUNC) &slim_sqrt_ladm_scr,     21},
    {"sugm_clime_ladm_scr",    (DL_FUNC) &sugm_clime_ladm_scr,    22},
    {"sugm_tiger_ladm_scr",    (DL_FUNC) &sugm_tiger_ladm_scr,    28},
    {NULL, NULL, 0}
};

void R_init_flare(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
