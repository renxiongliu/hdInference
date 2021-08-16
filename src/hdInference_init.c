#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void glasso_nonconvex_path(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *);
extern void glasso_nonconvex_path_new(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *);
extern void glasso_nonconvex_constrained_path(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *, void *,void *,void *);
extern void glasso_nonconvex_constrained_single(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *);
extern void glasso_nonconvex_path_v2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *);
extern void glasso_cMLE(void *, void *, void *, void *, void *, void *, void *,void *,void *,void *);
extern void glasso_nonconvex_path_inference(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *,void *,void *);
extern void glasso_nonconvex_path_inference_LR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *,void *,void *);
extern void glasso_nonconvex_path_inference_LR_v2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *,void *,void *);
extern void glasso_nonconvex_free(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *,void *,void *);
extern void glasso_nonconvex_zero(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *,void *,void *,void *,void *);

static const R_CMethodDef CEntries[] = {
  {"glasso_nonconvex_path", (DL_FUNC) &glasso_nonconvex_path, 14},
  {"glasso_nonconvex_path_new", (DL_FUNC) &glasso_nonconvex_path_new, 15},
  {"glasso_nonconvex_constrained_path", (DL_FUNC) &glasso_nonconvex_constrained_path, 18},
  {"glasso_nonconvex_constrained_single", (DL_FUNC) &glasso_nonconvex_constrained_single,15}, 
  {"glasso_nonconvex_path_v2", (DL_FUNC) &glasso_nonconvex_path_v2,15}, 
  {"glasso_cMLE", (DL_FUNC) &glasso_cMLE,10}, 
  {"glasso_nonconvex_path_inference", (DL_FUNC) &glasso_nonconvex_path_inference,17}, 
  {"glasso_nonconvex_path_inference_LR", (DL_FUNC) &glasso_nonconvex_path_inference_LR,17}, 
  {"glasso_nonconvex_path_inference_LR_v2", (DL_FUNC) &glasso_nonconvex_path_inference_LR_v2,17}, 
  {"glasso_nonconvex_free", (DL_FUNC) &glasso_nonconvex_free,17}, 
  {"glasso_nonconvex_zero", (DL_FUNC) &glasso_nonconvex_zero,17}, 
  {NULL, NULL, 0}
};

void R_init_hdInference(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
