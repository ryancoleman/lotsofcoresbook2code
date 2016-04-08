#ifndef INTRIN_SSE_MULT_ADJ_SU3_MAT_HWVEC_H
#define INTRIN_SSE_MULT_ADJ_SU3_MAT_HWVEC_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_adj_su3_mat_hwvec(su3_matrixf *aa, half_wilson_vectorf *bb, half_wilson_vectorf *cc);
#ifdef __cplusplus
}; 
#endif
#endif
