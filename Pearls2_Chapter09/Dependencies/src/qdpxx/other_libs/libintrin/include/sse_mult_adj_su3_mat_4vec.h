#ifndef INTRIN_SSE_MULT_ADJ_SU3_MAT_4VEC_H
#define INTRIN_SSE_MULT_ADJ_SU3_MAT_4VEC_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_adj_su3_mat_4vec(su3_matrixf aa[4], su3_vectorf *bb, su3_vectorf *cc0, su3_vectorf *cc1, su3_vectorf *cc2, su3_vectorf *cc3);
#ifdef __cplusplus
}; 
#endif
#endif
