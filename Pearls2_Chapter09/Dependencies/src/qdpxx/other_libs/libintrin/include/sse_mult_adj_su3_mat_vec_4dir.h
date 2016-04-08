#ifndef INTRIN_SSE_MULT_ADJ_SU3_MAT_VEC_4DIR_H
#define INTRIN_SSE_MULT_ADJ_SU3_MAT_VEC_4DIR_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_adj_su3_mat_vec_4dir(su3_matrixf aa[4], su3_vectorf *bb, su3_vectorf cc[4]);
#ifdef __cplusplus
}; 
#endif
#endif
