#ifndef INTRIN_SSE_MULT_SU3_MAT_VEC_SUM_4DIR_H
#define INTRIN_SSE_MULT_SU3_MAT_VEC_SUM_4DIR_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_su3_mat_vec_sum_4dir(su3_matrixf aa[4], su3_vectorf* bb0, su3_vectorf* bb1, su3_vectorf* bb2, su3_vectorf* bb3, su3_vectorf* cc);
#ifdef __cplusplus
}; 
#endif
#endif
