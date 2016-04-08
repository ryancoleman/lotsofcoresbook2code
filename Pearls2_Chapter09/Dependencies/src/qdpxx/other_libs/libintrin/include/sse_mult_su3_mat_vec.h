#ifndef INTRIN_SSE_MULT_SU3_MAT_VEC_H
#define INTRIN_SSE_MULT_SU3_MAT_VEC_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_su3_mat_vec(su3_matrixf *aa, su3_vectorf* bb, su3_vectorf* cc);
#ifdef __cplusplus
}; 
#endif
#endif
