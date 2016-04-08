#ifndef INTRIN_SSE_SCALAR_MULT_ADD_SU3_MATRIX_H
#define INTRIN_SSE_SCALAR_MULT_ADD_SU3_MATRIX_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_scalar_mult_add_su3_matrix(su3_matrixf* aa, su3_matrixf* bb, float cc, su3_matrixf* dd);
#ifdef __cplusplus
}; 
#endif
#endif
