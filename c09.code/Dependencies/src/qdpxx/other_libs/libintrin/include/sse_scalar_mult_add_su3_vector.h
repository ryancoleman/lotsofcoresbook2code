#ifndef INTRIN_SSE_SCALAR_MULT_ADD_SU3_VECTOR_H
#define INTRIN_SSE_SCALAR_MULT_ADD_SU3_VECTOR_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_scalar_mult_add_su3_vector(su3_vectorf* aa, su3_vectorf* bb, float cc, su3_vectorf* dd);
#ifdef __cplusplus
}; 
#endif
#endif
