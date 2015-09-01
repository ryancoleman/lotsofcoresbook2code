#ifndef INTRIN_SSE_MULT_SU3_NN_H
#define INTRIN_SSE_MULT_SU3_NN_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_mult_su3_nn(su3_matrixf* aa, su3_matrixf* bb, su3_matrixf* cc);
#ifdef __cplusplus
}; 
#endif
#endif
