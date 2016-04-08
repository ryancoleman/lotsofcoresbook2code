#ifndef INTRIN_SSE_SU3_PROJECTOR_H
#define INTRIN_SSE_SU3_PROJECTOR_H

#include "inline_sse.h" 
#ifdef __cplusplus
extern "C" { 
#endif

void
intrin_sse_su3_projector(su3_vectorf* aa, su3_vectorf* bb, su3_matrixf* cc);
#ifdef __cplusplus
}; 
#endif
#endif
