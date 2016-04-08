// $Id: sse_blas_vscal3_g5.h,v 1.4 2007-08-20 17:08:14 uid4709 Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VSCAL_G5
#define QDP_SSE_BLAS_VSCAL_G5

#if defined(__GNUC__)

#include "qdp_config.h"
#include "qdp_sse_intrin.h"

namespace QDP {

#if BASE_PRECISION==32



// (Vector) out = (*scalep)* P_{+} X
inline
void scal_g5ProjPlus(REAL32 *Out, REAL32* scalep, REAL32 *X, int n_4vec)
{
  // GNUC vector type  
  v4sf vscalep = _mm_load_ss(scalep); // High bytes 
  vscalep =_mm_shuffle_ps(vscalep, vscalep, 0 );

  REAL32 rzero=0;

  // A zero vector
  v4sf vzero = _mm_load_ss(&rzero);
  vzero = _mm_shuffle_ps(vzero, vzero, 0);

  for(int i=0; i < n_4vec; i++) {

    // Spin Component 0: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+ 0, _mm_mul_ps(vscalep, _mm_load_ps(X+ 0)));
     
    // Spin Component 0: z2r, z2i, SpinComponent 1: z0r, z0i
    _mm_store_ps(Out+ 4, _mm_mul_ps(vscalep, _mm_load_ps(X+ 4)));


    // Spin Component 1: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+ 8, _mm_mul_ps(vscalep, _mm_load_ps(X+ 8)));

    
    // Spin Component 2: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+12, vzero);

    // Spin Component 2: z2r, z2i, z0r, z0r
    _mm_store_ps(Out+16, vzero);

    // Spin Component 3: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+20, vzero);

    // Update offsets
    Out += 24;  X += 24;
  }

}


// (Vector) out = (*scalep)* P_{-} X
inline
void scal_g5ProjMinus(REAL32 *Out, REAL32* scalep, REAL32 *X, int n_4vec)
{
  // GNUC vector type
  

  v4sf vscalep = _mm_load_ss(scalep);
  vscalep = _mm_shuffle_ps(vscalep, vscalep,0);
  REAL32 rzero = (REAL32)0;
  
  // A zero vector
  v4sf vzero = _mm_load_ss(&rzero);
  vzero = _mm_shuffle_ps(vzero, vzero, 0);
  
  for(int i=0; i < n_4vec; i++) {
    // Spin Component 0: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+0, vzero);

    // Spin Component 0: z2r, z2i, Spin Component 1: z0r, z0r
    _mm_store_ps(Out+4, vzero);

    // Spin Component 1: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+8, vzero);


    // Spin Component 2: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+ 12, _mm_mul_ps(vscalep, _mm_load_ps(X+ 12)));
     
    // Spin Component 2: z2r, z2i, SpinComponent 3: z0r, z0i
    _mm_store_ps(Out+ 16, _mm_mul_ps(vscalep, _mm_load_ps(X+ 16)));


    // Spin Component 3: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+ 20, _mm_mul_ps(vscalep, _mm_load_ps(X+ 20)));

    

    // Update offsets
    Out += 24;  X += 24;
  }

}


#endif

} // namespace QDP;


#endif // GNUC

#endif // guard
