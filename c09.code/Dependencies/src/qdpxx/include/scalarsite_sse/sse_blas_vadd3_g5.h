// $Id: sse_blas_vadd3_g5.h,v 1.3 2007-06-10 14:32:11 edwards Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VADD_G5
#define QDP_SSE_BLAS_VADD_G5

#if defined(__GNUC__)

#include "qdp_config.h"
#include "qdp_sse_intrin.h"
namespace QDP {

#if BASE_PRECISION==32



// (Vector) out = X + P_{+} (Vector) Y 
inline
void add_g5ProjPlus(REAL32 *Out, REAL32 *X, REAL32 *Y,int n_4vec)
{

  for(int i=0; i < n_4vec; i++) {

    // Spin Component 0: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+ 0, _mm_add_ps(_mm_load_ps(Y+ 0), _mm_load_ps(X+ 0)));
     
    // Spin Component 0: z2r, z2i, SpinComponent 1: z0r, z0i
    _mm_store_ps(Out+ 4, _mm_add_ps(_mm_load_ps(Y+ 4), _mm_load_ps(X+ 4)));


    // Spin Component 1: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+ 8, _mm_add_ps(_mm_load_ps(Y+ 8), _mm_load_ps(X+ 8)));

    // Spin Component 2: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+12, _mm_load_ps(X+ 12));

    // Spin Component 2: z2r, z2i, z0r, z0r
    _mm_store_ps(Out+16, _mm_load_ps(X+ 16));

    // Spin Component 3: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+20, _mm_load_ps(X+ 20));

    // Update offsets
    Out += 24; Y += 24; X += 24;
  }

}

// (Vector) out = X + (Vector) P{-} Y
inline
void add_g5ProjMinus(REAL32 *Out,REAL32 *X, REAL32 *Y,int n_4vec)
{

  for(int i=0; i < n_4vec; i++) {
    // Spin Component 0: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+0,_mm_load_ps(X+ 0));

    // Spin Component 0: z2r, z2i, Spin Component 1: z0r, z0r
    _mm_store_ps(Out+4,_mm_load_ps(X+ 4));

    // Spin Component 1: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+8, _mm_load_ps(X+ 8));


    // Spin Component 2: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+12, _mm_add_ps(_mm_load_ps(X+ 12), _mm_load_ps(Y+ 12)));
     
    // Spin Component 2: z2r, z2i, SpinComponent 3: z0r, z0i
    _mm_store_ps(Out+16, _mm_add_ps(_mm_load_ps(X+ 16), _mm_load_ps(Y+ 16)));

    // Spin Component 3: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+20, _mm_add_ps(_mm_load_ps(X+ 20), _mm_load_ps(Y+ 20)));

    // Update offsets
    Out += 24; Y += 24; X += 24;
  }

}



// AXMY  versions
// (Vector) out = (Scalar) (*scalep) * (Vector) Y - (Vector) P{+} X
inline
void sub_g5ProjPlus(REAL32 *Out, REAL32 *X, REAL32 *Y,int n_4vec)
{
  for(int i=0; i < n_4vec; i++) {
    // Spin Component 0: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+ 0, _mm_sub_ps(_mm_load_ps(X+ 0), _mm_load_ps(Y + 0)));
     
    // Spin Component 0: z2r, z2i, SpinComponent 1: z0r, z0i
    _mm_store_ps(Out+ 4, _mm_sub_ps(_mm_load_ps(X+ 4), _mm_load_ps(Y + 4)));


    // Spin Component 1: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+ 8, _mm_sub_ps(_mm_load_ps(X+ 8), _mm_load_ps(Y+ 8)));
    
  
    // Spin Component 2: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+12, _mm_load_ps(X + 12));
        
    
    // Spin Component 2: z2r, z2i, z0r, z0r
    _mm_store_ps(Out+16, _mm_load_ps(X+ 16));


    // Spin Component 3: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+20, _mm_load_ps(X+20));

    // Update offsets
    Out += 24; Y += 24; X += 24;
  }


}

// (Vector) out = Y - (Vector) P{-} X
inline
void sub_g5ProjMinus(REAL32 *Out,REAL32 *X, REAL32 *Y,int n_4vec)
{

  for(int i=0; i < n_4vec; i++) {
    // Spin Component 0: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+0, _mm_load_ps(X+0));

    // Spin Component 0: z2r, z2i, Spin Component 1: z0r, z0r
    _mm_store_ps(Out+4, _mm_load_ps(X+4));

    // Spin Component 1: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+8, _mm_load_ps(X+8));

    // Spin Component 2: z0r, z0i, z1r, z1i
    _mm_store_ps(Out+12, _mm_sub_ps(_mm_load_ps(X+ 12), _mm_load_ps(Y+ 12)));
     
    // Spin Component 2: z2r, z2i, SpinComponent 3: z0r, z0i
    _mm_store_ps(Out+16, _mm_sub_ps(_mm_load_ps(X+ 16), _mm_load_ps(Y+ 16)));

    // Spin Component 3: z1r, z1i, z2r, z2i
    _mm_store_ps(Out+20, _mm_sub_ps(_mm_load_ps(X+ 20), _mm_load_ps(Y+ 20)));

    // Update offsets
    Out += 24; Y += 24; X += 24;
  }

}


#endif

} // namespace QDP;


#endif // GNUC

#endif // guard
