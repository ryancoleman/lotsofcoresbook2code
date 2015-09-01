// $Id: qdp_scalarvecsite_sse_blas.h,v 1.5 2007-06-10 14:32:11 edwards Exp $
/*! @file
 * @brief Blas optimizations
 * 
 * Generic and maybe SSE optimizations of basic operations
 */


#ifndef QDP_SCALARVECSITE_SSE_BLAS_H
#define QDP_SCALARVECSITE_SSE_BLAS_H

namespace QDP {

// #define QDP_SCALARVECSITE_BLAS_DEBUG

#define QDP_SCALARVECSITE_USE_EVALUATE
#include "qdp_sse_intrin.h"

typedef PScalar<PColorMatrix<RComplex<ILattice<REAL32,4> >, 3> >       TCMat;
typedef PScalar<PColorVector<RComplex<ILattice<REAL32,4> >, 3> >       TCVec;
typedef PSpinVector<PColorVector<RComplex<ILattice<REAL32,4> >, 3>, 4> TDirac;
typedef PScalar<PScalar<RScalar<IScalar<REAL32> > > >                  TScal;


// d = Scalar*ColorMatrix
template<>
inline BinaryReturn<TScal, TCMat, OpMultiply>::Type_t
operator*(const TScal& l, const TCMat& r)
{
  BinaryReturn<TScal, TCMat, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARVECSITE_BLAS_DEBUG)
  cout << "M*S" << endl;
#endif

  vReal32 vscale = vmk1(l.elem().elem().elem().elem());

  d.elem().elem(0,0).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(0,0).real().elem_v());
  d.elem().elem(0,0).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(0,0).imag().elem_v());
  d.elem().elem(0,1).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(0,1).real().elem_v());
  d.elem().elem(0,1).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(0,1).imag().elem_v());
  d.elem().elem(0,2).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(0,2).real().elem_v());
  d.elem().elem(0,2).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(0,2).imag().elem_v());
  d.elem().elem(1,0).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(1,0).real().elem_v());
  d.elem().elem(1,0).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(1,0).imag().elem_v());
  d.elem().elem(1,1).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(1,1).real().elem_v());
  d.elem().elem(1,1).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(1,1).imag().elem_v());
  d.elem().elem(1,2).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(1,2).real().elem_v());
  d.elem().elem(1,2).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(1,2).imag().elem_v());
  d.elem().elem(2,0).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(2,0).real().elem_v());
  d.elem().elem(2,0).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(2,0).imag().elem_v());
  d.elem().elem(2,1).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(2,1).real().elem_v());
  d.elem().elem(2,1).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(2,1).imag().elem_v());
  d.elem().elem(2,2).real().elem_v() = _mm_mul_ps(vscale, r.elem().elem(2,2).real().elem_v());
  d.elem().elem(2,2).imag().elem_v() = _mm_mul_ps(vscale, r.elem().elem(2,2).imag().elem_v());

  return d;
}


#if defined(QDP_SCALARVECSITE_BLAS_DEBUG)
#undef QDP_SCALARVECSITE_BLAS_DEBUG
#endif

#if defined(QDP_SCALARVECSITE_USE_EVALUATE)
#undef QDP_SCALARVECSITE_USE_EVALUATE
#endif


} // namespace QDP;

#endif  // guard
 
