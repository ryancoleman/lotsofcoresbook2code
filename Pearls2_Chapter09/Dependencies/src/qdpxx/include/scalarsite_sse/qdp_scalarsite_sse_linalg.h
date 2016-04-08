// -*- C++ -*-
// $Id: qdp_scalarsite_sse_linalg.h,v 1.17 2009-04-17 00:46:36 bjoo Exp $

/*! @file
 * @brief Blas optimizations
 *
 * Blas optimizations of basic operations
 */

#ifndef QDP_SCALARSITE_SSE_LINALG_H
#define QDP_SCALARSITE_SSE_LINALG_H

// These SSE asm instructions are only supported under GCC/G++
#if defined(__GNUC__)

namespace QDP {


// #define QDP_SCALARSITE_DEBUG

#define QDP_SCALARSITE_USE_EVALUATE


/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */

// Use this def just to safe some typing later on in the file
typedef RComplex<REAL32>  RComplexFloat;

// Types needed for the expression templates. 
typedef PScalar<PColorMatrix<RComplexFloat, 3> > TCol;
typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 2> TVec2;
typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 4> TVec4;



#include "sse_mult_su3_nn.h"
#include "sse_mult_su3_an.h"

#include "sse_mult_su3_na.h"
#include "sse_mult_su3_mat_vec.h"
#include "sse_mult_adj_su3_mat_vec.h"
#include "sse_add_su3_vector.h"
#include "sse_mult_su3_mat_hwvec.h"
#include "sse_mult_adj_su3_mat_hwvec.h"


// Optimized version of  
//    PColorMatrix<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * PColorMatrix<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiply>::Type_t
operator*(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	  const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*M" << endl;
#endif
  
  su3_matrixf* lm = (su3_matrixf *) &(l.elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *) &(r.elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *) &(d.elem(0,0).real());

  intrin_sse_mult_su3_nn(lm,rm,dm);

  return d;
}



// Optimized version of  
//    PScalar<PColorMatrix<RComplexFloat,3>> <- PScalar<PColorMatrix<RComplexFloat,3>> * 
//                         PScalar<PColorMatrix<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >,
  PScalar<PColorMatrix<RComplexFloat,3> >, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	  const PScalar<PColorMatrix<RComplexFloat,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PScalar<PColorMatrix<RComplexFloat,3> >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*PSc<M>" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *)&(r.elem().elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *)&(d.elem().elem(0,0).real());

  intrin_sse_mult_su3_nn(lm, rm, dm);

  return d;
}


// Optimized version of  
//    PScalar<PScalar<PColorMatrix<RComplexFloat,3>>> <- PScalar<PScalar<PColorMatrix<RComplexFloat,3>>> * 
//                         PScalar<PScalar<PColorMatrix<RComplexFloat,3>>>
template<>
inline BinaryReturn<PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >,
  PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, OpMultiply>::Type_t
operator*(const PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >& l, 
	  const PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >& r)
{
  BinaryReturn<PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, 
    PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<PSc<M>>*PSc<PSc<M>>" << endl;
#endif

  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem().elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *)&(r.elem().elem().elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *)&(d.elem().elem().elem(0,0).real());

  intrin_sse_mult_su3_nn(lm, rm, dm);

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- adj(PColorMatrix<RComplexFloat,3>) * PColorMatrix<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	    const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*M" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *) &(l.elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *) &(r.elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *) &(d.elem(0,0).real());

  intrin_sse_mult_su3_an(lm,rm,dm);

  return d;
}


// Optimized version of  
//   PScalar<PColorMatrix<RComplexFloat,3>> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PScalar<PColorMatrix<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PScalar<PColorMatrix<RComplexFloat,3> >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PScalar<PColorMatrix<RComplexFloat,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PScalar<PColorMatrix<RComplexFloat,3> >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*PSc<M>" << endl;
#endif

  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *)&(r.elem().elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *)&(d.elem().elem(0,0).real());

  intrin_sse_mult_su3_an(lm, rm, dm);


  return d;
}


// Optimized version of  
//   PScalar<PScalar<PColorMatrix<RComplexFloat,3>>> <- adj(PScalar<PScalar<PColorMatrix<RComplexFloat,3>>>) * 
//        PScalar<PScalar<PColorMatrix<RComplexFloat,3>>>
template<>
inline BinaryReturn<PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, 
  PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >& l, 
	    const PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >& r)
{
  BinaryReturn<PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, 
    PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<PSc<M>>)*PSc<PSc<M>>" << endl;
#endif

  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem().elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *)&(r.elem().elem().elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *)&(d.elem().elem().elem(0,0).real());

  intrin_sse_mult_su3_an(lm, rm, dm);

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * adj(PColorMatrix<RComplexFloat,3>)
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiplyAdj>::Type_t
multiplyAdj(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	    const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*adj(M)" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *) &(l.elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *) &(r.elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *) &(d.elem(0,0).real());

  intrin_sse_mult_su3_na(lm,rm,dm);


  return d;
}


// Optimized version of  
//   PScalar<PColorMatrix<RComplexFloat,3>> <- PScalar<PColorMatrix<RComplexFloat,3>> * 
//          adj(PScalar<PColorMatrix<RComplexFloat,3>>)
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PScalar<PColorMatrix<RComplexFloat,3> >, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PScalar<PColorMatrix<RComplexFloat,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PScalar<PColorMatrix<RComplexFloat,3> >, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*adj(PSc<M>)" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *)&(r.elem().elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *)&(d.elem().elem(0,0).real());

  intrin_sse_mult_su3_na(lm, rm, dm);

  return d;
}


// Optimized version of  
//   PScalar<Pscalar<PColorMatrix<RComplexFloat,3>>> <- PScalar<PScalar<PColorMatrix<RComplexFloat,3>>> * 
//           adj(PScalar<PScalar<PColorMatrix<RComplexFloat,3>>>)
template<>
inline BinaryReturn<PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, 
  PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >& l, 
	    const PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >& r)
{
  BinaryReturn<PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, 
    PScalar<PScalar<PColorMatrix<RComplexFloat,3> > >, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<PSc<M>>*adj(PSc<PSc<M>>)" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem().elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *)&(r.elem().elem().elem(0,0).real());
  su3_matrixf* dm = (su3_matrixf *)&(d.elem().elem().elem(0,0).real());

  intrin_sse_mult_su3_na(lm, rm, dm);

  return d;
}


// Ooops, this macro does not exist!!

// Optimized version of  
//   PColorMatrix<RComplexFloat,3> <- adj(PColorMatrix<RComplexFloat,3>) * adj(PColorMatrix<RComplexFloat,3>)
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	       const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpAdjMultiplyAdj>::Type_t  d;

  PColorMatrix<RComplexFloat,3> tmp;

//

  su3_matrixf* lm = (su3_matrixf *) &(l.elem(0,0).real());
  su3_matrixf* rm = (su3_matrixf *) &(r.elem(0,0).real());
  su3_matrixf* tm = (su3_matrixf *) &(tmp.elem(0,0).real());

  // Do the adj*adj the hard way

  intrin_sse_mult_su3_nn(rm,lm,tm);

  // Take the adj(r*l) = adj(l)*adj(r)
  d.elem(0,0).real() =  tmp.elem(0,0).real();
  d.elem(0,0).imag() = -tmp.elem(0,0).imag();
  d.elem(0,1).real() =  tmp.elem(1,0).real();
  d.elem(0,1).imag() = -tmp.elem(1,0).imag();
  d.elem(0,2).real() =  tmp.elem(2,0).real();
  d.elem(0,2).imag() = -tmp.elem(2,0).imag();

  d.elem(1,0).real() =  tmp.elem(0,1).real();
  d.elem(1,0).imag() = -tmp.elem(0,1).imag();
  d.elem(1,1).real() =  tmp.elem(1,1).real();
  d.elem(1,1).imag() = -tmp.elem(1,1).imag();
  d.elem(1,2).real() =  tmp.elem(2,1).real();
  d.elem(1,2).imag() = -tmp.elem(2,1).imag();

  d.elem(2,0).real() =  tmp.elem(0,2).real();
  d.elem(2,0).imag() = -tmp.elem(0,2).imag();
  d.elem(2,1).real() =  tmp.elem(1,2).real();
  d.elem(2,1).imag() = -tmp.elem(1,2).imag();
  d.elem(2,2).real() =  tmp.elem(2,2).real();
  d.elem(2,2).imag() = -tmp.elem(2,2).imag();

  return d;
}


// Optimized version of  
//    PColorVector<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * PColorVector<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PVector<RComplexFloat,3,PColorVector>, OpMultiply>::Type_t
operator*(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	  const PVector<RComplexFloat,3,PColorVector>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PVector<RComplexFloat,3,PColorVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*V" << endl;
#endif

  su3_matrixf* lm = (su3_matrixf *)&(l.elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of  
//    PScalar<PColorVector<RComplexFloat,3>> <- PScalar<PColorMatrix<RComplexFloat,3>> * PScalar<PColorVector<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PScalar<PColorVector<RComplexFloat,3> >, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	  const PScalar<PColorVector<RComplexFloat,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PScalar<PColorVector<RComplexFloat,3> >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*PSc<V>" << endl;
#endif

  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem().elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem().elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of  
//    PColorVector<RComplexFloat,3> <- adj(PColorMatrix<RComplexFloat,3>) * PColorVector<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PVector<RComplexFloat,3,PColorVector>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	    const PVector<RComplexFloat,3,PColorVector>& r)
{
  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PVector<RComplexFloat,3,PColorVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*V" << endl;
#endif

  su3_matrixf* lm = (su3_matrixf *)&(l.elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of   StaggeredFermion <- ColorMatrix*StaggeredFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,1> <- PScalar<PColorMatrix<RComplexFloat,3>> * PSpinVector<PColorVector<RComplexFloat,3>,1>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PSpinVector<PColorVector<RComplexFloat,3>,1>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	  const PSpinVector<PColorVector<RComplexFloat,3>,1>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PSpinVector<PColorVector<RComplexFloat,3>,1>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*S" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of  
//    PScalar<PColorVector<RComplexFloat,3>> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PScalar<PColorVector<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PScalar<PColorVector<RComplexFloat,3> >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PScalar<PColorVector<RComplexFloat,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PScalar<PColorVector<RComplexFloat,3> >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*PSc<V>" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem().elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem().elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of   StaggeredFermion <- adj(ColorMatrix)*StaggeredFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,1> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PSpinVector<PColorVector<RComplexFloat,3>,1>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PSpinVector<PColorVector<RComplexFloat,3>,1>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PSpinVector<PColorVector<RComplexFloat,3>,1>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PSpinVector<PColorVector<RComplexFloat,3>,1>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*S" << endl;
#endif
  su3_matrixf* lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,2> <- PScalar<PColorMatrix<RComplexFloat,3>> * 
//                     PSpinVector<ColorVector<RComplexFloat,3>,2>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PSpinVector<PColorVector<RComplexFloat,3>,2>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
          const PSpinVector<PColorVector<RComplexFloat,3>,2>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PSpinVector<PColorVector<RComplexFloat,3>,2>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*H" << endl;
#endif
  su3_matrixf* lm =(su3_matrixf *)&( l.elem().elem(0,0).real() );
  half_wilson_vectorf* rh = (half_wilson_vectorf *)&(r.elem(0).elem(0).real()  );
  half_wilson_vectorf* dh = (half_wilson_vectorf *)&(d.elem(0).elem(0).real() );

  intrin_sse_mult_su3_mat_hwvec(lm,rh,dh);

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,2> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * 
//                     PSpinVector<ColorVector<RComplexFloat,3>,2>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PSpinVector<PColorVector<RComplexFloat,3>,2>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
            const PSpinVector<PColorVector<RComplexFloat,3>,2>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PSpinVector<PColorVector<RComplexFloat,3>,2>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*H" << endl;
#endif

  su3_matrixf* lm =(su3_matrixf *)&( l.elem().elem(0,0).real() );
  half_wilson_vectorf* rh = (half_wilson_vectorf *)&(r.elem(0).elem(0).real()  );
  half_wilson_vectorf* dh = (half_wilson_vectorf *)&(d.elem(0).elem(0).real() );

  intrin_sse_mult_adj_su3_mat_hwvec(lm,rh,dh);

  return d;
}


// Optimized version of  
//    PColorVector<RComplexFloat,3> <- PColorVector<RComplexFloat,3> + PColorVector<RComplexFloat,3>
template<>
inline BinaryReturn<PVector<RComplexFloat,3,PColorVector>, 
  PVector<RComplexFloat,3,PColorVector>, OpAdd>::Type_t
operator+(const PVector<RComplexFloat,3,PColorVector>& l, 
	  const PVector<RComplexFloat,3,PColorVector>& r)
{
  BinaryReturn<PVector<RComplexFloat,3,PColorVector>, 
    PVector<RComplexFloat,3,PColorVector>, OpAdd>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "V+V" << endl;
#endif

  su3_vectorf *lv = (su3_vectorf *)&( l.elem(0).real() );
  su3_vectorf *rv = (su3_vectorf *)&( r.elem(0).real() );
  su3_vectorf *dv = (su3_vectorf *)&( d.elem(0).real() );

  intrin_sse_add_su3_vector(lv,rv,dv);

  return d;
}


// Optimized version of   DiracFermion <- ColorMatrix*DiracFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,4> <- PScalar<PColorMatrix<RComplexFloat,3>> 
//                           * PSpinVector<PColorVector<RComplexFloat,3>,4>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PSpinVector<PColorVector<RComplexFloat,3>,4>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	  const PSpinVector<PColorVector<RComplexFloat,3>,4>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PSpinVector<PColorVector<RComplexFloat,3>,4>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*D" << endl;
#endif
  su3_matrixf *lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_vectorf *rv0 = (su3_vectorf *)&(r.elem(0).elem(0).real());
  su3_vectorf *rv1 = (su3_vectorf *)&(r.elem(1).elem(0).real());
  su3_vectorf *rv2 = (su3_vectorf *)&(r.elem(2).elem(0).real());
  su3_vectorf *rv3 = (su3_vectorf *)&(r.elem(3).elem(0).real());

  su3_vectorf *dv0 = (su3_vectorf *)&(d.elem(0).elem(0).real());
  su3_vectorf *dv1 = (su3_vectorf *)&(d.elem(1).elem(0).real());
  su3_vectorf *dv2 = (su3_vectorf *)&(d.elem(2).elem(0).real());
  su3_vectorf *dv3 = (su3_vectorf *)&(d.elem(3).elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv0,dv0);
  intrin_sse_mult_su3_mat_vec(lm,rv1,dv1);
  intrin_sse_mult_su3_mat_vec(lm,rv2,dv2);
  intrin_sse_mult_su3_mat_vec(lm,rv3,dv3);

  return d;
}


// Optimized version of   DiracFermion <- adj(ColorMatrix)*DiracFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,4> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>)
//                           * PSpinVector<PColorVector<RComplexFloat,3>,4>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PSpinVector<PColorVector<RComplexFloat,3>,4>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PSpinVector<PColorVector<RComplexFloat,3>,4>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PSpinVector<PColorVector<RComplexFloat,3>,4>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*D" << endl;
#endif
  su3_matrixf *lm = (su3_matrixf *)&(l.elem().elem(0,0).real());
  su3_vectorf *rv0 = (su3_vectorf *)&(r.elem(0).elem(0).real());
  su3_vectorf *rv1 = (su3_vectorf *)&(r.elem(1).elem(0).real());
  su3_vectorf *rv2 = (su3_vectorf *)&(r.elem(2).elem(0).real());
  su3_vectorf *rv3 = (su3_vectorf *)&(r.elem(3).elem(0).real());

  su3_vectorf *dv0 = (su3_vectorf *)&(d.elem(0).elem(0).real());
  su3_vectorf *dv1 = (su3_vectorf *)&(d.elem(1).elem(0).real());
  su3_vectorf *dv2 = (su3_vectorf *)&(d.elem(2).elem(0).real());
  su3_vectorf *dv3 = (su3_vectorf *)&(d.elem(3).elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv0,dv0);
  intrin_sse_mult_adj_su3_mat_vec(lm,rv1,dv1);
  intrin_sse_mult_adj_su3_mat_vec(lm,rv2,dv2);
  intrin_sse_mult_adj_su3_mat_vec(lm,rv3,dv3);

  return d;
}


// Optimized version of  
//    PScalar<PColorVector<RComplexFloat,3>> <- PScalar<PColorVector<RComplexFloat,3>> 
//                                            + PScalar<PColorVector<RComplexFloat,3>>
template<>
inline BinaryReturn<PScalar<PColorVector<RComplexFloat,3> >, 
  PScalar<PColorVector<RComplexFloat,3> >, OpAdd>::Type_t
operator+(const PScalar<PColorVector<RComplexFloat,3> >& l, 
	  const PScalar<PColorVector<RComplexFloat,3> >& r)
{
  BinaryReturn<PScalar<PColorVector<RComplexFloat,3> >, 
    PScalar<PColorVector<RComplexFloat,3> >, OpAdd>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<V>+PSc<V>" << endl;
#endif

  su3_vectorf *lv = (su3_vectorf *)&(l.elem().elem(0).real());
  su3_vectorf *rv = (su3_vectorf *)&(r.elem().elem(0).real());
  su3_vectorf *dv = (su3_vectorf *)&(d.elem().elem(0).real());

  intrin_sse_add_su3_vector(lv,rv,dv);

  return d;
}


#if defined(QDP_SCALARSITE_USE_EVALUATE)
// NOTE: let these be subroutines to save space

//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix = adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix = adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix += LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix += adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix += LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix += adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix -= LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix -= adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix -= LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix -= adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s);



// Specialization to optimize the case
//   LatticeColorMatrix = LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< TCol, OLattice< TCol > > > >,
                 OLattice< TCol > >& rhs, 
	      const Subset& s);


//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix += LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< TCol, OLattice< TCol > > > >,
                 OLattice< TCol > >& rhs, 
	      const Subset& s);

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix -= LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< TCol, OLattice< TCol > > > >,
                 OLattice< TCol > >& rhs, 
	      const Subset& s);

//-------------------------------------------------------------------

// Specialization to optimize the case   
//    LatticeHalfFermion = LatticeColorMatrix * LatticeHalfFermion
// NOTE: let this be a subroutine to save space
template<>
void evaluate(OLattice< TVec2 >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TVec2, OLattice< TVec2 > > > >,
	                    OLattice< TVec2 > >& rhs,
	      const Subset& s);

#endif


/*! @} */   // end of group optimizations

#if defined(QDP_SCALARSITE_DEBUG)
#undef QDP_SCALARSITE_DEBUG
#endif

#if defined(QDP_SCALARSITE_USE_EVALUATE)
#undef QDP_SCALARSITE_USE_EVALUATE
#endif

} // namespace QDP;

#endif  // defined(__GNUC__)

#endif
