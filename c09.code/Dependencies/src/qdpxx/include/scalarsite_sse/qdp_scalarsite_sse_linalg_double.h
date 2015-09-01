// -*- C++ -*-
// $Id: qdp_scalarsite_sse_linalg_double.h,v 1.2 2009-02-10 17:06:37 bjoo Exp $

/*! @file
 * @brief Blas optimizations
 *
 * Blas optimizations of basic operations
 */

#ifndef QDP_SCALARSITE_SSE_LINALG_DOUBLE_H
#define QDP_SCALARSITE_SSE_LINALG_DOUBLE_H

#include "sse_linalg_mm_su3_double.h"


namespace QDP {

// #define QDP_SCALARSITE_DEBUG

#define QDP_SCALARSITE_USE_EVALUATE


/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */


// Types needed for the expression templates. 
  typedef PScalar<PColorMatrix<RComplex<REAL64>, 3> > DCol;

  // Vector and half vector not yet ready.
#if 0
  typedef PSpinVector<PColorVector<RComplex<REAL64>, 3>, 2> DVec2;
  typedef PSpinVector<PColorVector<RComplex<REAL64>, 3>, 4> DVec4;
#endif



// Optimized version of  
//    PColorMatrix<RComplex<REAL64>,3> <- PColorMatrix<RComplex<REAL64>,3> * PColorMatrix<RComplex<REAL64>,3>
template<>
inline BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
  PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpMultiply>::Type_t
operator*(const PMatrix<RComplex<REAL64>,3,PColorMatrix>& l, 
	  const PMatrix<RComplex<REAL64>,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
    PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*M" << endl;
#endif
  
  REAL64* lm = (REAL64 *) &(l.elem(0,0).real());
  REAL64* rm = (REAL64 *) &(r.elem(0,0).real());
  REAL64* dm = (REAL64 *) &(d.elem(0,0).real());

  // Use Unaligned here
  ssed_m_eq_mm_u(dm,lm,rm,1);

  return d;
}



// Optimized version of  
//    PScalar<PColorMatrix<RComplex<REAL64>,3>> <- PScalar<PColorMatrix<RComplex<REAL64>,3>> * 
//                         PScalar<PColorMatrix<RComplex<REAL64>,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >,
  PScalar<PColorMatrix<RComplex<REAL64>,3> >, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	  const PScalar<PColorMatrix<RComplex<REAL64>,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PScalar<PColorMatrix<RComplex<REAL64>,3> >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*PSc<M>" << endl;
#endif
  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  REAL64* rm = (REAL64 *)&(r.elem().elem(0,0).real());
  REAL64* dm = (REAL64 *)&(d.elem().elem(0,0).real());
  
  // Use unaligned here
  ssed_m_eq_mm_u(dm, lm, rm, 1);

  return d;
}


// Optimized version of  
//    PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>> <- PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>> * 
//                         PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>>
template<>
inline BinaryReturn<PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >,
  PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, OpMultiply>::Type_t
operator*(const PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >& l, 
	  const PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >& r)
{
  BinaryReturn<PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, 
    PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<PSc<M>>*PSc<PSc<M>>" << endl;
#endif

  REAL64* lm = (REAL64 *)&(l.elem().elem().elem(0,0).real());
  REAL64* rm = (REAL64 *)&(r.elem().elem().elem(0,0).real());
  REAL64* dm = (REAL64 *)&(d.elem().elem().elem(0,0).real());

  // Use unaligned here
  ssed_m_eq_mm_u(dm,lm, rm,1);

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplex<REAL64>,3> <- adj(PColorMatrix<RComplex<REAL64>,3>) * PColorMatrix<RComplex<REAL64>,3>
template<>
inline BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
  PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<RComplex<REAL64>,3,PColorMatrix>& l, 
	    const PMatrix<RComplex<REAL64>,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
    PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*M" << endl;
#endif
  REAL64* lm = (REAL64 *) &(l.elem(0,0).real());
  REAL64* rm = (REAL64 *) &(r.elem(0,0).real());
  REAL64* dm = (REAL64 *) &(d.elem(0,0).real());

  ssed_m_eq_hm_u(dm,lm,rm,1);

  return d;
}


// Optimized version of  
//   PScalar<PColorMatrix<RComplex<REAL64>,3>> <- adj(PScalar<PColorMatrix<RComplex<REAL64>,3>>) * PScalar<PColorMatrix<RComplex<REAL64>,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PScalar<PColorMatrix<RComplex<REAL64>,3> >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	    const PScalar<PColorMatrix<RComplex<REAL64>,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PScalar<PColorMatrix<RComplex<REAL64>,3> >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*PSc<M>" << endl;
#endif

  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  REAL64* rm = (REAL64 *)&(r.elem().elem(0,0).real());
  REAL64* dm = (REAL64 *)&(d.elem().elem(0,0).real());

  ssed_m_eq_hm_u(dm, lm, rm, 1);


  return d;
}


// Optimized version of  
//   PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>> <- adj(PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>>) * 
//        PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>>
template<>
inline BinaryReturn<PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, 
  PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >& l, 
	    const PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >& r)
{
  BinaryReturn<PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, 
    PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<PSc<M>>)*PSc<PSc<M>>" << endl;
#endif

  REAL64* lm = (REAL64 *)&(l.elem().elem().elem(0,0).real());
  REAL64* rm = (REAL64 *)&(r.elem().elem().elem(0,0).real());
  REAL64* dm = (REAL64 *)&(d.elem().elem().elem(0,0).real());

  ssed_m_eq_hm_u(dm,lm, rm,1);

  return d;
}


// Optimized version of  
//   PColorMatrix<RComplex<REAL64>,3> <- PColorMatrix<RComplex<REAL64>,3> * adj(PColorMatrix<RComplex<REAL64>,3>)
template<>
inline BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
  PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpMultiplyAdj>::Type_t
multiplyAdj(const PMatrix<RComplex<REAL64>,3,PColorMatrix>& l, 
	    const PMatrix<RComplex<REAL64>,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
    PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*adj(M)" << endl;
#endif
  REAL64* lm = (REAL64 *) &(l.elem(0,0).real());
  REAL64* rm = (REAL64 *) &(r.elem(0,0).real());
  REAL64* dm = (REAL64 *) &(d.elem(0,0).real());

  ssed_m_eq_mh_u(dm,lm,rm,1);


  return d;
}


// Optimized version of  
//   PScalar<PColorMatrix<RComplex<REAL64>,3>> <- PScalar<PColorMatrix<RComplex<REAL64>,3>> * 
//          adj(PScalar<PColorMatrix<RComplex<REAL64>,3>>)
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PScalar<PColorMatrix<RComplex<REAL64>,3> >, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	    const PScalar<PColorMatrix<RComplex<REAL64>,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PScalar<PColorMatrix<RComplex<REAL64>,3> >, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*adj(PSc<M>)" << endl;
#endif
  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  REAL64* rm = (REAL64 *)&(r.elem().elem(0,0).real());
  REAL64* dm = (REAL64 *)&(d.elem().elem(0,0).real());

  ssed_m_eq_mh_u(dm,lm,rm,1);
  return d;
}


// Optimized version of  
//   PScalar<Pscalar<PColorMatrix<RComplex<REAL64>,3>>> <- PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>> * 
//           adj(PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3>>>)
template<>
inline BinaryReturn<PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, 
  PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >& l, 
	    const PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >& r)
{
  BinaryReturn<PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, 
    PScalar<PScalar<PColorMatrix<RComplex<REAL64>,3> > >, OpMultiplyAdj>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<PSc<M>>*adj(PSc<PSc<M>>)" << endl;
#endif
  REAL64* lm = (REAL64 *)&(l.elem().elem().elem(0,0).real());
  REAL64* rm = (REAL64 *)&(r.elem().elem().elem(0,0).real());
  REAL64* dm = (REAL64 *)&(d.elem().elem().elem(0,0).real());

  ssed_m_eq_mh_u(dm,lm,rm,1);
  return d;
}


// Ooops, this macro does not exist!!

// Optimized version of  
//   PColorMatrix<RComplex<REAL64>,3> <- adj(PColorMatrix<RComplex<REAL64>,3>) * adj(PColorMatrix<RComplex<REAL64>,3>)
template<>
inline BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
  PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const PMatrix<RComplex<REAL64>,3,PColorMatrix>& l, 
	       const PMatrix<RComplex<REAL64>,3,PColorMatrix>& r)
{
  BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
    PMatrix<RComplex<REAL64>,3,PColorMatrix>, OpAdjMultiplyAdj>::Type_t  d;

  PColorMatrix<RComplex<REAL64>,3> tmp;

//

  REAL64* lm = (REAL64 *) &(l.elem(0,0).real());
  REAL64* rm = (REAL64 *) &(r.elem(0,0).real());
  REAL64* dm = (REAL64 *) &(d.elem(0,0).real());
  ssed_m_eq_hh_u(dm,lm,rm,1);

  return d;
}


#if 0
// Optimized version of  
//    PColorVector<RComplex<REAL64>,3> <- PColorMatrix<RComplex<REAL64>,3> * PColorVector<RComplex<REAL64>,3>
template<>
inline BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
  PVector<RComplex<REAL64>,3,PColorVector>, OpMultiply>::Type_t
operator*(const PMatrix<RComplex<REAL64>,3,PColorMatrix>& l, 
	  const PVector<RComplex<REAL64>,3,PColorVector>& r)
{
  BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
    PVector<RComplex<REAL64>,3,PColorVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*V" << endl;
#endif

  REAL64* lm = (REAL64 *)&(l.elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv,dv);

  return d;
}



// Optimized version of  
//    PScalar<PColorVector<RComplex<REAL64>,3>> <- PScalar<PColorMatrix<RComplex<REAL64>,3>> * PScalar<PColorVector<RComplex<REAL64>,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PScalar<PColorVector<RComplex<REAL64>,3> >, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	  const PScalar<PColorVector<RComplex<REAL64>,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PScalar<PColorVector<RComplex<REAL64>,3> >, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*PSc<V>" << endl;
#endif

  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem().elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem().elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of  
//    PColorVector<RComplex<REAL64>,3> <- adj(PColorMatrix<RComplex<REAL64>,3>) * PColorVector<RComplex<REAL64>,3>
template<>
inline BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
  PVector<RComplex<REAL64>,3,PColorVector>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<RComplex<REAL64>,3,PColorMatrix>& l, 
	    const PVector<RComplex<REAL64>,3,PColorVector>& r)
{
  BinaryReturn<PMatrix<RComplex<REAL64>,3,PColorMatrix>, 
    PVector<RComplex<REAL64>,3,PColorVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(M)*V" << endl;
#endif

  REAL64* lm = (REAL64 *)&(l.elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of   StaggeredFermion <- ColorMatrix*StaggeredFermion
//    PSpinVector<PColorVector<RComplex<REAL64>,3>,1> <- PScalar<PColorMatrix<RComplex<REAL64>,3>> * PSpinVector<PColorVector<RComplex<REAL64>,3>,1>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PSpinVector<PColorVector<RComplex<REAL64>,3>,1>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	  const PSpinVector<PColorVector<RComplex<REAL64>,3>,1>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PSpinVector<PColorVector<RComplex<REAL64>,3>,1>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "M*S" << endl;
#endif
  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of  
//    PScalar<PColorVector<RComplex<REAL64>,3>> <- adj(PScalar<PColorMatrix<RComplex<REAL64>,3>>) * PScalar<PColorVector<RComplex<REAL64>,3>>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PScalar<PColorVector<RComplex<REAL64>,3> >, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	    const PScalar<PColorVector<RComplex<REAL64>,3> >& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PScalar<PColorVector<RComplex<REAL64>,3> >, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*PSc<V>" << endl;
#endif
  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem().elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem().elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of   StaggeredFermion <- adj(ColorMatrix)*StaggeredFermion
//    PSpinVector<PColorVector<RComplex<REAL64>,3>,1> <- adj(PScalar<PColorMatrix<RComplex<REAL64>,3>>) * PSpinVector<PColorVector<RComplex<REAL64>,3>,1>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PSpinVector<PColorVector<RComplex<REAL64>,3>,1>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	    const PSpinVector<PColorVector<RComplex<REAL64>,3>,1>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PSpinVector<PColorVector<RComplex<REAL64>,3>,1>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*S" << endl;
#endif
  REAL64* lm = (REAL64 *)&(l.elem().elem(0,0).real());
  su3_vectorf* rv = (su3_vectorf *)&(r.elem(0).elem(0).real());
  su3_vectorf* dv = (su3_vectorf *)&(d.elem(0).elem(0).real());

  intrin_sse_mult_adj_su3_mat_vec(lm,rv,dv);

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplex<REAL64>,3>,2> <- PScalar<PColorMatrix<RComplex<REAL64>,3>> * 
//                     PSpinVector<ColorVector<RComplex<REAL64>,3>,2>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PSpinVector<PColorVector<RComplex<REAL64>,3>,2>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
          const PSpinVector<PColorVector<RComplex<REAL64>,3>,2>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PSpinVector<PColorVector<RComplex<REAL64>,3>,2>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*H" << endl;
#endif
  REAL64* lm =(REAL64 *)&( l.elem().elem(0,0).real() );
  half_wilson_vectorf* rh = (half_wilson_vectorf *)&(r.elem(0).elem(0).real()  );
  half_wilson_vectorf* dh = (half_wilson_vectorf *)&(d.elem(0).elem(0).real() );

  intrin_sse_mult_su3_mat_hwvec(lm,rh,dh);

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplex<REAL64>,3>,2> <- adj(PScalar<PColorMatrix<RComplex<REAL64>,3>>) * 
//                     PSpinVector<ColorVector<RComplex<REAL64>,3>,2>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PSpinVector<PColorVector<RComplex<REAL64>,3>,2>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
            const PSpinVector<PColorVector<RComplex<REAL64>,3>,2>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PSpinVector<PColorVector<RComplex<REAL64>,3>,2>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*H" << endl;
#endif

  REAL64* lm =(REAL64 *)&( l.elem().elem(0,0).real() );
  half_wilson_vectorf* rh = (half_wilson_vectorf *)&(r.elem(0).elem(0).real()  );
  half_wilson_vectorf* dh = (half_wilson_vectorf *)&(d.elem(0).elem(0).real() );

  intrin_sse_mult_adj_su3_mat_hwvec(lm,rh,dh);

  return d;
}


// Optimized version of  
//    PColorVector<RComplex<REAL64>,3> <- PColorVector<RComplex<REAL64>,3> + PColorVector<RComplex<REAL64>,3>
template<>
inline BinaryReturn<PVector<RComplex<REAL64>,3,PColorVector>, 
  PVector<RComplex<REAL64>,3,PColorVector>, OpAdd>::Type_t
operator+(const PVector<RComplex<REAL64>,3,PColorVector>& l, 
	  const PVector<RComplex<REAL64>,3,PColorVector>& r)
{
  BinaryReturn<PVector<RComplex<REAL64>,3,PColorVector>, 
    PVector<RComplex<REAL64>,3,PColorVector>, OpAdd>::Type_t  d;

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
//    PSpinVector<PColorVector<RComplex<REAL64>,3>,4> <- PScalar<PColorMatrix<RComplex<REAL64>,3>> 
//                           * PSpinVector<PColorVector<RComplex<REAL64>,3>,4>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PSpinVector<PColorVector<RComplex<REAL64>,3>,4>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	  const PSpinVector<PColorVector<RComplex<REAL64>,3>,4>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PSpinVector<PColorVector<RComplex<REAL64>,3>,4>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*D" << endl;
#endif
  REAL64 *lm = (REAL64 *)&(l.elem().elem(0,0).real());
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
//    PSpinVector<PColorVector<RComplex<REAL64>,3>,4> <- adj(PScalar<PColorMatrix<RComplex<REAL64>,3>>)
//                           * PSpinVector<PColorVector<RComplex<REAL64>,3>,4>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
  PSpinVector<PColorVector<RComplex<REAL64>,3>,4>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplex<REAL64>,3> >& l, 
	    const PSpinVector<PColorVector<RComplex<REAL64>,3>,4>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplex<REAL64>,3> >, 
    PSpinVector<PColorVector<RComplex<REAL64>,3>,4>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*D" << endl;
#endif
  REAL64 *lm = (REAL64 *)&(l.elem().elem(0,0).real());
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
//    PScalar<PColorVector<RComplex<REAL64>,3>> <- PScalar<PColorVector<RComplex<REAL64>,3>> 
//                                            + PScalar<PColorVector<RComplex<REAL64>,3>>
template<>
inline BinaryReturn<PScalar<PColorVector<RComplex<REAL64>,3> >, 
  PScalar<PColorVector<RComplex<REAL64>,3> >, OpAdd>::Type_t
operator+(const PScalar<PColorVector<RComplex<REAL64>,3> >& l, 
	  const PScalar<PColorVector<RComplex<REAL64>,3> >& r)
{
  BinaryReturn<PScalar<PColorVector<RComplex<REAL64>,3> >, 
    PScalar<PColorVector<RComplex<REAL64>,3> >, OpAdd>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<V>+PSc<V>" << endl;
#endif

  su3_vectorf *lv = (su3_vectorf *)&(l.elem().elem(0).real());
  su3_vectorf *rv = (su3_vectorf *)&(r.elem().elem(0).real());
  su3_vectorf *dv = (su3_vectorf *)&(d.elem().elem(0).real());

  intrin_sse_add_su3_vector(lv,rv,dv);

  return d;
}

#endif


#if defined(QDP_SCALARSITE_USE_EVALUATE)
// NOTE: let these be subroutines to save space

//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix = adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix = adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix += LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix += adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix += LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix += adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix -= LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix -= adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >, 
	                    Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix -= LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);

// Specialization to optimize the case   
//    LatticeColorMatrix -= adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< DCol, OLattice< DCol > > > > >,
	                    OLattice< DCol > >& rhs,
	      const Subset& s);



// Specialization to optimize the case
//   LatticeColorMatrix = LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< DCol, OLattice< DCol > > > >,
                 OLattice< DCol > >& rhs, 
	      const Subset& s);


//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix += LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< DCol, OLattice< DCol > > > >,
                 OLattice< DCol > >& rhs, 
	      const Subset& s);

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix -= LatticeColorMatrix
template<>
void evaluate(OLattice< DCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< DCol, OLattice< DCol > > > >,
                 OLattice< DCol > >& rhs, 
	      const Subset& s);


#if 0
//-------------------------------------------------------------------

// Specialization to optimize the case   
//    LatticeHalfFermion = LatticeColorMatrix * LatticeHalfFermion
// NOTE: let this be a subroutine to save space
template<>
void evaluate(OLattice< DVec2 >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< DCol, OLattice< DCol > > >, 
	                    Reference<QDPType< DVec2, OLattice< DVec2 > > > >,
	                    OLattice< DVec2 > >& rhs,
	      const Subset& s);

#endif


/*! @} */   // end of group optimizations

#if defined(QDP_SCALARSITE_DEBUG)
#undef QDP_SCALARSITE_DEBUG
#endif

#if defined(QDP_SCALARSITE_USE_EVALUATE)
#undef QDP_SCALARSITE_USE_EVALUATE
#endif

#endif

} // namespace QDP;


#endif
