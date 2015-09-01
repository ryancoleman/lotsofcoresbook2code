// -*- C++ -*-
// $Id: qdp_scalarsite_generic_linalg.h,v 1.10 2009-07-14 20:08:41 bjoo Exp $

/*! @file
 * @brief Generic optimizations
 *
 * Generic version of optimized basic operations
 */

#ifndef QDP_SCALARSITE_GENERIC_LINALG_H
#define QDP_SCALARSITE_GENERIC_LINALG_H

namespace QDP {

/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */

// Use this def just to safe some typing later on in the file
typedef RComplex<REAL>  RComplexFloat;


#include "scalarsite_generic/generic_mult_nn.h"
#include "scalarsite_generic/generic_mult_na.h"
#include "scalarsite_generic/generic_mult_an.h"
#include "scalarsite_generic/generic_mult_aa.h"
#include "scalarsite_generic/generic_mat_vec.h"
// #include "scalarsite_generic/generic_adj_mat_vec.h" -- No longer used."
#include "scalarsite_generic/generic_addvec.h"


// #define QDP_SCALARSITE_DEBUG

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

  _inline_generic_mult_su3_nn(l,r,d);

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

  _inline_generic_mult_su3_nn(l.elem(),r.elem(),d.elem());

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

  _inline_generic_mult_su3_an(l,r,d);

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

  _inline_generic_mult_su3_an(l.elem(),r.elem(),d.elem());

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

  _inline_generic_mult_su3_na(l,r,d);

  return d;
}


// Optimized version of  
//   PScalar<PColorMatrix<RComplexFloat,3>> <- PScalar<PColorMatrix<RComplexFloat,3>> * adj(PScalar<PColorMatrix<RComplexFloat,3>>)
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

  _inline_generic_mult_su3_na(l.elem(),r.elem(),d.elem());

  return d;
}


#if 0
// Do not use this routine - for some unknown reason on P4 it is slower than
// the default code!!

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

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*adj(PSc<M>)" << endl;
#endif

  _inline_generic_mult_su3_aa(l,r,d);

  return d;
}
#endif


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

  _inline_generic_mult_su3_mat_vec(l,r,d);

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

  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(),d.elem());

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

  _inline_generic_mult_adj_su3_mat_vec(l,r,d);
  
  return d;
}


// Optimized version of   StaggeredFermion <- ColorMatrix*StaggeredFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,1> <- PScalar<PColorMatrix<RComplexFloat,3>> * PSpinVector<PColorVector<RComplexFloat,3>,1>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PVector<PColorVector<RComplexFloat,3>,1,PSpinVector>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	  const PVector<PColorVector<RComplexFloat,3>,1,PSpinVector>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PVector<PColorVector<RComplexFloat,3>,1,PSpinVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*S" << endl;
#endif

  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));

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

  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(),d.elem());

  return d;
}


// Optimized version of   StaggeredFermion <- adj(ColorMatrix)*StaggeredFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,1> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * PSpinVector<PColorVector<RComplexFloat,3>,1>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PVector<PColorVector<RComplexFloat,3>,1,PSpinVector>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PVector<PColorVector<RComplexFloat,3>,1,PSpinVector>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PVector<PColorVector<RComplexFloat,3>,1,PSpinVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*S" << endl;
#endif

  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,2> <- PScalar<PColorMatrix<RComplexFloat,3>> * 
//                     PSpinVector<ColorVector<RComplexFloat,3>,2>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PVector<PColorVector<RComplexFloat,3>,2,PSpinVector>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
          const PVector<PColorVector<RComplexFloat,3>,2,PSpinVector>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PVector<PColorVector<RComplexFloat,3>,2,PSpinVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*H" << endl;
#endif

  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));

  return d;
}


// Optimized version of    HalfFermion <- ColorMatrix*HalfFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,2> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>) * 
//                     PSpinVector<ColorVector<RComplexFloat,3>,2>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PVector<PColorVector<RComplexFloat,3>,2,PSpinVector>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
            const PVector<PColorVector<RComplexFloat,3>,2,PSpinVector>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PVector<PColorVector<RComplexFloat,3>,2,PSpinVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*H" << endl;
#endif

  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));

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

  _inline_generic_add_su3_vector(l,r,d);

  return d;
}


// Optimized version of   DiracFermion <- ColorMatrix*DiracFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,4> <- PScalar<PColorMatrix<RComplexFloat,3>> 
//                           * PSpinVector<PColorVector<RComplexFloat,3>,4>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PVector<PColorVector<RComplexFloat,3>,4,PSpinVector>, OpMultiply>::Type_t
operator*(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	  const PVector<PColorVector<RComplexFloat,3>,4,PSpinVector>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PVector<PColorVector<RComplexFloat,3>,4,PSpinVector>, OpMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "PSc<M>*D" << endl;
#endif

  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));
  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(2),d.elem(2));
  _inline_generic_mult_su3_mat_vec(l.elem(),r.elem(3),d.elem(3));

  return d;
}


// Optimized version of   DiracFermion <- adj(ColorMatrix)*DiracFermion
//    PSpinVector<PColorVector<RComplexFloat,3>,4> <- adj(PScalar<PColorMatrix<RComplexFloat,3>>)
//                           * PSpinVector<PColorVector<RComplexFloat,3>,4>
template<>
inline BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
  PVector<PColorVector<RComplexFloat,3>,4,PSpinVector>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<PColorMatrix<RComplexFloat,3> >& l, 
	    const PVector<PColorVector<RComplexFloat,3>,4,PSpinVector>& r)
{
  BinaryReturn<PScalar<PColorMatrix<RComplexFloat,3> >, 
    PVector<PColorVector<RComplexFloat,3>,4,PSpinVector>, OpAdjMultiply>::Type_t  d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(PSc<M>)*D" << endl;
#endif

  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(0),d.elem(0));
  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(1),d.elem(1));
  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(2),d.elem(2));
  _inline_generic_mult_adj_su3_mat_vec(l.elem(),r.elem(3),d.elem(3));

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

  _inline_generic_add_su3_vector(l.elem(),r.elem(),d.elem());

  return d;
}


#if 1

////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 07 August, 2008
////////////////////////////////

// the wrappers for the function to be threaded
#include "qdp_scalarsite_generic_linalg_wrapper.h"

// Specialization to optimize the case   
//    LatticeHalfFermion = LatticeColorMatrix * LatticeHalfFermion
// NOTE: let this be a subroutine to save space
template<>
inline 
void evaluate(OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PSpinVector<PColorVector<RComplexFloat, 3>, 2>, 
	      OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > > > >,
	      OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > >& rhs,
	      const Subset& s)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "specialized QDP_H_M_times_H" << endl;
#endif

  ////////////////////
  // Original code
  ////////////////////
  //typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  //typedef OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > H;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const H& r = static_cast<const H&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    
    int totalSize = s.end() - s.start() + 1;
    
    int base = s.start();
    
    ordered_linalg_user_arg a(d, l, r, base);
    
    dispatch_to_threads(totalSize, a, ordered_linalg_evaluate_userfunc);

    ////////////////////
    // Original code
    ////////////////////
    // Ordered Way - loop through sites and save a table lookup
    //for(int i=s.start(); i <= s.end(); i++) { 
      
    //_inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
    //			       r.elem(i).elem(0),
    //			       d.elem(i).elem(0));
    //_inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
    //			       r.elem(i).elem(1),
    //			       d.elem(i).elem(1));
    //}

  }
  else { 

    int totalSize = s.numSiteTable();

    const int *tab = s.siteTable().slice(); 

    unordered_linalg_user_arg arg(d, l, r, tab);

    dispatch_to_threads(totalSize, arg, unordered_linalg_evaluate_userfunc);
    
    ////////////////////
    // Original code
    ////////////////////
    // Unordered Way - do a site table lookup
    //for(int j=0; j < s.numSiteTable(); j++) { 
    //int i = tab[j];
      
    //_inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
    //			       r.elem(i).elem(0),
    //			       d.elem(i).elem(0));
    //_inline_generic_mult_su3_mat_vec(l.elem(i).elem(),
    //			       r.elem(i).elem(1),
    //			       d.elem(i).elem(1));
    //}
  }
}

#endif

/*! @} */   // end of group optimizations

#if defined(QDP_SCALARSITE_DEBUG)
#undef QDP_SCALARSITE_DEBUG
#endif

} // namespace QDP;

#endif
