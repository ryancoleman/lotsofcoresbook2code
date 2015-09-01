// $Id: qdp_scalarsite_sse_blas_double.h,v 1.8 2009-07-14 20:08:42 bjoo Exp $

/*! @file
 * @brief Generic Scalarsite  optimization hooks
 * 
 * 
 */


#ifndef QDP_SCALARSITE_SSE_BLAS_DOUBLE_H
#define QDP_SCALARSITE_SSE_BLAS_DOUBLE_H

#include "qdp_dispatch.h"

#include "scalarsite_sse/sse_blas_vaxpy4_double.h"
#include "scalarsite_sse/sse_blas_vaypx4_double.h"
#include "scalarsite_sse/sse_blas_vaxmyz4_double.h"
#include "scalarsite_sse/sse_blas_vaxpbyz4_double.h"
#include "scalarsite_sse/sse_blas_vaxmbyz4_double.h"
#include "scalarsite_sse/sse_blas_vscal4_double.h"
#include "scalarsite_sse/sse_blas_local_sumsq_double.h"
#include "scalarsite_sse/sse_blas_local_vcdot_real_double.h"
#include "scalarsite_sse/sse_blas_local_vcdot_double.h"

namespace QDP {

  namespace ThreadReductions{ 
    extern REAL64* norm2_results;
    extern REAL64* innerProd_results;

  }

// Types needed for the expression templates. 
// Bugger Staggered! For Wilson, Ns=4 is nice and cache line
// aligned on a 64byte line (4 vectors=>12 complex=24 doubles.

typedef PSpinVector<PColorVector<RComplex<REAL64>, 3>, 4> DVec;
typedef PScalar<PScalar<RScalar<REAL64> > >  DScal;


////////////////////////////////
// Threading evaluates
//
// by Xu Guo, EPCC, 6 October, 2008
////////////////////////////////

// the wrappers for the functions to be threaded

#include "qdp_dispatch.h"
#include "qdp_scalarsite_sse_blas_double_wrapper.h"


// #define DEBUG_BLAS
// TVec is the LatticeFermion from qdp_dwdefs.h with the OLattice<> stripped
// from around it

// TScalar is the usual Real, with the OScalar<> stripped from it
//
// THis is simply to make the code more readable, and reduces < < s and > >s
// in the template arguments

// d += Scalar*Vec

template<>
inline
void evaluate(OLattice< DVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< DScal, OScalar < DScal > > >,
	      Reference< QDPType< DVec, OLattice< DVec > > > >,
	      OLattice< DVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "y += a*x" << endl;
#endif

  const OLattice< DVec >& x = static_cast<const OLattice< DVec > &>(rhs.expression().right());
  const OScalar< DScal >& a = static_cast<const OScalar< DScal > &> (rhs.expression().left());
  
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;

  if( s.hasOrderedRep() ) { 
    REAL64* xptr = (REAL64 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL64* yptr = &(d.elem(s.start()).elem(0).elem(0).real());

    int total_n_4vec = s.end()-s.start()+1;
   
    ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end()-s.start()+1);
    //vaxpy4(yptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy4_double_user_arg arg(x, d, aptr, 1, tab, vaxpy4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(d.elem(i).elem(0).elem(0).real());
      vaxpy4(yptr, aptr, xptr, 1);
      }*/
  }
}

// d -= Scalar*Vec
template<>
inline
void evaluate(OLattice< DVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< DScal, OScalar < DScal > > >,
	      Reference< QDPType< DVec, OLattice< DVec > > > >,
	      OLattice< DVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "y -= a*x" << endl;
#endif

  const OLattice< DVec >& x = static_cast<const OLattice< DVec > &>(rhs.expression().right());
  const OScalar< DScal >& a = static_cast<const OScalar< DScal > &> (rhs.expression().left());

  // - sign as y -= ax <=> y = y-ax = -ax + y = axpy with -a 
  REAL64 ar = -( a.elem().elem().elem().elem());
  REAL64* aptr = &ar;
  if( s.hasOrderedRep() ) { 
    REAL64* xptr = (REAL64 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL64* yptr = &(d.elem(s.start()).elem(0).elem(0).real());

    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    
    //int n_4vec = (s.end()-s.start()+1);
    //vaxpy4(yptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy4_double_user_arg arg(x, d, aptr, 1, tab, vaxpy4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(d.elem(i).elem(0).elem(0).real());
      vaxpy4(yptr, aptr, xptr, 1);
      }*/
  }
	
}


// z = ax + y
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > >,
	        Reference< QDPType< DVec, OLattice< DVec > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x + y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().right());

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.left());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.right());
  // Set pointers 
  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;

  if( s.hasOrderedRep() ) { 
    
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    //    int n_4vec = (s.end()-s.start()+1);    
    
    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 
      // y = ax + y => AXPY

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,aptr,xptr, n_4vec);
    }
    else { 
      // z = ax + y => AXPYZ
      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxpyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
    }
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());

      vaxpyz4(zptr, aptr, xptr, yptr, 1);
      }*/

  }

}


// Vec = Vec + Scal*Vec
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        Reference< QDPType< DVec, OLattice< DVec > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "z = y + a*x" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.left());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.right());
  // Set pointers 
  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 
      // y = y + ax 
      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,aptr,xptr, n_4vec);
    }
    else {
    

	ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxpyz4};
	
	dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
	
	////////////////
	// Original code
	////////////////
	//vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
    }

  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*

    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
      }*/
  }

}

// Vec = Scalar*Vec - Vec
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > >,
	        Reference< QDPType< DVec, OLattice< DVec > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x - y" << endl;
#endif


  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().right());

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.left());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.right());
  // Set pointers 
  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;
  if( s.hasOrderedRep() ) {
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
    
    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxmy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxmy4(yptr,aptr,xptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxmyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      // vaxmyz4(zptr, aptr, xptr, yptr, n_4vec);
    }

  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxmyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxmyz4(zptr, aptr, xptr, yptr, 1);
      }*/
  }

}


template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        Reference< QDPType< DVec, OLattice< DVec > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "z = y - a*x" << endl;
#endif

  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.left());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.right());
  // Set pointers etc.

  // -ve sign as y - ax = -ax + y  = axpy with -a.
  REAL64 ar =  -a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;
  if( s.hasOrderedRep() ) { 

    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,aptr,xptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxpyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
    }
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
      }*/
  }

}

// Vec += Vec * Scalar (AXPY)
template<>
inline
void evaluate(OLattice< DVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< DVec, OLattice< DVec > > >,
	      Reference< QDPType< DScal, OScalar < DScal > > > >,
	      OLattice< DVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "y += x*a" << endl;
#endif

  const OLattice< DVec >& x = static_cast<const OLattice< DVec > &>(rhs.expression().left());
  const OScalar< DScal >& a = static_cast<const OScalar< DScal > &> (rhs.expression().right());
  
  REAL64 ar = a.elem().elem().elem().elem();
  REAL64* aptr = &ar;
  
  if( s.hasOrderedRep() ) { 
    REAL64* xptr = (REAL64 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL64* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    // cout << "Specialised axpy a ="<< ar << endl;
     
    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    //int n_4vec = (s.end()-s.start()+1);
    //vaxpy4(yptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy4_double_user_arg arg(x, d, aptr, 1, tab, vaxpy4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(d.elem(i).elem(0).elem(0).real());
      
      vaxpy4(yptr, aptr, xptr, 1);
      }*/
  }


}


// Vec -= Vec *Scalar 
template<>
inline
void evaluate(OLattice< DVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< DVec, OLattice< DVec > > >,
	      Reference< QDPType< DScal, OScalar < DScal > > > >,
	      OLattice< DVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "y -= x*a" << endl;
#endif

  const OLattice< DVec >& x = static_cast<const OLattice< DVec > &>(rhs.expression().left());
  const OScalar< DScal >& a = static_cast<const OScalar< DScal > &> (rhs.expression().right());

  // - sign as y -= ax <=> y = y-ax = -ax + y = axpy with -a 
  REAL64 ar = -( a.elem().elem().elem().elem());
  REAL64* aptr = &ar;

  if( s.hasOrderedRep() ) { 
    REAL64* xptr = (REAL64 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL64* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
        
    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end()-s.start()+1);
    //vaxpy4(yptr, aptr, xptr, n_4vec);
  } 
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy4_double_user_arg arg(x, d, aptr, 1, tab, vaxpy4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = &(d.elem(i).elem(0).elem(0).real());
      
      vaxpy4(yptr, aptr, xptr, 1);
      }*/
  }

}


// Vec = Vec *Scalar  + Vec (AXPY)
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > >,
	        Reference< QDPType< DVec, OLattice< DVec > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a + y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().right());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,    
    Reference< QDPType< DScal, OScalar< DScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.right());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.left());
  // Set pointers 
  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;
  if( s.hasOrderedRep() ) { 

    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
    
    if( zptr == yptr ) {

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,aptr,xptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxpyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
    }

  } 
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];

      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = (REAL64 *)  &(d.elem(i).elem(0).elem(0).real());
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
    }*/
  }

}


// Vec = Vec + Vec * Scalar (AXPY)
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        Reference< QDPType< DVec, OLattice< DVec > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "z = y + x*a" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,    
    Reference< QDPType< DScal, OScalar< DScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.right());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.left());
  // Set pointers 
  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;

  if( s.hasOrderedRep() ) {
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,aptr,xptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxpyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
    }

    
  } 
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr = (REAL64 *)&(d.elem(i).elem(0).elem(0).real());
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
      }*/
  }

}


// Vec = Vec*Scalar - Vec (AXMY)
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > >,
	        Reference< QDPType< DVec, OLattice< DVec > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a - y" << endl;
#endif

  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().right());


  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,    
    Reference< QDPType< DScal, OScalar< DScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.right());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.left());
  // Set pointers 
  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxmy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxmy4(yptr,aptr,xptr, n_4vec);
    }
    else {

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxmyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxmyz4(zptr, aptr, xptr, yptr, n_4vec);
    }


  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxmyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxmyz4(zptr, aptr, xptr, yptr, 1);
      }*/
  }

}


// Vec = Vec - Vec*Scalar (AXPY with -Scalar)
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        Reference< QDPType< DVec, OLattice< DVec > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "z = y - x*a" << endl;
#endif

  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&> (rhs.expression().left());

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,    
    Reference< QDPType< DScal, OScalar< DScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode.right());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode.left());
  // Set pointers etc.

  // -ve sign as y - ax = -ax + y  = axpy with -a.
  REAL64 ar =  -a.elem().elem().elem().elem();
  REAL64 *aptr = (REAL64 *)&ar;

  if( s.hasOrderedRep() ) { 

    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, aptr, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,aptr,xptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, aptr, xptr, yptr, vaxpyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpyz4(zptr, aptr, xptr, yptr, n_4vec);
    }

  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, aptr, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxpyz4(zptr, aptr, xptr, yptr, 1);
      }*/
  }

}


template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	       Reference< QDPType< DVec, OLattice< DVec > > >,
	       Reference< QDPType< DVec, OLattice< DVec > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "BJ: v+v " << endl;
#endif

  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(rhs.expression().left());
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(rhs.expression().right());

  REAL64 one = 1;

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, &one, xptr, vaxpy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpy4(yptr,&one,xptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, &one, xptr, yptr, vaxpyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxpyz4(zptr, &one, xptr, yptr, n_4vec);
    }

  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, &one, 1, tab, vaxpyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      
      REAL64* xptr = (REAL64 *)&(x.elem(i).elem(0).elem(0).real());
      REAL64* yptr = (REAL64 *)&(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxpyz4(zptr,&one, xptr, yptr, 1);
      }*/
  }


}

template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	       Reference< QDPType< DVec, OLattice< DVec > > >,
	       Reference< QDPType< DVec, OLattice< DVec > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "BJ: v-v " << endl;
#endif

  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(rhs.expression().left());
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(rhs.expression().right());
  REAL64 one=1;

  if( s.hasOrderedRep() ) { 

    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;

    if( zptr == yptr ) { 

      ordered_sse_vaxOpy4_double_user_arg arg = {yptr, &one, xptr, vaxmy4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxmy4(yptr,&one,xptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpyz4_double_user_arg arg = {zptr, &one, xptr, yptr, vaxmyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpyz4_double_evaluate_function);
      
      ////////////////
      // Original code
      ////////////////
      //vaxmyz4(zptr, &one, xptr, yptr, n_4vec);
    }

  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpyz4_double_user_arg arg(x, y, d, &one, 1, tab, vaxmyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      vaxmyz4(zptr,&one, xptr, yptr, 1);
     
      }*/
  }

}


// Vec = Scal * Vec
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpMultiply,
	       Reference< QDPType< DScal, OScalar< DScal > > >,
	       Reference< QDPType< DVec, OLattice< DVec > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "BJ: v = a*v " << endl;
#endif
  const OLattice< DVec > &x = static_cast<const OLattice< DVec >&>(rhs.expression().right());
  const OScalar< DScal > &a = static_cast<const OScalar< DScal >&>(rhs.expression().left());

  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = &ar;  
  
  if( s.hasOrderedRep() ) {

    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    /*int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vscal4_double_user_arg arg = {zptr,aptr,xptr};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vscal4_double_evaluate_function);
    */
    ////////////////
    // Original code
    ////////////////
    int n_4vec = (s.end()-s.start()+1);
    vscal4(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    /*
    int totalSize = s.numSiteTable();

    unordered_sse_vscal4_double_user_arg arg = {x, d, aptr, 1, tab};

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal4_double_evaluate_function);
    */
    ////////////////
    // Original code
    ////////////////
    
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *zptr =  &(d.elem(i).elem(0).elem(0).real());

      
      vscal4(zptr, aptr, xptr, 1);
      }
  }

}


template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpMultiply,
	       Reference< QDPType< DVec, OLattice< DVec > > >,
	       Reference< QDPType< DScal, OScalar< DScal > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "BJ: v = v*a " << endl;
#endif

  const OLattice< DVec > &x = static_cast<const OLattice< DVec >&>(rhs.expression().left());
  const OScalar< DScal > &a = static_cast<const OScalar< DScal >&>(rhs.expression().right());

  REAL64 ar =  a.elem().elem().elem().elem();
  REAL64 *aptr = &ar;  

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vscal4_double_user_arg arg = {zptr,aptr,xptr};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vscal4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end()-s.start()+1);
    //vscal4(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal4_double_user_arg arg(x, d, aptr, 1, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*  for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      vscal4(zptr, aptr, xptr, 1);
      }*/
  }
}


// v *= a
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpMultiplyAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< DScal, OScalar< DScal > > > >,
	       OScalar< DScal > > &rhs,
	       const Subset& s)
{
  const OScalar< DScal >& a = static_cast< const OScalar<DScal >&>(rhs.expression().child());


#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v *= a, a = " << a << endl;
#endif
  
  REAL64 ar = a.elem().elem().elem().elem();
  if( s.hasOrderedRep() ) { 

    REAL64* xptr = &(d.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr = xptr;
    
    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vscal4_double_user_arg arg = {zptr, &ar, xptr};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vscal4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end()-s.start()+1);
    //vscal4(zptr,&ar, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal4_double_user_arg arg(d, d, &ar, 1, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];

      REAL64* xptr = &(d.elem(i).elem(0).elem(0).real());
      REAL64* zptr = xptr;
      
      vscal4(zptr,&ar, xptr, 1);
      }*/
  }
}

// v /= a
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpDivideAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< DScal, OScalar< DScal > > > >,
	       OScalar< DScal > > &rhs,
	       const Subset& s)
{
  const OScalar< DScal >& a = static_cast< const OScalar<DScal >&>(rhs.expression().child());


#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v /= a, a = " << a << endl;
#endif
  
  REAL64 ar = (REAL64)1/a.elem().elem().elem().elem();
  if( s.hasOrderedRep() ) { 
    REAL64* xptr = &(d.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr = xptr;
    
    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vscal4_double_user_arg arg = {zptr, &ar, xptr};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vscal4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end()-s.start()+1);
    //vscal4(zptr,&ar, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal4_double_user_arg arg(d, d, &ar, 1, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal4_double_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];

      REAL64* xptr = &(d.elem(i).elem(0).elem(0).real());
      REAL64* zptr = xptr;
      
      vscal4(zptr,&ar, xptr, 1);
      }*/
  }
}


// v += v
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAddAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< DVec, OLattice< DVec > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
  const OLattice< DVec >& x = static_cast< const OLattice<DVec >&>(rhs.expression().child());

 

#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v += v" << endl;
#endif
  REAL64 one = 1;

  if( s.hasOrderedRep() ) {

    
    REAL64 *xptr = (REAL64 *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *)(&d.elem(s.start()).elem(0).elem(0).real());

    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vaxOpy4_double_user_arg arg = {yptr, &one, xptr, vaxpy4};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end() - s.start()+1);  
    //vaxpy4(yptr, &one, xptr,n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy4_double_user_arg arg(x, d, &one, 1, tab, vaxpy4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];

      REAL64 *xptr = (REAL64 *)(&x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *)(&d.elem(i).elem(0).elem(0).real());
      
      vaxpy4(yptr, &one, xptr,1);

      }*/
  }

}


// v -= v
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpSubtractAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< DVec, OLattice< DVec > > > >,
	       OLattice< DVec > > &rhs,
	       const Subset& s)
{
  const OLattice< DVec >& x = static_cast< const OLattice<DVec >&>(rhs.expression().child());

 

#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v -= v" << endl;
#endif
  REAL64 mone = (REAL64)-1;
    
  if( s.hasOrderedRep() ) { 

    REAL64 *xptr = (REAL64 *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *)(&d.elem(s.start()).elem(0).elem(0).real());

    int total_n_4vec = s.end()-s.start()+1;

    ordered_sse_vaxOpy4_double_user_arg arg = {yptr, &mone, xptr, vaxpy4};

    dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpy4_double_evaluate_function);
      
    ////////////////
    // Original code
    ////////////////
    //int n_4vec = (s.end() - s.start()+1);    
    //vaxpy4(yptr, &mone, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy4_double_user_arg arg(x, d, &mone, 1, tab, vaxpy4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64 *xptr = (REAL64 *)(&x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *)(&d.elem(i).elem(0).elem(0).real());
    
      vaxpy4(yptr, &mone, xptr, 1);
   
      }*/
  }

}


// z = ax + by
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x + b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.left());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.left());
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr,bptr, vaxpby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 

      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }


    
  }
  else { 
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxpbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxpbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
      }*/
  }

}


// z = xa + by
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a + b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,
    Reference< QDPType< DScal, OScalar< DScal > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN2;


  
  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.left());

  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.left());
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());
  
  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr,bptr, vaxpby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxpbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxpbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }*/
  }
}

// z = ax + yb
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x + y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // type of a*x
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN1;

  // type of y*b
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,
    Reference< QDPType< DScal, OScalar< DScal > > > > BN2;


  
  // get the binary nodes
  // a*x node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());

  // y*b node
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.left());

  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.right());

  
  // get b and y out of the binary node
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.left());

  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.right());

  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr, bptr, vaxpby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxpbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*

    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1

      vaxpbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }
    */
  
  }
}

// z = xa + yb
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a + y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,
    Reference< QDPType< DScal, OScalar< DScal > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.left());
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.left());

  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.right());
  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);
    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr, bptr, vaxpby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxpbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxpbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
      for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxpbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }
   
    */

  }
}


// z = ax - by
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x - b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.left());
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.left());
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) {
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr, bptr, vaxmby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxmbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }
   
    */
  }
}


// z = xa - by
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a - b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,
    Reference< QDPType< DScal, OScalar< DScal > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN2;


  
  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.left());

  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.left());
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());
  
  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr, bptr, vaxmby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxmbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
      for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }
   
    */
  }
}

// z = ax - yb
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DScal, OScalar< DScal > > >,
	         Reference< QDPType< DVec, OLattice< DVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x - y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // type of a*x
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DScal, OScalar< DScal > > >,
    Reference< QDPType< DVec, OLattice< DVec > > > > BN1;

  // type of y*b
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,
    Reference< QDPType< DScal, OScalar< DScal > > > > BN2;


  
  // get the binary nodes
  // a*x node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());

  // y*b node
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.left());

  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.right());

  
  // get b and y out of the binary node
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.left());

  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.right());

  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr, bptr, vaxmby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxmbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*

    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }
   
    */
  
  }
}

// z = xa - yb
template<>
inline
void evaluate( OLattice< DVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< DVec, OLattice< DVec > > >,
	         Reference< QDPType< DScal, OScalar< DScal > > > > >,
	        OLattice< DVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a - y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< DVec, OLattice< DVec > > >,
    Reference< QDPType< DScal, OScalar< DScal > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< DVec >& x = static_cast<const OLattice< DVec >&>(mulNode1.left());
  const OScalar< DScal >& a = static_cast<const OScalar< DScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OLattice< DVec >& y = static_cast<const OLattice< DVec >&>(mulNode2.left());

  const OScalar< DScal >& b = static_cast<const OScalar< DScal >&>(mulNode2.right());
  
  // Set pointers 
  REAL64 *aptr = (REAL64 *)&(a.elem().elem().elem().elem());
  REAL64 *bptr = (REAL64 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 

    REAL64 *xptr = (REAL64 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL64 *yptr = (REAL64 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL64* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_4vec = (s.end()-s.start()+1);

    int total_n_4vec = s.end()-s.start()+1;
   
    if( zptr == yptr ) { 

      ordered_sse_vaxOpby4_double_user_arg arg = {yptr, aptr, xptr, bptr, vaxmby4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpby4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmby4(yptr, aptr, xptr, bptr, n_4vec);
    }
    else { 
      ordered_sse_vaxOpbyz4_double_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmbyz4};

      dispatch_to_threads(total_n_4vec, arg, ordered_sse_vaxOpbyz4_double_evaluate_function);

      ////////////////
      // Original code
      ////////////////
      //vaxmbyz4(zptr, aptr, xptr, bptr, yptr, n_4vec);
    }

    
  }
  else { 
    const int* tab = s.siteTable().slice();    
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpbyz4_double_user_arg arg(x, y, d, aptr, bptr,  1, tab, vaxmbyz4);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpbyz4_double_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
      for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
   
      REAL64 *xptr = (REAL64 *) &(x.elem(i).elem(0).elem(0).real());
      REAL64 *yptr = (REAL64 *) &(y.elem(i).elem(0).elem(0).real());
      REAL64* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmbyz4(zptr, aptr, xptr, bptr, yptr, 1);
   
    }
   
    */
  
  }
}



// Global norm squared of a vector...
template<>
inline UnaryReturn<OLattice< DVec >, FnNorm2>::Type_t
norm2(const QDPType<DVec ,OLattice< DVec > >& s1, const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using BJ sumsq" << endl;
#endif

  if ( s.hasOrderedRep() ) {

#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ sumsq " << endl;
#endif

    ordered_norm_double_user_arg arg;
    int n4vec = s.end()-s.start()+1;

    arg.vptr = (REAL64*)&(s1.elem(s.start()).elem(0).elem(0).real());
    arg.results = ThreadReductions::norm2_results;
    arg.func = local_sumsq4;
    dispatch_to_threads(n4vec, arg, ordered_norm_double_func);
    REAL64 lsum=arg.results[0];
    for(int i=1; i < qdpNumThreads(); i++) { 
	lsum += arg.results[i];
    }
    UnaryReturn< OLattice< DVec >, FnNorm2>::Type_t  gsum(lsum);
    QDPInternal::globalSum(gsum);
    return gsum;

#if 0 
    const REAL64 *s1ptr =  &(s1.elem(s.start()).elem(0).elem(0).real());
    
    // Has Type OScalar< PScalar < PScalar < RScalar < REAL64 > > > >

    REAL64 lsum; // local_sumsq4 zeros this
    int n_4vec = (s.end() - s.start() + 1);    
    local_sumsq4(&lsum,(REAL64 *)s1ptr, n_4vec); 
    UnaryReturn< OLattice< DVec >, FnNorm2>::Type_t  gsum(lsum);
    QDPInternal::globalSum(gsum);
    return gsum;
#endif
  }
  else {

    // Has Type OScalar< PScalar < PScalar < RScalar < REAL64 > > > >
    REAL64 lsum =(REAL64)0;
    REAL64 ltmp =(REAL64)0;

    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64* s1ptr = (REAL64 *)&(s1.elem(i).elem(0).elem(0).real());
      local_sumsq4(&ltmp,s1ptr,1); 
      lsum +=ltmp;
    }

    UnaryReturn< OLattice< DVec >, FnNorm2>::Type_t  gsum(lsum);
    QDPInternal::globalSum(gsum);
    return gsum;
  }
}


template<>
inline UnaryReturn<OLattice< DVec >, FnNorm2>::Type_t
norm2(const QDPType<DVec ,OLattice< DVec > >& s1)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using BJ sumsq all" << endl;
#endif
  ordered_norm_double_user_arg arg;
  int n4vec = all.end()-all.start()+1;

  arg.vptr = (REAL64*)&(s1.elem(all.start()).elem(0).elem(0).real());
  arg.func = local_sumsq4;
  arg.results = ThreadReductions::norm2_results;


  dispatch_to_threads(n4vec, arg, ordered_norm_double_func);

  // Sum partial results
  REAL64 lsum=arg.results[0];
  for(int i=1; i < qdpNumThreads(); i++) {
      lsum += arg.results[i];
  }
  UnaryReturn< OLattice< DVec >, FnNorm2>::Type_t  gsum(lsum);
  QDPInternal::globalSum(gsum);
  return gsum;

#if 0
  int n_4vec = (all.end() - all.start() + 1);
  const REAL64 *s1ptr =  &(s1.elem(all.start()).elem(0).elem(0).real());
    


  REAL64 lsum = (REAL64)0;
  local_sumsq4(&lsum, (REAL64 *)s1ptr, n_4vec); 
  UnaryReturn< OLattice< DVec >, FnNorm2>::Type_t  gsum(lsum);
  QDPInternal::globalSum(gsum);
  return gsum;
#endif
}



template<>
inline  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t
innerProduct(const QDPType< DVec, OLattice<DVec> > &v1,
	     const QDPType< DVec, OLattice<DVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProduct all" << endl;
#endif
    ordered_inner_product_double_user_arg arg;
    int n4vec = all.end()-all.start()+1;
    arg.xptr = (REAL64*)&(v1.elem(all.start()).elem(0).elem(0).real());
    arg.yptr = (REAL64*)&(v2.elem(all.start()).elem(0).elem(0).real());
    arg.func = local_vcdot4;
    arg.results = ThreadReductions::innerProd_results;

    dispatch_to_threads(n4vec, arg, ordered_inner_product_double_func);
    REAL64 ip[2] = { arg.results[0],arg.results[1] };

    for( int i=1; i < qdpNumThreads(); i++) {
      ip[0] += arg.results[2*i];
      ip[1] += arg.results[2*i+1];
    }

    QDPInternal::globalSumArray(ip,2);     

    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RComplex<PScalar<REAL64> > > > >
    BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t lprod;

    lprod.elem().elem().elem().real() = ip[0];
    lprod.elem().elem().elem().imag() = ip[1];
    

    return lprod;
#if 0
  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RComplex<PScalar<REAL64> > > > >
  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t lprod;
  // Inner product is accumulated internally in DOUBLE
  REAL64 ip[2];
  ip[0]=0;
  ip[1]=0;

  // Length of subset 
  unsigned long n_4vec = (all.end() - all.start() + 1);
    
  // Call My CDOT
  local_vcdot4(&(ip[0]),
	      (REAL64 *)&(v1.elem(all.start()).elem(0).elem(0).real()),
	      (REAL64 *)&(v2.elem(all.start()).elem(0).elem(0).real()),
	      n_4vec);


  // Global sum -- still on a vector of doubles
  QDPInternal::globalSumArray(ip,2);

  // Downcast (and possibly lose precision) here 
  lprod.elem().elem().elem().real() = ip[0];
  lprod.elem().elem().elem().imag() = ip[1];

  // Return
  return lprod;
#endif

}

template<>
inline  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t
innerProduct(const QDPType< DVec, OLattice<DVec> > &v1,
	     const QDPType< DVec, OLattice<DVec> > &v2, 
	     const Subset& s)
{
  
  if( s.hasOrderedRep() ) {
#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ: innerProduct s" << endl;
#endif
    ordered_inner_product_double_user_arg arg;
    int n4vec = s.end()-s.start()+1;
    arg.xptr = (REAL64*)&(v1.elem(s.start()).elem(0).elem(0).real());
    arg.yptr = (REAL64*)&(v2.elem(s.start()).elem(0).elem(0).real());
    arg.func = local_vcdot4;
    arg.results = ThreadReductions::innerProd_results;

    dispatch_to_threads(n4vec, arg, ordered_inner_product_double_func);
    REAL64 ip[2] = { arg.results[0],arg.results[1] };

    for( int i=1; i < qdpNumThreads(); i++) {
      ip[0] += arg.results[2*i];
      ip[1] += arg.results[2*i+1];
    }

    QDPInternal::globalSumArray(ip,2);     

    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RComplex<PScalar<REAL64> > > > >
    BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t lprod;

    lprod.elem().elem().elem().real() = ip[0];
    lprod.elem().elem().elem().imag() = ip[1];
    

    return lprod;

#if 0
    ip[0] = 0;
    ip[1] = 0;

    unsigned long n_4vec = (s.end() - s.start() + 1);
    local_vcdot4(&(ip[0]),
		(REAL64 *)&(v1.elem(s.start()).elem(0).elem(0).real()),
		(REAL64 *)&(v2.elem(s.start()).elem(0).elem(0).real()),
		n_4vec);


#endif


  }
  else {

    BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t lprod;
    REAL64 ip[2], ip_tmp[2];
    ip[0] = 0;
    ip[1] = 0;

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 

      int i=tab[j];
      
      local_vcdot4(&(ip_tmp[0]),
		  (REAL64 *)&(v1.elem(i).elem(0).elem(0).real()),
		  (REAL64 *)&(v2.elem(i).elem(0).elem(0).real()),
		  1);
      
      ip[0] += ip_tmp[0];
      ip[1] += ip_tmp[1];
    }

    QDPInternal::globalSumArray(ip,2);

    lprod.elem().elem().elem().real() = ip[0];
    lprod.elem().elem().elem().imag() = ip[1];
    

    return lprod;

  }
}


// Inner Product Real
template<>
inline  
BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t
innerProductReal(const QDPType< DVec, OLattice<DVec> > &v1,
		 const QDPType< DVec, OLattice<DVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL64> > > > >
  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t lprod;
  // Inner product is accumulated internally in DOUBLE
  REAL64 ip_re=0;

  // Length of subset 
  unsigned long n_4vec = (all.end() - all.start() + 1);

  // Call My CDOT
  local_vcdot_real4(&ip_re,
		   (REAL64 *)&(v1.elem(all.start()).elem(0).elem(0).real()),
		   (REAL64 *)&(v2.elem(all.start()).elem(0).elem(0).real()),
		   n_4vec);

  // Global sum
  QDPInternal::globalSum(ip_re);

  // Whether CDOT did anything or not ip_re and ip_im should 
  // now be right. Assign them to the ReturnType
  lprod.elem().elem().elem().elem() = ip_re;


  // Return
  return lprod;
}


template<>
inline  
BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t
innerProductReal(const QDPType< DVec, OLattice<DVec> > &v1,
		 const QDPType< DVec, OLattice<DVec> > &v2, 
		 const Subset& s)
{
  if( s.hasOrderedRep() ) {
#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ: innerProductReal s" << endl;
#endif

    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RScalar<PScalar<REAL64> > > > >
    BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t lprod;
    REAL64 ip_re=0;

    unsigned long n_4vec = (s.end() - s.start() + 1);
    local_vcdot_real4(&ip_re,
		     (REAL64 *)&(v1.elem(s.start()).elem(0).elem(0).real()),
		     (REAL64 *)&(v2.elem(s.start()).elem(0).elem(0).real()),
		     n_4vec);

    QDPInternal::globalSum(ip_re);
    lprod.elem().elem().elem().elem() = ip_re;


    return lprod;
  }
  else {


    BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t lprod;
    REAL64 ip_re=0, ip_re_tmp;


    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 

      int i=tab[j];
      
      local_vcdot_real4(&ip_re_tmp,
		  (REAL64 *)&(v1.elem(i).elem(0).elem(0).real()),
		  (REAL64 *)&(v2.elem(i).elem(0).elem(0).real()),
		  1);
      
      ip_re += ip_re_tmp;
    }
    QDPInternal::globalSum(ip_re);
    lprod.elem().elem().elem().elem() = ip_re;
    return lprod;
  }
}


template<>
inline UnaryReturn<OLattice< DVec >, FnNorm2>::Type_t
norm2(const multi1d< OLattice< DVec > >& s1)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using SSE multi1d sumsq all" << endl;
#endif

  int n_4vec = (all.end() - all.start() + 1);
  REAL64 ltmp = 0;
  for(int n=0; n < s1.size(); ++n)
  {
    const REAL64* s1ptr =  &(s1[n].elem(all.start()).elem(0).elem(0).real());
    
    // I am relying on this being a Double here 
    REAL64 lltmp=0;
    local_sumsq4(&lltmp, (REAL64*)s1ptr, n_4vec); 

    ltmp += lltmp;
  }

  UnaryReturn< OLattice< DVec >, FnNorm2>::Type_t  lsum(ltmp);
  QDPInternal::globalSum(lsum);
  return lsum;
}


template<>
inline  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OLattice<DVec> > &v1,
	     const multi1d< OLattice<DVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: multi1d innerProduct all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RComplex<PScalar<REAL64> > > > >
  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProduct>::Type_t lprod;

  // Inner product is accumulated internally in DOUBLE
  REAL64 ip[2];
  ip[0]=0;
  ip[1]=0;

  // Length of subset 
  unsigned long n_4vec = (all.end() - all.start() + 1);
    
  for(int n=0; n < v1.size(); ++n)
  {
    REAL64 iip[2];
    iip[0]=0;
    iip[1]=0;

    // Call My CDOT
    local_vcdot4(&(iip[0]), 
		(REAL64 *)&(v1[n].elem(all.start()).elem(0).elem(0).real()),
		(REAL64 *)&(v2[n].elem(all.start()).elem(0).elem(0).real()),
		n_4vec);
    
    ip[0] += iip[0];
    ip[1] += iip[1];
  }

  // Global sum -- still on a vector of doubles
  QDPInternal::globalSumArray(ip,2);

  // Downcast (and possibly lose precision) here 
  lprod.elem().elem().elem().real() = ip[0];
  lprod.elem().elem().elem().imag() = ip[1];

  // Return
  return lprod;
}



// Inner Product Real
template<>
inline  
BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<DVec> > &v1,
		 const multi1d< OLattice<DVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal(multi1d) all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL64> > > > >
  BinaryReturn< OLattice<DVec>, OLattice<DVec>, FnInnerProductReal>::Type_t lprod;
  // Inner product is accumulated internally in DOUBLE
  REAL64 ip_re=0;

  // Length of subset 
  unsigned long n_4vec = (all.end() - all.start() + 1);

  for(int n=0; n < v1.size(); ++n)
  {
    REAL64 iip_re=0;

    // Call My CDOT
    local_vcdot_real4(&iip_re,
		     (REAL64 *)&(v1[n].elem(all.start()).elem(0).elem(0).real()),
		     (REAL64 *)&(v2[n].elem(all.start()).elem(0).elem(0).real()),
		     n_4vec);

    ip_re += iip_re;
  }

  // Global sum
  QDPInternal::globalSum(ip_re);

  // Whether CDOT did anything or not ip_re and ip_im should 
  // now be right. Assign them to the ReturnType
  lprod.elem().elem().elem().elem() = ip_re;


  // Return
  return lprod;
}


} // namespace QDP;


#endif  // guard
 
