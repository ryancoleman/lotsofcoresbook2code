// $Id: qdp_scalarsite_sse_blas.h,v 1.30 2009-09-15 20:48:51 bjoo Exp $
/*! @file
 * @brief Blas optimizations
 * 
 * Generic and maybe SSE optimizations of basic operations
 */

#ifndef QDP_SCALARSITE_SSE_BLAS_H
#define QDP_SCALARSITE_SSE_BLAS_H

#include "qdp_config.h"

namespace QDP {

  namespace ThreadReductions{ 
    extern REAL64* norm2_results;
  };

  // Forward declarations of BLAS routines
  void vaxpy3(REAL32 *Out, REAL32 *scalep,REAL32 *InScale, REAL32 *Add,int n_3vec);
  void vaxmy3(REAL32 *Out, REAL32 *scalep,REAL32 *InScale, REAL32 *Sub,int n_3vec);
  void vadd(REAL32 *Out, REAL32 *In1, REAL32 *In2, int n_3vec);
  void vsub(REAL32 *Out, REAL32 *In1, REAL32 *In2, int n_3vec);
  void vscal(REAL32 *Out, REAL32 *scalep, REAL32 *In, int n_3vec);
  void vaxpby3(REAL32* Out, REAL32* a, REAL32* x, REAL32* b, REAL32* y, int n_3vec);
  void vaxmby3(REAL32* Out, REAL32* a, REAL32* x, REAL32* b, REAL32* y, int n_3vec);
  
  void local_sumsq_24_48(REAL64 *Out, REAL32 *In, int n_3vec);
  
  void local_vcdot(REAL64 *Out_re, REAL64 *Out_im, REAL32 *V1, REAL32 *V2, int n_3vec);
  void local_vcdot_real(REAL64 *Out_re,  REAL32 *V1, REAL32 *V2, int n_3vec);
  

  typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 4> TVec;
  typedef PScalar<PScalar<RScalar<REAL32> > >  TScal;


  ////////////////////////////////
  // Threading evaluates
  //
  // by Xu Guo, EPCC, 6 October, 2008
  ////////////////////////////////
  // the wrappers for the functions to be threaded
#include "qdp_scalarsite_sse_blas_wrapper.h"
  


  /* #define DEBUG_BLAS_VAXMBY */
  /* #define DEBUG_BLAS_VAXPBY */
  
#define QDP_SCALARSITE_USE_EVALUATE


// TVec is the LatticeFermion from qdp_dwdefs.h with the OLattice<> stripped
// from around it

// TScalar is the usual Real, with the OScalar<> stripped from it
//
// THis is simply to make the code more readable, and reduces < < s and > >s
// in the template arguments


#if defined(QDP_SCALARSITE_USE_EVALUATE)

// d += Scalar*Vec
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< TScal, OScalar < TScal > > >,
	      Reference< QDPType< TVec, OLattice< TVec > > > >,
	      OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: y += a*x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());
  
  REAL32 ar = a.elem().elem().elem().elem();
  REAL32* aptr = &ar;
  if( s.hasOrderedRep() ) { 
    REAL32* xptr = (REAL32 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL32* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    // cout << "Specialised axpy a ="<< ar << endl;

    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {yptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
   
    const int* tab = s.siteTable().slice();

    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_y_user_arg arg(x, d, aptr, tab, 1, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_y_evaluate_function);

    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL32* xptr = (REAL32 *)&(x.elem(i).elem(0).elem(0).real());
      REAL32* yptr = &(d.elem(i).elem(0).elem(0).real());
    
      vaxpy3(yptr, aptr, xptr, yptr, 24);
      }
    */
  }


 }

// d -= Scalar*Vec
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< TScal, OScalar < TScal > > >,
	      Reference< QDPType< TVec, OLattice< TVec > > > >,
	      OLattice< TVec > > &rhs,
	      const Subset& s)
  {

#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: y -= a*x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());

  // - sign as y -= ax <=> y = y-ax = -ax + y = axpy with -a 
  REAL32 ar = -( a.elem().elem().elem().elem());
  REAL32* aptr = &ar;
  if( s.hasOrderedRep() ) { 
    REAL32* xptr = (REAL32 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL32* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {yptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////    
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
   
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_y_user_arg arg(x, d, aptr, tab, 1, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_y_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      REAL32* xptr = (REAL32 *)&(x.elem(i).elem(0).elem(0).real());
      REAL32* yptr = &(d.elem(i).elem(0).elem(0).real());
      vaxpy3(yptr, aptr, xptr, yptr, 24);
   
      }*/
  }
	
  }

// z = ax + y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = a*x + y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right());

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.right());
  // Set pointers 
  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////  
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxpy3(zptr, aptr, xptr, yptr, 24);
   
      }*/
  }

 }

// Vec = Vec + Scal*Vec
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = y + a*x" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.right());
  // Set pointers 
  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
    REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    vaxpy3(zptr, aptr, xptr, yptr, 24);
    
   
    }*/
  }
 }

// Vec = Scalar*Vec - Vec
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = a*x - y" << endl;
#endif


  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right());

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.right());
  // Set pointers 
  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep() ) { 
    
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxmy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxmy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxmy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      
      vaxmy3(zptr, aptr, xptr, yptr,24);
      }*/
  }
 }

template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = y - a*x" << endl;
#endif

  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.right());
  // Set pointers etc.

  // -ve sign as y - ax = -ax + y  = axpy with -a.
  REAL32 ar =  -a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep()) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxpy3(zptr, aptr, xptr, yptr, 24);
      }*/
  }
 }

// Vec += Vec * Scalar (AXPY)
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< TVec, OLattice< TVec > > >,
	      Reference< QDPType< TScal, OScalar < TScal > > > >,
	      OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: y += x*a" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().left());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().right());
  
  REAL32 ar = a.elem().elem().elem().elem();
  REAL32* aptr = &ar;
  if( s.hasOrderedRep() ) { 

    REAL32* xptr = (REAL32 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL32* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {yptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    //////////////// 
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_y_user_arg arg(x, d, aptr, tab, 1, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_y_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(d.elem(i).elem(0).elem(0).real());
      vaxpy3(yptr, aptr, xptr, yptr,24);   
      }*/
  }

 }


// Vec -= Vec *Scalar 
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< TVec, OLattice< TVec > > >,
	      Reference< QDPType< TScal, OScalar < TScal > > > >,
	      OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: y -= x*a" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().left());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().right());

  // - sign as y -= ax <=> y = y-ax = -ax + y = axpy with -a 
  REAL32 ar = -( a.elem().elem().elem().elem());
  REAL32* aptr = &ar;
  if( s.hasOrderedRep() ) {

    REAL32* xptr = (REAL32 *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL32* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {yptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////     
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_y_user_arg arg(x, d, aptr, tab, 1, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_y_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(d.elem(i).elem(0).elem(0).real());
   
      vaxpy3(yptr, aptr, xptr, yptr, 24);
      }*/
  }
	
 }


// Vec = Vec *Scalar  + Vec (AXPY)
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TScal, OScalar< TScal > > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = x*a + y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,    
    Reference< QDPType< TScal, OScalar< TScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.right());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.left());
  // Set pointers 
  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
      vaxpy3(zptr, aptr, xptr, yptr, 24);
      }*/
  }
 }


// Vec = Vec + Vec * Scalar (AXPY)
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TScal, OScalar< TScal > > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = y + x*a" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,    
    Reference< QDPType< TScal, OScalar< TScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.right());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.left());
  // Set pointers 
  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep() ) {
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
    REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    vaxpy3(zptr, aptr, xptr, yptr, 24);
    }*/
  }
}


// Vec = Vec*Scalar - Vec (AXMY)
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TScal, OScalar< TScal > > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = x*a - y" << endl;
#endif

  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right());


  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,    
    Reference< QDPType< TScal, OScalar< TScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().left());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.right());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.left());
  // Set pointers 
  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if ( s.hasOrderedRep() ) {
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxmy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxmy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxmy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmy3(zptr, aptr, xptr, yptr, 24);
      }*/
  }
}


// Vec = Vec - Vec*Scalar (AXPY with -Scalar)
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       BinaryNode<OpMultiply, 
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TScal, OScalar< TScal > > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE: z = y - x*a" << endl;
#endif

  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,    
    Reference< QDPType< TScal, OScalar< TScal > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.right());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode.left());
  // Set pointers etc.

  // -ve sign as y - ax = -ax + y  = axpy with -a.
  REAL32 ar =  -a.elem().elem().elem().elem();
  REAL32 *aptr = (REAL32 *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpy3_user_arg arg = {zptr, aptr, xptr, yptr, vaxpy3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpy3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpy3_z_user_arg arg(x, y, d, aptr, tab, vaxpy3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpy3_z_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxpy3(zptr, aptr, xptr, yptr, 24);
      }*/
  }
}


template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "SSE: v+v " << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(rhs.expression().left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vOp_user_arg arg = {zptr, xptr, yptr, vadd};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vOp_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vadd(zptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vOp_z_user_arg arg(x, y, d, tab, vadd);

    dispatch_to_threads(totalSize, arg, unordered_vOp_z_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vadd(zptr, xptr, yptr, 24);
      }*/
  }
}

template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "SSE: v-v " << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(rhs.expression().left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vOp_user_arg arg = {zptr, xptr, yptr, vsub};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vOp_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vsub(zptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vOp_z_user_arg arg(x, y, d, tab, vsub);

    dispatch_to_threads(totalSize, arg, unordered_vOp_z_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vsub(zptr, xptr, yptr, 24);
      }*/
  }
}

// Vec = Scal * Vec
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpMultiply,
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "SSE: v = a*v " << endl;
#endif
  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(rhs.expression().right());
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().left());

  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = &ar;  
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vscal_user_arg arg = {zptr, aptr, xptr};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vscal_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////  
    //int n_3vec = (s.end()-s.start()+1)*24;
    
    //vscal(zptr, aptr, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal_user_arg arg(x, d, aptr, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      vscal(zptr, aptr, xptr, 24);
      }*/
  }
}

template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpMultiply,
	       Reference< QDPType< TVec, OLattice< TVec > > >,
	       Reference< QDPType< TScal, OScalar< TScal > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS
  cout << "SSE: v = v*a " << endl;
#endif

  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(rhs.expression().left());
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().right());

  REAL32 ar =  a.elem().elem().elem().elem();
  REAL32 *aptr = &ar;  
  if( s.hasOrderedRep() ) {
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vscal_user_arg arg = {zptr, aptr, xptr};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vscal_evaluate_function);
    
    ////////////////
    // Original code
    //////////////// 
    //int n_3vec = (s.end()-s.start()+1)*24;

    //vscal(zptr, aptr, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal_user_arg arg(x, d, aptr, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      vscal(zptr, aptr, xptr, 24);
      }*/
  }
}


//-----------------------------------------------------------------------------
// v *= a
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpMultiplyAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< TScal, OScalar< TScal > > > >,
	       OScalar< TScal > > &rhs,
	       const Subset& s)
{
  const OScalar< TScal >& a = static_cast< const OScalar<TScal >&>(rhs.expression().child());


#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v *= a, a = " << a << endl;
#endif
  
  REAL32 ar = a.elem().elem().elem().elem();
  if( s.hasOrderedRep() ) { 
    REAL32 * xptr = &(d.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr = xptr;
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vscal_user_arg arg = {zptr, &ar, xptr};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vscal_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vscal(zptr,&ar, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal_user_arg arg(d, d, &ar, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 * xptr = &(d.elem(i).elem(0).elem(0).real());
      REAL32 * zptr = xptr;
      
      vscal(zptr, &ar, xptr, 24);
      }*/
  }
}

// v /= a
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpDivideAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< TScal, OScalar< TScal > > > >,
	       OScalar< TScal > > &rhs,
	       const Subset& s)
{
  const OScalar< TScal >& a = static_cast< const OScalar<TScal >&>(rhs.expression().child());


#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v /= a, a = " << a << endl;
#endif
  
  REAL32 ar = (REAL)1/a.elem().elem().elem().elem();
  if( s.hasOrderedRep() ) {
    REAL32 * xptr = &(d.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr = xptr;
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vscal_user_arg arg = {zptr, &ar, xptr};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vscal_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vscal(zptr,&ar, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vscal_user_arg arg(d, d, &ar, tab);

    dispatch_to_threads(totalSize, arg, unordered_sse_vscal_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 * xptr = &(d.elem(i).elem(0).elem(0).real());
      REAL32 * zptr = xptr;
      
      vscal(zptr,&ar, xptr, 24);
      }*/
  }
}

// v += v
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAddAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
  const OLattice< TVec >& x = static_cast< const OLattice<TVec >&>(rhs.expression().child());

 

#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v += v" << endl;
#endif
  if(s.hasOrderedRep() ) {
    REAL32 *xptr = (REAL32 *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *)(&d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vOp_user_arg arg = {yptr, yptr, xptr, vadd};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vOp_evaluate_function);
    
    ////////////////
    // Original code
    //////////////// 
    /*int n_3vec = (s.end() - s.start()+1)*24;
    REAL32 *xptr = (REAL32 *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *)(&d.elem(s.start()).elem(0).elem(0).real());
    REAL one = 1;
    vadd(yptr, yptr, xptr,n_3vec);*/
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vOp_y_user_arg arg(x, d, tab, vadd);

    dispatch_to_threads(totalSize, arg, unordered_sse_vOp_y_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *)(&x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *)(&d.elem(i).elem(0).elem(0).real());
      vadd(yptr, yptr, xptr,24);
      }*/
  }
}

// v -= v
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpSubtractAssign &op,
	       const QDPExpr< 
	       UnaryNode<OpIdentity,
	       Reference< QDPType< TVec, OLattice< TVec > > > >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{
  const OLattice< TVec >& x = static_cast< const OLattice<TVec >&>(rhs.expression().child());

 

#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: v -= v" << endl;
#endif
  if( s.hasOrderedRep() ) { 
    //int n_3vec = (s.end() - s.start()+1)*24;
    REAL32 *xptr = (REAL32 *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *)(&d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vOp_user_arg arg = {yptr, yptr, xptr, vsub};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vOp_evaluate_function);
    
    ////////////////
    // Original code
    //////////////// 
    //vsub(yptr, yptr, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vOp_y_user_arg arg(x, d, tab, vsub);

    dispatch_to_threads(totalSize, arg, unordered_sse_vOp_y_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *)(&x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *)(&d.elem(i).elem(0).elem(0).real());
      
      vsub(yptr, yptr, xptr, 24);
      }*/
  }
}


// z = ax + by
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x + b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
     
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    //////////////// 
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxpby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      vaxpby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}


// z = xa + by
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a + b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,
    Reference< QDPType< TScal, OScalar< TScal > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN2;


  
  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.left());

  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////  
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxpby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
     
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      vaxpby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}

// z = ax + yb
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x + y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // type of a*x
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN1;

  // type of y*b
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,
    Reference< QDPType< TScal, OScalar< TScal > > > > BN2;


  
  // get the binary nodes
  // a*x node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());

  // y*b node
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());

  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());

  
  // get b and y out of the binary node
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.left());

  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.right());

  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) {

    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    //////////////// 
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxpby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());


      vaxpby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}

// z = xa + yb
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a + y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,
    Reference< QDPType< TScal, OScalar< TScal > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.left());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.left());

  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.right());
  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxpby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////     
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxpby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. i and s.end() are inclusive so add +1
      vaxpby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}

// z = ax - by
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x - b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////       
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxmby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmby3(zptr, aptr, xptr, bptr, yptr, 24);

      }*/
  }
}


// z = xa - by
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a - b*y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,
    Reference< QDPType< TScal, OScalar< TScal > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN2;


  
  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.left());

  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right());

  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxmby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}

// z = ax - yb
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TScal, OScalar< TScal > > >,
	         Reference< QDPType< TVec, OLattice< TVec > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = a*x - y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // type of a*x
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN1;

  // type of y*b
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,
    Reference< QDPType< TScal, OScalar< TScal > > > > BN2;


  
  // get the binary nodes
  // a*x node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());

  // y*b node
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());

  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());

  
  // get b and y out of the binary node
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.left());

  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.right());

  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 * zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxmby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 * zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}

// z = xa - yb
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > >,
	        BinaryNode<OpMultiply, 
	         Reference< QDPType< TVec, OLattice< TVec > > >,
	         Reference< QDPType< TScal, OScalar< TScal > > > > >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "z = x*a - y*b" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TVec, OLattice< TVec > > >,
    Reference< QDPType< TScal, OScalar< TScal > > > > BN;

  // get the binary node
  const BN &mulNode1 = static_cast<const BN&> (rhs.expression().left());
  const BN &mulNode2 = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the binary node
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.left());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.left());

  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.right());
  
  // Set pointers 
  REAL32 *aptr = (REAL32 *)&(a.elem().elem().elem().elem());
  REAL32 *bptr = (REAL32 *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL32 *xptr = (REAL32 *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL32 *yptr = (REAL32 *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL32 *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    int total_n_3vec = (s.end()-s.start()+1);

    ordered_sse_vaxOpby3_user_arg arg = {zptr, aptr, xptr, bptr, yptr, vaxmby3};

    dispatch_to_threads(total_n_3vec, arg, ordered_sse_vaxOpby3_evaluate_function);
    
    ////////////////
    // Original code
    ////////////////   
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    //int n_3vec = (s.end()-s.start()+1)*24;
    //vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    
    int totalSize = s.numSiteTable();

    unordered_sse_vaxOpby3_user_arg arg(aptr, x, bptr, y, d, tab, vaxmby3);

    dispatch_to_threads(totalSize, arg, unordered_sse_vaxOpby3_evaluate_function );
    
    ////////////////
    // Original code
    ////////////////
    /*
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];
      REAL32 *xptr = (REAL32 *) &(x.elem(i).elem(0).elem(0).real());
      REAL32 *yptr = (REAL32 *) &(y.elem(i).elem(0).elem(0).real());
      REAL32 *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      vaxmby3(zptr, aptr, xptr, bptr, yptr, 24);
      }*/
  }
}



//-----------------------------------------------------------------------------

#endif     // if defined(QDP_SCALARSITE_USE_EVALUATE)



#if 1
// Global norm squared of a vector...
template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
norm2(const QDPType<TVec ,OLattice< TVec > >& s1, const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using SSE sumsq" << endl;
#endif

  if ( s.hasOrderedRep() ) {

#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ sumsq " << endl;
#endif
    int n_real = (s.end() - s.start() + 1);
    ordered_sse_norm_single_user_arg arg;
    arg.vptr = (REAL32*) &(s1.elem(s.start()).elem(0).elem(0).real());
    arg.results = ThreadReductions::norm2_results;
    arg.func = local_sumsq_24_48;
    dispatch_to_threads(n_real, arg, ordered_norm_single_func);
    
    REAL64 ltmp=arg.results[0];
    for(int i=1; i < qdpNumThreads(); i++) { 
	ltmp += arg.results[i];
    }
    // Use specialized sum for REAL64
    QDPInternal::globalSum(ltmp);
    UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  lsum(ltmp);
    return lsum;
  }
  else {
    REAL64 ltmp1=0;


    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL64 ltmp2=0;
      REAL32* s1ptr = (REAL32 *)&(s1.elem(i).elem(0).elem(0).real()); 
      local_sumsq_24_48(&ltmp2, s1ptr, 1);
      ltmp1 += ltmp2;
    }

    // Use specialized sum for REAL64
    QDPInternal::globalSum(ltmp1);
    UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  lsum(ltmp1);

    return lsum;
    
  }
}

template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
norm2(const QDPType<TVec ,OLattice< TVec > >& s1)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using SSE sumsq all" << endl;
#endif

  int n_real = (all.end() - all.start() + 1);
  ordered_sse_norm_single_user_arg arg;
  arg.vptr = (REAL32*) &(s1.elem(all.start()).elem(0).elem(0).real());
  arg.results = ThreadReductions::norm2_results;
  arg.func = local_sumsq_24_48;
  dispatch_to_threads(n_real, arg, ordered_norm_single_func);
  
    
  // I am relying on this being a Double here 
  REAL64 ltmp=arg.results[0];
  for(int i=1; i < qdpNumThreads(); i++) { 
	ltmp += arg.results[i];
  }
  // Do the sum with specialized QMP_sum_double
  QDPInternal::globalSum(ltmp);
  UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  lsum(ltmp);
  return lsum;
}
#endif


template<>
inline  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t
innerProduct(const QDPType< TVec, OLattice<TVec> > &v1,
	     const QDPType< TVec, OLattice<TVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProduct all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t lprod;
  // Inner product is accumulated internally in REAL64
  REAL64 ip[2];
  ip[0]=0;
  ip[1]=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*Ns;
    
  // Call My CDOT
  local_vcdot((REAL64 *)&(ip[0]), (REAL64 *)&(ip[1]),
	      (REAL32 *)&(v1.elem(all.start()).elem(0).elem(0).real()),
	      (REAL32 *)&(v2.elem(all.start()).elem(0).elem(0).real()),
	      (int)n_3vec);


  // Global sum -- still on a vector of doubles
  QDPInternal::globalSumArray(ip,2);

  // Downcast (and possibly lose precision) here 
  lprod.elem().elem().elem().real() = ip[0];
  lprod.elem().elem().elem().imag() = ip[1];

  // Return
  return lprod;
}

template<>
inline  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t
innerProduct(const QDPType< TVec, OLattice<TVec> > &v1,
	     const QDPType< TVec, OLattice<TVec> > &v2, 
	     const Subset& s)
{

  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t lprod;

  REAL64 ip[2];
  if( s.hasOrderedRep() ) {
#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ: innerProduct s" << endl;
#endif

    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >

   
    ip[0] = 0;
    ip[1] = 0;

    unsigned long n_3vec = (s.end() - s.start() + 1)*Ns;
    local_vcdot((REAL64 *)&(ip[0]), (REAL64 *)&(ip[1]),
		(REAL32 *)&(v1.elem(s.start()).elem(0).elem(0).real()),
		(REAL32 *)&(v2.elem(s.start()).elem(0).elem(0).real()),
		(int)n_3vec);

    
  }
  else {
    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >

    REAL64 ip_tmp[2];
    ip[0] = 0;
    ip[1] = 0;

    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];

      local_vcdot((REAL64 *)&(ip_tmp[0]), (REAL64 *)&(ip_tmp[1]),
		  (REAL32 *)&(v1.elem(i).elem(0).elem(0).real()),
		  (REAL32 *)&(v2.elem(i).elem(0).elem(0).real()),
		  (int)Ns);
      
      ip[0] += ip_tmp[0];
      ip[1] += ip_tmp[1];
    }
  }

  QDPInternal::globalSumArray(ip,2);
  
  lprod.elem().elem().elem().real() = ip[0];
  lprod.elem().elem().elem().imag() = ip[1];

  return lprod;
}



// Inner Product Real
template<>
inline  
BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t
innerProductReal(const QDPType< TVec, OLattice<TVec> > &v1,
		 const QDPType< TVec, OLattice<TVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t lprod;
  // Inner product is accumulated internally in REAL64
  REAL64 ip_re=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*Ns;

  // Call My CDOT
  local_vcdot_real((DOUBLE*)&ip_re,
		   (REAL32 *)&(v1.elem(all.start()).elem(0).elem(0).real()),
		   (REAL32 *)&(v2.elem(all.start()).elem(0).elem(0).real()),
		   (int)n_3vec);

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
BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t
innerProductReal(const QDPType< TVec, OLattice<TVec> > &v1,
		 const QDPType< TVec, OLattice<TVec> > &v2, 
		 const Subset& s)
{
  
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t lprod;

  REAL64 ip_re=0;

  if( s.hasOrderedRep() ) {
#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ: innerProductReal s" << endl;
#endif

    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RScalar<PScalar<REAL> > > > >
    unsigned long n_3vec = (s.end() - s.start() + 1)*Ns;
    local_vcdot_real((REAL64 *)&ip_re,
		     (REAL32 *)&(v1.elem(s.start()).elem(0).elem(0).real()),
		     (REAL32 *)&(v2.elem(s.start()).elem(0).elem(0).real()),
		     (int)n_3vec);




  }
  else {
    REAL64 ip_re_tmp=0;
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];

      local_vcdot_real((REAL64 *)&ip_re_tmp,
		       (REAL32 *)&(v1.elem(i).elem(0).elem(0).real()),
		       (REAL32 *)&(v2.elem(i).elem(0).elem(0).real()),
		       (int)Ns);

      ip_re += ip_re_tmp;

    }


  }

  QDPInternal::globalSum(ip_re);
  lprod.elem().elem().elem().elem() = ip_re;
  return lprod;
}

#if 1
template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
norm2(const multi1d< OLattice< TVec > >& s1)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using SSE multi1d sumsq all" << endl;
#endif

  int n_real = (all.end() - all.start() + 1);
  REAL64 ltmp = 0;
  for(int n=0; n < s1.size(); ++n)
  {
    const REAL32 *s1ptr =  &(s1[n].elem(all.start()).elem(0).elem(0).real());
    
    // I am relying on this being a Double here 
    REAL64 lltmp=0;
    local_sumsq_24_48(&lltmp, (REAL32 *)s1ptr, n_real); 

    ltmp += lltmp;
  }

  QDPInternal::globalSum(ltmp);
  UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  lsum(ltmp);
  
  return lsum;
}


template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
  norm2(const multi1d< OLattice< TVec > >& s1, const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using SSE multi1d sumsq all" << endl;
#endif

  REAL64 ltmp = 0;

  if( s.hasOrderedRep() ) { 
    int n_real = (s.end() - s.start() + 1);

    for(int n=0; n < s1.size(); ++n) {
  
      const REAL32 *s1ptr =  &(s1[n].elem(s.start()).elem(0).elem(0).real());
    
      // I am relying on this being a Double here 
      REAL64 lltmp=0;
      local_sumsq_24_48(&lltmp, (REAL32 *)s1ptr, n_real); 
      
      ltmp += lltmp;
    }
  }
  else  {

    const int* tab = s.siteTable().slice();
    for(int n=0; n < s1.size(); ++n) { 
      for(int j=0; j < s.numSiteTable(); j++) { 
	int i=tab[j];
	REAL64 lltmp=0;
	const REAL32 *s1ptr =  &(s1[n].elem(i).elem(0).elem(0).real());
	local_sumsq_24_48(&lltmp,(REAL32 *)s1ptr,1);
	ltmp += lltmp;
      }
    }   

  }
  QDPInternal::globalSum(ltmp);
  UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  lsum(ltmp);

  return lsum;
}
#endif

template<>
inline  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OLattice<TVec> > &v1,
	     const multi1d< OLattice<TVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: multi1d innerProduct all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t lprod;
  // Inner product is accumulated internally in REAL64
  REAL64 ip[2];
  ip[0]=0;
  ip[1]=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*Ns;
    
  for(int n=0; n < v1.size(); ++n)
  {
    REAL64 iip[2];
    iip[0]=0;
    iip[1]=0;

    // Call My CDOT
    local_vcdot((REAL64 *)&(iip[0]), (REAL64 *)&(iip[1]),
		(REAL32 *)&(v1[n].elem(all.start()).elem(0).elem(0).real()),
		(REAL32 *)&(v2[n].elem(all.start()).elem(0).elem(0).real()),
		(int)n_3vec);
    
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

template<>
inline  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OLattice<TVec> > &v1,
	     const multi1d< OLattice<TVec> > &v2,
	     const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: multi1d innerProduct subset" << endl;
#endif

  
  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProduct>::Type_t lprod;
  // Inner product is accumulated internally in REAL64
  REAL64 ip[2];
  ip[0]=0;
  ip[1]=0;

  if( s.hasOrderedRep() ) { 
    // Length of subset 
    unsigned long n_3vec = (s.end() - s.start() + 1)*Ns;
    
    for(int n=0; n < v1.size(); ++n) {
      
      REAL64 iip[2];
      iip[0]=0;
      iip[1]=0;
      
      // Call My CDOT
      local_vcdot((REAL64 *)&(iip[0]), 
		  (REAL64 *)&(iip[1]),
		  (REAL32 *)&(v1[n].elem(s.start()).elem(0).elem(0).real()),
		  (REAL32 *)&(v2[n].elem(s.start()).elem(0).elem(0).real()),
		  (int)n_3vec);
    
      ip[0] += iip[0];
      ip[1] += iip[1];
    }
  }
  else { 
     // Length of an atom
    unsigned long n_3vec = Ns;

    // Site table
    const int* tab = s.siteTable().slice();

    // Loop through N5
    for(int n=0; n < v1.size(); ++n) { 

      // Loop through site table 
      for(int j=0; j < s.numSiteTable(); j++) { 
	int i=tab[j];    
	REAL64 iip[2];
	iip[0]=0;
	iip[1]=0;
	
	// Call My CDOT
	local_vcdot((REAL64 *)&(iip[0]), (REAL64 *)&(iip[1]),
		    (REAL32 *)&(v1[n].elem(i).elem(0).elem(0).real()),
		    (REAL32 *)&(v2[n].elem(i).elem(0).elem(0).real()),
		    (int)n_3vec);
	
	ip[0] += iip[0];
	ip[1] += iip[1];
      }
    }
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
BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<TVec> > &v1,
		 const multi1d< OLattice<TVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal(multi1d) all" << endl;
#endif

  // This BinaryReturn hasType_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t lprod;
  // Inner product is accumulated internally in REAL64
  REAL64 ip_re=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*Ns;

  for(int n=0; n < v1.size(); ++n)
  {
    REAL64 iip_re=0;

    // Call My CDOT
    local_vcdot_real((REAL64 *)&iip_re,
		     (REAL32 *)&(v1[n].elem(all.start()).elem(0).elem(0).real()),
		     (REAL32 *)&(v2[n].elem(all.start()).elem(0).elem(0).real()),
		     (int)n_3vec);

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

// Inner Product Real
template<>
inline  
BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<TVec> > &v1,
		 const multi1d< OLattice<TVec> > &v2,
		 const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal(multi1d) all" << endl;
#endif

  // This BinaryReturn hasType_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t lprod;
  // Inner product is accumulated internally in REAL64
  REAL64 ip_re=0;

  // Length of subset 
  if ( s.hasOrderedRep() ) {
    unsigned long n_3vec = (s.end() - s.start() + 1)*Ns;

    for(int n=0; n < v1.size(); ++n) {
      
      REAL64 iip_re=0;
      
      // Call My CDOT
      local_vcdot_real((REAL64 *)&iip_re,
		       (REAL32 *)&(v1[n].elem(s.start()).elem(0).elem(0).real()),
		       (REAL32 *)&(v2[n].elem(s.start()).elem(0).elem(0).real()),
		       (int)n_3vec);
      
      ip_re += iip_re;
    }
  }
  else { 
    unsigned long n_3vec = Ns;
    const int* tab = s.siteTable().slice();

    for(int n=0; n < v1.size(); ++n) {
      for(int j=0; j< s.numSiteTable(); j++) { 
	const int i=tab[j];
	 REAL64 iip_re=0;
      
	 // Call My CDOT
	 local_vcdot_real((REAL64 *)&iip_re,
		       (REAL32 *)&(v1[n].elem(i).elem(0).elem(0).real()),
		       (REAL32 *)&(v2[n].elem(i).elem(0).elem(0).real()),
		       (int)n_3vec);
      
	 ip_re += iip_re;
      }
    }
  }
  // Global sum
  QDPInternal::globalSum(ip_re);

  // Whether CDOT did anything or not ip_re and ip_im should 
  // now be right. Assign them to the ReturnType
  lprod.elem().elem().elem().elem() = ip_re;


  // Return
  return lprod;
}




#if defined(QDP_SCALARSITE_DEBUG)
#undef QDP_SCALARSITE_DEBUG
#endif

#if defined(QDP_SCALARSITE_USE_EVALUATE)
#undef QDP_SCALARSITE_USE_EVALUATE
#endif

  
} // namespace QDP;

#endif  // guard
 
