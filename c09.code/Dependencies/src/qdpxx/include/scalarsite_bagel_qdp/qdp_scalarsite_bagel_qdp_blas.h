// $Id: qdp_scalarsite_bagel_qdp_blas.h,v 1.4 2007-06-10 14:32:10 edwards Exp $

/*! @file
 * @brief Generic Scalarsite  optimization hooks
 * 
 * 
 */


#ifndef QDP_SCALARSITE_BAGEL_QDP_BLAS_H
#define QDP_SCALARSITE_BAGEL_QDP_BLAS_H

#include "bagel_qdp.h"

namespace QDP {

// Types needed for the expression templates. 
// TVec has outer Ns template so it ought to work for staggered as well
typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, Ns> TVec;
typedef PScalar<PScalar<RScalar<REAL> > >  TScal;

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
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference< QDPType< TScal, OScalar < TScal > > >,
	      Reference< QDPType< TVec, OLattice< TVec > > > >,
	      OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "y += a*x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());
  
  REAL ar = a.elem().elem().elem().elem();
  REAL* aptr = &ar;

  if( s.hasOrderedRep() ) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    // cout << "Specialised axpy a ="<< ar << endl;
    
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];

      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = &(d.elem(i).elem(0).elem(0).real());
    
      qdp_vaxpy3(yptr, aptr, xptr, yptr, 4);
    }
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
  QDPIO::cout << "y -= a*x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());

  // - sign as y -= ax <=> y = y-ax = -ax + y = axpy with -a 
  REAL ar = -( a.elem().elem().elem().elem());
  REAL* aptr = &ar;

  if( s.hasOrderedRep() ) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];

      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = &(d.elem(i).elem(0).elem(0).real());
      qdp_vaxpy3(yptr, aptr, xptr, yptr, 4);
   
    }
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
  QDPIO::cout << "z = a*x + y" << endl;
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
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      qdp_vaxpy3(zptr, aptr, xptr, yptr, 4);
   
    }
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
  QDPIO::cout << "z = y + a*x" << endl;
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
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());

      qdp_vaxpy3(zptr, aptr, xptr, yptr, 4);
    }
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
  QDPIO::cout << "z = a*x - y" << endl;
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
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxmy3(zptr, aptr, xptr, yptr, n_3vec);

  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());

      qdp_vaxmy3(zptr, aptr, xptr, yptr, 4);
    }
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
  QDPIO::cout << "z = y - a*x" << endl;
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
  REAL ar =  -a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr, aptr, xptr, yptr, n_3vec);

  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());

      qdp_vaxpy3(zptr, aptr, xptr, yptr, 4);
    }
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
  QDPIO::cout << "y += x*a" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().left());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().right());
  
  REAL ar = a.elem().elem().elem().elem();
  REAL* aptr = &ar;

  if( s.hasOrderedRep() ) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    // cout << "Specialised axpy a ="<< ar << endl;
    
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(d.elem(i).elem(0).elem(0).real());
      qdp_vaxpy3(yptr, aptr, xptr, yptr, 4);
    }
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
  QDPIO::cout << "y -= x*a" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().left());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().right());

  // - sign as y -= ax <=> y = y-ax = -ax + y = axpy with -a 
  REAL ar = -( a.elem().elem().elem().elem());
  REAL* aptr = &ar;

  if( s.hasOrderedRep() ) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = &(d.elem(s.start()).elem(0).elem(0).real());
    
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(yptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(d.elem(i).elem(0).elem(0).real());
      qdp_vaxpy3(yptr, aptr, xptr, yptr, 4);
    }
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
  QDPIO::cout << "z = x*a + y" << endl;
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
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) { 
     REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
  
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxpy3(zptr, aptr, xptr, yptr, 4);

    }
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
  QDPIO::cout << "z = y + x*a" << endl;
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
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxpy3(zptr, aptr, xptr, yptr, 4);

    }
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
  QDPIO::cout << "z = x*a - y" << endl;
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
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxmy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxmy3(zptr, aptr, xptr, yptr, 4);

    }
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
  QDPIO::cout << "z = y - x*a" << endl;
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
  REAL ar =  -a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr, aptr, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxpy3(zptr, aptr, xptr, yptr, 4);

    }
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
  cout << "BJ: v+v " << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(rhs.expression().left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(rhs.expression().right());
  REAL one = 1;

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpy3(zptr,&one, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxpy3(zptr,&one, xptr, yptr, 4);
    }
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
  cout << "BJ: v-v " << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(rhs.expression().left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(rhs.expression().right());
  REAL one=1;

  if( s.hasOrderedRep()) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;

    qdp_vaxmy3(zptr,&one, xptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxmy3(zptr,&one, xptr, yptr, 4);
    }
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
  cout << "BJ: v = a*v " << endl;
#endif
  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(rhs.expression().right());
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().left());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = &ar;  

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vscal3(zptr, aptr, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      qdp_vscal3(zptr, aptr, xptr, 4);
      
    }
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
  cout << "BJ: v = v*a " << endl;
#endif

  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(rhs.expression().left());
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().right());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = &ar;  

  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    int n_3vec = (s.end()-s.start()+1)*4;

    qdp_vscal3(zptr, aptr, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      qdp_vscal3(zptr, aptr, xptr, 4);
    }
  }        

}

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
  
  REAL ar = a.elem().elem().elem().elem();

  if( s.hasOrderedRep() ) { 
    REAL* xptr = &(d.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr = xptr;
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vscal3(zptr,&ar, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(d.elem(i).elem(0).elem(0).real());
      REAL* zptr = xptr;
      qdp_vscal3(zptr, &ar, xptr, 4);
    }
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
  
  REAL ar = (REAL)1/a.elem().elem().elem().elem();
  if(s.hasOrderedRep()) { 
    REAL* xptr = &(d.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr = xptr;
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vscal3(zptr,&ar, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(d.elem(i).elem(0).elem(0).real());
      REAL* zptr = xptr;
      qdp_vscal3(zptr, &ar, xptr, 4);
    }
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



  if( s.hasOrderedRep() ) {
    int n_3vec = (s.end() - s.start()+1)*4;
    REAL *xptr = (REAL *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *)(&d.elem(s.start()).elem(0).elem(0).real());

    qdp_vadd3(yptr, yptr, xptr,n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *)(&d.elem(i).elem(0).elem(0).real());

      qdp_vadd3(yptr, yptr, xptr,4);
    }
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
    int n_3vec = (s.end() - s.start()+1)*4;
    REAL *xptr = (REAL *)(&x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *)(&d.elem(s.start()).elem(0).elem(0).real());
  
    qdp_vsub3(yptr, yptr, xptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *)(&d.elem(i).elem(0).elem(0).real());

      qdp_vsub3(yptr, yptr, xptr, 4);
    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      qdp_vaxpby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
          
      qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep()) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
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
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_3vec = (s.end()-s.start()+1)*4;
    qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, n_3vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      qdp_vaxmby3(zptr, aptr, xptr, bptr, yptr, 4);

    }
  }        

}


// Global norm squared of a vector...
template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
norm2(const QDPType<TVec ,OLattice< TVec > >& s1, const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using BJ sumsq" << endl;
#endif

  DOUBLE lsum=(DOUBLE)0;

  if ( s.hasOrderedRep() ) {

#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ sumsq " << endl;
#endif
    int n_3vec = (s.end() - s.start() + 1)*4;
    const REAL *s1ptr =  &(s1.elem(s.start()).elem(0).elem(0).real());

    qdp_lsum2(&lsum,(REAL *)s1ptr, n_3vec); 

  }
  else {

    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i = tab[j];
      const REAL *s1ptr =  &(s1.elem(i).elem(0).elem(0).real());
      DOUBLE ltmp=(DOUBLE)0;
      qdp_lsum2(&ltmp,(REAL *)s1ptr, 4); 
      lsum += ltmp;
    }


  }

  // Global sum in double for sure
  QDPInternal::globalSum(lsum);
  
  // Possibly downcast 
  UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  gsum(lsum);
  return gsum;
  
}


template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
norm2(const QDPType<TVec ,OLattice< TVec > >& s1)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "Using BJ sumsq all" << endl;
#endif


  int n_3vec = (all.end() - all.start() + 1)*4;
  const REAL *s1ptr =  &(s1.elem(all.start()).elem(0).elem(0).real());

  DOUBLE lsum = 0;
  qdp_lsum2(&lsum, (REAL *)s1ptr, n_3vec); 
  QDPInternal::globalSum(lsum);

  UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  gsum(lsum);
  return gsum;
}



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
  // Inner product is accumulated internally in DOUBLE
  DOUBLE ip[2] = {0,0};

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*4;
    
  // Call My CDOT
  qdp_lcdot(&(ip[0]), &(ip[1]),
	      (REAL *)&(v1.elem(all.start()).elem(0).elem(0).real()),
	      (REAL *)&(v2.elem(all.start()).elem(0).elem(0).real()),
	      n_3vec);


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
  DOUBLE ip[2]={0,0};

  if( s.hasOrderedRep() ) {
#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ: innerProduct s" << endl;
#endif

    // This BinaryReturn has Type_t
    // OScalar<OScalar<OScalar<RComplex<PScalar<REAL> > > > >

    unsigned long n_3vec = (s.end() - s.start() + 1)*4;
    qdp_lcdot(&(ip[0]), &(ip[1]),
		(REAL *)&(v1.elem(s.start()).elem(0).elem(0).real()),
		(REAL *)&(v2.elem(s.start()).elem(0).elem(0).real()),
		n_3vec);

  }
  else {

    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++){
      int i=tab[j];
      DOUBLE ip_tmp[2]={0,0};       

      qdp_lcdot(&(ip_tmp[0]), &(ip_tmp[1]),
		(REAL *)&(v1.elem(i).elem(0).elem(0).real()),
		(REAL *)&(v2.elem(i).elem(0).elem(0).real()),
		4);

      ip[0]+=ip_tmp[0];
      ip[1]+=ip_tmp[1];
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
  // Inner product is accumulated internally in DOUBLE
  DOUBLE ip_re=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*4;

  // Call My CDOT
  qdp_lcdotr(&ip_re,
		   (REAL *)&(v1.elem(all.start()).elem(0).elem(0).real()),
		   (REAL *)&(v2.elem(all.start()).elem(0).elem(0).real()),
		   n_3vec);

  // Global sum
  QDPInternal::globalSum(ip_re);

  // Whether CDOT did anything or not ip_re and ip_im should 
  // now be right. Assign them to the ReturnType
  lprod.elem().elem().elem().elem() = ip_re;
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

  DOUBLE ip_re=0;


  if( s.hasOrderedRep() ) {
#ifdef DEBUG_BLAS
    QDPIO::cout << "BJ: innerProductReal s" << endl;
#endif

    unsigned long n_3vec = (s.end() - s.start() + 1)*4;
    qdp_lcdotr(&ip_re,
	       (REAL *)&(v1.elem(s.start()).elem(0).elem(0).real()),
	       (REAL *)&(v2.elem(s.start()).elem(0).elem(0).real()),
	       n_3vec);

  }
  else {

    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++){
      int i=tab[j];
      DOUBLE ip_re_tmp=0;
      qdp_lcdotr(&ip_re_tmp,
		 (REAL *)&(v1.elem(i).elem(0).elem(0).real()),
		 (REAL *)&(v2.elem(i).elem(0).elem(0).real()),
		 4);
     
      ip_re +=ip_re_tmp;
    }
  }

  QDPInternal::globalSum(ip_re);
  lprod.elem().elem().elem().elem() = ip_re;

  return lprod;

}



template<>
inline UnaryReturn<OLattice< TVec >, FnNorm2>::Type_t
norm2(const multi1d< OLattice< TVec > >& s1)
{

#ifdef DEBUG_BLAS
  QDPIO::cout << "Using SSE multi1d sumsq all" << endl;
#endif

  int n_3vec = (all.end() - all.start() + 1)*4;
  DOUBLE ltmp = 0;

  for(int n=0; n < s1.size(); ++n)
  {
    const REAL* s1ptr =  &(s1[n].elem(all.start()).elem(0).elem(0).real());
    
    // I am relying on this being a Double here 
    DOUBLE lltmp=0;
    qdp_lsum2(&lltmp, (REAL*)s1ptr, n_3vec); 

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
  
  DOUBLE ltmp = 0;
  for(int n=0; n < s1.size(); ++n) {
    
    
    if( s.hasOrderedRep() ) { 

      int n_3vec = (s.end() - s.start() + 1)*4;
      
      const REAL* s1ptr =  &(s1[n].elem(s.start()).elem(0).elem(0).real());
      
      // I am relying on this being a Double here 
      DOUBLE lltmp=0;
      qdp_lsum2(&lltmp, (REAL*)s1ptr, n_3vec); 
      
      ltmp += lltmp;
    }
    else { 
      const int* tab= s.siteTable().slice();
      for(int j=0; j < s.numSiteTable(); j++) {
	int i = tab[j];

	const REAL* s1ptr =  &(s1[n].elem(i).elem(0).elem(0).real());
    
	// I am relying on this being a Double here 
	DOUBLE lltmp=0;
	
	qdp_lsum2(&lltmp, (REAL*)s1ptr, 4); 
	
	ltmp += lltmp;
      }
    }
  }

  QDPInternal::globalSum(ltmp);
  UnaryReturn< OLattice< TVec >, FnNorm2>::Type_t  lsum(ltmp);
  return lsum;

}


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
  // Inner product is accumulated internally in DOUBLE
  DOUBLE ip[2];
  ip[0]=0;
  ip[1]=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*4;
    
  for(int n=0; n < v1.size(); ++n)
  {
    DOUBLE iip[2]={0,0};

  
    // Call My CDOT
    qdp_lcdot(&(iip[0]), &(iip[1]),
	      (REAL *)&(v1[n].elem(all.start()).elem(0).elem(0).real()),
	      (REAL *)&(v2[n].elem(all.start()).elem(0).elem(0).real()),
	      n_3vec);
    
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
BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<TVec> > &v1,
		 const multi1d< OLattice<TVec> > &v2)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal(multi1d) all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t lprod;

  // Inner product is accumulated internally in DOUBLE
  DOUBLE ip_re=0;

  // Length of subset 
  unsigned long n_3vec = (all.end() - all.start() + 1)*4;

  for(int n=0; n < v1.size(); ++n)
  {
    DOUBLE iip_re=0;

    // Call My CDOT
    qdp_lcdotr(&iip_re,
	       (REAL *)&(v1[n].elem(all.start()).elem(0).elem(0).real()),
	       (REAL *)&(v2[n].elem(all.start()).elem(0).elem(0).real()),
	       n_3vec);

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

template<>
inline  
BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<TVec> > &v1,
		 const multi1d< OLattice<TVec> > &v2, const Subset& s)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "BJ: innerProductReal(multi1d) all" << endl;
#endif

  // This BinaryReturn has Type_t
  // OScalar<OScalar<OScalar<RScalar<PScalar<REAL> > > > >
  BinaryReturn< OLattice<TVec>, OLattice<TVec>, FnInnerProductReal>::Type_t lprod;

  // Inner product is accumulated internally in DOUBLE
  DOUBLE ip_re=0;

  for(int n=0; n < v1.size(); ++n) {
    

    if( s.hasOrderedRep() ) {
      // Length of subset 
      unsigned long n_3vec = (s.end() - s.start() + 1)*4;

      DOUBLE iip_re=0;

      // Call My CDOT
      qdp_lcdotr(&iip_re,
		 (REAL *)&(v1[n].elem(s.start()).elem(0).elem(0).real()),
		 (REAL *)&(v2[n].elem(s.start()).elem(0).elem(0).real()),
		 n_3vec);

      ip_re += iip_re;
    }
    else { 
      const int* tab=s.siteTable().slice();
      for(int j=0; j < s.numSiteTable(); j++) { 
	int i=tab[j];

	DOUBLE iip_re=0;

	// Call My CDOT
	qdp_lcdotr(&iip_re,
		   (REAL *)&(v1[n].elem(i).elem(0).elem(0).real()),
		   (REAL *)&(v2[n].elem(i).elem(0).elem(0).real()),
		   4);

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


} // namespace QDP;

#endif  // guard
 
