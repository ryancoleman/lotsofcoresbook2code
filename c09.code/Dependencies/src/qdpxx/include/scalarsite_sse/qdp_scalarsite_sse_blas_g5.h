// $Id: qdp_scalarsite_sse_blas_g5.h,v 1.8 2009-02-11 20:50:45 bjoo Exp $

/*! @file
 * @brief Generic Scalarsite  optimization hooks
 * 
 * 
 */


#ifndef QDP_SCALARSITE_SSE_BLAS_G5_H
#define QDP_SCALARSITE_SSE_BLAS_G5_H

#include "scalarsite_sse/qdp_scalarsite_sse_blas_g5_includes.h"

using namespace QDP;

namespace QDP {

// Types needed for the expression templates. 
// TVec has outer Ns template so it ought to work for staggered as well
typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 4> TVec;
typedef PScalar<PScalar<RScalar<REAL32> > >  TScal;

// #define DEBUG_BLAS_G6
// TVec is the LatticeFermion from qdp_dwdefs.h with the OLattice<> stripped
// from around it

// TScalar is the usual Real, with the OScalar<> stripped from it
//
// THis is simply to make the code more readable, and reduces < < s and > >s
// in the template arguments

// d += Scalar*ChiralProjPlus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                     Reference< QDPType< TScal, OScalar < TScal > > >,
	                     UnaryNode< FnChiralProjectPlus, Reference< QDPType<TVec,OLattice<TVec> > > >
>, 
	                     OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y += a*P{+}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right().child());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());

  REAL ar = a.elem().elem().elem().elem();
  REAL* aptr = &ar;
  
  if( s.hasOrderedRep()) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
        
    int n_4vec = (s.end()-s.start()+1);
    xpayz_g5ProjPlus(yptr, aptr,yptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
     
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
      xpayz_g5ProjPlus(yptr, aptr,yptr, xptr, 1);
      
    }
  }
  
}

// d += Scalar*ChiralProjMinus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                     Reference< QDPType< TScal, OScalar < TScal > > >,
	                     UnaryNode< FnChiralProjectMinus, Reference< QDPType<TVec,OLattice<TVec> > > >
>, 
	                     OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y += a*P{-}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right().child());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());

  
  REAL ar = a.elem().elem().elem().elem();
  REAL* aptr = &ar;

  if( s.hasOrderedRep()) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
  
  
    int n_4vec = (s.end()-s.start()+1);
    xpayz_g5ProjMinus(yptr, aptr,yptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
     
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
      xpayz_g5ProjMinus(yptr, aptr,yptr, xptr, 1);
      
    }
  }

}


// d -= Scalar*ChiralProjPlus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                     Reference< QDPType< TScal, OScalar < TScal > > >,
	                     UnaryNode< FnChiralProjectPlus, Reference< QDPType<TVec,OLattice<TVec> > > >
>, 
	                     OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y -= a*P{+}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right().child());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());

  
  REAL ar = a.elem().elem().elem().elem();
  REAL* aptr = &ar;

  if( s.hasOrderedRep()) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
  
  
    int n_4vec = (s.end()-s.start()+1);
    xmayz_g5ProjPlus(yptr, aptr,yptr, xptr, n_4vec);
  }
  else {
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
      xmayz_g5ProjPlus(yptr, aptr,yptr, xptr, 1);
      
    }
  }
  
}

// d -= Scalar*ChiralProjMinus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                     Reference< QDPType< TScal, OScalar < TScal > > >,
	                     UnaryNode< FnChiralProjectMinus, Reference< QDPType<TVec,OLattice<TVec> > > >
>, 
	                     OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y -= a*P{-}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().right().child());
  const OScalar< TScal >& a = static_cast<const OScalar< TScal > &> (rhs.expression().left());

  
  REAL ar = a.elem().elem().elem().elem();
  REAL* aptr = &ar;

  if( s.hasOrderedRep()) {

    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
  
  
    int n_4vec = (s.end()-s.start()+1);
    xmayz_g5ProjMinus(yptr, aptr,yptr, xptr, n_4vec);
  }
  else {
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
      xmayz_g5ProjMinus(yptr, aptr,yptr, xptr, 1);
      
    }
  }

}


// d += ChiralProjPlus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<
	                    UnaryNode< FnChiralProjectPlus, Reference< QDPType<TVec,OLattice<TVec> > > >, 
	                    OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y += P{+}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().child());

  if( s.hasOrderedRep()) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
  
  
    int n_4vec = (s.end()-s.start()+1);
    add_g5ProjPlus(yptr, yptr, xptr, n_4vec);
  }
  else {
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
  
    
      add_g5ProjPlus(yptr, yptr, xptr, 1);
    }
  }


  
}


// d += ChiralProjMinus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<
	                    UnaryNode< FnChiralProjectMinus, Reference< QDPType<TVec,OLattice<TVec> > > >, 
	                    OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y += P{-}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().child());

  if( s.hasOrderedRep()) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
  
  
    int n_4vec = (s.end()-s.start()+1);
    add_g5ProjMinus(yptr, yptr, xptr, n_4vec);
  }
  else {
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
  
    
      add_g5ProjMinus(yptr, yptr, xptr, 1);
    }
  }

}


// d -= ChiralProjPlus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<
	                    UnaryNode< FnChiralProjectPlus, Reference< QDPType<TVec,OLattice<TVec> > > >, 
	                    OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y -= P{+}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().child());

  if( s.hasOrderedRep()) { 
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
    
    int n_4vec = (s.end()-s.start()+1);
    sub_g5ProjPlus(yptr, yptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i=tab[j];
      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
      
      sub_g5ProjPlus(yptr, yptr, xptr, 1); 
    }
  }
}


// d += ChiralProjMinus(Vec);
template<>
inline
void evaluate(OLattice< TVec >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<
	                    UnaryNode< FnChiralProjectMinus, Reference< QDPType<TVec,OLattice<TVec> > > >, 
	                    OLattice< TVec > > &rhs,
	      const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "y -= P{-}x" << endl;
#endif

  const OLattice< TVec >& x = static_cast<const OLattice< TVec > &>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
  
    REAL* xptr = (REAL *)&(x.elem(s.start()).elem(0).elem(0).real());
    REAL* yptr = (REAL *)&(d.elem(s.start()).elem(0).elem(0).real());
  
  
    int n_4vec = (s.end()-s.start()+1);
    sub_g5ProjMinus(yptr, yptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) {
      int i=tab[j];

      REAL* xptr = (REAL *)&(x.elem(i).elem(0).elem(0).real());
      REAL* yptr = (REAL *)&(d.elem(i).elem(0).elem(0).real());
      
      sub_g5ProjMinus(yptr, yptr, xptr, 1);
    }

  }
  
}

// d = x + a P_{+} y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	                  Reference< QDPType< TVec, OLattice< TVec > > >,
	                  BinaryNode<OpMultiply, 
	                             Reference< QDPType< TScal, OScalar< TScal > > >,
	                  UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice<TVec> > > >
                                    > 
                          >,
	                  OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x + a*P{+} y" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,    
    UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode.right().child());
  // Set pointers 
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;


  if ( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xpayz_g5ProjPlus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      xpayz_g5ProjPlus(zptr, aptr, xptr, yptr, 1);

    }
  }

}

// d = x + a P_{-} y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpAdd,
	                  Reference< QDPType< TVec, OLattice< TVec > > >,
	                  BinaryNode<OpMultiply, 
	                             Reference< QDPType< TScal, OScalar< TScal > > >,
	                  UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice<TVec> > > >
                                    > 
                          >,
	                  OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x + a*P{-} y" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,    
    UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode.right().child());
  // Set pointers 
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xpayz_g5ProjMinus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      xpayz_g5ProjMinus(zptr, aptr, xptr, yptr, 1);

    }
  }

}

// d = x - a P_{+} y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	                  Reference< QDPType< TVec, OLattice< TVec > > >,
	                  BinaryNode<OpMultiply, 
	                             Reference< QDPType< TScal, OScalar< TScal > > >,
	                  UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice<TVec> > > >
                                    > 
                          >,
	                  OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x - a*P{+} y" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,    
    UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode.right().child());
  // Set pointers 
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if( s.hasOrderedRep()) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xmayz_g5ProjPlus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      xmayz_g5ProjPlus(zptr, aptr, xptr, yptr, 1);

    }
  }

}

// d = x - a P_{-} y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpSubtract,
	                  Reference< QDPType< TVec, OLattice< TVec > > >,
	                  BinaryNode<OpMultiply, 
	                             Reference< QDPType< TScal, OScalar< TScal > > >,
	                  UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice<TVec> > > >
                                    > 
                          >,
	                  OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x - a*P{-} y" << endl;
#endif


  // Peel the stuff out of the expression

  // y is the left side of rhs
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&> (rhs.expression().left());

  // ax is the right side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,    
    UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN;

  // get the binary node
  const BN &mulNode = static_cast<const BN&> (rhs.expression().right());

  // get a and x out of the bynary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode.right().child());
  // Set pointers 
  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;

  if ( s.hasOrderedRep() ) { 

    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xmayz_g5ProjMinus(zptr, aptr, xptr, yptr, n_4vec);

  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      xmayz_g5ProjMinus(zptr, aptr, xptr, yptr, 1);

    }
  }

}

// d = ax + P+ y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpAdd,
	           BinaryNode<OpMultiply, 
	              Reference< QDPType< TScal, OScalar< TScal > > >,
	              Reference< QDPType< TVec, OLattice< TVec  > > > 
                   >,
	           UnaryNode<FnChiralProjectPlus, Reference< QDPType<TVec, OLattice<TVec> > > >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + P{+} y" << endl;
#endif


  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right().child());

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
    int n_4vec = (s.end()-s.start()+1);
    axpyz_g5ProjPlus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axpyz_g5ProjPlus(zptr, aptr, xptr, yptr, 1);

    }
  }

}

// d = ax + P- y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpAdd,
	           BinaryNode<OpMultiply, 
	              Reference< QDPType< TScal, OScalar< TScal > > >,
	              Reference< QDPType< TVec, OLattice< TVec  > > > 
                   >,
	           UnaryNode<FnChiralProjectMinus, Reference< QDPType<TVec, OLattice<TVec> > > >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + P{-} y" << endl;
#endif


  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right().child());

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
    int n_4vec = (s.end()-s.start()+1);
    axpyz_g5ProjMinus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axpyz_g5ProjMinus(zptr, aptr, xptr, yptr, 1);

    }
  }
}


// d = ax - P+ y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpSubtract,
	           BinaryNode<OpMultiply, 
	              Reference< QDPType< TScal, OScalar< TScal > > >,
	              Reference< QDPType< TVec, OLattice< TVec  > > > 
                   >,
	           UnaryNode<FnChiralProjectPlus, Reference< QDPType<TVec, OLattice<TVec> > > >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + P{+} y" << endl;
#endif


  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right().child());

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

  if( s.hasOrderedRep()) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axmyz_g5ProjPlus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axmyz_g5ProjPlus(zptr, aptr, xptr, yptr, 1);

    }
  }
}

// d = ax - P- y
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpSubtract,
	           BinaryNode<OpMultiply, 
	              Reference< QDPType< TScal, OScalar< TScal > > >,
	              Reference< QDPType< TVec, OLattice< TVec  > > > 
                   >,
	           UnaryNode<FnChiralProjectMinus, Reference< QDPType<TVec, OLattice<TVec> > > >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{
#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + P{-} y" << endl;
#endif


  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&> (rhs.expression().right().child());

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
    int n_4vec = (s.end()-s.start()+1);
    axmyz_g5ProjMinus(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axmyz_g5ProjMinus(zptr, aptr, xptr, yptr, 1);

    }
  }
}

// Vec = Scal * P_{+} Vec
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpMultiply,
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > >
               >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  cout << "BJ: v = a*P{+}v " << endl;
#endif

  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(rhs.expression().right().child());
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().left());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = &ar;  

  if( s.hasOrderedRep()) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    int n_4vec = (s.end()-s.start()+1);
    
    scal_g5ProjPlus(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      scal_g5ProjPlus(zptr, aptr, xptr, 1);
    }
  }

}

// Vec = Scal * P_{-} Vec
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	       BinaryNode<OpMultiply,
	       Reference< QDPType< TScal, OScalar< TScal > > >,
	       UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > >
               >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  cout << "BJ: v = a*P{-}v " << endl;
#endif

  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(rhs.expression().right().child());
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().left());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = &ar;  
  
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    int n_4vec = (s.end()-s.start()+1);
    
    scal_g5ProjMinus(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());

      scal_g5ProjMinus(zptr, aptr, xptr, 1);
    }
  }

}

// z = ax + bP+ y
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
	         UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > >
                > 
               >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + b*P+y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN2;

  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right().child());

  
  // Set pointers 
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if ( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axpbyz_g5ProjPlus(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      axpbyz_g5ProjPlus(zptr, aptr, xptr, bptr, yptr, 1);
    }
  }

}

// z = ax + bP- y
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
	         UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > >
                > 
               >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + b*P-y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN2;

  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right().child());

  
  // Set pointers 
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());

  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axpbyz_g5ProjMinus(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      axpbyz_g5ProjMinus(zptr, aptr, xptr, bptr, yptr, 1);
    }
  }

}

// z = ax - bP+ y
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
	         UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > >
                > 
               >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x - b*P+y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    UnaryNode< FnChiralProjectPlus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN2;

  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right().child());

  
  // Set pointers 
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep()) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axmbyz_g5ProjPlus(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      axmbyz_g5ProjPlus(zptr, aptr, xptr, bptr, yptr, 1);
    }
  }

}

// z = ax - bP- y
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
	         UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > >
                > 
               >,
	       OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x - b*P-y" << endl;
#endif

  // Peel the stuff out of the expression
  // y is the right side of rhs

  // ax is the left side of rhs and is in a binary node
  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    Reference< QDPType< TVec, OLattice< TVec > > > > BN1;

  typedef BinaryNode<OpMultiply, 
    Reference< QDPType< TScal, OScalar< TScal > > >,
    UnaryNode< FnChiralProjectMinus, Reference< QDPType< TVec, OLattice< TVec > > > > > BN2;

  // get the binary node
  const BN1 &mulNode1 = static_cast<const BN1&> (rhs.expression().left());
  const BN2 &mulNode2 = static_cast<const BN2&> (rhs.expression().right());

  // get a and x out of the binary node
  const OScalar< TScal >& a = static_cast<const OScalar< TScal >&>(mulNode1.left());
  const OLattice< TVec >& x = static_cast<const OLattice< TVec >&>(mulNode1.right());
  
  // get b and y out of the binary node
  const OScalar< TScal >& b = static_cast<const OScalar< TScal >&>(mulNode2.left());
  const OLattice< TVec >& y = static_cast<const OLattice< TVec >&>(mulNode2.right().child());

  
  // Set pointers 
  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());
  if ( s.hasOrderedRep() ) {
    
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axmbyz_g5ProjMinus(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
    
    
      axmbyz_g5ProjMinus(zptr, aptr, xptr, bptr, yptr, 1);
    }
  }

}

// Vec = Scal * GammaConst<Ns,Ns-1>* Vec
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	          BinaryNode<OpMultiply,
	          Reference< QDPType< TScal, OScalar< TScal > > >,
	            BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                Reference< QDPType< TVec, OLattice< TVec > > >
                    >	        
                  >,
	          OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*(GammaConst<Ns,Ns*Ns-1>()*x)" << endl;
#endif


   typedef BinaryNode< 
            OpGammaConstMultiply,
            GammaConst<Ns,Ns*Ns-1>,
	    Reference< QDPType< TVec, OLattice< TVec > > >
    > BN1;
        
  const OScalar< TScal > &a = static_cast<const OScalar< TScal >&>(rhs.expression().left());
  const BN1 &node = static_cast<const BN1&>(rhs.expression().right());
  
  const OLattice< TVec > &x = static_cast<const OLattice< TVec >&>(node.right());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = &ar;  
  if ( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    int n_4vec = (s.end()-s.start()+1);
   
    scal_g5(zptr, aptr, xptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      scal_g5(zptr, aptr, xptr, 1);
    }
  }
 
}


// Vec = Vec - a*Gamma5*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpSubtract,
	           Reference< QDPType< TVec, OLattice< TVec > > >, 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                Reference< QDPType< TVec, OLattice< TVec > > >
                     >	        
                   >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x - a GammaConst<Ns,Ns*Ns-1>()*y" << endl;
#endif
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(rhs.expression().left());

  typedef BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                Reference< QDPType< TVec, OLattice< TVec > > >
                     >	        
    > MN;
  const MN& mul_node = static_cast<const MN&>(rhs.expression().right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>& >(mul_node.left());

  typedef BinaryNode< 
    OpGammaConstMultiply,
    GammaConst<Ns,Ns*Ns-1>,
    Reference< QDPType< TVec, OLattice< TVec > > >
    > GN;

  const GN& gamma_node = static_cast<const GN&>(mul_node.right());
  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(gamma_node.right());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xmayz_g5(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      xmayz_g5(zptr, aptr, xptr, yptr, 1);
      
    }
  }

}


// Vec = a*Vec + b*Gamma5*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpAdd,
	           BinaryNode<OpMultiply,
	             Reference< QDPType<TScal, OScalar< TScal > > >,  
	             Reference< QDPType<TVec, OLattice< TVec > > >
	           >, 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                Reference< QDPType< TVec, OLattice< TVec > > >
                     >	        
                   >
                 >,
	         OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + b*GammaConst<Ns,Ns*Ns-1>()*y" << endl;
#endif

  typedef BinaryNode<OpMultiply,
            Reference< QDPType<TScal, OScalar< TScal > > >,  
            Reference< QDPType<TVec, OLattice< TVec > > >
    > MN1;

  typedef BinaryNode<OpMultiply,
            Reference< QDPType< TScal, OScalar< TScal > > >,
	      BinaryNode< 
                OpGammaConstMultiply,
 	        GammaConst<Ns,Ns*Ns-1>,
	        Reference< QDPType< TVec, OLattice< TVec > > >
            >	        
    > MN2;

  typedef BinaryNode< 
            OpGammaConstMultiply,
            GammaConst<Ns,Ns*Ns-1>,
            Reference< QDPType< TVec, OLattice< TVec > > >
    > GN;

  const MN1& mulNode1 = static_cast< const MN1& >(rhs.expression().left());
  const MN2& mulNode2 = static_cast< const MN2& >(rhs.expression().right());
  const GN& gammaNode = static_cast< const GN& >(mulNode2.right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>&>(mulNode1.left());
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(mulNode1.right());

  const OScalar<TScal>& b = static_cast<const OScalar<TScal>&>(mulNode2.left());
  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(gammaNode.right());


  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axpbyz_g5(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL* zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axpbyz_g5(zptr, aptr, xptr, bptr, yptr, 1);
      
    }
  }
  
}

// Vec = Gamma_5 *( a*Vec - b*Vec )
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpGammaConstMultiply,
	                    GammaConst<Ns,Ns*Ns-1>, 
	                    BinaryNode<OpSubtract,
	                      BinaryNode<OpMultiply,
	                        Reference< QDPType<TScal, OScalar< TScal > > >,  
	                        Reference< QDPType<TVec, OLattice< TVec > > >
	                      >, 
 	                      BinaryNode<OpMultiply,
	                        Reference< QDPType< TScal, OScalar< TScal > > >,
	                        Reference< QDPType< TVec, OLattice< TVec > > >
                              >	        
                            >
                 >,
	         OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = GammaConst<Ns,Ns*Ns-1>()*(ax - by)" << endl;
#endif

  typedef BinaryNode<OpSubtract,
	             BinaryNode<OpMultiply,
	               Reference< QDPType<TScal, OScalar< TScal > > >,  
	               Reference< QDPType<TVec, OLattice< TVec > > >
	             >, 
 	             BinaryNode<OpMultiply,
	               Reference< QDPType< TScal, OScalar< TScal > > >,
	               Reference< QDPType< TVec, OLattice< TVec > > >
                     >	        
    > AXMBY;

  const AXMBY& axmby_node = static_cast<const AXMBY&>(rhs.expression().right());


  typedef BinaryNode<OpMultiply,
            Reference< QDPType<TScal, OScalar< TScal > > >,  
            Reference< QDPType<TVec, OLattice< TVec > > >
    > MN;

  const MN& mulNode1 = static_cast<const MN& >( axmby_node.left());
  const MN& mulNode2 = static_cast<const MN& >( axmby_node.right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>&>(mulNode1.left());
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(mulNode1.right());

  const OScalar<TScal>& b = static_cast<const OScalar<TScal>&>(mulNode2.left());
  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(mulNode2.right());


  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    g5_axmbyz(zptr, aptr, xptr, bptr, yptr, n_4vec); 
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      g5_axmbyz(zptr, aptr, xptr, bptr, yptr, 1); 
      
    }
  }

 
}


// Vec = a*Vec + b*Gamma5*i*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpAdd,
	           BinaryNode<OpMultiply,
	             Reference< QDPType<TScal, OScalar< TScal > > >,  
	             Reference< QDPType<TVec, OLattice< TVec > > >
	           >, 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType<TVec, OLattice< TVec > > >
	                >
                     >	        
                   >
                 >,
	         OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + b*GammaConst<Ns,Ns*Ns-1>()*timesI(y)" << endl;
#endif

  typedef BinaryNode<OpMultiply,
            Reference< QDPType<TScal, OScalar< TScal > > >,  
            Reference< QDPType<TVec, OLattice< TVec > > >
    > MN1;

  typedef BinaryNode<OpMultiply,
            Reference< QDPType< TScal, OScalar< TScal > > >,
	      BinaryNode< 
                OpGammaConstMultiply,
 	        GammaConst<Ns,Ns*Ns-1>,
	        UnaryNode< FnTimesI, 
                           Reference< QDPType< TVec, OLattice< TVec > > >
                >
            >	        
    > MN2;

  typedef BinaryNode< 
            OpGammaConstMultiply,
            GammaConst<Ns,Ns*Ns-1>,
            UnaryNode< FnTimesI,
                       Reference< QDPType< TVec, OLattice< TVec > > >
            >
    > GN;

  typedef UnaryNode< FnTimesI, 
                     Reference< QDPType< TVec, OLattice< TVec > > > 
    > IN;

  const MN1& mulNode1 = static_cast< const MN1& >(rhs.expression().left());
  const MN2& mulNode2 = static_cast< const MN2& >(rhs.expression().right());
  const GN& gammaNode = static_cast< const GN& >(mulNode2.right());
  const IN& mulINode = static_cast< const IN& >(gammaNode.right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>&>(mulNode1.left());
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(mulNode1.right());

  const OScalar<TScal>& b = static_cast<const OScalar<TScal>&>(mulNode2.left());
  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(mulINode.child());


  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axpbyz_ig5(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axpbyz_ig5(zptr, aptr, xptr, bptr, yptr, 1);
      
    }
  }

}

// Vec = a*Vec - b*Gamma5*i*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpSubtract,
	           BinaryNode<OpMultiply,
	             Reference< QDPType<TScal, OScalar< TScal > > >,  
	             Reference< QDPType<TVec, OLattice< TVec > > >
	           >, 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType<TVec, OLattice< TVec > > >
	                >
                     >	        
                   >
                 >,
	         OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = a*x + b*GammaConst<Ns,Ns*Ns-1>()*timesI(y)" << endl;
#endif

  typedef BinaryNode<OpMultiply,
            Reference< QDPType<TScal, OScalar< TScal > > >,  
            Reference< QDPType<TVec, OLattice< TVec > > >
    > MN1;

  typedef BinaryNode<OpMultiply,
            Reference< QDPType< TScal, OScalar< TScal > > >,
	      BinaryNode< 
                OpGammaConstMultiply,
 	        GammaConst<Ns,Ns*Ns-1>,
	        UnaryNode< FnTimesI, 
                           Reference< QDPType< TVec, OLattice< TVec > > >
                >
            >	        
    > MN2;

  typedef BinaryNode< 
            OpGammaConstMultiply,
            GammaConst<Ns,Ns*Ns-1>,
            UnaryNode< FnTimesI,
                       Reference< QDPType< TVec, OLattice< TVec > > >
            >
    > GN;

  typedef UnaryNode< FnTimesI, 
                     Reference< QDPType< TVec, OLattice< TVec > > > 
    > IN;

  const MN1& mulNode1 = static_cast< const MN1& >(rhs.expression().left());
  const MN2& mulNode2 = static_cast< const MN2& >(rhs.expression().right());
  const GN& gammaNode = static_cast< const GN& >(mulNode2.right());
  const IN& mulINode = static_cast< const IN& >(gammaNode.right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>&>(mulNode1.left());
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(mulNode1.right());

  const OScalar<TScal>& b = static_cast<const OScalar<TScal>&>(mulNode2.left());
  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(mulINode.child());


  REAL *aptr = (REAL *)&(a.elem().elem().elem().elem());
  REAL *bptr = (REAL *)&(b.elem().elem().elem().elem());
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());
    
  
    
    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    axmbyz_ig5(zptr, aptr, xptr, bptr, yptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      axmbyz_ig5(zptr, aptr, xptr, bptr, yptr, 1);
      
    }
  }

}

// Vec = Vec + a*Gamma5*i*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpAdd,
	           Reference< QDPType< TVec, OLattice< TVec > > >, 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
                   >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x + a GammaConst<Ns,Ns*Ns-1>()*i*y" << endl;
#endif
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(rhs.expression().left());

  typedef BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
    	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
    > MN;
  const MN& mul_node = static_cast<const MN&>(rhs.expression().right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>& >(mul_node.left());

  typedef BinaryNode< 
    OpGammaConstMultiply,
    GammaConst<Ns,Ns*Ns-1>,
    UnaryNode<FnTimesI, 
              Reference< QDPType< TVec, OLattice< TVec > > >
    >
   > GN;

  typedef UnaryNode<FnTimesI,
              Reference< QDPType< TVec, OLattice< TVec > > >
    > IN;

  const GN& gamma_node = static_cast<const GN&>(mul_node.right());
  const IN& timesI_node = static_cast<const IN&>(gamma_node.right());

  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(timesI_node.child());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());

    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xpayz_ig5(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      xpayz_ig5(zptr, aptr, xptr, yptr, 1);
      
    }
  }

}

// Vec = Vec - a*Gamma5*i*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAssign &op,
	       const QDPExpr< 
	         BinaryNode<OpSubtract,
	           Reference< QDPType< TVec, OLattice< TVec > > >, 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
                   >
                 >,
	        OLattice< TVec > > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z = x - a GammaConst<Ns,Ns*Ns-1>()*i*y" << endl;
#endif
  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(rhs.expression().left());

  typedef BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
    > MN;
  const MN& mul_node = static_cast<const MN&>(rhs.expression().right());

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>& >(mul_node.left());

  typedef BinaryNode< 
    OpGammaConstMultiply,
    GammaConst<Ns,Ns*Ns-1>,
    UnaryNode<FnTimesI, 
              Reference< QDPType< TVec, OLattice< TVec > > >
    >
   > GN;

  typedef UnaryNode<FnTimesI,
              Reference< QDPType< TVec, OLattice< TVec > > >
    > IN;

  const GN& gamma_node = static_cast<const GN&>(mul_node.right());
  const IN& timesI_node = static_cast<const IN&>(gamma_node.right());

  const OLattice<TVec>& y = static_cast<const OLattice<TVec>&>(timesI_node.child());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;
  if( s.hasOrderedRep()) {

    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL *yptr = (REAL *) &(y.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xmayz_ig5(zptr, aptr, xptr, yptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *yptr = (REAL *) &(y.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      
      // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
      xmayz_ig5(zptr, aptr, xptr, yptr, 1);
      
    }
  }

}



// Vec += a*Gamma5*i*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpAddAssign &op,
	       const QDPExpr< 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
	           >,
	           OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z += a GammaConst<Ns,Ns*Ns-1>()*i*y" << endl;
#endif

#if 0
  typedef BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
    	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
    > MN;

  const MN& mul_node = static_cast<const MN&>(rhs.expression().right());
#endif

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>& >(rhs.expression().left());

  typedef BinaryNode< 
    OpGammaConstMultiply,
    GammaConst<Ns,Ns*Ns-1>,
    UnaryNode<FnTimesI, 
              Reference< QDPType< TVec, OLattice< TVec > > >
    >
   > GN;

  typedef UnaryNode<FnTimesI,
              Reference< QDPType< TVec, OLattice< TVec > > >
    > IN;

  const GN& gamma_node = static_cast<const GN&>(rhs.expression().right());
  const IN& timesI_node = static_cast<const IN&>(gamma_node.right());

  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(timesI_node.child());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;
  if( s.hasOrderedRep() ) {
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xpayz_ig5(zptr, aptr, zptr, xptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      xpayz_ig5(zptr, aptr, zptr, xptr, 1);
      
    }
  }

}


// Vec -= a*Gamma5*i*Vec
//
template<>
inline
void evaluate( OLattice< TVec > &d,
	       const OpSubtractAssign &op,
	       const QDPExpr< 
 	           BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
	           >,
	           OLattice< TVec > 
               > &rhs,
	       const Subset& s)
{

#ifdef DEBUG_BLAS_G5
  QDPIO::cout << "z -= a GammaConst<Ns,Ns*Ns-1>()*i*y" << endl;
#endif

#if 0 
  typedef BinaryNode<OpMultiply,
	             Reference< QDPType< TScal, OScalar< TScal > > >,
	             BinaryNode< 
	                OpGammaConstMultiply,
 	                GammaConst<Ns,Ns*Ns-1>,
    	                UnaryNode<FnTimesI,
	                          Reference< QDPType< TVec, OLattice< TVec > > >
	                >
                     >	        
    > MN;

  const MN& mul_node = static_cast<const MN&>(rhs.expression().right());
#endif

  const OScalar<TScal>& a = static_cast<const OScalar<TScal>& >(rhs.expression().left());

  typedef BinaryNode< 
    OpGammaConstMultiply,
    GammaConst<Ns,Ns*Ns-1>,
    UnaryNode<FnTimesI, 
              Reference< QDPType< TVec, OLattice< TVec > > >
    >
   > GN;

  typedef UnaryNode<FnTimesI,
              Reference< QDPType< TVec, OLattice< TVec > > >
    > IN;

  const GN& gamma_node = static_cast<const GN&>(rhs.expression().right());
  const IN& timesI_node = static_cast<const IN&>(gamma_node.right());

  const OLattice<TVec>& x = static_cast<const OLattice<TVec>&>(timesI_node.child());

  REAL ar =  a.elem().elem().elem().elem();
  REAL *aptr = (REAL *)&ar;
  if( s.hasOrderedRep() ) { 
    REAL *xptr = (REAL *) &(x.elem(s.start()).elem(0).elem(0).real());
    REAL* zptr =  &(d.elem(s.start()).elem(0).elem(0).real());


    // Get the no of 3vecs. s.start() and s.end() are inclusive so add +1
    int n_4vec = (s.end()-s.start()+1);
    xmayz_ig5(zptr, aptr, zptr, xptr, n_4vec);
  }
  else { 
    const int* tab=s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      REAL *xptr = (REAL *) &(x.elem(i).elem(0).elem(0).real());
      REAL *zptr =  &(d.elem(i).elem(0).elem(0).real());
      
      xmayz_ig5(zptr, aptr, zptr, xptr, 1);
      
    }
  }

}


} // namespace QDP;

#endif  // guard
 
