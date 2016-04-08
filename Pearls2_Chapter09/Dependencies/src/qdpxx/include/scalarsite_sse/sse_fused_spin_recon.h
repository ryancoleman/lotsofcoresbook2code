#ifndef SSE_FUSED_SPIN_RECON_H
#define SSE_FUSED_SPIN_RECON_H

#include "sse_mult_su3_mat_hwvec.h"

namespace QDP {
 
// Convenience Types
typedef PColorVector<RComplex<REAL32>, 3> ColVec32;
typedef PSpinVector< ColVec32, 4> Spin4_32;
typedef PSpinVector< ColVec32, 2> Spin2_32;
typedef PColorMatrix<RComplex<REAL32>, 3> ColMat32;


// Tell the compiler what the result of the multiply 
// half spinor with gauge field operations are
// OScalar
template<>
struct BinaryReturn< QDPType< PScalar<ColMat32>, OScalar< PScalar<ColMat32> > >,
		     QDPType< Spin2_32,           OScalar< Spin2_32 > >,
		     OpMultiply > 
{
  typedef  OScalar<Spin2_32> Type_t; 
};

// Tell the compiler what the result of opreations are
// For OLattice
template<>
struct BinaryReturn< QDPType< PScalar<ColMat32>, OLattice< PScalar<ColMat32> > >,
		     QDPType< Spin2_32,           OLattice< Spin2_32 > >,
		     OpMultiply >
{
  typedef  OLattice<Spin2_32> Type_t;
};

/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir0Minus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir0MinusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir0MinusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir0MinusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir0Minus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir0MinusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir0Minus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir0MinusProd>::Type_t
sreconDir0Minus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir0MinusProd>::Type_t ret;
  ret = spinReconstructDir0Minus(tmp);
  return ret;
}

// This rewrites spinReconstructDir0Minus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir0MinusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir0MinusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir0MinusProd>::Type_t >::Expression_t
spinReconstructDir0Minus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir0MinusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;


  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir0MinusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}

/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir0Plus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir0PlusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir0PlusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir0PlusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir0Plus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir0PlusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir0Plus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir0PlusProd>::Type_t
sreconDir0Plus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir0PlusProd>::Type_t ret;
  ret = spinReconstructDir0Plus(tmp);
  return ret;
}

// This rewrites spinReconstructDir0Plus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir0PlusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir0PlusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir0PlusProd>::Type_t >::Expression_t
spinReconstructDir0Plus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir0PlusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;


  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir0PlusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}

/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir1Minus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir1MinusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir1MinusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir1MinusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir1Minus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir1MinusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir1Minus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir1MinusProd>::Type_t
sreconDir1Minus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir1MinusProd>::Type_t ret;
  ret = spinReconstructDir1Minus(tmp);
  return ret;
}

// This rewrites spinReconstructDir1Minus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir1MinusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir1MinusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir1MinusProd>::Type_t >::Expression_t
spinReconstructDir1Minus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir1MinusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;


  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir1MinusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}

/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir1Plus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir1PlusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir1PlusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir1PlusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir1Plus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir1PlusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir1Plus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir1PlusProd>::Type_t
sreconDir1Plus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir1PlusProd>::Type_t ret;
  ret = spinReconstructDir1Plus(tmp);
  return ret;
}

// This rewrites spinReconstructDir1Plus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir1PlusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir1PlusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir1PlusProd>::Type_t >::Expression_t
spinReconstructDir1Plus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir1PlusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;


  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir1PlusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}


/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir2Minus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir2MinusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir2MinusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir2MinusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir2Minus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir2MinusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir2Minus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir2MinusProd>::Type_t
sreconDir2Minus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir2MinusProd>::Type_t ret;
  ret = spinReconstructDir2Minus(tmp);
  return ret;
}

// This rewrites spinReconstructDir2Minus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir2MinusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir2MinusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir2MinusProd>::Type_t >::Expression_t
spinReconstructDir2Minus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir2MinusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;


  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir2MinusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}

/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir2Plus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir2PlusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir2PlusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir2PlusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir2Plus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir2PlusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir2Plus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir2PlusProd>::Type_t
sreconDir2Plus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir2PlusProd>::Type_t ret;
  ret = spinReconstructDir2Plus(tmp);
  return ret;
}

// This rewrites spinReconstructDir2Plus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir2PlusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir2PlusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir2PlusProd>::Type_t >::Expression_t
spinReconstructDir2Plus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir2PlusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;


  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir2PlusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}


/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir3Minus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir3MinusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir3MinusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir3MinusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir3Minus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir3MinusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir3Minus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir3MinusProd>::Type_t
sreconDir3Minus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir3MinusProd>::Type_t ret;
  ret = spinReconstructDir3Minus(tmp);
  return ret;
}

// This rewrites spinReconstructDir3Minus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir3MinusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir3MinusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir3MinusProd>::Type_t >::Expression_t
spinReconstructDir3Minus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir3MinusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;

  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir3MinusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}

/////////////////////////////////////////////////////////////////////////
/// Reconstruct Dir3Plus Fused 
/////////////////////////////////////////////////////////////////////////
struct FnSReconDir3PlusProd 
{
  PETE_EMPTY_CONSTRUCTORS(FnSReconDir3PlusProd)
  
  template<class T1, class T2, int N>
  inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir3PlusProd>::Type_t
  operator()(const PScalar<T1>& a, const PSpinVector<T2,N>& b) const 
  {
    return sreconDir3Plus(a, b);
  }
};


// Proclaim that the return type of the fused op is the same as the return type
// of the SpinReconstruction on a normal half vector
template<class T1, class T2, int N>
struct BinaryReturn< PScalar<T1>, PSpinVector<T2, N>, FnSReconDir3PlusProd>
{
  typedef typename UnaryReturn< PSpinVector<T2,N>, FnSpinReconstructDir3Plus >::Type_t  Type_t;
};

// For Generic Subtypes of OLattice<> and OScalar<>
// Do the generic thing
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N> , FnSReconDir3PlusProd>::Type_t
sreconDir3Plus(const PScalar<T1>& a, const PSpinVector<T2,N>& b) {
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t tmp ;
  tmp=a*b;
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, FnSReconDir3PlusProd>::Type_t ret;
  ret = spinReconstructDir3Plus(tmp);
  return ret;
}

// This rewrites spinReconstructDir3Plus( Matrix * Vec ) 
// as a fused reconstruct: FnSReconDir3PlusProd< Matrix, Vec >
template<class T1, class T2, template<class> class C, int N  >
inline typename MakeReturn<BinaryNode< FnSReconDir3PlusProd,
				       typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
				       typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>,
				       typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir3PlusProd>::Type_t >::Expression_t
spinReconstructDir3Plus(const QDPExpr<
                                  BinaryNode< OpMultiply, 
		                      Reference< QDPType< PScalar<T1> , C<PScalar<T1> > > >, 
			              Reference< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >
		               >, 
                               typename BinaryReturn< QDPType< PScalar<T1> , C<PScalar<T1> > >,
		                                      QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > >,
		                                      OpMultiply >::Type_t >& l)
{

  
  typedef BinaryNode< FnSReconDir3PlusProd,
    typename CreateLeaf< QDPType< PScalar<T1> , C<PScalar<T1> > > >::Leaf_t,
    typename CreateLeaf< QDPType< PSpinVector<T2,N>, C<PSpinVector<T2,N> > > >::Leaf_t>  Tree_t;

  return MakeReturn<Tree_t,typename BinaryReturn< C<PScalar<T1> >, C<PSpinVector<T2,N> >, FnSReconDir3PlusProd>::Type_t>::make(
	      Tree_t(
                CreateLeaf<QDPType<PScalar<T1> ,C<PScalar<T1> > > >::make(l.expression().left()),
                CreateLeaf<QDPType<PSpinVector<T2,N>,C<PSpinVector<T2,N> > > >::make(l.expression().right())
              )
         );
    

}



//////////////////////////////////////////////////////////////////////////
/// The Specialist Ops themesleves
//////////////////////////////////////////////////////////////////////////

// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t
sreconDir0Minus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir0Minus(&(d.elem(0).elem(0).real()),
			  &(ret.elem(0).elem(0).real()),
			  1);

  return ret;
}


// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0PlusProd>::Type_t
sreconDir0Plus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir0Plus(&(d.elem(0).elem(0).real()),
			 &(ret.elem(0).elem(0).real()),
			 1);


  return ret;
}

// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir1MinusProd>::Type_t
sreconDir1Minus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir1Minus(&(d.elem(0).elem(0).real()),
			 &(ret.elem(0).elem(0).real()),
			 1);

  return ret;
}

// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir1PlusProd>::Type_t
sreconDir1Plus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);


  inlineSpinReconDir1Plus(&(d.elem(0).elem(0).real()),
			 &(ret.elem(0).elem(0).real()),
			 1);

  return ret;
}

// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir2MinusProd>::Type_t
sreconDir2Minus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir2Minus(&(d.elem(0).elem(0).real()),
			  &(ret.elem(0).elem(0).real()),
			  1);

  return ret;
}

// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir2PlusProd>::Type_t
sreconDir2Plus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir2Plus(&(d.elem(0).elem(0).real()),
			 &(ret.elem(0).elem(0).real()),
			 1);

  return ret;

}

// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir3MinusProd>::Type_t
sreconDir3Minus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0MinusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir3Minus(&(d.elem(0).elem(0).real()),
			 &(ret.elem(0).elem(0).real()),
			 1);

  return ret;
}


// The actual fused op - Replace this with specialist code
template<>
inline BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir3PlusProd>::Type_t
sreconDir3Plus(const PScalar<ColMat32>& a, const PSpinVector<ColVec32,2>& b) 
{
  PSpinVector<ColVec32,2> d ;
  BinaryReturn< PScalar<ColMat32>, PSpinVector<ColVec32,2>, FnSReconDir0PlusProd>::Type_t  ret;

  su3_matrixf *am = (su3_matrixf *)&( a.elem().elem(0,0).real());
  half_wilson_vectorf *bh = (half_wilson_vectorf *)&( b.elem(0).elem(0).real());
  half_wilson_vectorf *dh = (half_wilson_vectorf *)&( d.elem(0).elem(0).real());

  intrin_sse_mult_su3_mat_hwvec(am,bh,dh);

  inlineSpinReconDir3Plus(&(d.elem(0).elem(0).real()),
			 &(ret.elem(0).elem(0).real()),
			 1);

  return ret;
}


} // namespace QDP

#endif
