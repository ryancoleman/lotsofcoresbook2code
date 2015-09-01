// -*- C++ -*-

/*! \file
 * \brief Primitive Color Vector
 */


#ifndef QDP_PRIMCOLORVEC_H
#define QDP_PRIMCOLORVEC_H

namespace QDP {

//-------------------------------------------------------------------------------------
/*! \addtogroup primcolorvector Color vector primitive
 * \ingroup primvector
 *
 * Primitive type that transforms like a Color vector
 *
 * @{
 */

//! Primitive color Vector class
template <class T, int N> class PColorVector : public PVector<T, N, PColorVector>
{
public:
  //! PColorVector = PColorVector
  /*! Set equal to another PColorVector */
  template<class T1>
  inline
  PColorVector& operator=(const PColorVector<T1,N>& rhs) 
    {
      this->assign(rhs);
      return *this;
    }

};

/*! @} */  // end of group primcolorvector

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1, int N>
struct WordType<PColorVector<T1,N> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};



template<class T1, int N>
struct SinglePrecType< PColorVector<T1,N> >
{
  typedef PColorVector< typename SinglePrecType<T1>::Type_t, N> Type_t;
};


template<class T1, int N>
struct DoublePrecType< PColorVector<T1,N> >
{
  typedef PColorVector< typename DoublePrecType<T1>::Type_t, N> Type_t;
};


// Internally used scalars
template<class T, int N>
struct InternalScalar<PColorVector<T,N> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving other indices along
template<class T, int N>
struct PrimitiveScalar<PColorVector<T,N> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T, int N>
struct LatticeScalar<PColorVector<T,N> > {
  typedef PColorVector<typename LatticeScalar<T>::Type_t, N>  Type_t;
};

//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PColorVector) -> PColorVector
template<class T1, int N, class Op>
struct UnaryReturn<PColorVector<T1,N>, Op> {
  typedef PColorVector<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};
// Default binary(PScalar,PColorVector) -> PColorVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PScalar<T1>, PColorVector<T2,N>, Op> {
  typedef PColorVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PColorMatrix,PColorVector) -> PColorVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PColorMatrix<T1,N>, PColorVector<T2,N>, Op> {
  typedef PColorVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PColorVector,PScalar) -> PColorVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PColorVector<T1,N>, PScalar<T2>, Op> {
  typedef PColorVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PColorVector,PColorVector) -> PColorVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, Op> {
  typedef PColorVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<PScalar<T2>, OpCast<T1> > {
  typedef PScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif


// Assignment is different
template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, OpAssign > {
  typedef PColorVector<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, OpAddAssign > {
  typedef PColorVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, OpSubtractAssign > {
  typedef PColorVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PScalar<T2>, OpMultiplyAssign > {
  typedef PColorVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PScalar<T2>, OpDivideAssign > {
  typedef PColorVector<T1,N> &Type_t;
};
 

// ColorVector
template<class T, int N>
struct UnaryReturn<PColorVector<T,N>, FnNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PColorVector<T,N>, FnLocalNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnLocalInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnLocalInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};




//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

// Peeking and poking
//! Extract color vector components 
template<class T, int N>
struct UnaryReturn<PColorVector<T,N>, FnPeekColorVector > {
  typedef PScalar<typename UnaryReturn<T, FnPeekColorVector>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PColorVector<T,N>, FnPeekColorVector>::Type_t
peekColor(const PColorVector<T,N>& l, int row)
{
  typename UnaryReturn<PColorVector<T,N>, FnPeekColorVector>::Type_t  d;

  // Note, do not need to propagate down since the function is eaten at this level
  d.elem() = l.elem(row);
  return d;
}

//! Insert color vector components
template<class T1, class T2, int N>
inline PColorVector<T1,N>&
pokeColor(PColorVector<T1,N>& l, const PScalar<T2>& r, int row)
{
  // Note, do not need to propagate down since the function is eaten at this level
  l.elem(row) = r.elem();
  return l;
}


//-----------------------------------------------------------------------------
// Contraction for color vectors
// colorContract 
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, PColorVector<T3,N>, FnColorContract> {
  typedef PScalar<typename TrinaryReturn<T1, T2, T3, FnColorContract>::Type_t>  Type_t;
};

//! dest  = colorContract(Qvec1,Qvec2,Qvec3)
/*!
 * Performs:
 *  \f$dest = \sum_{i,j,k} \epsilon^{i,j,k} V1^{i} V2^{j} V3^{k}\f$
 *
 * This routine is completely unrolled for 3 colors
 */
template<class T1, class T2, class T3>
inline typename TrinaryReturn<PColorVector<T1,3>, PColorVector<T2,3>, PColorVector<T3,3>, FnColorContract>::Type_t
colorContract(const PColorVector<T1,3>& s1, const PColorVector<T2,3>& s2, const PColorVector<T3,3>& s3)
{
  typename TrinaryReturn<PColorVector<T1,3>, PColorVector<T2,3>, PColorVector<T3,3>, FnColorContract>::Type_t  d;

  // Permutations: +(0,1,2)+(1,2,0)+(2,0,1)-(1,0,2)-(0,2,1)-(2,1,0)

  // d = \epsilon^{i,j,k} V1^{i} V2^{j} V3^{k}
  d.elem() = (s1.elem(0)*s2.elem(1)
           -  s1.elem(1)*s2.elem(0))*s3.elem(2)
           + (s1.elem(1)*s2.elem(2)
           -  s1.elem(2)*s2.elem(1))*s3.elem(0)
           + (s1.elem(2)*s2.elem(0)
           -  s1.elem(0)*s2.elem(2))*s3.elem(1);

  return d;
}


//-----------------------------------------------------------------------------
// Contraction for color vectors
// colorContract 
template<class T1, class T2, int N>
struct BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnColorVectorContract> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnColorVectorContract>::Type_t>  Type_t;
};

//! dest  = colorVectorContract(Qvec1,Qvec2)
/*!
 * Performs:
 *  \f$dest = \sum_{i} V1^{i} V2^{i}\f$
 */
template<class T1, class T2, int N>
inline typename BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnColorVectorContract>::Type_t
colorVectorContract(const PColorVector<T1,N>& s1, const PColorVector<T2,N>& s2)
{
  typename BinaryReturn<PColorVector<T1,N>, PColorVector<T2,N>, FnColorVectorContract>::Type_t  d;

  // d = V1^{i} V2^{i}
  d.elem() = s1.elem(0)*s2.elem(0);
  for(int i=1; i < N; ++i)
    d.elem() += s1.elem(i)*s2.elem(i);

  return d;
}


//-----------------------------------------------------------------------------
// diquark color cross product   s1 X s2
//! Contraction for color vectors
template<class T1, class T2>
struct BinaryReturn<PColorVector<T1,3>, PColorVector<T2,3>, FnColorCrossProduct> {
  typedef PColorVector<typename BinaryReturn<T1, T2, FnColorCrossProduct>::Type_t, 3>  Type_t;
};

//! dest  = colorCrossProduct(Qvec1,Qvec2)
/*!
 * Performs:
 *  \f$dest^{i} = \sum_{j,k} \epsilon^{i,j,k} V1^{j} V2^{k}\f$
 *
 * This routine is completely unrolled for 3 colors
 */
template<class T1, class T2>
inline typename BinaryReturn<PColorVector<T1,3>, PColorVector<T2,3>, FnColorCrossProduct>::Type_t
colorCrossProduct(const PColorVector<T1,3>& s1, const PColorVector<T2,3>& s2)
{
  typename BinaryReturn<PColorVector<T1,3>, PColorVector<T2,3>, FnColorCrossProduct>::Type_t  d;
  
  d.elem(0) = s1.elem(1)*s2.elem(2) - s1.elem(2)*s2.elem(1);
  d.elem(1) = s1.elem(2)*s2.elem(0) - s1.elem(0)*s2.elem(2);
  d.elem(2) = s1.elem(0)*s2.elem(1) - s1.elem(1)*s2.elem(0);

 return d;
}



} // namespace QDP

#endif

