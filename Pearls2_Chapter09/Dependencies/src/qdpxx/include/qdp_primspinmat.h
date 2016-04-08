// -*- C++ -*-

/*! \file
 * \brief Primitive Spin Matrix
 */

#ifndef QDP_PRIMSPINMAT_H
#define QDP_PRIMSPINMAT_H

namespace QDP {


//-------------------------------------------------------------------------------------
/*! \addtogroup primspinmatrix Spin matrix primitive
 * \ingroup primmatrix
 *
 * Primitive type that transforms like a Spin Matrix
 *
 * @{
 */


//! Primitive Spin Matrix class
/*! 
   * Spin matrix class support gamma matrix algebra 
   *
   * NOTE: the class is mostly empty - it is the specialized versions below
   * that know for a fixed size how gamma matrices (constants) should act
   * on the spin vectors.
   */
template <class T, int N> class PSpinMatrix : public PMatrix<T, N, PSpinMatrix>
{
public:
  //! PSpinMatrix = PScalar
  /*! Fill with primitive scalar */
  template<class T1>
  inline
  PSpinMatrix& operator=(const PScalar<T1>& rhs)
    {
      this->assign(rhs);
      return *this;
    }

  //! PSpinMatrix = PSpinMatrix
  /*! Set equal to another PSpinMatrix */
  template<class T1>
  inline
  PSpinMatrix& operator=(const PSpinMatrix<T1,N>& rhs) 
    {
      this->assign(rhs);
      return *this;
    }

};

/*! @} */   // end of group primspinmatrix

#if 0
//! Specialization of primitive spin Matrix class for 4 spin components
/*! 
 * Spin matrix class support gamma matrix algebra for 4 spin components
 */
template <class T> class PSpinMatrix<T,4> : public PMatrix<T,4, PSpinMatrix>
{
public:
};
#endif


//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1, int N>
struct WordType<PSpinMatrix<T1,N> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

template<class T1, int N>
struct SinglePrecType<PSpinMatrix<T1, N> >
{
  typedef PSpinMatrix< typename SinglePrecType<T1>::Type_t , N > Type_t;
};

template<class T1, int N>
struct DoublePrecType<PSpinMatrix<T1, N> >
{
  typedef PSpinMatrix< typename DoublePrecType<T1>::Type_t , N > Type_t;
};

// Internally used scalars
template<class T, int N>
struct InternalScalar<PSpinMatrix<T,N> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive into a scalar leaving grid alone
template<class T, int N>
struct PrimitiveScalar<PSpinMatrix<T,N> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T, int N>
struct LatticeScalar<PSpinMatrix<T,N> > {
  typedef PSpinMatrix<typename LatticeScalar<T>::Type_t, N>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PSpinMatrix) -> PSpinMatrix
template<class T1, int N, class Op>
struct UnaryReturn<PSpinMatrix<T1,N>, Op> {
  typedef PSpinMatrix<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};

// Default binary(PScalar,PSpinMatrix) -> PSpinMatrix
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, Op> {
  typedef PSpinMatrix<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinMatrix,PSpinMatrix) -> PSpinMatrix
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, Op> {
  typedef PSpinMatrix<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinMatrix,PScalar) -> PSpinMatrix
template<class T1, int N, class T2, class Op>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, Op> {
  typedef PSpinMatrix<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<PSpinMatrix<T2,N>, OpCast<T1> > {
  typedef PScalar<typename UnaryReturn<T, OpCast>::Type_t, N>  Type_t;
//  typedef T1 Type_t;
};
#endif


// Assignment is different
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, OpAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, OpAddAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, OpSubtractAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, OpMultiplyAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, OpAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, OpAddAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, OpSubtractAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, OpMultiplyAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, OpDivideAssign > {
  typedef PSpinMatrix<T1,N> &Type_t;
};
 


// SpinMatrix
template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnTrace > {
  typedef PScalar<typename UnaryReturn<T, FnTrace>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnRealTrace > {
  typedef PScalar<typename UnaryReturn<T, FnRealTrace>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnImagTrace > {
  typedef PScalar<typename UnaryReturn<T, FnImagTrace>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnLocalNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnTraceMultiply> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnTraceMultiply> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnTraceMultiply> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnLocalInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnLocalInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnLocalInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnLocalInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnLocalInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnLocalInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};



// Gamma algebra
template<int m, class T2, int N, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, PSpinMatrix<T2,N>, OpGammaConstMultiply> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<PSpinMatrix<T2,N>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, PSpinMatrix<T2,N>, OpGammaTypeMultiply> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaType>
struct BinaryReturn<PSpinMatrix<T2,N>, GammaType<N>, OpMultiplyGammaType> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};


// Gamma algebra
template<int m, class T2, int N, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, PSpinMatrix<T2,N>, OpGammaConstDPMultiply> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<PSpinMatrix<T2,N>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, PSpinMatrix<T2,N>, OpGammaTypeDPMultiply> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<PSpinMatrix<T2,N>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef PSpinMatrix<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};




//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------
/*! \addtogroup primspinmatrix */
/*! @{ */

// SpinMatrix class primitive operations

// trace = traceSpin(source1)
/*! This only acts on spin indices and is diagonal in all other indices */
template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnTraceSpin > {
  typedef PScalar<typename UnaryReturn<T, FnTraceSpin>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PSpinMatrix<T,N>, FnTraceSpin>::Type_t
traceSpin(const PSpinMatrix<T,N>& s1)
{
  typename UnaryReturn<PSpinMatrix<T,N>, FnTraceSpin>::Type_t  d;

  // Since the spin index is eaten, do not need to pass on function by
  // calling trace(...) again
  d.elem() = s1.elem(0,0);
  for(int i=1; i < N; ++i)
    d.elem() += s1.elem(i,i);

  return d;
}

//! traceSpinMultiply(source1,source2)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnTraceSpinMultiply> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PSpinMatrix<T1,N>& l, const PSpinMatrix<T2,N>& r)
{
  typename BinaryReturn<PSpinMatrix<T1,N>, PSpinMatrix<T2,N>, FnTraceSpinMultiply>::Type_t  d;

  // The traceSpin is eaten here
  d.elem() = l.elem(0,0) * r.elem(0,0);
  for(int k=1; k < N; ++k)
    d.elem() += l.elem(0,k) * r.elem(k,0);

  for(int j=1; j < N; ++j)
    for(int k=0; k < N; ++k)
      d.elem() += l.elem(j,k) * r.elem(k,j);

  return d;
}

//! PScalar = traceSpinMultiply(PSpinMatrix,PScalar)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnTraceSpinMultiply> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PSpinMatrix<T1,N>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PSpinMatrix<T1,N>, PScalar<T2>, FnTraceSpinMultiply>::Type_t  d;

  // The traceSpin is eaten here
  d.elem() = l.elem(0,0) * r.elem();
  for(int k=1; k < N; ++k)
    d.elem() += l.elem(k,k) * r.elem();

  return d;
}

// PScalar = traceSpinMultiply(PScalar,PSpinMatrix)
template<class T1, class T2, int N>
struct BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnTraceSpinMultiply> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceSpinMultiply>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PScalar<T1>& l, const PSpinMatrix<T2,N>& r)
{
  typename BinaryReturn<PScalar<T1>, PSpinMatrix<T2,N>, FnTraceSpinMultiply>::Type_t  d;

  // The traceSpin is eaten here
  d.elem() = l.elem() * r.elem(0,0);
  for(int k=1; k < N; ++k)
    d.elem() += l.elem() * r.elem(k,k);

  return d;
}



/*! Specialise the return type */
template <class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnTransposeSpin > {
  typedef PSpinMatrix<typename UnaryReturn<T, FnTransposeSpin>::Type_t, N> Type_t;
};

//! PSpinMatrix = transposeSpin(PSpinMatrix) 
/*! t = transposeSpin(source1) - SpinMatrix specialization -- where the work is actually done */
template<class T, int N>
inline typename UnaryReturn<PSpinMatrix<T,N>, FnTransposeSpin >::Type_t
transposeSpin(const PSpinMatrix<T,N>& s1)
{
  typename UnaryReturn<PSpinMatrix<T,N>, FnTransposeSpin>::Type_t d;
 
  for(int i=0; i < N; i++) { 
    for(int j=0; j < N; j++) { 
      // Transpose, so flip indices
      d.elem(i,j) = s1.elem(j,i);
    }
  }
  return d;
}


//-----------------------------------------------
// OuterProduct must be handled specially for each color and spin
// The problem is the traits class - I have no way to say to PVector's
//  transform into a PMatrix but downcast the trait to a PColorMatrix or
//  PSpinMatrix

//! PSpinMatrix = outerProduct(PSpinVector, PSpinVector)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnOuterProduct> {
  typedef PSpinMatrix<typename BinaryReturn<T1, T2, FnOuterProduct>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnOuterProduct>::Type_t
outerProduct(const PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnOuterProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    for(int j=0; j < N; ++j)
      d.elem(i,j) = outerProduct(l.elem(i),r.elem(j));

  return d;
}


//-----------------------------------------------
// Optimization of traceSpin(outerProduct(PSpinVector, PSpinVector))

//! PScalar = traceSpinOuterProduct(PSpinVector, PSpinVector)
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnTraceSpinOuterProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnOuterProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnTraceSpinOuterProduct>::Type_t
traceSpinOuterProduct(const PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnTraceSpinOuterProduct>::Type_t  d;

  d.elem() = outerProduct(l.elem(0),r.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += outerProduct(l.elem(i),r.elem(i));

  return d;
}


//-----------------------------------------------
// Peeking and poking
//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T, int N>
struct UnaryReturn<PSpinMatrix<T,N>, FnPeekSpinMatrix > {
  typedef PScalar<typename UnaryReturn<T, FnPeekSpinMatrix>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PSpinMatrix<T,N>, FnPeekSpinMatrix>::Type_t
peekSpin(const PSpinMatrix<T,N>& l, int row, int col)
{
  typename UnaryReturn<PSpinMatrix<T,N>, FnPeekSpinMatrix>::Type_t  d;

  // Note, do not need to propagate down since the function is eaten at this level
  d.elem() = l.elem(row,col);
  return d;
}

//! Insert spin matrix components
template<class T1, class T2, int N>
inline PSpinMatrix<T1,N>&
pokeSpin(PSpinMatrix<T1,N>& l, const PScalar<T2>& r, int row, int col)
{
  // Note, do not need to propagate down since the function is eaten at this level
  l.elem(row,col) = r.elem();
  return l;
}



//-----------------------------------------------

// SpinMatrix<4> = Gamma<4,m> * SpinMatrix<4>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<GammaConst<4,0>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,0>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,0>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;
  
  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = r.elem(0,i);
    d.elem(1,i) = r.elem(1,i);
    d.elem(2,i) = r.elem(2,i);
    d.elem(3,i) = r.elem(3,i);
  }

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,1>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,1>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,1>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;
  
  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(3,i));
    d.elem(1,i) = timesI(r.elem(2,i));
    d.elem(2,i) = timesMinusI(r.elem(1,i));
    d.elem(3,i) = timesMinusI(r.elem(0,i));
  }

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,2>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,2>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,2>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(3,i);
    d.elem(1,i) = r.elem(2,i);
    d.elem(2,i) = r.elem(1,i);
    d.elem(3,i) = -r.elem(0,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,3>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,3>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,3>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesMinusI(r.elem(0,i));
    d.elem(1,i) = timesI(r.elem(1,i));
    d.elem(2,i) = timesMinusI(r.elem(2,i));
    d.elem(3,i) = timesI(r.elem(3,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,4>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,4>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,4>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(2,i));
    d.elem(1,i) = timesMinusI(r.elem(3,i));
    d.elem(2,i) = timesMinusI(r.elem(0,i));
    d.elem(3,i) = timesI(r.elem(1,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,5>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,5>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,5>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(1,i);
    d.elem(1,i) = r.elem(0,i);
    d.elem(2,i) = -r.elem(3,i);
    d.elem(3,i) = r.elem(2,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,6>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,6>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,6>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesMinusI(r.elem(1,i));
    d.elem(1,i) = timesMinusI(r.elem(0,i));
    d.elem(2,i) = timesMinusI(r.elem(3,i));
    d.elem(3,i) = timesMinusI(r.elem(2,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,7>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,7>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,7>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = r.elem(2,i);
    d.elem(1,i) = r.elem(3,i);
    d.elem(2,i) = -r.elem(0,i);
    d.elem(3,i) = -r.elem(1,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,8>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,8>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,8>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = r.elem(2,i);
    d.elem(1,i) = r.elem(3,i);
    d.elem(2,i) = r.elem(0,i);
    d.elem(3,i) = r.elem(1,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,9>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,9>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,9>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(1,i));
    d.elem(1,i) = timesI(r.elem(0,i));
    d.elem(2,i) = timesMinusI(r.elem(3,i));
    d.elem(3,i) = timesMinusI(r.elem(2,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,10>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,10>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,10>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(1,i);
    d.elem(1,i) = r.elem(0,i);
    d.elem(2,i) = r.elem(3,i);
    d.elem(3,i) = -r.elem(2,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,11>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,11>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,11>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesMinusI(r.elem(2,i));
    d.elem(1,i) = timesI(r.elem(3,i));
    d.elem(2,i) = timesMinusI(r.elem(0,i));
    d.elem(3,i) = timesI(r.elem(1,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,12>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,12>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,12>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(0,i));
    d.elem(1,i) = timesMinusI(r.elem(1,i));
    d.elem(2,i) = timesMinusI(r.elem(2,i));
    d.elem(3,i) = timesI(r.elem(3,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,13>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,13>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,13>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(3,i);
    d.elem(1,i) = r.elem(2,i);
    d.elem(2,i) = -r.elem(1,i);
    d.elem(3,i) = r.elem(0,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,14>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,14>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,14>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesMinusI(r.elem(3,i));
    d.elem(1,i) = timesMinusI(r.elem(2,i));
    d.elem(2,i) = timesMinusI(r.elem(1,i));
    d.elem(3,i) = timesMinusI(r.elem(0,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,15>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,15>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,15>, PSpinMatrix<T2,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = r.elem(0,i);
    d.elem(1,i) = r.elem(1,i);
    d.elem(2,i) = -r.elem(2,i);
    d.elem(3,i) = -r.elem(3,i);
  }
  
  return d;
}


// SpinMatrix<4> = SpinMatrix<4> * Gamma<4,m>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,0>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,0>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,0>, OpGammaConstMultiply>::Type_t  d; 

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,0);
    d.elem(i,1) =  l.elem(i,1);
    d.elem(i,2) =  l.elem(i,2);
    d.elem(i,3) =  l.elem(i,3);
  }
 
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,1>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,1>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,1>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesMinusI(l.elem(i,3));
    d.elem(i,1) = timesMinusI(l.elem(i,2));
    d.elem(i,2) = timesI(l.elem(i,1));
    d.elem(i,3) = timesI(l.elem(i,0));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,2>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,2>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,2>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = -l.elem(i,3);
    d.elem(i,1) =  l.elem(i,2);
    d.elem(i,2) =  l.elem(i,1);
    d.elem(i,3) = -l.elem(i,0);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,3>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,3>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,3>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesMinusI(l.elem(i,0));
    d.elem(i,1) = timesI(l.elem(i,1));
    d.elem(i,2) = timesMinusI(l.elem(i,2));
    d.elem(i,3) = timesI(l.elem(i,3));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,4>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,4>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,4>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesMinusI(l.elem(i,2));
    d.elem(i,1) = timesI(l.elem(i,3));
    d.elem(i,2) = timesI(l.elem(i,0));
    d.elem(i,3) = timesMinusI(l.elem(i,1));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,5>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,5>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,5>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,1);
    d.elem(i,1) = -l.elem(i,0);
    d.elem(i,2) =  l.elem(i,3);
    d.elem(i,3) = -l.elem(i,2);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,6>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,6>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,6>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesMinusI(l.elem(i,1));
    d.elem(i,1) = timesMinusI(l.elem(i,0));
    d.elem(i,2) = timesMinusI(l.elem(i,3));
    d.elem(i,3) = timesMinusI(l.elem(i,2));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,7>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,7>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,7>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = -l.elem(i,2);
    d.elem(i,1) = -l.elem(i,3);
    d.elem(i,2) =  l.elem(i,0);
    d.elem(i,3) =  l.elem(i,1);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,8>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,8>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,8>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,2);
    d.elem(i,1) =  l.elem(i,3);
    d.elem(i,2) =  l.elem(i,0);
    d.elem(i,3) =  l.elem(i,1);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,9>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,9>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,9>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,1));
    d.elem(i,1) = timesI(l.elem(i,0));
    d.elem(i,2) = timesMinusI(l.elem(i,3));
    d.elem(i,3) = timesMinusI(l.elem(i,2));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,10>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,10>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,10>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,1);
    d.elem(i,1) = -l.elem(i,0);
    d.elem(i,2) = -l.elem(i,3);
    d.elem(i,3) =  l.elem(i,2);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,11>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,11>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,11>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesMinusI(l.elem(i,2));
    d.elem(i,1) = timesI(l.elem(i,3));
    d.elem(i,2) = timesMinusI(l.elem(i,0));
    d.elem(i,3) = timesI(l.elem(i,1));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,12>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,12>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,12>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,0));
    d.elem(i,1) = timesMinusI(l.elem(i,1));
    d.elem(i,2) = timesMinusI(l.elem(i,2));
    d.elem(i,3) = timesI(l.elem(i,3));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,13>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,13>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,13>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,3);
    d.elem(i,1) = -l.elem(i,2);
    d.elem(i,2) =  l.elem(i,1);
    d.elem(i,3) = -l.elem(i,0);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,14>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,14>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,14>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesMinusI(l.elem(i,3));
    d.elem(i,1) = timesMinusI(l.elem(i,2));
    d.elem(i,2) = timesMinusI(l.elem(i,1));
    d.elem(i,3) = timesMinusI(l.elem(i,0));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,15>, OpGammaConstMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConst<4,15>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConst<4,15>, OpGammaConstMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,0);
    d.elem(i,1) =  l.elem(i,1);
    d.elem(i,2) = -l.elem(i,2);
    d.elem(i,3) = -l.elem(i,3);
  }
  
  return d;
}


//-----------------------------------------------

// SpinMatrix<4> = GammaDP<4,m> * SpinMatrix<4>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<GammaConstDP<4,0>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,0>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,0>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;
  
  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = r.elem(0,i);
    d.elem(1,i) = r.elem(1,i);
    d.elem(2,i) = r.elem(2,i);
    d.elem(3,i) = r.elem(3,i);
  }

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,1>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,1>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,1>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;
  
  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesMinusI(r.elem(3,i));
    d.elem(1,i) = timesMinusI(r.elem(2,i));
    d.elem(2,i) = timesI(r.elem(1,i));
    d.elem(3,i) = timesI(r.elem(0,i));
  }

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,2>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,2>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,2>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(3,i);
    d.elem(1,i) =  r.elem(2,i);
    d.elem(2,i) =  r.elem(1,i);
    d.elem(3,i) = -r.elem(0,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,3>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,3>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,3>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(0,i));
    d.elem(1,i) = timesMinusI(r.elem(1,i));
    d.elem(2,i) = timesI(r.elem(2,i));
    d.elem(3,i) = timesMinusI(r.elem(3,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,4>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,4>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,4>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesMinusI(r.elem(2,i));
    d.elem(1,i) = timesI(r.elem(3,i));
    d.elem(2,i) = timesI(r.elem(0,i));
    d.elem(3,i) = timesMinusI(r.elem(1,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,5>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,5>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,5>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(1,i);
    d.elem(1,i) =  r.elem(0,i);
    d.elem(2,i) = -r.elem(3,i);
    d.elem(3,i) =  r.elem(2,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,6>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,6>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,6>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(1,i));
    d.elem(1,i) = timesI(r.elem(0,i));
    d.elem(2,i) = timesI(r.elem(3,i));
    d.elem(3,i) = timesI(r.elem(2,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,7>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,7>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,7>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) =  r.elem(2,i);
    d.elem(1,i) =  r.elem(3,i);
    d.elem(2,i) = -r.elem(0,i);
    d.elem(3,i) = -r.elem(1,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,8>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,8>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,8>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) =  r.elem(0,i);
    d.elem(1,i) =  r.elem(1,i);
    d.elem(2,i) = -r.elem(2,i);
    d.elem(3,i) = -r.elem(3,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,9>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,9>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,9>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(3,i));
    d.elem(1,i) = timesI(r.elem(2,i));
    d.elem(2,i) = timesI(r.elem(1,i));
    d.elem(3,i) = timesI(r.elem(0,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,10>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,10>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,10>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) =  r.elem(3,i);
    d.elem(1,i) = -r.elem(2,i);
    d.elem(2,i) =  r.elem(1,i);
    d.elem(3,i) = -r.elem(0,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,11>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,11>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,11>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(0,i));
    d.elem(1,i) = timesMinusI(r.elem(1,i));
    d.elem(2,i) = timesMinusI(r.elem(2,i));
    d.elem(3,i) = timesI(r.elem(3,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,12>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,12>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,12>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(2,i));
    d.elem(1,i) = timesMinusI(r.elem(3,i));
    d.elem(2,i) = timesI(r.elem(0,i));
    d.elem(3,i) = timesMinusI(r.elem(1,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,13>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,13>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,13>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(1,i);
    d.elem(1,i) =  r.elem(0,i);
    d.elem(2,i) =  r.elem(3,i);
    d.elem(3,i) = -r.elem(2,i);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,14>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,14>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,14>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = timesI(r.elem(1,i));
    d.elem(1,i) = timesI(r.elem(0,i));
    d.elem(2,i) = timesMinusI(r.elem(3,i));
    d.elem(3,i) = timesMinusI(r.elem(2,i));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,15>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,15>&, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,15>, PSpinMatrix<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = -r.elem(2,i);
    d.elem(1,i) = -r.elem(3,i);
    d.elem(2,i) = -r.elem(0,i);
    d.elem(3,i) = -r.elem(1,i);
  }
  
  return d;
}


// SpinMatrix<4> = SpinMatrix<4> * GammaDP<4,m>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,0>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,0>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,0>, OpGammaConstDPMultiply>::Type_t  d; 

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,0);
    d.elem(i,1) =  l.elem(i,1);
    d.elem(i,2) =  l.elem(i,2);
    d.elem(i,3) =  l.elem(i,3);
  }
 
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,1>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,1>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,1>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,3));
    d.elem(i,1) = timesI(l.elem(i,2));
    d.elem(i,2) = timesMinusI(l.elem(i,1));
    d.elem(i,3) = timesMinusI(l.elem(i,0));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,2>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,2>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,2>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = -l.elem(i,3);
    d.elem(i,1) =  l.elem(i,2);
    d.elem(i,2) =  l.elem(i,1);
    d.elem(i,3) = -l.elem(i,0);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,3>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,3>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,3>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,0));
    d.elem(i,1) = timesMinusI(l.elem(i,1));
    d.elem(i,2) = timesI(l.elem(i,2));
    d.elem(i,3) = timesMinusI(l.elem(i,3));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,4>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,4>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,4>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,2));
    d.elem(i,1) = timesMinusI(l.elem(i,3));
    d.elem(i,2) = timesMinusI(l.elem(i,0));
    d.elem(i,3) = timesI(l.elem(i,1));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,5>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,5>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,5>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,1);
    d.elem(i,1) = -l.elem(i,0);
    d.elem(i,2) =  l.elem(i,3);
    d.elem(i,3) = -l.elem(i,2);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,6>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,6>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,6>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,1));
    d.elem(i,1) = timesI(l.elem(i,0));
    d.elem(i,2) = timesI(l.elem(i,3));
    d.elem(i,3) = timesI(l.elem(i,2));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,7>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,7>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,7>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = -l.elem(i,2);
    d.elem(i,1) = -l.elem(i,3);
    d.elem(i,2) =  l.elem(i,0);
    d.elem(i,3) =  l.elem(i,1);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,8>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,8>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,8>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,0);
    d.elem(i,1) =  l.elem(i,1);
    d.elem(i,2) = -l.elem(i,2);
    d.elem(i,3) = -l.elem(i,3);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,9>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,9>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,9>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,3));
    d.elem(i,1) = timesI(l.elem(i,2));
    d.elem(i,2) = timesI(l.elem(i,1));
    d.elem(i,3) = timesI(l.elem(i,0));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,10>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,10>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,10>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = -l.elem(i,3);
    d.elem(i,1) =  l.elem(i,2);
    d.elem(i,2) = -l.elem(i,1);
    d.elem(i,3) =  l.elem(i,0);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,11>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,11>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,11>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,0));
    d.elem(i,1) = timesMinusI(l.elem(i,1));
    d.elem(i,2) = timesMinusI(l.elem(i,2));
    d.elem(i,3) = timesI(l.elem(i,3));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,12>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,12>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,12>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,2));
    d.elem(i,1) = timesMinusI(l.elem(i,3));
    d.elem(i,2) = timesI(l.elem(i,0));
    d.elem(i,3) = timesMinusI(l.elem(i,1));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,13>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,13>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,13>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) =  l.elem(i,1);
    d.elem(i,1) = -l.elem(i,0);
    d.elem(i,2) = -l.elem(i,3);
    d.elem(i,3) =  l.elem(i,2);
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,14>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,14>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,14>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = timesI(l.elem(i,1));
    d.elem(i,1) = timesI(l.elem(i,0));
    d.elem(i,2) = timesMinusI(l.elem(i,3));
    d.elem(i,3) = timesMinusI(l.elem(i,2));
  }
  
  return d;
}

template<class T2>
inline typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,15>, OpGammaConstDPMultiply>::Type_t
operator*(const PSpinMatrix<T2,4>& l, const GammaConstDP<4,15>&)
{
  typename BinaryReturn<PSpinMatrix<T2,4>, GammaConstDP<4,15>, OpGammaConstDPMultiply>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(i,0) = -l.elem(i,2);
    d.elem(i,1) = -l.elem(i,3);
    d.elem(i,2) = -l.elem(i,0);
    d.elem(i,3) = -l.elem(i,1);
  }
  
  return d;
}


//-----------------------------------------------------------------------------
//! PSpinVector<T,4> = P_+ * PSpinVector<T,4>
template<class T>
inline typename UnaryReturn<PSpinMatrix<T,4>, FnChiralProjectPlus>::Type_t
chiralProjectPlus(const PSpinMatrix<T,4>& s1)
{
  typename UnaryReturn<PSpinMatrix<T,4>, FnChiralProjectPlus>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    d.elem(0,i) = s1.elem(0,i);
    d.elem(1,i) = s1.elem(1,i);
    zero_rep(d.elem(2,i));
    zero_rep(d.elem(3,i));
  }

  return d;
}

//! PSpinVector<T,4> = P_- * PSpinVector<T,4>
template<class T>
inline typename UnaryReturn<PSpinMatrix<T,4>, FnChiralProjectMinus>::Type_t
chiralProjectMinus(const PSpinMatrix<T,4>& s1)
{
  typename UnaryReturn<PSpinMatrix<T,4>, FnChiralProjectMinus>::Type_t  d;

  for(int i=0; i < 4; ++i)
  {
    zero_rep(d.elem(0,i));
    zero_rep(d.elem(1,i));
    d.elem(2,i) = s1.elem(2,i);
    d.elem(3,i) = s1.elem(3,i);
  }

  return d;
}

//------------------------------------------
// PScalar = traceSpinQuarkContract13(PSpinMatrix,PSpinMatrix)
template<class T1, class T2>
struct BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnTraceSpinQuarkContract13> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnTraceSpinQuarkContract13>::Type_t>  Type_t;
};

//! PScalar = traceSpinQuarkContract13(PSpinMatrix,PSpinMatrix)
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnTraceSpinQuarkContract13>::Type_t
traceSpinQuarkContract13(const PSpinMatrix<T1,4>& l, const PSpinMatrix<T2,4>& r)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnTraceSpinQuarkContract13>::Type_t  d;

  d.elem() = quarkContractXX(l.elem(0,0), r.elem(0,0));
  for(int k=1; k < 4; ++k)
    d.elem() += quarkContractXX(l.elem(k,0), r.elem(k,0));

  for(int j=1; j < 4; ++j)
    for(int k=0; k < 4; ++k)
      d.elem() += quarkContractXX(l.elem(k,j), r.elem(k,j));

  return d;
}


// quark propagator contraction
template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract13>::Type_t
quarkContract13(const PSpinMatrix<T1,4>& s1, const PSpinMatrix<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract13>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(0,i), s2.elem(0,j));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(k,i), s2.elem(k,j));
    }

  return d;
}

template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract14>::Type_t
quarkContract14(const PSpinMatrix<T1,4>& s1, const PSpinMatrix<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract14>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(0,i), s2.elem(j,0));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(k,i), s2.elem(j,k));
    }

  return d;
}

template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract23>::Type_t
quarkContract23(const PSpinMatrix<T1,4>& s1, const PSpinMatrix<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract23>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(i,0), s2.elem(0,j));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(i,k), s2.elem(k,j));
    }

  return d;
}

template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract24>::Type_t
quarkContract24(const PSpinMatrix<T1,4>& s1, const PSpinMatrix<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract24>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(i,0), s2.elem(j,0));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(i,k), s2.elem(j,k));
    }

  return d;
}

template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract12>::Type_t
quarkContract12(const PSpinMatrix<T1,4>& s1, const PSpinMatrix<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract12>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(0,0), s2.elem(i,j));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(k,k), s2.elem(i,j));
    }

  return d;
}

template<class T1, class T2>
inline typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract34>::Type_t
quarkContract34(const PSpinMatrix<T1,4>& s1, const PSpinMatrix<T2,4>& s2)
{
  typename BinaryReturn<PSpinMatrix<T1,4>, PSpinMatrix<T2,4>, FnQuarkContract34>::Type_t  d;

  for(int j=0; j < 4; ++j)
    for(int i=0; i < 4; ++i)
    {
      d.elem(i,j) = quarkContractXX(s1.elem(i,j), s2.elem(0,0));
      for(int k=1; k < 4; ++k)
	d.elem(i,j) += quarkContractXX(s1.elem(i,j), s2.elem(k,k));
    }

  return d;
}

/*! @} */   // end of group primspinmatrix

} // namespace QDP

#endif
