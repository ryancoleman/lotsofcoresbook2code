// -*- C++ -*-

/*! \file
 * \brief Primitive Vector
 */


#ifndef QDP_PRIMVECTOR_H
#define QDP_PRIMVECTOR_H

namespace QDP {


//-------------------------------------------------------------------------------------
/*! \addtogroup primvector Vector primitive
 * \ingroup fiber
 *
 * Primitive type that transforms like a vector
 *
 * @{
 */

//! Primitive Vector class
/*!
 * All vector classes inherit this class
 * NOTE: For efficiency, there can be no virtual methods, so the data
 * portion is a part of the generic class, hence it is called a domain
 * and not a category
 */
template <class T, int N, template<class,int> class C> class PVector
{
public:
  PVector() { }
  ~PVector() { }

  typedef C<T,N>  CC;

  //! PVector = PVector
  /*! Set equal to another PVector */
  template<class T1>
  inline
  CC& assign(const C<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);

      return static_cast<CC&>(*this);
    }

  //! PVector = PVector
  /*! Set equal to another PVector */
  template<class T1>
  inline
  CC& operator=(const C<T1,N>& rhs) 
    {
      return assign(rhs);
    }

  //! PVector += PVector
  template<class T1>
  inline
  CC& operator+=(const C<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem(i);

      return static_cast<CC&>(*this);
    }

  //! PVector -= PVector
  template<class T1>
  inline
  CC& operator-=(const C<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem(i);

      return static_cast<CC&>(*this);
    }

  //! PVector *= PScalar
  template<class T1>
  inline
  CC& operator*=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) *= rhs.elem();

      return static_cast<CC&>(*this);
    }

  //! PVector /= PScalar
  template<class T1>
  inline
  CC& operator/=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) /= rhs.elem();

      return static_cast<CC&>(*this);
    }


#if 0
  // NOTE: intentially avoid defining a copy constructor - let the compiler
  // generate one via the bit copy mechanism. This effectively achieves
  // the first form of the if below (QDP_USE_ARRAY_INITIALIZER) without having
  // to use that syntax which is not strictly legal in C++.

  //! Deep copy constructor
#if defined(QDP_USE_ARRAY_INITIALIZER)
  PVector(const PVector& a) : F(a.F) {}
#else
  /*! This is a copy form - legal but not necessarily efficient */
  PVector(const PVector& a)
    {
     
      for(int i=0; i < N; ++i)
	F[i] = a.F[i];
    }
#endif
#endif


public:
  T& elem(int i) {return F[i];}
  const T& elem(int i) const {return F[i];}

private:
  T F[N];
 
};


//! Stream input
template<class T, int N, template<class,int> class C>  
inline
std::istream& operator>>(std::istream& s, PVector<T,N,C>& d)
{
  for(int i=0; i < N; ++i)
    s >> d.elem(i);

  return s;
}

template<class T, int N, template<class,int> class C>  
inline
StandardInputStream& operator>>(StandardInputStream& s, PVector<T,N,C>& d)
{
  for(int i=0; i < N; ++i)
    s >> d.elem(i);

  return s;
}

//! Stream output
template<class T, int N, template<class,int> class C>  
inline
std::ostream& operator<<(std::ostream& s, const PVector<T,N,C>& d)
{
  for(int i=0; i < N; ++i)
    s << d.elem(i);

  return s;
}

//! Stream output
template<class T, int N, template<class,int> class C>  
inline
StandardOutputStream& operator<<(StandardOutputStream& s, const PVector<T,N,C>& d)
{
  for(int i=0; i < N; ++i)
    s << d.elem(i);

  return s;
}


//! Text input
template<class T, int N, template<class,int> class C>  
inline
TextReader& operator>>(TextReader& txt, PVector<T,N,C>& d)
{
  for(int i=0; i < N; ++i)
    txt >> d.elem(i);

  return txt;
}

//! Text output
template<class T, int N, template<class,int> class C>  
inline
TextWriter& operator<<(TextWriter& txt, const PVector<T,N,C>& d)
{
  for(int i=0; i < N; ++i)
    txt << d.elem(i);

  return txt;
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T, int N, template<class,int> class C> 
inline
XMLWriter& operator<<(XMLWriter& xml, const PVector<T,N,C>& d)
{
  xml.openTag("Vector");

  XMLWriterAPI::AttributeList alist;

  // Copy into another array first
  for(int i=0; i < N; ++i)
  {
    alist.clear();
    alist.push_back(XMLWriterAPI::Attribute("row", i));

    xml.openTag("elem", alist);
    xml << d.elem(i);
    xml.closeTag();
  }

  xml.closeTag();  // Vector
  return xml;
}
#endif
/*! @} */  // end of group primvector


//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1, int N, template<class,int> class C>
struct WordType<PVector<T1,N,C> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

template<class T1, int N, template<class, int> class C> 
struct SinglePrecType< PVector<T1,N,C> >
{
  typedef PVector< typename SinglePrecType<T1>::Type_t, N, C> Type_t;
};

template<class T1, int N, template<class, int> class C> 
struct DoublePrecType< PVector<T1,N,C> >
{
  typedef PVector< typename DoublePrecType<T1>::Type_t, N, C> Type_t;
};

// Internally used scalars
template<class T, int N, template<class,int> class C>
struct InternalScalar<PVector<T,N,C> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving grid alone
template<class T, int N, template<class,int> class C>
struct PrimitiveScalar<PVector<T,N,C> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T, int N, template<class,int> class C>
struct LatticeScalar<PVector<T,N,C> > {
  typedef C<typename LatticeScalar<T>::Type_t, N>  Type_t;
};

//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PVector) -> PVector
template<class T1, int N, template<class,int> class C, class Op>
struct UnaryReturn<PVector<T1,N,C>, Op> {
  typedef C<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};
// Default binary(PScalar,PVector) -> PVector
template<class T1, class T2, int N, template<class,int> class C, class Op>
struct BinaryReturn<PScalar<T1>, PVector<T2,N,C>, Op> {
  typedef C<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PMatrix,PVector) -> PVector
template<class T1, class T2, int N, template<class,int> class C1, 
  template<class,int> class C2, class Op>
struct BinaryReturn<PMatrix<T1,N,C1>, PVector<T2,N,C2>, Op> {
  typedef C2<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PVector,PScalar) -> PVector
template<class T1, class T2, int N, template<class,int> class C, class Op>
struct BinaryReturn<PVector<T1,N,C>, PScalar<T2>, Op> {
  typedef C<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PVector,PVector) -> PVector
template<class T1, class T2, int N, template<class,int> class C, class Op>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, Op> {
  typedef C<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<PScalar<T2>, OpCast<T1> > {
  typedef PScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif


// Assignment is different
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpAssign > {
  typedef C<T1,N> &Type_t;
};
 
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpAddAssign > {
  typedef C<T1,N> &Type_t;
};
 
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpSubtractAssign > {
  typedef C<T1,N> &Type_t;
};
 
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpMultiplyAssign > {
  typedef C<T1,N> &Type_t;
};
 
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpDivideAssign > {
  typedef C<T1,N> &Type_t;
};
 



//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup primvector */
/*! @{ */

// Primitive Vectors

template<class T1, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, OpUnaryPlus>::Type_t
operator+(const PVector<T1,N,C>& l)
{
  typename UnaryReturn<PVector<T1,N,C>, OpUnaryPlus>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = +l.elem(i);
  return d;
}


template<class T1, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, OpUnaryMinus>::Type_t
operator-(const PVector<T1,N,C>& l)
{
  typename UnaryReturn<PVector<T1,N,C>, OpUnaryMinus>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = -l.elem(i);
  return d;
}


template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpAdd>::Type_t
operator+(const PVector<T1,N,C>& l, const PVector<T2,N,C>& r)
{
  typename BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) + r.elem(i);
  return d;
}


template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpSubtract>::Type_t
operator-(const PVector<T1,N,C>& l, const PVector<T2,N,C>& r)
{
  typename BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) - r.elem(i);
  return d;
}


// PVector * PScalar
template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpMultiply>::Type_t
operator*(const PVector<T1,N,C>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

// Optimized  PVector * adj(PScalar)
template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const PVector<T1,N,C>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpMultiplyAdj>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = multiplyAdj(l.elem(i), r.elem());
  return d;
}


// PScalar * PVector
template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PScalar<T1>, PVector<T2,N,C>, OpMultiply>::Type_t
operator*(const PScalar<T1>& l, const PVector<T2,N,C>& r)
{
  typename BinaryReturn<PScalar<T1>, PVector<T2,N,C>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}

// Optimized  adj(PScalar) * PVector
template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PScalar<T1>, PVector<T2,N,C>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<T1>& l, const PVector<T2,N,C>& r)
{
  typename BinaryReturn<PScalar<T1>, PVector<T2,N,C>, OpAdjMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = adjMultiply(l.elem(), r.elem(i));
  return d;
}


// PMatrix * PVector
template<class T1, class T2, int N, template<class,int> class C1, template<class,int> class C2>
inline typename BinaryReturn<PMatrix<T1,N,C1>, PVector<T2,N,C2>, OpMultiply>::Type_t
operator*(const PMatrix<T1,N,C1>& l, const PVector<T2,N,C2>& r)
{
  typename BinaryReturn<PMatrix<T1,N,C1>, PVector<T2,N,C2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
  {
    d.elem(i) = l.elem(i,0) * r.elem(0);
    for(int j=1; j < N; ++j)
      d.elem(i) += l.elem(i,j) * r.elem(j);
  }

  return d;
}

// Optimized  adj(PMatrix)*PVector
template<class T1, class T2, int N, template<class,int> class C1, template<class,int> class C2>
inline typename BinaryReturn<PMatrix<T1,N,C1>, PVector<T2,N,C2>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<T1,N,C1>& l, const PVector<T2,N,C2>& r)
{
  typename BinaryReturn<PMatrix<T1,N,C1>, PVector<T2,N,C2>, OpAdjMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
  {
    d.elem(i) = adjMultiply(l.elem(0,i), r.elem(0));
    for(int j=1; j < N; ++j)
      d.elem(i) += adjMultiply(l.elem(j,i), r.elem(j));
  }

  return d;
}


template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpDivide>::Type_t
operator/(const PVector<T1,N,C>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PVector<T1,N,C>, PScalar<T2>, OpDivide>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) / r.elem();
  return d;
}



//! PVector = Re(PVector)
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnReal>::Type_t
real(const PVector<T,N,C>& s1)
{
  typename UnaryReturn<PVector<T,N,C>, FnReal>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = real(s1.elem(i));

  return d;
}


//! PVector = Im(PVector)
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnImag>::Type_t
imag(const PVector<T,N,C>& s1)
{
  typename UnaryReturn<PVector<T,N,C>, FnImag>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = imag(s1.elem(i));

  return d;
}


//! PVector<T> = (PVector<T> , PVector<T>)
template<class T1, class T2, int N, template<class,int> class C>
inline typename BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, FnCmplx>::Type_t
cmplx(const PVector<T1,N,C>& s1, const PVector<T2,N,C>& s2)
{
  typename BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, FnCmplx>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = cmplx(s1.elem(i), s2.elem(i));

  return d;
}


//-----------------------------------------------------------------------------
// These functions always return bool
//! isnan
template<class T1, int N, template<class,int> class C>
struct UnaryReturn<PVector<T1,N,C>, FnIsNan> {
  bool Type_t;
};

template<class T1, int N, template<class,int> class C>
inline bool
isnan(const PVector<T1,N,C>& l)
{
  bool d = false;

  for(int i=0; i < N; ++i)
    d |= isnan(l.elem(i));

  return d;
}

//! isinf
template<class T1, int N, template<class,int> class C>
struct UnaryReturn<PVector<T1,N,C>, FnIsInf> {
  bool Type_t;
};

template<class T1, int N, template<class,int> class C>
inline bool
isinf(const PVector<T1,N,C>& l)
{
  bool d = false;

  for(int i=0; i < N; ++i)
    d |= isinf(l.elem(i));

  return d;
}

//! isnormal
template<class T1, int N, template<class,int> class C>
struct UnaryReturn<PVector<T1,N,C>, FnIsNormal> {
  bool Type_t;
};

template<class T1, int N, template<class,int> class C>
inline bool
isnormal(const PVector<T1,N,C>& l)
{
  bool d = true;

  for(int i=0; i < N; ++i)
    d &= isnormal(l.elem(i));

  return d;
}

//! isfinite
template<class T1, int N, template<class,int> class C>
struct UnaryReturn<PVector<T1,N,C>, FnIsFinite> {
  bool Type_t;
};

template<class T1, int N, template<class,int> class C>
inline bool
isfinite(const PVector<T1,N,C>& l)
{
  bool d = true;

  for(int i=0; i < N; ++i)
    d &= isfinite(l.elem(i));

  return d;
}


//-----------------------------------------------------------------------------
// Functions
// Conjugate
template<class T1, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, FnConjugate>::Type_t
conj(const PVector<T1,N,C>& l)
{
  typename UnaryReturn<PVector<T1,N,C>, FnConjugate>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = conj(l.elem(i));

  return d;
}

//! PVector = i * PVector
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnTimesI>::Type_t
timesI(const PVector<T,N,C>& s1)
{
  typename UnaryReturn<PVector<T,N,C>, FnTimesI>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = timesI(s1.elem(i));

  return d;
}

//! PVector = -i * PVector
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnTimesMinusI>::Type_t
timesMinusI(const PVector<T,N,C>& s1)
{
  typename UnaryReturn<PVector<T,N,C>, FnTimesMinusI>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = timesMinusI(s1.elem(i));

  return d;
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnGetSite>::Type_t
getSite(const PVector<T,N,C>& s1, int innersite)
{ 
  typename UnaryReturn<PVector<T,N,C>, FnGetSite>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = getSite(s1.elem(i), innersite);

  return d;
}

//! Extract color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnPeekColorVector>::Type_t
peekColor(const PVector<T,N,C>& l, int row)
{
  typename UnaryReturn<PVector<T,N,C>, FnPeekColorVector>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekColor(l.elem(i),row);
  return d;
}

//! Extract color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnPeekColorMatrix>::Type_t
peekColor(const PVector<T,N,C>& l, int row, int col)
{
  typename UnaryReturn<PVector<T,N,C>, FnPeekColorMatrix>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekColor(l.elem(i),row,col);
  return d;
}

//! Extract spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnPeekSpinVector>::Type_t
peekSpin(const PVector<T,N,C>& l, int row)
{
  typename UnaryReturn<PVector<T,N,C>, FnPeekSpinVector>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekSpin(l.elem(i),row);
  return d;
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnPeekSpinMatrix>::Type_t
peekSpin(const PVector<T,N,C>& l, int row, int col)
{
  typename UnaryReturn<PVector<T,N,C>, FnPeekSpinMatrix>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekSpin(l.elem(i),row,col);
  return d;
}

//! Insert color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, FnPokeColorVector>::Type_t&
pokeColor(PVector<T1,N,C>& l, const PVector<T2,N,C>& r, int row)
{
  typedef typename UnaryReturn<PVector<T1,N,C>, FnPokeColorVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeColor(l.elem(i),r.elem(i),row);
  return static_cast<Return_t&>(l);
}

//! Insert color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, FnPokeColorVector>::Type_t&
pokeColor(PVector<T1,N,C>& l, const PVector<T2,N,C>& r, int row, int col)
{
  typedef typename UnaryReturn<PVector<T1,N,C>, FnPokeColorVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeColor(l.elem(i),r.elem(i),row,col);
  return static_cast<Return_t&>(l);
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, FnPokeSpinVector>::Type_t&
pokeSpin(PVector<T1,N,C>& l, const PVector<T2,N,C>& r, int row)
{
  typedef typename UnaryReturn<PVector<T1,N,C>, FnPokeSpinVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeSpin(l.elem(i),r.elem(i),row);
  return static_cast<Return_t&>(l);
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T1,N,C>, FnPokeSpinVector>::Type_t&
pokeSpin(PVector<T1,N,C>& l, const PVector<T2,N,C>& r, int row, int col)
{
  typedef typename UnaryReturn<PVector<T1,N,C>, FnPokeSpinVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeSpin(l.elem(i),r.elem(i),row,col);
  return static_cast<Return_t&>(l);
}


//! dest = 0
template<class T, int N, template<class,int> class C> 
inline void 
zero_rep(PVector<T,N,C>& dest) 
{
  for(int i=0; i < N; ++i)
    zero_rep(dest.elem(i));
}

//! dest = (mask) ? s1 : dest
template<class T, class T1, int N, template<class,int> class C> 
inline void 
copymask(PVector<T,N,C>& d, const PScalar<T1>& mask, const PVector<T,N,C>& s1) 
{
  for(int i=0; i < N; ++i)
    copymask(d.elem(i),mask.elem(),s1.elem(i));
}


//! dest [some type] = source [some type]
template<class T, class T1, int N, template<class,int> class C>
inline void 
copy_site(PVector<T,N,C>& d, int isite, const PVector<T1,N,C>& s1)
{
  for(int i=0; i < N; ++i)
    copy_site(d.elem(i), isite, s1.elem(i));
}

//! dest [some type] = source [some type]
template<class T, class T1, int N, template<class,int> class C>
inline void 
copy_site(PVector<T,N,C>& d, int isite, const PScalar<T1>& s1)
{
  for(int i=0; i < N; ++i)
    copy_site(d.elem(i), isite, s1.elem());
}


//! gather several inner sites together
template<class T, class T1, int N, template<class,int> class C>
inline void 
gather_sites(PVector<T,N,C>& d, 
	     const PVector<T1,N,C>& s0, int i0, 
	     const PVector<T1,N,C>& s1, int i1,
	     const PVector<T1,N,C>& s2, int i2,
	     const PVector<T1,N,C>& s3, int i3)
{
  for(int i=0; i < N; ++i)
    gather_sites(d.elem(i), 
		 s0.elem(i), i0, 
		 s1.elem(i), i1, 
		 s2.elem(i), i2, 
		 s3.elem(i), i3);
}


//! dest  = random  
template<class T, int N, template<class,int> class C, class T1, class T2>
inline void
fill_random(PVector<T,N,C>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  // Loop over rows the slowest
  for(int i=0; i < N; ++i)
    fill_random(d.elem(i), seed, skewed_seed, seed_mult);
}


//! dest  = gaussian
template<class T, int N, template<class,int> class C>
inline void
fill_gaussian(PVector<T,N,C>& d, PVector<T,N,C>& r1, PVector<T,N,C>& r2)
{
  for(int i=0; i < N; ++i)
    fill_gaussian(d.elem(i), r1.elem(i), r2.elem(i));
}


#if 0
// Global sum over site indices only
template<class T, int N, template<class,int> class C>
struct UnaryReturn<PVector<T,N,C>, FnSum > {
  typedef C<typename UnaryReturn<T, FnSum>::Type_t, N>  Type_t;
};

template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnSum>::Type_t
sum(const PVector<T,N,C>& s1)
{
  typename UnaryReturn<PVector<T,N,C>, FnSum>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = sum(s1.elem(i));

  return d;
}
#endif


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T, int N, template<class,int> class C>
struct UnaryReturn<PVector<T,N,C>, FnNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N, template<class,int> class C>
struct UnaryReturn<PVector<T,N,C>, FnLocalNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T, int N, template<class,int> class C>
inline typename UnaryReturn<PVector<T,N,C>, FnLocalNorm2>::Type_t
localNorm2(const PVector<T,N,C>& s1)
{
  typename UnaryReturn<PVector<T,N,C>, FnLocalNorm2>::Type_t  d;

  d.elem() = localNorm2(s1.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += localNorm2(s1.elem(i));

  return d;
}


//! PScalar<T> = InnerProduct(adj(PVector<T1>)*PVector<T1>)
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, FnInnerProduct > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, FnLocalInnerProduct > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N, template<class,int> class C>
inline PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>
localInnerProduct(const PVector<T1,N,C>& s1, const PVector<T2,N,C>& s2)
{
  PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  d;

  d.elem() = localInnerProduct(s1.elem(0), s2.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += localInnerProduct(s1.elem(i), s2.elem(i));

  return d;
}


//! PScalar<T> = InnerProductReal(adj(PVector<T1>)*PVector<T1>)
/*!
 * return  realpart of InnerProduct(adj(s1)*s2)
 */
template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, FnInnerProductReal > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N, template<class,int> class C>
struct BinaryReturn<PVector<T1,N,C>, PVector<T2,N,C>, FnLocalInnerProductReal > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N, template<class,int> class C>
inline PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>
localInnerProductReal(const PVector<T1,N,C>& s1, const PVector<T2,N,C>& s2)
{
  PScalar<typename BinaryReturn<T1,T2, FnLocalInnerProductReal>::Type_t>  d;

  d.elem() = localInnerProductReal(s1.elem(0), s2.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += localInnerProductReal(s1.elem(i), s2.elem(i));

  return d;
}


// This PVector<T1,N,C> stuff versus PSpinVector<T1,N> is causing problems. 
// When searching for type matching functions, the language does not allow
// for varying template arguments to match a function. We should just move
// away from PVector to use PSpinVector and PColorVector. However, have to
// replicate all the functions. Uggh - another day...

//
////! PVector<T> = localInnerProduct(adj(PScalar<T1>)*PVector<T1>)
//template<class T1, class T2, int N, template<class,int> class C>
//struct BinaryReturn<PScalar<T1>, PVector<T2,N,C>, FnLocalInnerProduct> {
//  typedef PVector<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t, N, C>  Type_t;
//};
//
//template<class T1, class T2, int N, template<class,int> class C>
//inline PVector<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t,N,C>
//localInnerProduct(const PScalar<T1>& s1, const PVector<T2,N,C>& s2)
//{
//  typename BinaryReturn<PScalar<T1>, PVector<T2,N,C>, FnLocalInnerProduct>::Type_t  d;
//
//  for(int i=0; i < N; ++i)
//    d.elem(i) = localInnerProduct(s1.elem(0), s2.elem(i));
//
//  return d;
//}
//
//
////! PScalar<T> = InnerProductReal(adj(PScalar<T1>)*PVector<T1>)
///*!
// * return  realpart of InnerProduct(adj(s1)*s2)
// */
//template<class T1, class T2, int N, template<class,int> class C>
//struct BinaryReturn<PScalar<T1>, PVector<T2,N,C>, FnLocalInnerProductReal > {
//  typedef PVector<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t, N,C>  Type_t;
//};
//
//template<class T1, class T2, int N, template<class,int> class C>
//inline PVector<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t,N,C>
//localInnerProductReal(const PScalar<T1>& s1, const PVector<T2,N,C>& s2)
//{
//  typename BinaryReturn<PScalar<T1>, PVector<T2,N,C>, FnLocalInnerProductReal>::Type_t  d;
//
//  for(int i=0; i < N; ++i)
//    d.elem(i) = localInnerProductReal(s1.elem(), s2.elem(i));
//
//  return d;
//}


//! PVector<T> = where(PScalar, PVector, PVector)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N, template<class,int> class C>
struct TrinaryReturn<PScalar<T1>, PVector<T2,N,C>, PVector<T3,N,C>, FnWhere> {
  typedef C<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N, template<class,int> class C>
inline typename TrinaryReturn<PScalar<T1>, PVector<T2,N,C>, PVector<T3,N,C>, FnWhere>::Type_t
where(const PScalar<T1>& a, const PVector<T2,N,C>& b, const PVector<T3,N,C>& c)
{
  typename TrinaryReturn<PScalar<T1>, PVector<T2,N,C>, PVector<T3,N,C>, FnWhere>::Type_t  d;

  // Not optimal - want to have where outside assignment
  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem(i));

  return d;
}

/*! @} */  // end of group primvector

} // namespace QDP

#endif
