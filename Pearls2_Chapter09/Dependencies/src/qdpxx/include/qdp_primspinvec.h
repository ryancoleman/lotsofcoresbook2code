// -*- C++ -*-

/*! \file
 * \brief Primitive Spin Vector
 */


#ifndef QDP_PRIMSPINVEC_H
#define QDP_PRIMSPINVEC_H

namespace QDP {


//-------------------------------------------------------------------------------------
/*! \addtogroup primspinvector Spin vector primitive
 * \ingroup primvector
 *
 * Primitive type that transforms like a Spin vector
 *
 * @{
 */

//! Primitive spin Vector class
/*! 
 * Spin vector class supports gamma matrix algebra 
 *
 * NOTE: the class is mostly empty - it is the specialized versions below
 * that know for a fixed size how gamma matrices (constants) should act
 * on the spin vectors.
 */

template <class T, int N> class PSpinVector
{
public:
  PSpinVector() { }
  ~PSpinVector() {} 

  template<class T1>
  inline
  PSpinVector& assign(const PSpinVector<T1, N>& rhs)
  {
    for(int i=0; i < N; i++) 
      elem(i) = rhs.elem(i);

    return *this;
  }

  template<class T1> 
  inline 
  PSpinVector& operator=(const PSpinVector<T1,N>& rhs)
  {
    return assign(rhs);
  }

  //! PSpinVector += PSpinVector
  template<class T1>
  inline
  PSpinVector& operator+=(const PSpinVector<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem(i);

      return *this;
    }


  template<class T1>
  inline
  PSpinVector& operator*=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) *= rhs.elem();

      return *this;
    }

  template<class T1>
  inline
  PSpinVector& operator-=(const PSpinVector<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem(i);

      return *this;
    }

  //! PSpinVector /= PScalar
  template<class T1>
  inline
  PSpinVector& operator/=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) /= rhs.elem();

      return *this;
    }

  T& elem(int i) {return F[i];}
  const T& elem(int i) const {return F[i];}
private:
  T F[N] QDP_ALIGN16; 
};

//! Stream input
template<class T, int N>  
inline
std::istream& operator>>(std::istream& s, PSpinVector<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s >> d.elem(i);

  return s;
}

template<class T, int N>  
inline
StandardInputStream& operator>>(StandardInputStream& s, PSpinVector<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s >> d.elem(i);

  return s;
}
//! Stream output
template<class T, int N>
inline
std::ostream& operator<<(std::ostream& s, const PSpinVector<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s << d.elem(i);

  return s;
}

//! Stream output
template<class T, int N>
inline
StandardOutputStream& operator<<(StandardOutputStream& s, const PSpinVector<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s << d.elem(i);

  return s;
}


//! Text input
template<class T, int N>
inline
TextReader& operator>>(TextReader& txt, PSpinVector<T,N>& d)
{
  for(int i=0; i < N; ++i)
    txt >> d.elem(i);

  return txt;
}

//! Text output
template<class T, int N>  
inline
TextWriter& operator<<(TextWriter& txt, const PSpinVector<T,N>& d)
{
  for(int i=0; i < N; ++i)
    txt << d.elem(i);

  return txt;
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T, int N>
inline
XMLWriter& operator<<(XMLWriter& xml, const PSpinVector<T,N>& d)
{
  xml.openTag("SpinVector");

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



 


// Primitive Vectors

template<class T1, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, OpUnaryPlus>::Type_t
operator+(const PSpinVector<T1,N>& l)
{
  typename UnaryReturn<PSpinVector<T1,N>, OpUnaryPlus>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = +l.elem(i);
  return d;
}


template<class T1, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, OpUnaryMinus>::Type_t
operator-(const PSpinVector<T1,N>& l)
{
  typename UnaryReturn<PSpinVector<T1,N>, OpUnaryMinus>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = -l.elem(i);
  return d;
}


template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpAdd>::Type_t
operator+(const PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) + r.elem(i);
  return d;
}


template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpSubtract>::Type_t
operator-(const PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) - r.elem(i);
  return d;
}


// PSpinVector * PScalar
template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpMultiply>::Type_t
operator*(const PSpinVector<T1,N>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

// Optimized  PSpinVector * adj(PScalar)
template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const PSpinVector<T1,N>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpMultiplyAdj>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = multiplyAdj(l.elem(i), r.elem());
  return d;
}


// PScalar * PSpinVector
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t
operator*(const PScalar<T1>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}

// Optimized  adj(PScalar) * PSpinVector
template<class T1, class T2, int N>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<T1>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, OpAdjMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = adjMultiply(l.elem(), r.elem(i));
  return d;
}


// PMatrix * PSpinVector
template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinMatrix<T1,N>, PSpinVector<T2,N>, OpMultiply>::Type_t
operator*(const PSpinMatrix<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinMatrix<T1,N>, PSpinVector<T2,N>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
  {
    d.elem(i) = l.elem(i,0) * r.elem(0);
    for(int j=1; j < N; ++j)
      d.elem(i) += l.elem(i,j) * r.elem(j);
  }

  return d;
}


// PMatrix * PSpinVector
template<class T1, class T2,  template<class,int> class C, int N>
inline typename BinaryReturn<PMatrix<T1,N,C>, PSpinVector<T2,N>, OpMultiply>::Type_t
operator*(const PSpinMatrix<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinMatrix<T1,N>, PSpinVector<T2,N>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
  {
    d.elem(i) = l.elem(i,0) * r.elem(0);
    for(int j=1; j < N; ++j)
      d.elem(i) += l.elem(i,j) * r.elem(j);
  }

  return d;
}

// Optimized  adj(PMatrix)*PSpinVector
template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinMatrix<T1,N>, PSpinVector<T2,N>, OpAdjMultiply>::Type_t
adjMultiply(const PSpinMatrix<T1,N>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PSpinMatrix<T1,N>, PSpinVector<T2,N>, OpAdjMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
  {
    d.elem(i) = adjMultiply(l.elem(0,i), r.elem(0));
    for(int j=1; j < N; ++j)
      d.elem(i) += adjMultiply(l.elem(j,i), r.elem(j));
  }

  return d;
}

// Optimized  adj(PMatrix)*PVector
template<class T1, class T2, int N, template<class,int> class C1>
inline typename BinaryReturn<PMatrix<T1,N,C1>, PSpinVector<T2,N>, OpAdjMultiply>::Type_t
adjMultiply(const PMatrix<T1,N,C1>& l, const PSpinVector<T2,N>& r)
{
  typename BinaryReturn<PMatrix<T1,N,C1>, PSpinVector<T2,N>, OpAdjMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
  {
    d.elem(i) = adjMultiply(l.elem(0,i), r.elem(0));
    for(int j=1; j < N; ++j)
      d.elem(i) += adjMultiply(l.elem(j,i), r.elem(j));
  }

  return d;
}

template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpDivide>::Type_t
operator/(const PSpinVector<T1,N>& l, const PScalar<T2>& r)
{
  typename BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpDivide>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) / r.elem();
  return d;
}



//! PSpinVector = Re(PSpinVector)
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnReal>::Type_t
real(const PSpinVector<T,N>& s1)
{
  typename UnaryReturn<PSpinVector<T,N>, FnReal>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = real(s1.elem(i));

  return d;
}


//! PSpinVector = Im(PSpinVector)
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnImag>::Type_t
imag(const PSpinVector<T,N>& s1)
{
  typename UnaryReturn<PSpinVector<T,N>, FnImag>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = imag(s1.elem(i));

  return d;
}


//! PSpinVector<T> = (PSpinVector<T> , PSpinVector<T>)
template<class T1, class T2, int N>
inline typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnCmplx>::Type_t
cmplx(const PSpinVector<T1,N>& s1, const PSpinVector<T2,N>& s2)
{
  typename BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnCmplx>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = cmplx(s1.elem(i), s2.elem(i));

  return d;
}


//-----------------------------------------------------------------------------
// These functions always return bool
//! isnan
template<class T1, int N>
struct UnaryReturn<PSpinVector<T1,N>, FnIsNan> {
  bool Type_t;
};

template<class T1, int N>
inline bool
isnan(const PSpinVector<T1,N>& l)
{
  bool d = false;

  for(int i=0; i < N; ++i)
    d |= isnan(l.elem(i));

  return d;
}

//! isinf
template<class T1, int N>
struct UnaryReturn<PSpinVector<T1,N>, FnIsInf> {
  bool Type_t;
};

template<class T1, int N>
inline bool
isinf(const PSpinVector<T1,N>& l)
{
  bool d = false;

  for(int i=0; i < N; ++i)
    d |= isinf(l.elem(i));

  return d;
}

//! isnormal
template<class T1, int N>
struct UnaryReturn<PSpinVector<T1,N>, FnIsNormal> {
  bool Type_t;
};

template<class T1, int N>
inline bool
isnormal(const PSpinVector<T1,N>& l)
{
  bool d = true;

  for(int i=0; i < N; ++i)
    d &= isnormal(l.elem(i));

  return d;
}

//! isfinite
template<class T1, int N>
struct UnaryReturn<PSpinVector<T1,N>, FnIsFinite> {
  bool Type_t;
};

template<class T1, int N>
inline bool
isfinite(const PSpinVector<T1,N>& l)
{
  bool d = true;

  for(int i=0; i < N; ++i)
    d &= isfinite(l.elem(i));

  return d;
}


//-----------------------------------------------------------------------------
// Functions
// Conjugate
template<class T1, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, FnConjugate>::Type_t
conj(const PSpinVector<T1,N>& l)
{
  typename UnaryReturn<PSpinVector<T1,N>, FnConjugate>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = conj(l.elem(i));

  return d;
}

//! PSpinVector = i * PSpinVector
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnTimesI>::Type_t
timesI(const PSpinVector<T,N>& s1)
{
  typename UnaryReturn<PSpinVector<T,N>, FnTimesI>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = timesI(s1.elem(i));

  return d;
}

//! PSpinVector = -i * PSpinVector
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnTimesMinusI>::Type_t
timesMinusI(const PSpinVector<T,N>& s1)
{
  typename UnaryReturn<PSpinVector<T,N>, FnTimesMinusI>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = timesMinusI(s1.elem(i));

  return d;
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnGetSite>::Type_t
getSite(const PSpinVector<T,N>& s1, int innersite)
{ 
  typename UnaryReturn<PSpinVector<T,N>, FnGetSite>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = getSite(s1.elem(i), innersite);

  return d;
}

//! Extract color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnPeekColorVector>::Type_t
peekColor(const PSpinVector<T,N>& l, int row)
{
  typename UnaryReturn<PSpinVector<T,N>, FnPeekColorVector>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekColor(l.elem(i),row);
  return d;
}

//! Extract color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnPeekColorMatrix>::Type_t
peekColor(const PSpinVector<T,N>& l, int row, int col)
{
  typename UnaryReturn<PSpinVector<T,N>, FnPeekColorMatrix>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekColor(l.elem(i),row,col);
  return d;
}


//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnPeekSpinMatrix>::Type_t
peekSpin(const PSpinVector<T,N>& l, int row, int col)
{
  typename UnaryReturn<PSpinVector<T,N>, FnPeekSpinMatrix>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = peekSpin(l.elem(i),row,col);
  return d;
}

//! Insert color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, FnPokeColorVector>::Type_t&
pokeColor(PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r, int row)
{
  typedef typename UnaryReturn<PSpinVector<T1,N>, FnPokeColorVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeColor(l.elem(i),r.elem(i),row);
  return static_cast<Return_t&>(l);
}

//! Insert color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, FnPokeColorVector>::Type_t&
pokeColor(PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r, int row, int col)
{
  typedef typename UnaryReturn<PSpinVector<T1,N>, FnPokeColorVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeColor(l.elem(i),r.elem(i),row,col);
  return static_cast<Return_t&>(l);
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, FnPokeSpinVector>::Type_t&
pokeSpin(PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r, int row)
{
  typedef typename UnaryReturn<PSpinVector<T1,N>, FnPokeSpinVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeSpin(l.elem(i),r.elem(i),row);
  return static_cast<Return_t&>(l);
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2, int N>
inline typename UnaryReturn<PSpinVector<T1,N>, FnPokeSpinVector>::Type_t&
pokeSpin(PSpinVector<T1,N>& l, const PSpinVector<T2,N>& r, int row, int col)
{
  typedef typename UnaryReturn<PSpinVector<T1,N>, FnPokeSpinVector>::Type_t  Return_t;

  for(int i=0; i < N; ++i)
    pokeSpin(l.elem(i),r.elem(i),row,col);
  return static_cast<Return_t&>(l);
}


//! dest = 0
template<class T, int N> 
inline void 
zero_rep(PSpinVector<T,N>& dest) 
{
  for(int i=0; i < N; ++i)
    zero_rep(dest.elem(i));
}

//! dest = (mask) ? s1 : dest
template<class T, class T1, int N> 
inline void 
copymask(PSpinVector<T,N>& d, const PScalar<T1>& mask, const PSpinVector<T,N>& s1) 
{
  for(int i=0; i < N; ++i)
    copymask(d.elem(i),mask.elem(),s1.elem(i));
}


//! dest [some type] = source [some type]
template<class T, class T1, int N>
inline void 
copy_site(PSpinVector<T,N>& d, int isite, const PSpinVector<T1,N>& s1)
{
  for(int i=0; i < N; ++i)
    copy_site(d.elem(i), isite, s1.elem(i));
}

//! dest [some type] = source [some type]
template<class T, class T1, int N>
inline void 
copy_site(PSpinVector<T,N>& d, int isite, const PScalar<T1>& s1)
{
  for(int i=0; i < N; ++i)
    copy_site(d.elem(i), isite, s1.elem());
}


//! gather several inner sites together
template<class T, class T1, int N>
inline void 
gather_sites(PSpinVector<T,N>& d, 
	     const PSpinVector<T1,N>& s0, int i0, 
	     const PSpinVector<T1,N>& s1, int i1,
	     const PSpinVector<T1,N>& s2, int i2,
	     const PSpinVector<T1,N>& s3, int i3)
{
  for(int i=0; i < N; ++i)
    gather_sites(d.elem(i), 
		 s0.elem(i), i0, 
		 s1.elem(i), i1, 
		 s2.elem(i), i2, 
		 s3.elem(i), i3);
}



#if 0
// Global sum over site indices only
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSum > {
  typedef PSpinVectortypename UnaryReturn<T, FnSum>::Type_t, N>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnSum>::Type_t
sum(const PSpinVector<T,N>& s1)
{
  typename UnaryReturn<PSpinVector<T,N>, FnSum>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = sum(s1.elem(i));

  return d;
}
#endif




template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnLocalNorm2>::Type_t
localNorm2(const PSpinVector<T,N>& s1)
{
  typename UnaryReturn<PSpinVector<T,N>, FnLocalNorm2>::Type_t  d;

  d.elem() = localNorm2(s1.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += localNorm2(s1.elem(i));

  return d;
}




//! PSpinVector<T> = where(PScalar, PSpinVector, PSpinVector)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<PScalar<T1>, PSpinVector<T2,N>, PSpinVector<T3,N>, FnWhere> {
  typedef PSpinVector<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<PScalar<T1>, PSpinVector<T2,N>, PSpinVector<T3,N>, FnWhere>::Type_t
where(const PScalar<T1>& a, const PSpinVector<T2,N>& b, const PSpinVector<T3,N>& c)
{
  typename TrinaryReturn<PScalar<T1>, PSpinVector<T2,N>, PSpinVector<T3,N>, FnWhere>::Type_t  d;

  // Not optimal - want to have where outside assignment
  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem(i));

  return d;
}


//! Specialization of primitive spin Vector class for 4 spin components
/*! 
 * Spin vector class supports gamma matrix algebra for 4 spin components
 */


//! Specialization of primitive spin Vector class for 2 spin components
/*! 
 * Spin vector class supports gamma matrix algebra for 2 spin components
 * NOTE: this can be used for spin projection tricks of a 4 component spinor
 * to 2 spin components, or a 2 spin component Dirac fermion in 2 dimensions
 */


/*! @} */   // end of group primspinvec

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1, int N>
struct WordType<PSpinVector<T1,N> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

// Fixed Precision
template<class T1, int N>
struct SinglePrecType< PSpinVector<T1, N> > 
{
  typedef PSpinVector< typename SinglePrecType<T1>::Type_t, N> Type_t;
};

template<class T1, int N>
struct DoublePrecType< PSpinVector<T1, N> > 
{
  typedef PSpinVector< typename DoublePrecType<T1>::Type_t, N> Type_t;
};

// Internally used scalars
template<class T, int N>
struct InternalScalar<PSpinVector<T,N> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive into a scalar leaving grid alone
template<class T, int N>
struct PrimitiveScalar<PSpinVector<T,N> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T, int N>
struct LatticeScalar<PSpinVector<T,N> > {
  typedef PSpinVector<typename LatticeScalar<T>::Type_t, N>  Type_t;
};

//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PSpinVector) -> PSpinVector
template<class T1, int N, class Op>
struct UnaryReturn<PSpinVector<T1,N>, Op> {
  typedef PSpinVector<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};
// Default binary(PScalar,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PScalar<T1>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinMatrix,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn< PSpinMatrix<T1,N>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


// Default binary(PMatrix,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, template <class,int> class C1, class Op>
struct BinaryReturn< PMatrix<T1,N,C1>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


// Default binary(PSpinVector,PScalar) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(PSpinVector,PSpinVector) -> PSpinVector
template<class T1, class T2, int N, class Op>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, Op> {
  typedef PSpinVector<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
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
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpAddAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, OpSubtractAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpMultiplyAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PScalar<T2>, OpDivideAssign > {
  typedef PSpinVector<T1,N> &Type_t;
};
 


// SpinVector
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnLocalNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnLocalInnerProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<PSpinVector<T1,N>, PSpinVector<T2,N>, FnLocalInnerProductReal> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};


template<class T1, class T2, int N>
inline PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>
localInnerProduct(const PSpinVector<T1,N>& s1, const PSpinVector<T2,N>& s2)
{
  PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  d;

  d.elem() = localInnerProduct(s1.elem(0), s2.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += localInnerProduct(s1.elem(i), s2.elem(i));

  return d;
}

template<class T1, class T2, int N>
inline PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>
localInnerProductReal(const PSpinVector<T1,N>& s1, const PSpinVector<T2,N>& s2)
{
  PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  d;

  d.elem() = localInnerProductReal(s1.elem(0), s2.elem(0));
  for(int i=1; i < N; ++i)
    d.elem() += localInnerProductReal(s1.elem(i), s2.elem(i));

  return d;
}

// Gamma algebra
template<int m, class T2, int N>
struct BinaryReturn<GammaConst<N,m>, PSpinVector<T2,N>, OpGammaConstMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N>
struct BinaryReturn<GammaType<N>, PSpinVector<T2,N>, OpGammaTypeMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

// Gamma algebra
template<int m, class T2, int N>
struct BinaryReturn<GammaConstDP<N,m>, PSpinVector<T2,N>, OpGammaConstDPMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N>
struct BinaryReturn<GammaTypeDP<N>, PSpinVector<T2,N>, OpGammaTypeDPMultiply> {
  typedef PSpinVector<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

// Generic Spin projection
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProject > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProject>::Type_t, (N>>1) >  Type_t;
};

// spin projection for each direction
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir0Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir0Plus>::Type_t, (N>>1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir1Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir1Plus>::Type_t, (N>>1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir2Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir2Plus>::Type_t, (N>>1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir3Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir3Plus>::Type_t, (N>>1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir0Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir0Minus>::Type_t, (N>>1) > Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir1Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir1Minus>::Type_t, (N>>1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir2Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir2Minus>::Type_t, (N>>1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinProjectDir3Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinProjectDir3Minus>::Type_t, (N>>1) >  Type_t;
};


// Generic Spin reconstruction
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstruct > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstruct>::Type_t, (N<<1) >  Type_t;
};

// spin reconstruction for each direction
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir0Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir0Plus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir1Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir1Plus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir2Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir2Plus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir3Plus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir3Plus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir0Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir0Minus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir1Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir1Minus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir2Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir2Minus>::Type_t, (N<<1) >  Type_t;
};

template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnSpinReconstructDir3Minus > {
  typedef PSpinVector<typename UnaryReturn<T, FnSpinReconstructDir3Minus>::Type_t, (N<<1) >  Type_t;
};




//! dest  = random  
template<class T, int N,  class T1, class T2>
inline void
fill_random(PSpinVector<T,N>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  // Loop over rows the slowest
  for(int i=0; i < N; ++i)
    fill_random(d.elem(i), seed, skewed_seed, seed_mult);
}


//! dest  = gaussian
template<class T, int N>
inline void
fill_gaussian(PSpinVector<T,N>& d, PSpinVector<T,N>& r1, PSpinVector<T,N>& r2)
{
  for(int i=0; i < N; ++i)
    fill_gaussian(d.elem(i), r1.elem(i), r2.elem(i));
}

//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup primspinvector */
/*! @{ */

// Peeking and poking
//! Extract spin vector components 
template<class T, int N>
struct UnaryReturn<PSpinVector<T,N>, FnPeekSpinVector > {
  typedef PScalar<typename UnaryReturn<T, FnPeekSpinVector>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<PSpinVector<T,N>, FnPeekSpinVector>::Type_t
peekSpin(const PSpinVector<T,N>& l, int row)
{
  typename UnaryReturn<PSpinVector<T,N>, FnPeekSpinVector>::Type_t  d;

  // Note, do not need to propagate down since the function is eaten at this level
  d.elem() = l.elem(row);
  return d;
}

//! Insert spin vector components
template<class T1, class T2, int N>
inline PSpinVector<T1,N>&
pokeSpin(PSpinVector<T1,N>& l, const PScalar<T2>& r, int row)
{
  // Note, do not need to propagate down since the function is eaten at this level
  l.elem(row) = r.elem();
  return l;
}



// SpinVector<4> = Gamma<4,m> * SpinVector<4>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<GammaConst<4,0>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,0>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,0>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;
  
  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) =  r.elem(2);
  d.elem(3) =  r.elem(3);

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,1>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,1>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,1>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;
  
  d.elem(0) = timesI(r.elem(3));
  d.elem(1) = timesI(r.elem(2));
  d.elem(2) = timesMinusI(r.elem(1));
  d.elem(3) = timesMinusI(r.elem(0));

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,2>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,2>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,2>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(3);
  d.elem(1) =  r.elem(2);
  d.elem(2) =  r.elem(1);
  d.elem(3) = -r.elem(0);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,3>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,3>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,3>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(0));
  d.elem(1) = timesI(r.elem(1));
  d.elem(2) = timesMinusI(r.elem(2));
  d.elem(3) = timesI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,4>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,4>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,4>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(2));
  d.elem(1) = timesMinusI(r.elem(3));
  d.elem(2) = timesMinusI(r.elem(0));
  d.elem(3) = timesI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,5>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,5>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,5>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) = -r.elem(3);
  d.elem(3) =  r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,6>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,6>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,6>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(1));
  d.elem(1) = timesMinusI(r.elem(0));
  d.elem(2) = timesMinusI(r.elem(3));
  d.elem(3) = timesMinusI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,7>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,7>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,7>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) =  r.elem(2);
  d.elem(1) =  r.elem(3);
  d.elem(2) = -r.elem(0);
  d.elem(3) = -r.elem(1);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,8>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,8>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,8>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) =  r.elem(2);
  d.elem(1) =  r.elem(3);
  d.elem(2) =  r.elem(0);
  d.elem(3) =  r.elem(1);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,9>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,9>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,9>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(1));
  d.elem(1) = timesI(r.elem(0));
  d.elem(2) = timesMinusI(r.elem(3));
  d.elem(3) = timesMinusI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,10>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,10>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,10>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) =  r.elem(3);
  d.elem(3) = -r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,11>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,11>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,11>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(2));
  d.elem(1) = timesI(r.elem(3));
  d.elem(2) = timesMinusI(r.elem(0));
  d.elem(3) = timesI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,12>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,12>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,12>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(0));
  d.elem(1) = timesMinusI(r.elem(1));
  d.elem(2) = timesMinusI(r.elem(2));
  d.elem(3) = timesI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,13>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,13>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,13>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = -r.elem(3);
  d.elem(1) =  r.elem(2);
  d.elem(2) = -r.elem(1);
  d.elem(3) =  r.elem(0);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,14>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,14>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,14>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(3));
  d.elem(1) = timesMinusI(r.elem(2));
  d.elem(2) = timesMinusI(r.elem(1));
  d.elem(3) = timesMinusI(r.elem(0));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConst<4,15>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<4,15>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConst<4,15>, PSpinVector<T2,4>, OpGammaConstMultiply>::Type_t  d;

  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) = -r.elem(2);
  d.elem(3) = -r.elem(3);
  
  return d;
}


// SpinVector<2> = SpinProject(SpinVector<4>)
// There are 4 cases here for Nd=4 for each forward/backward direction
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir0Minus>::Type_t
spinProjectDir0Minus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir0Minus>::Type_t  d;

  /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
   *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
   *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r + a3i} + i{a0i - a3r} )
   *      ( b1r + i b1i )     ( {a1r + a2i} + i{a1i - a2r} )

   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
   *      ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
   */
  d.elem(0) = s1.elem(0) - timesI(s1.elem(3));
  d.elem(1) = s1.elem(1) - timesI(s1.elem(2));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir1Minus>::Type_t
spinProjectDir1Minus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir1Minus>::Type_t  d;

  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
	 
   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a3r} + i{a0i + a3i} )
   *      ( b1r + i b1i )     ( {a1r - a2r} + i{a1i - a2i} )
      
   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
   *      ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
   */
  d.elem(0) = s1.elem(0) + s1.elem(3);
  d.elem(1) = s1.elem(1) - s1.elem(2);

  return d;
}
    
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir2Minus>::Type_t
spinProjectDir2Minus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir2Minus>::Type_t  d;

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
      
   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *      ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */
  d.elem(0) = s1.elem(0) - timesI(s1.elem(2));
  d.elem(1) = s1.elem(1) + timesI(s1.elem(3));

  return d;
}
    
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir3Minus>::Type_t
spinProjectDir3Minus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir3Minus>::Type_t  d;

  /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
   *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
   *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
      
   * Therefore the top components are
      
   *      ( b0r + i b0i )  =  ( {a0r - a2r} + i{a0i - a2i} )
   *      ( b1r + i b1i )     ( {a1r - a3r} + i{a1i - a3i} )

   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   *      ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */
  d.elem(0) = s1.elem(0) - s1.elem(2);
  d.elem(1) = s1.elem(1) - s1.elem(3);

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir0Plus>::Type_t
spinProjectDir0Plus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir0Plus>::Type_t  d;

  /*                              ( 1  0  0 +i)  ( a0 )    ( a0 + i a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1 +i  0)  ( a1 )  = ( a1 + i a2 )
   *                    0         ( 0 -i  1  0)  ( a2 )    ( a2 - i a1 )
   *                              (-i  0  0  1)  ( a3 )    ( a3 - i a0 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r - a3i} + i{a0i + a3r} )
   *      ( b1r + i b1i )     ( {a1r - a2i} + i{a1i + a2r} )

   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r + a1i} + i{a2i - a1r} )  =  ( b1i - i b1r )
   *      ( b3r + i b3i )     ( {a3r + a0i} + i{a3i - a0r} )     ( b0i - i b0r ) 
   */
  d.elem(0) = s1.elem(0) + timesI(s1.elem(3));
  d.elem(1) = s1.elem(1) + timesI(s1.elem(2));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir1Plus>::Type_t
spinProjectDir1Plus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir1Plus>::Type_t  d;

  /*                              ( 1  0  0 -1)  ( a0 )    ( a0 - a3 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  1  0)  ( a1 )  = ( a1 + a2 )
   *                    1         ( 0  1  1  0)  ( a2 )    ( a2 + a1 )
   *                              (-1  0  0  1)  ( a3 )    ( a3 - a0 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r - a3r} + i{a0i - a3i} )
   *      ( b1r + i b1i )     ( {a1r + a2r} + i{a1i + a2i} )

   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
   *      ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
   */
  d.elem(0) = s1.elem(0) - s1.elem(3);
  d.elem(1) = s1.elem(1) + s1.elem(2);

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir2Plus>::Type_t
spinProjectDir2Plus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir2Plus>::Type_t  d;

  /*                              ( 1  0  i  0)  ( a0 )    ( a0 + i a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0 -i)  ( a1 )  = ( a1 - i a3 )
   *                    2         (-i  0  1  0)  ( a2 )    ( a2 - i a0 )
   *                              ( 0  i  0  1)  ( a3 )    ( a3 + i a1 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r - a2i} + i{a0i + a2r} )
   *      ( b1r + i b1i )     ( {a1r + a3i} + i{a1i - a3r} )

   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *      ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
   */
  d.elem(0) = s1.elem(0) + timesI(s1.elem(2));
  d.elem(1) = s1.elem(1) - timesI(s1.elem(3));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir3Plus>::Type_t
spinProjectDir3Plus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnSpinProjectDir3Plus>::Type_t  d;

  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )

   * Therefore the top components are

   *      ( b0r + i b0i )  =  ( {a0r + a2r} + i{a0i + a2i} )
   *      ( b1r + i b1i )     ( {a1r + a3r} + i{a1i + a3i} )

   * The bottom components of be may be reconstructed using the formula

   *      ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *      ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */
  d.elem(0) = s1.elem(0) + s1.elem(2);
  d.elem(1) = s1.elem(1) + s1.elem(3);

  return d;
}


// SpinVector<4> = SpinReconstruct(SpinVector<2>)
// There are 4 cases here for Nd=4 for each forward/backward direction
template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir0Minus>::Type_t
spinReconstructDir0Minus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir0Minus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = timesI(s1.elem(1));
  d.elem(3) = timesI(s1.elem(0));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir1Minus>::Type_t
spinReconstructDir1Minus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir1Minus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = -s1.elem(1);
  d.elem(3) = s1.elem(0);

  return d;
}


template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir2Minus>::Type_t
spinReconstructDir2Minus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir2Minus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = timesI(s1.elem(0));
  d.elem(3) = timesMinusI(s1.elem(1));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir3Minus>::Type_t
spinReconstructDir3Minus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir3Minus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = -s1.elem(0);
  d.elem(3) = -s1.elem(1);

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir0Plus>::Type_t
spinReconstructDir0Plus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir0Plus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = timesMinusI(s1.elem(1));
  d.elem(3) = timesMinusI(s1.elem(0));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir1Plus>::Type_t
spinReconstructDir1Plus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir1Plus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = s1.elem(1);
  d.elem(3) = -s1.elem(0);

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir2Plus>::Type_t
spinReconstructDir2Plus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir2Plus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = timesMinusI(s1.elem(0));
  d.elem(3) = timesI(s1.elem(1));

  return d;
}

template<class T>
inline typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir3Plus>::Type_t
spinReconstructDir3Plus(const PSpinVector<T,2>& s1)
{
  typename UnaryReturn<PSpinVector<T,2>, FnSpinReconstructDir3Plus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  d.elem(2) = s1.elem(0);
  d.elem(3) = s1.elem(1);

  return d;
}

//-----------------------------------------------

// SpinVector<4> = GammaDP<4,m> * SpinVector<4>
// There are 16 cases here for Nd=4
template<class T2>
inline typename BinaryReturn<GammaConstDP<4,0>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,0>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,0>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;
  
  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) =  r.elem(2);
  d.elem(3) =  r.elem(3);

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,1>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,1>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,1>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;
  
  d.elem(0) = timesMinusI(r.elem(3));
  d.elem(1) = timesMinusI(r.elem(2));
  d.elem(2) = timesI(r.elem(1));
  d.elem(3) = timesI(r.elem(0));

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,2>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,2>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,2>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(3);
  d.elem(1) =  r.elem(2);
  d.elem(2) =  r.elem(1);
  d.elem(3) = -r.elem(0);

  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,3>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,3>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,3>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(0));
  d.elem(1) = timesMinusI(r.elem(1));
  d.elem(2) = timesI(r.elem(2));
  d.elem(3) = timesMinusI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,4>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,4>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,4>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesMinusI(r.elem(2));
  d.elem(1) = timesI(r.elem(3));
  d.elem(2) = timesI(r.elem(0));
  d.elem(3) = timesMinusI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,5>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,5>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,5>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) = -r.elem(3);
  d.elem(3) =  r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,6>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,6>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,6>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(1));
  d.elem(1) = timesI(r.elem(0));
  d.elem(2) = timesI(r.elem(3));
  d.elem(3) = timesI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,7>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,7>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,7>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) =  r.elem(2);
  d.elem(1) =  r.elem(3);
  d.elem(2) = -r.elem(0);
  d.elem(3) = -r.elem(1);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,8>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,8>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,8>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) =  r.elem(0);
  d.elem(1) =  r.elem(1);
  d.elem(2) = -r.elem(2);
  d.elem(3) = -r.elem(3);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,9>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,9>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,9>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(3));
  d.elem(1) = timesI(r.elem(2));
  d.elem(2) = timesI(r.elem(1));
  d.elem(3) = timesI(r.elem(0));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,10>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,10>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,10>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) =  r.elem(3);
  d.elem(1) = -r.elem(2);
  d.elem(2) =  r.elem(1);
  d.elem(3) = -r.elem(0);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,11>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,11>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,11>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(0));
  d.elem(1) = timesMinusI(r.elem(1));
  d.elem(2) = timesMinusI(r.elem(2));
  d.elem(3) = timesI(r.elem(3));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,12>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,12>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,12>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(2));
  d.elem(1) = timesMinusI(r.elem(3));
  d.elem(2) = timesI(r.elem(0));
  d.elem(3) = timesMinusI(r.elem(1));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,13>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,13>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,13>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(1);
  d.elem(1) =  r.elem(0);
  d.elem(2) =  r.elem(3);
  d.elem(3) = -r.elem(2);
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,14>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,14>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,14>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = timesI(r.elem(1));
  d.elem(1) = timesI(r.elem(0));
  d.elem(2) = timesMinusI(r.elem(3));
  d.elem(3) = timesMinusI(r.elem(2));
  
  return d;
}

template<class T2>
inline typename BinaryReturn<GammaConstDP<4,15>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t
operator*(const GammaConstDP<4,15>&, const PSpinVector<T2,4>& r)
{
  typename BinaryReturn<GammaConstDP<4,15>, PSpinVector<T2,4>, OpGammaConstDPMultiply>::Type_t  d;

  d.elem(0) = -r.elem(2);
  d.elem(1) = -r.elem(3);
  d.elem(2) = -r.elem(0);
  d.elem(3) = -r.elem(1);
  
  return d;
}


//-----------------------------------------------------------------------------
//! PSpinVector<T,4> = P_+ * PSpinVector<T,4>
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectPlus>::Type_t
chiralProjectPlus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectPlus>::Type_t  d;

  d.elem(0) = s1.elem(0);
  d.elem(1) = s1.elem(1);
  zero_rep(d.elem(2));
  zero_rep(d.elem(3));

  return d;
}

//! PSpinVector<T,4> = P_- * PSpinVector<T,4>
template<class T>
inline typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectMinus>::Type_t
chiralProjectMinus(const PSpinVector<T,4>& s1)
{
  typename UnaryReturn<PSpinVector<T,4>, FnChiralProjectMinus>::Type_t  d;

  zero_rep(d.elem(0));
  zero_rep(d.elem(1));
  d.elem(2) = s1.elem(2);
  d.elem(3) = s1.elem(3);

  return d;
}


/*! @} */   // end of group primspinvector

} // namespace QDP

#endif
