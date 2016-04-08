// -*- C++ -*-

/*! \file
 * \brief Primitive Scalar
 */

#ifndef QDP_PRIMSCALAR_H
#define QDP_PRIMSCALAR_H

namespace QDP {


//-------------------------------------------------------------------------------------
/*! \addtogroup primscalar Scalar primitive
 * \ingroup fiber
 *
 * Primitive Scalar is a placeholder for no primitive structure
 *
 * @{
 */

//! Primitive Scalar
/*! Placeholder for no primitive structure */
template<class T> class PScalar
{
public:
  PScalar() {}
  ~PScalar() {}

  //---------------------------------------------------------
  //! construct dest = const
  PScalar(const typename WordType<T>::Type_t& rhs) : F(rhs) {}

  //! construct dest = rhs
  template<class T1>
  PScalar(const PScalar<T1>& rhs) : F(rhs.elem()) {}

  //! construct dest = rhs
  template<class T1>
  PScalar(const T1& rhs) : F(rhs) {}

  //---------------------------------------------------------
#if 0
  //! dest = const
  /*! Fill with an integer constant. Will be promoted to underlying word type */
  inline
  PScalar& operator=(const typename WordType<T>::Type_t& rhs)
    {
      elem() = rhs;
      return *this;
    }
#endif

  //! PScalar = PScalar
  /*! Set equal to another PScalar */
  template<class T1>
  inline
  PScalar& operator=(const PScalar<T1>& rhs) 
    {
      elem() = rhs.elem();
      return *this;
    }

  //! PScalar += PScalar
  template<class T1>
  inline
  PScalar& operator+=(const PScalar<T1>& rhs) 
    {
      elem() += rhs.elem();
      return *this;
    }

  //! PScalar -= PScalar
  template<class T1>
  inline
  PScalar& operator-=(const PScalar<T1>& rhs) 
    {
      elem() -= rhs.elem();
      return *this;
    }

  //! PScalar *= PScalar
  template<class T1>
  inline
  PScalar& operator*=(const PScalar<T1>& rhs) 
    {
      elem() *= rhs.elem();
      return *this;
    }

  //! PScalar /= PScalar
  template<class T1>
  inline
  PScalar& operator/=(const PScalar<T1>& rhs) 
    {
      elem() /= rhs.elem();
      return *this;
    }

  //! PScalar %= PScalar
  template<class T1>
  inline
  PScalar& operator%=(const PScalar<T1>& rhs) 
    {
      elem() %= rhs.elem();
      return *this;
    }

  //! PScalar |= PScalar
  template<class T1>
  inline
  PScalar& operator|=(const PScalar<T1>& rhs) 
    {
      elem() |= rhs.elem();
      return *this;
    }

  //! PScalar &= PScalar
  template<class T1>
  inline
  PScalar& operator&=(const PScalar<T1>& rhs) 
    {
      elem() &= rhs.elem();
      return *this;
    }

  //! PScalar ^= PScalar
  template<class T1>
  inline
  PScalar& operator^=(const PScalar<T1>& rhs) 
    {
      elem() ^= rhs.elem();
      return *this;
    }

  //! PScalar <<= PScalar
  template<class T1>
  inline
  PScalar& operator<<=(const PScalar<T1>& rhs) 
    {
      elem() <<= rhs.elem();
      return *this;
    }

  //! PScalar >>= PScalar
  template<class T1>
  inline
  PScalar& operator>>=(const PScalar<T1>& rhs) 
    {
      elem() >>= rhs.elem();
      return *this;
    }

  //! Deep copies here
  PScalar(const PScalar& a): F(a.F) {/* fprintf(stderr,"copy PScalar\n"); */}

public:
  T& elem() {return F;}
  const T& elem() const {return F;}

private:
  T F;
};


// Input
//! Ascii input
template<class T>
inline
std::istream& operator>>(std::istream& s, PScalar<T>& d)
{
  return s >> d.elem();
}

//! Ascii input
template<class T>
inline
StandardInputStream& operator>>(StandardInputStream& s, PScalar<T>& d)
{
  return s >> d.elem();
}

// Output
//! Ascii output
template<class T>
inline
std::ostream& operator<<(std::ostream& s, const PScalar<T>& d)
{
  return s << d.elem();
}

//! Ascii output
template<class T>
inline
StandardOutputStream& operator<<(StandardOutputStream& s, const PScalar<T>& d)
{
  return s << d.elem();
}

//! Text input
template<class T>
inline
TextReader& operator>>(TextReader& txt, PScalar<T>& d)
{
  return txt >> d.elem();
}

//! Text output
template<class T>
inline
TextWriter& operator<<(TextWriter& txt, const PScalar<T>& d)
{
  return txt << d.elem();
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T>
inline
XMLWriter& operator<<(XMLWriter& xml, const PScalar<T>& d)
{
  return xml << d.elem();
}

//! XML input
template<class T>
inline
void read(XMLReader& xml, const std::string& path, PScalar<T>& d)
{
  read(xml, path, d.elem());
}
#endif

/*! @} */  // end of group primscalar


//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T>
struct WordType<PScalar<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};

// Fixed Precision Types 
template<class T>
struct SinglePrecType<PScalar<T> >
{
  typedef PScalar< typename SinglePrecType<T>::Type_t > Type_t;
};

template<class T>
struct DoublePrecType<PScalar<T> >
{
  typedef PScalar< typename DoublePrecType<T>::Type_t > Type_t;
};

// Internally used scalars
template<class T>
struct InternalScalar<PScalar<T> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Internally used real scalars
template<class T>
struct RealScalar<PScalar<T> > {
  typedef PScalar<typename RealScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar<PScalar<T> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T>
struct LatticeScalar<PScalar<T> > {
  typedef PScalar<typename LatticeScalar<T>::Type_t>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(PScalar) -> PScalar
template<class T1, class Op>
struct UnaryReturn<PScalar<T1>, Op> {
  typedef PScalar<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default binary(PScalar,PScalar) -> PScalar
template<class T1, class T2, class Op>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, Op> {
  typedef PScalar<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};


#if 0
template<class T1, class T2>
struct UnaryReturn<PScalar<T2>, OpCast<T1> > {
  typedef PScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif

// Assignment is different
template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpAssign > {
  typedef PScalar<T1> &Type_t;
};

template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpAddAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpSubtractAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpMultiplyAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpDivideAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpModAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpBitwiseOrAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpBitwiseAndAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpBitwiseXorAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpLeftShiftAssign > {
  typedef PScalar<T1> &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpRightShiftAssign > {
  typedef PScalar<T1> &Type_t;
};
 



//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup primscalar */
/*! @{ */

// Primitive Scalars

// ! PScalar
template<class T>
struct UnaryReturn<PScalar<T>, OpNot > {
  typedef PScalar<typename UnaryReturn<T, OpNot>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<PScalar<T1>, OpNot>::Type_t
operator!(const PScalar<T1>& l)
{
  return ! l.elem();
}

// + PScalar
template<class T1>
inline typename UnaryReturn<PScalar<T1>, OpUnaryPlus>::Type_t
operator+(const PScalar<T1>& l)
{
  return +l.elem();
}

// - PScalar
template<class T1>
inline typename UnaryReturn<PScalar<T1>, OpUnaryMinus>::Type_t
operator-(const PScalar<T1>& l)
{
  return -l.elem();
}

// PScalar + PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpAdd>::Type_t
operator+(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() + r.elem();
}

// PScalar - PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpSubtract>::Type_t
operator-(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() - r.elem();
}

// PScalar * PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpMultiply>::Type_t
operator*(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() * r.elem();
}

// Optimized  adj(PMatrix)*PMatrix
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return adjMultiply(l.elem(), r.elem());
}

// Optimized  PMatrix*adj(PMatrix)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return multiplyAdj(l.elem(), r.elem());
}

// Optimized  PMatrix*adj(PMatrix)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PSpinVector<T2,4>, OpMultiplyAdj>::Type_t
multiplyAdj(const PScalar<T1>& l, const PSpinVector<T2,4>& r)
{
  return multiplyAdj(l.elem(), r.elem());
}

// Optimized  adj(PMatrix)*adj(PMatrix)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return adjMultiplyAdj(l.elem(), r.elem());
}

// PScalar / PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpDivide>::Type_t
operator/(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() / r.elem();
}


// PScalar << PScalar
template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpLeftShift > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpLeftShift>::Type_t
operator<<(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() << r.elem();
}

// PScalar >> PScalar
template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpRightShift > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpRightShift>::Type_t
operator>>(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() >> r.elem();
}

// PScalar % PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpMod>::Type_t
operator%(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() % r.elem();
}

// PScalar ^ PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpBitwiseXor>::Type_t
operator^(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() ^ r.elem();
}

// PScalar & PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() & r.elem();
}

// PScalar | PScalar
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpBitwiseOr>::Type_t
operator|(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() | r.elem();
}


// Comparisons
template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpLT > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpLT>::Type_t
operator<(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() < r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpLE > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpLE>::Type_t
operator<=(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() <= r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpGT > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpGT>::Type_t
operator>(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() > r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpGE > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpGE>::Type_t
operator>=(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() >= r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpEQ > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpEQ>::Type_t
operator==(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() == r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpNE > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpNE>::Type_t
operator!=(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() != r.elem();
}


template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpAnd > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpAnd>::Type_t
operator&&(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() && r.elem();
}


template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, OpOr > {
  typedef PScalar<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, OpOr>::Type_t
operator||(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return l.elem() || r.elem();
}


//-----------------------------------------------------------------------------
// Functions

// Adjoint
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnAdjoint>::Type_t
adj(const PScalar<T1>& s1)
{
  return adj(s1.elem());
}


// Conjugate
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnConjugate>::Type_t
conj(const PScalar<T1>& s1)
{
  return conj(s1.elem());
}


// Transpose
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnTranspose>::Type_t
transpose(const PScalar<T1>& s1)
{
  return transpose(s1.elem());
}


// TRACE
// trace = Trace(source1)
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnTrace>::Type_t
trace(const PScalar<T1>& s1)
{
  return trace(s1.elem());
}


// trace = Re(Trace(source1))
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnRealTrace>::Type_t
realTrace(const PScalar<T1>& s1)
{
  return realTrace(s1.elem());
}


// trace = Im(Trace(source1))
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnImagTrace>::Type_t
imagTrace(const PScalar<T1>& s1)
{
  return imagTrace(s1.elem());
}


// trace = colorTrace(source1)
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnTraceColor>::Type_t
traceColor(const PScalar<T1>& s1)
{
  return traceColor(s1.elem());
}


//! PScalar = traceSpin(PScalar)
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnTraceSpin>::Type_t
traceSpin(const PScalar<T1>& s1)
{
  return traceSpin(s1.elem());
}

//! PScalar = transposeSpin(PScalar)
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnTransposeSpin>::Type_t
transposeSpin(const PScalar<T1>& s1)
{
  return transposeSpin(s1.elem());
}

//! PScalar = trace(PScalar * PScalar)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnTraceMultiply>::Type_t
traceMultiply(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return traceMultiply(l.elem(), r.elem());
}

//! PScalar = traceColor(PScalar * PScalar)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnTraceColorMultiply>::Type_t
traceColorMultiply(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return traceMultiply(l.elem(), r.elem());
}

//! PScalar = traceSpin(PScalar * PScalar)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnTraceSpinMultiply>::Type_t
traceSpinMultiply(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return traceMultiply(l.elem(), r.elem());
}

//! PScalar = traceSpin(outerProduct(PScalar, PScalar))
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnTraceSpinOuterProduct>::Type_t
traceSpinOuterProduct(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return traceSpinOuterProduct(l.elem(), r.elem());
}

//! PScalar = outerProduct(PScalar, PScalar)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const PScalar<T1>& l, const PScalar<T2>& r)
{
  return outerProduct(l.elem(),r.elem());
}


//! PScalar = Re(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnReal>::Type_t
real(const PScalar<T>& s1)
{
  return real(s1.elem());
}


// PScalar = Im(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnImag>::Type_t
imag(const PScalar<T>& s1)
{
  return imag(s1.elem());
}


// ArcCos
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnArcCos>::Type_t
acos(const PScalar<T1>& s1)
{
  return acos(s1.elem());
}

// ArcSin
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnArcSin>::Type_t
asin(const PScalar<T1>& s1)
{
  return asin(s1.elem());
}

// ArcTan
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnArcTan>::Type_t
atan(const PScalar<T1>& s1)
{
  return atan(s1.elem());
}

// Ceil(ing)
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnCeil>::Type_t
ceil(const PScalar<T1>& s1)
{
  return ceil(s1.elem());
}

// Cos
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnCos>::Type_t
cos(const PScalar<T1>& s1)
{
  return cos(s1.elem());
}

// Cosh
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnHypCos>::Type_t
cosh(const PScalar<T1>& s1)
{
  return cosh(s1.elem());
}

// Exp
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnExp>::Type_t
exp(const PScalar<T1>& s1)
{
  return exp(s1.elem());
}

// Fabs
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnFabs>::Type_t
fabs(const PScalar<T1>& s1)
{
  return fabs(s1.elem());
}

// Floor
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnFloor>::Type_t
floor(const PScalar<T1>& s1)
{
  return floor(s1.elem());
}

// Log
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnLog>::Type_t
log(const PScalar<T1>& s1)
{
  return log(s1.elem());
}

// Log10
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnLog10>::Type_t
log10(const PScalar<T1>& s1)
{
  return log10(s1.elem());
}

// Sin
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnSin>::Type_t
sin(const PScalar<T1>& s1)
{
  return sin(s1.elem());
}

// Sinh
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnHypSin>::Type_t
sinh(const PScalar<T1>& s1)
{
  return sinh(s1.elem());
}

// Sqrt
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnSqrt>::Type_t
sqrt(const PScalar<T1>& s1)
{
  return sqrt(s1.elem());
}

// Tan
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnTan>::Type_t
tan(const PScalar<T1>& s1)
{
  return tan(s1.elem());
}

// Tanh
template<class T1>
inline typename UnaryReturn<PScalar<T1>, FnHypTan>::Type_t
tanh(const PScalar<T1>& s1)
{
  return tanh(s1.elem());
}


//-----------------------------------------------------------------------------
// These functions always return bool
//! isnan
template<class T1>
struct UnaryReturn<PScalar<T1>, FnIsNan> {
  bool Type_t;
};

template<class T1>
inline bool
isnan(const PScalar<T1>& s1)
{
  return isnan(s1.elem());
}

//! isinf
template<class T1>
struct UnaryReturn<PScalar<T1>, FnIsInf> {
  bool Type_t;
};

template<class T1>
inline bool
isinf(const PScalar<T1>& s1)
{
  return isinf(s1.elem());
}

//! isnormal
template<class T1>
struct UnaryReturn<PScalar<T1>, FnIsNormal> {
  bool Type_t;
};

template<class T1>
inline bool
isnormal(const PScalar<T1>& s1)
{
  return isnormal(s1.elem());
}

//! isfinite
template<class T1>
struct UnaryReturn<PScalar<T1>, FnIsFinite> {
  bool Type_t;
};

template<class T1>
inline bool
isfinite(const PScalar<T1>& s1)
{
  return isfinite(s1.elem());
}


//-----------------------------------------------------------------------------
//! PScalar<T> = pow(PScalar<T> , PScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnPow>::Type_t
pow(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return pow(s1.elem(), s2.elem());
}

//! PScalar<T> = atan2(PScalar<T> , PScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnArcTan2>::Type_t
atan2(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return atan2(s1.elem(), s2.elem());
}


//! PScalar<T> = (PScalar<T> , PScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnCmplx>::Type_t
cmplx(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return cmplx(s1.elem(), s2.elem());
}



// Global Functions
// PScalar = i * PScalar
template<class T>
inline typename UnaryReturn<PScalar<T>, FnTimesI>::Type_t
timesI(const PScalar<T>& s1)
{
  return timesI(s1.elem());
}

// PScalar = -i * PScalar
template<class T>
inline typename UnaryReturn<PScalar<T>, FnTimesMinusI>::Type_t
timesMinusI(const PScalar<T>& s1)
{
  return timesMinusI(s1.elem());
}


//! dest [float type] = source [seed type]
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSeedToFloat>::Type_t
seedToFloat(const PScalar<T>& s1)
{
  return seedToFloat(s1.elem());
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T>
inline typename UnaryReturn<PScalar<T>, FnGetSite>::Type_t
getSite(const PScalar<T>& s1, int innersite)
{
  return getSite(s1.elem(), innersite);
}

//! Extract color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T>
inline typename UnaryReturn<PScalar<T>, FnPeekColorVector>::Type_t
peekColor(const PScalar<T>& l, int row)
{
  return peekColor(l.elem(),row);
}

//! Extract color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T>
inline typename UnaryReturn<PScalar<T>, FnPeekColorMatrix>::Type_t
peekColor(const PScalar<T>& l, int row, int col)
{
  return peekColor(l.elem(),row,col);
}

//! Extract spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<PScalar<T>, FnPeekSpinVector>::Type_t
peekSpin(const PScalar<T>& l, int row)
{
  return peekSpin(l.elem(),row);
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<PScalar<T>, FnPeekSpinMatrix>::Type_t
peekSpin(const PScalar<T>& l, int row, int col)
{
  return peekSpin(l.elem(),row,col);
}


//! Insert color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2>
inline PScalar<T1>&
pokeColor(PScalar<T1>& l, const PScalar<T2>& r, int row)
{
  pokeColor(l.elem(),r.elem(),row);
  return l;
}

//! Insert color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T1, class T2>
inline PScalar<T1>&
pokeColor(PScalar<T1>& l, const PScalar<T2>& r, int row, int col)
{
  pokeColor(l.elem(),r.elem(),row,col);
  return l;
}

//! Insert spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline PScalar<T1>&
pokeSpin(PScalar<T1>& l, const PScalar<T2>& r, int row)
{
  pokeSpin(l.elem(),r.elem(),row);
  return l;
}

//! Insert spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T1, class T2>
inline PScalar<T1>&
pokeSpin(PScalar<T1>& l, const PScalar<T2>& r, int row, int col)
{
  pokeSpin(l.elem(),r.elem(),row,col);
  return l;
}


//-----------------------------------------------------------------------------
//! PScalar = Gamma<N,m> * PScalar
template<class T2, int N, int m>
inline typename BinaryReturn<GammaConst<N,m>, PScalar<T2>, OpGammaConstMultiply>::Type_t
operator*(const GammaConst<N,m>& l, const PScalar<T2>& r)
{
  return l * r.elem();
}

//! PScalar = PScalar * Gamma<N,m>
template<class T2, int N, int m>
inline typename BinaryReturn<PScalar<T2>, GammaConst<N,m>, OpGammaConstMultiply>::Type_t
operator*(const PScalar<T2>& l, const GammaConst<N,m>& r)
{
  return l.elem() * r;
}

//-----------------------------------------------------------------------------
//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir0Minus>::Type_t
spinProjectDir0Minus(const PScalar<T>& s1)
{
  return spinProjectDir0Minus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir0Minus>::Type_t
spinReconstructDir0Minus(const PScalar<T>& s1)
{
  return spinReconstructDir0Minus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir1Minus>::Type_t
spinProjectDir1Minus(const PScalar<T>& s1)
{
  return spinProjectDir1Minus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir1Minus>::Type_t
spinReconstructDir1Minus(const PScalar<T>& s1)
{
  return spinReconstructDir1Minus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir2Minus>::Type_t
spinProjectDir2Minus(const PScalar<T>& s1)
{
  return spinProjectDir2Minus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir2Minus>::Type_t
spinReconstructDir2Minus(const PScalar<T>& s1)
{
  return spinReconstructDir2Minus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir3Minus>::Type_t
spinProjectDir3Minus(const PScalar<T>& s1)
{
  return spinProjectDir3Minus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir3Minus>::Type_t
spinReconstructDir3Minus(const PScalar<T>& s1)
{
  return spinReconstructDir3Minus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir0Plus>::Type_t
spinProjectDir0Plus(const PScalar<T>& s1)
{
  return spinProjectDir0Plus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir0Plus>::Type_t
spinReconstructDir0Plus(const PScalar<T>& s1)
{
  return spinReconstructDir0Plus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir1Plus>::Type_t
spinProjectDir1Plus(const PScalar<T>& s1)
{
  return spinProjectDir1Plus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir1Plus>::Type_t
spinReconstructDir1Plus(const PScalar<T>& s1)
{
  return spinReconstructDir1Plus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir2Plus>::Type_t
spinProjectDir2Plus(const PScalar<T>& s1)
{
  return spinProjectDir2Plus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir2Plus>::Type_t
spinReconstructDir2Plus(const PScalar<T>& s1)
{
  return spinReconstructDir2Plus(s1.elem());
}


//! PScalar = SpinProject(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinProjectDir3Plus>::Type_t
spinProjectDir3Plus(const PScalar<T>& s1)
{
  return spinProjectDir3Plus(s1.elem());
}

//! PScalar = SpinReconstruct(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnSpinReconstructDir3Plus>::Type_t
spinReconstructDir3Plus(const PScalar<T>& s1)
{
  return spinReconstructDir3Plus(s1.elem());
}

//-----------------------------------------------------------------------------
//! PScalar = chiralProjectPlus(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnChiralProjectPlus>::Type_t
chiralProjectPlus(const PScalar<T>& s1)
{
  return chiralProjectPlus(s1.elem());
}

//! PScalar = chiralProjectMinus(PScalar)
template<class T>
inline typename UnaryReturn<PScalar<T>, FnChiralProjectMinus>::Type_t
chiralProjectMinus(const PScalar<T>& s1)
{
  return chiralProjectMinus(s1.elem());
}


//-----------------------------------------------------------------------------
// quark propagator contraction
template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnQuarkContract13>::Type_t
quarkContract13(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return quarkContract13(s1.elem(), s2.elem());
}

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnQuarkContract14>::Type_t
quarkContract14(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return quarkContract14(s1.elem(), s2.elem());
}

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnQuarkContract23>::Type_t
quarkContract23(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return quarkContract23(s1.elem(), s2.elem());
}

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnQuarkContract24>::Type_t
quarkContract24(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return quarkContract24(s1.elem(), s2.elem());
}

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnQuarkContract12>::Type_t
quarkContract12(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return quarkContract12(s1.elem(), s2.elem());
}

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnQuarkContract34>::Type_t
quarkContract34(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return quarkContract34(s1.elem(), s2.elem());
}


//-----------------------------------------------------------------------------
// Contraction for color matrices
// colorContract 
//! dest  = colorContract(Qprop1,Qprop2,Qprop3)
/*!
 * This routine is completely unrolled for 3 colors
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<PScalar<T1>, PScalar<T2>, PScalar<T3>, FnColorContract> {
  typedef PScalar<typename TrinaryReturn<T1, T2, T3, FnColorContract>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<PScalar<T1>, PScalar<T2>, PScalar<T3>, FnColorContract>::Type_t
colorContract(const PScalar<T1>& s1, const PScalar<T2>& s2, const PScalar<T3>& s3)
{
  return colorContract(s1.elem(), s2.elem(), s3.elem());
}


//-----------------------------------------------------------------------------
// Contraction of two colorvectors
//! dest  = colorVectorContract(Qvec1,Qvec2)
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, FnColorVectorContract> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnColorVectorContract>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnColorVectorContract>::Type_t
colorVectorContract(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return colorVectorContract(s1.elem(), s2.elem());
}



//-----------------------------------------------------------------------------
// Cross product for color vectors
//! dest  = colorCrossProduct(Qvec1,Qvec2)
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, FnColorCrossProduct> {
  typedef PScalar<typename BinaryReturn<T1, T2, FnColorCrossProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnColorCrossProduct>::Type_t
colorCrossProduct(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return colorCrossProduct(s1.elem(), s2.elem());
}



//-----------------------------------------------------------------------------
//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline void 
copymask(PScalar<T>& d, const PScalar<T1>& mask, const PScalar<T>& s1) 
{
  copymask(d.elem(),mask.elem(),s1.elem());
}

//! dest  = random  
template<class T, class T1, class T2>
inline void
fill_random(PScalar<T>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  fill_random(d.elem(), seed, skewed_seed, seed_mult);
}


//! dest  = gaussian  
template<class T>
inline void
fill_gaussian(PScalar<T>& d, PScalar<T>& r1, PScalar<T>& r2)
{
  fill_gaussian(d.elem(), r1.elem(), r2.elem());
}


#if 1
// Global sum over site indices only
template<class T>
struct UnaryReturn<PScalar<T>, FnSum > {
  typedef PScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<PScalar<T>, FnSum>::Type_t
sum(const PScalar<T>& s1)
{
  return sum(s1.elem());
}
#endif


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<PScalar<T>, FnNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<PScalar<T>, FnLocalNorm2 > {
  typedef PScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<PScalar<T>, FnLocalNorm2>::Type_t
localNorm2(const PScalar<T>& s1)
{
  return localNorm2(s1.elem());
}

// Global max
template<class T>
struct UnaryReturn<PScalar<T>, FnGlobalMax> {
  typedef PScalar<typename UnaryReturn<T, FnGlobalMax>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<PScalar<T>, FnGlobalMax>::Type_t
globalMax(const PScalar<T>& s1)
{
  return globalMax(s1.elem());
}


// Global min
template<class T>
struct UnaryReturn<PScalar<T>, FnGlobalMin> {
  typedef PScalar<typename UnaryReturn<T, FnGlobalMin>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<PScalar<T>, FnGlobalMin>::Type_t
globalMin(const PScalar<T>& s1)
{
  return globalMin(s1.elem());
}


//! PScalar<T> = InnerProduct(adj(PScalar<T1>)*PScalar<T2>)
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, FnInnerProduct > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, FnLocalInnerProduct > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return localInnerProduct(s1.elem(), s2.elem());
}


//! PScalar<T> = InnerProductReal(adj(PMatrix<T1>)*PMatrix<T1>)
template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, FnInnerProductReal > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<PScalar<T1>, PScalar<T2>, FnLocalInnerProductReal > {
  typedef PScalar<typename BinaryReturn<T1, T2, FnLocalInnerProductReal>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PScalar<T1>, PScalar<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const PScalar<T1>& s1, const PScalar<T2>& s2)
{
  return localInnerProductReal(s1.elem(), s2.elem());
}


//! PScalar<T> = where(PScalar, PScalar, PScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<PScalar<T1>, PScalar<T2>, PScalar<T3>, FnWhere> {
  typedef PScalar<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<PScalar<T1>, PScalar<T2>, PScalar<T3>, FnWhere>::Type_t
where(const PScalar<T1>& a, const PScalar<T2>& b, const PScalar<T3>& c)
{
  return where(a.elem(), b.elem(), c.elem());
}


//-----------------------------------------------------------------------------
//! QDP Int to int primitive in conversion routine
template<class T> 
inline int 
toInt(const PScalar<T>& s) 
{
  return toInt(s.elem());
}

//! QDP Real to float primitive in conversion routine
template<class T> 
inline float
toFloat(const PScalar<T>& s) 
{
  return toFloat(s.elem());
}

//! QDP Double to double primitive in conversion routine
template<class T> 
inline double
toDouble(const PScalar<T>& s) 
{
  return toDouble(s.elem());
}

//! QDP Boolean to bool primitive in conversion routine
template<class T> 
inline bool
toBool(const PScalar<T>& s) 
{
  return toBool(s.elem());
}

//! QDP Wordtype to primitive wordtype
template<class T> 
inline typename WordType< PScalar<T> >::Type_t
toWordType(const PScalar<T>& s) 
{
  return toWordType(s.elem());
}


//-----------------------------------------------------------------------------
// Other operations
//! dest = 0
template<class T> 
inline void 
zero_rep(PScalar<T>& dest) 
{
  zero_rep(dest.elem());
}

//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
cast_rep(T& d, const PScalar<T1>& s1)
{
  cast_rep(d, s1.elem());
}

//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
cast_rep(PScalar<T>& d, const PScalar<T1>& s1)
{
  cast_rep(d.elem(), s1.elem());
}

//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
copy_site(PScalar<T>& d, int isite, const PScalar<T1>& s1)
{
  copy_site(d.elem(), isite, s1.elem());
}

//! gather several inner sites together
template<class T, class T1>
inline void 
gather_sites(PScalar<T>& d, 
	     const PScalar<T1>& s0, int i0, 
	     const PScalar<T1>& s1, int i1,
	     const PScalar<T1>& s2, int i2,
	     const PScalar<T1>& s3, int i3)
{
  gather_sites(d.elem(), s0.elem(), i0, s1.elem(), i1, s2.elem(), i2, s3.elem(), i3);
}

/*! @} */  // end of group primscalar

} // namespace QDP

#endif
