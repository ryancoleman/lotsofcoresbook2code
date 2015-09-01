// -*- C++ -*-

/*! \file
 * \brief Reality
 */


#ifndef QDP_REALITY_H
#define QDP_REALITY_H

#include <sstream>

namespace QDP {


//-------------------------------------------------------------------------------------
/*! \addtogroup rscalar Scalar reality
 * \ingroup fiber
 *
 * Reality Scalar is a type for objects that are only real - no imaginary part
 *
 * @{
 */

//! Scalar reality (not complex)
template<class T> class RScalar
{
public:
  RScalar() {}
  ~RScalar() {}

  //---------------------------------------------------------
  //! construct dest = const
  RScalar(const typename WordType<T>::Type_t& rhs) : F(rhs) {}

  //! construct dest = rhs
  template<class T1>
  RScalar(const RScalar<T1>& rhs) : F(rhs.elem()) {}

  //! construct dest = rhs
  template<class T1>
  RScalar(const T1& rhs) : F(rhs) {}

  //---------------------------------------------------------
#if 0
  //! dest = const
  /*! Fill with a constant. Will be promoted to underlying word type */
  inline
  RScalar& operator=(const typename WordType<T>::Type_t& rhs)
    {
      elem() = rhs;
      return *this;
    }
#endif

  //! RScalar = RScalar
  /*! Set equal to another RScalar */
  template<class T1>
  inline
  RScalar& operator=(const RScalar<T1>& rhs) 
    {
      elem() = rhs.elem();
      return *this;
    }

  //! RScalar += RScalar
  template<class T1>
  inline
  RScalar& operator+=(const RScalar<T1>& rhs) 
    {
      elem() += rhs.elem();
      return *this;
    }

  //! RScalar -= RScalar
  template<class T1>
  inline
  RScalar& operator-=(const RScalar<T1>& rhs) 
    {
      elem() -= rhs.elem();
      return *this;
    }

  //! RScalar *= RScalar
  template<class T1>
  inline
  RScalar& operator*=(const RScalar<T1>& rhs) 
    {
      elem() *= rhs.elem();
      return *this;
    }

  //! RScalar /= RScalar
  template<class T1>
  inline
  RScalar& operator/=(const RScalar<T1>& rhs) 
    {
      elem() /= rhs.elem();
      return *this;
    }

  //! RScalar %= RScalar
  template<class T1>
  inline
  RScalar& operator%=(const RScalar<T1>& rhs) 
    {
      elem() %= rhs.elem();
      return *this;
    }

  //! RScalar |= RScalar
  template<class T1>
  inline
  RScalar& operator|=(const RScalar<T1>& rhs) 
    {
      elem() |= rhs.elem();
      return *this;
    }

  //! RScalar &= RScalar
  template<class T1>
  inline
  RScalar& operator&=(const RScalar<T1>& rhs) 
    {
      elem() &= rhs.elem();
      return *this;
    }

  //! RScalar ^= RScalar
  template<class T1>
  inline
  RScalar& operator^=(const RScalar<T1>& rhs) 
    {
      elem() ^= rhs.elem();
      return *this;
    }

  //! RScalar <<= RScalar
  template<class T1>
  inline
  RScalar& operator<<=(const RScalar<T1>& rhs) 
    {
      elem() <<= rhs.elem();
      return *this;
    }

  //! RScalar >>= RScalar
  template<class T1>
  inline
  RScalar& operator>>=(const RScalar<T1>& rhs) 
    {
      elem() >>= rhs.elem();
      return *this;
    }


  //! Do deep copies here
  RScalar(const RScalar& a): F(a.F) {}

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
std::istream& operator>>(std::istream& s, RScalar<T>& d)
{
  return s >> d.elem();
}

//! Ascii input
template<class T>
inline
StandardInputStream& operator>>(StandardInputStream& s, RScalar<T>& d)
{
  return s >> d.elem();
}

//! Ascii output
template<class T> 
inline  
std::ostream& operator<<(std::ostream& s, const RScalar<T>& d)
{
  return s << d.elem();
}

//! Ascii output
template<class T> 
inline  
StandardOutputStream& operator<<(StandardOutputStream& s, const RScalar<T>& d)
{
  return s << d.elem();
}


//! Text input
template<class T>
inline
TextReader& operator>>(TextReader& s, RScalar<T>& d)
{
  return s >> d.elem();
}

//! Text output
template<class T> 
inline  
TextWriter& operator<<(TextWriter& s, const RScalar<T>& d)
{
  return s << d.elem();
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T>
inline
XMLWriter& operator<<(XMLWriter& xml, const RScalar<T>& d)
{
  return xml << d.elem();
}

//! XML input
template<class T>
inline
void read(XMLReader& xml, const std::string& path, RScalar<T>& d)
{
  read(xml, path, d.elem());
}
#endif

/*! @} */  // end of group rscalar


//-------------------------------------------------------------------------------------
/*! \addtogroup rcomplex Complex reality
 * \ingroup fiber
 *
 * Reality Complex is a type for objects that hold a real and imaginary part
 *
 * @{
 */

//! Reality complex
/*! All fields are either complex or scalar reality */
template<class T> class RComplex
{
public:
  RComplex() {}
  ~RComplex() {}

  //! Construct from two reality scalars
  template<class T1, class T2>
  RComplex(const RScalar<T1>& _re, const RScalar<T2>& _im): re(_re.elem()), im(_im.elem()) {}

  //! Construct from two scalars
  template<class T1, class T2>
  RComplex(const T1& _re, const T2& _im): re(_re), im(_im) {}

  //---------------------------------------------------------
  //! RComplex = RScalar
  /*! Set the real part and zero the imag part */
  template<class T1>
  inline
  RComplex& operator=(const RScalar<T1>& rhs) 
    {
      real() = rhs.elem();
      zero_rep(imag());
      return *this;
    }

  //! RComplex += RScalar
  template<class T1>
  inline
  RComplex& operator+=(const RScalar<T1>& rhs) 
    {
      real() += rhs.elem();
      return *this;
    }

  //! RComplex -= RScalar
  template<class T1>
  inline
  RComplex& operator-=(const RScalar<T1>& rhs) 
    {
      real() -= rhs.elem();
      return *this;
    }

  //! RComplex *= RScalar
  template<class T1>
  inline
  RComplex& operator*=(const RScalar<T1>& rhs) 
    {
      real() *= rhs.elem();
      imag() *= rhs.elem();
      return *this;
    }

  //! RComplex /= RScalar
  template<class T1>
  inline
  RComplex& operator/=(const RScalar<T1>& rhs) 
    {
      real() /= rhs.elem();
      imag() /= rhs.elem();
      return *this;
    }



  //! RComplex = RComplex
  /*! Set equal to another RComplex */
  template<class T1>
  inline
  RComplex& operator=(const RComplex<T1>& rhs) 
    {
      real() = rhs.real();
      imag() = rhs.imag();
      return *this;
    }

  //! RComplex += RComplex
  template<class T1>
  inline
  RComplex& operator+=(const RComplex<T1>& rhs) 
    {
      real() += rhs.real();
      imag() += rhs.imag();
      return *this;
    }

  //! RComplex -= RComplex
  template<class T1>
  inline
  RComplex& operator-=(const RComplex<T1>& rhs) 
    {
      real() -= rhs.real();
      imag() -= rhs.imag();
      return *this;
    }

  //! RComplex *= RComplex
  template<class T1>
  inline
  RComplex& operator*=(const RComplex<T1>& rhs) 
    {
      RComplex<T> d;
      d = *this * rhs;

      real() = d.real();
      imag() = d.imag();
      return *this;
    }

  //! RComplex /= RComplex
  template<class T1>
  inline
  RComplex& operator/=(const RComplex<T1>& rhs) 
    {
      RComplex<T> d;
      d = *this / rhs;

      real() = d.real();
      imag() = d.imag();
      return *this;
    }


  //! Deep copy constructor
  RComplex(const RComplex& a): re(a.re), im(a.im) {}

public:
  T& real() {return re;}
  const T& real() const {return re;}

  T& imag() {return im;}
  const T& imag() const {return im;}

private:
  T re;
  T im;
} QDP_ALIGN8;   // possibly force alignment


//! Stream output
template<class T>
inline
std::ostream& operator<<(std::ostream& s, const RComplex<T>& d)
{
  s << "( " << d.real() << " , " << d.imag() << " )";
  return s;
}

//! Stream output
template<class T>
inline
StandardOutputStream& operator<<(StandardOutputStream& s, const RComplex<T>& d)
{
  s << "( " << d.real() << " , " << d.imag() << " )";
  return s;
}

//! Text input
template<class T>
inline
TextReader& operator>>(TextReader& s, RComplex<T>& d)
{
  return s >> d.real() >> d.imag();
}

//! Text output
template<class T> 
inline  
TextWriter& operator<<(TextWriter& s, const RComplex<T>& d)
{
  return s << d.real() << d.imag();
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T>
inline
XMLWriter& operator<<(XMLWriter& xml, const RComplex<T>& d)
{
  xml.openTag("re");
  xml << d.real();
  xml.closeTag();
  xml.openTag("im");
  xml << d.imag();
  xml.closeTag();

  return xml;
}

//! XML input
template<class T>
inline
void read(XMLReader& xml, const std::string& xpath, RComplex<T>& d)
{
  std::ostringstream error_message;
  
  // XPath for the real part 
  std::string path_real = xpath + "/re";
	
  // XPath for the imaginary part.
  std::string path_imag = xpath + "/im";
	
  // Try and recursively get the real part
  try { 
    read(xml, path_real, d.real());
  }
  catch(const std::string &e) {
    error_message << "XPath Query: " << xpath << " Error: "
		  << "Failed to match real part of RComplex Object with self constructed path: " << path_real;
    
    throw error_message.str();
  }
	
  // Try and recursively get the imaginary part
  try {
    read(xml, path_imag, d.imag());
  }
  catch(const std::string &e) {
    error_message << "XPath Query: " << xpath <<" Error:"
		  <<"Failed to match imaginary part of RComplex Object with self constructed path: " << path_imag;
    
    throw error_message.str();
  }
}
#endif

/*! @} */   // end of group rcomplex

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T>
struct WordType<RScalar<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};

template<class T>
struct WordType<RComplex<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};

// Fixed types
template<class T> 
struct SinglePrecType<RScalar<T> >
{
  typedef RScalar<typename SinglePrecType<T>::Type_t>  Type_t;
};

template<class T> 
struct SinglePrecType<RComplex<T> >
{
  typedef RComplex<typename SinglePrecType<T>::Type_t>  Type_t;
};

template<class T> 
struct DoublePrecType<RScalar<T> >
{
  typedef RScalar<typename DoublePrecType<T>::Type_t>  Type_t;
};

template<class T> 
struct DoublePrecType<RComplex<T> >
{
  typedef RComplex<typename DoublePrecType<T>::Type_t>  Type_t;
};


// Internally used scalars
template<class T>
struct InternalScalar<RScalar<T> > {
  typedef RScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

template<class T>
struct InternalScalar<RComplex<T> > {
  typedef RScalar<typename InternalScalar<T>::Type_t>  Type_t;
};


// Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar<RScalar<T> > {
  typedef RScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

template<class T>
struct PrimitiveScalar<RComplex<T> > {
  typedef RScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T>
struct LatticeScalar<RScalar<T> > {
  typedef RScalar<typename LatticeScalar<T>::Type_t>  Type_t;
};

template<class T>
struct LatticeScalar<RComplex<T> > {
  typedef RComplex<typename LatticeScalar<T>::Type_t>  Type_t;
};


// Internally used real scalars
template<class T>
struct RealScalar<RScalar<T> > {
  typedef RScalar<typename RealScalar<T>::Type_t>  Type_t;
};

template<class T>
struct RealScalar<RComplex<T> > {
  typedef RScalar<typename RealScalar<T>::Type_t>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(RScalar) -> RScalar
template<class T1, class Op>
struct UnaryReturn<RScalar<T1>, Op> {
  typedef RScalar<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default unary(RComplex) -> RComplex
template<class T1, class Op>
struct UnaryReturn<RComplex<T1>, Op> {
  typedef RComplex<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default binary(RScalar,RScalar) -> RScalar
template<class T1, class T2, class Op>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, Op> {
  typedef RScalar<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(RComplex,RComplex) -> RComplex
template<class T1, class T2, class Op>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, Op> {
  typedef RComplex<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(RScalar,RComplex) -> RComplex
template<class T1, class T2, class Op>
struct BinaryReturn<RScalar<T1>, RComplex<T2>, Op> {
  typedef RComplex<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Default binary(RComplex,RScalar) -> RComplex
template<class T1, class T2, class Op>
struct BinaryReturn<RComplex<T1>, RScalar<T2>, Op> {
  typedef RComplex<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};




// RScalar
#if 0
template<class T1, class T2>
struct UnaryReturn<RScalar<T2>, OpCast<T1> > {
  typedef RScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif


template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpAddAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpSubtractAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpMultiplyAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpDivideAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpModAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpModAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpBitwiseOrAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpBitwiseAndAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpBitwiseXorAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpLeftShiftAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpRightShiftAssign > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2, class T3>
struct TrinaryReturn<RScalar<T1>, RScalar<T2>, RScalar<T3>, FnColorContract> {
  typedef RScalar<typename TrinaryReturn<T1, T2, T3, FnColorContract>::Type_t>  Type_t;
};

// RScalar
// Gamma algebra
template<int N, int m, class T2, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, RScalar<T2>, OpGammaConstMultiply> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<RScalar<T2>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, RScalar<T2>, OpGammaTypeMultiply> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaType>
struct BinaryReturn<RScalar<T2>, GammaType<N>, OpMultiplyGammaType> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};


// RScalar
// Gamma algebra
template<int N, int m, class T2, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, RScalar<T2>, OpGammaConstDPMultiply> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<RScalar<T2>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, RScalar<T2>, OpGammaTypeDPMultiply> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<RScalar<T2>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef RScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};



// RComplex
// Gamma algebra
template<int N, int m, class T2, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, RComplex<T2>, OpGammaConstMultiply> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<RComplex<T2>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, RComplex<T2>, OpGammaTypeMultiply> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaType>
struct BinaryReturn<RComplex<T2>, GammaType<N>, OpMultiplyGammaType> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};


// Gamma algebra
template<int N, int m, class T2, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, RComplex<T2>, OpGammaConstDPMultiply> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<RComplex<T2>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, RComplex<T2>, OpGammaTypeDPMultiply> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<RComplex<T2>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef RComplex<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};


// Assignment is different
template<class T1, class T2 >
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpAssign > {
//  typedef RComplex<T1> &Type_t;
  typedef RComplex<typename BinaryReturn<T1, T2, OpAssign>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpAddAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpSubtractAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiplyAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpDivideAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpModAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpModAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpBitwiseOrAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpBitwiseAndAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpBitwiseXorAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpLeftShiftAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, OpRightShiftAssign > {
  typedef RComplex<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t>  Type_t;
};
 
template<class T1, class T2, class T3>
struct TrinaryReturn<RComplex<T1>, RComplex<T2>, RComplex<T3>, FnColorContract> {
  typedef RComplex<typename TrinaryReturn<T1, T2, T3, FnColorContract>::Type_t>  Type_t;
};






//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

/*! \addtogroup rscalar
 * @{ 
 */

// Scalar Reality
template<class T>
struct UnaryReturn<RScalar<T>, OpNot > {
  typedef RScalar<typename UnaryReturn<T, OpNot>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RScalar<T1>, OpNot>::Type_t
operator!(const RScalar<T1>& l)
{
  return ! l.elem();
}


template<class T1>
inline typename UnaryReturn<RScalar<T1>, OpUnaryPlus>::Type_t
operator+(const RScalar<T1>& l)
{
  return +l.elem();
}


template<class T1>
inline typename UnaryReturn<RScalar<T1>, OpUnaryMinus>::Type_t
operator-(const RScalar<T1>& l)
{
  return -l.elem();
}


template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpAdd>::Type_t
operator+(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem()+r.elem();
}


template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpSubtract>::Type_t
operator-(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() - r.elem();
}


template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpMultiply>::Type_t
operator*(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() * r.elem();
}

// Optimized  adj(RScalar)*RScalar
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const RScalar<T1>& l, const RScalar<T2>& r)
{
  /*! NOTE: removed transpose here !!!!!  */

//  return transpose(l.elem()) * r.elem();
  return l.elem() * r.elem();
}

// Optimized  RScalar*adj(RScalar)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const RScalar<T1>& l, const RScalar<T2>& r)
{
  /*! NOTE: removed transpose here !!!!!  */

//  return l.elem() * transpose(r.elem());
  return l.elem() * r.elem();
}

// Optimized  adj(RScalar)*adj(RScalar)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const RScalar<T1>& l, const RScalar<T2>& r)
{
  /*! NOTE: removed transpose here !!!!!  */

//  return transpose(l.elem()) * transpose(r.elem());
  return l.elem() * r.elem();
}


template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpDivide>::Type_t
operator/(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() / r.elem();
}



template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpLeftShift > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpLeftShift>::Type_t
operator<<(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() << r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpRightShift > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpRightShift>::Type_t
operator>>(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() >> r.elem();
}


template<class T1, class T2 >
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpMod>::Type_t
operator%(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() % r.elem();
}

template<class T1, class T2 >
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpBitwiseXor>::Type_t
operator^(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() ^ r.elem();
}

template<class T1, class T2 >
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() & r.elem();
}

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpBitwiseOr>::Type_t
operator|(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() | r.elem();
}



// Comparisons
template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpLT > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpLT>::Type_t
operator<(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() < r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpLE > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpLE>::Type_t
operator<=(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() <= r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpGT > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpGT>::Type_t
operator>(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() > r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpGE > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpGE>::Type_t
operator>=(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() >= r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpEQ > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpEQ>::Type_t
operator==(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() == r.elem();
}


template<class T1, class T2 >
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpNE > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpNE>::Type_t
operator!=(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() != r.elem();
}


template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpAnd > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpAnd>::Type_t
operator&&(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() && r.elem();
}


template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, OpOr > {
  typedef RScalar<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, OpOr>::Type_t
operator||(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() || r.elem();
}



//-----------------------------------------------------------------------------
// Functions

// Adjoint
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnAdjoint>::Type_t
adj(const RScalar<T1>& s1)
{
  /*! NOTE: removed transpose here !!!!!  */

//  return transpose(s1.elem()); // The complex nature has been eaten here
  return s1.elem(); // The complex nature has been eaten here
}


// Conjugate
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnConjugate>::Type_t
conj(const RScalar<T1>& s1)
{
  return s1.elem();  // The complex nature has been eaten here
}


// Transpose
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnTranspose>::Type_t
transpose(const RScalar<T1>& s1)
{
  /*! NOTE: removed transpose here !!!!!  */

//  return transpose(s1.elem());
  return s1.elem();
}



// TRACE
// trace = Trace(source1)
template<class T>
struct UnaryReturn<RScalar<T>, FnTrace > {
  typedef RScalar<typename UnaryReturn<T, FnTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnTrace>::Type_t
trace(const RScalar<T1>& s1)
{
//  return trace(s1.elem());

  /*! NOTE: removed trace here !!!!!  */
  return s1.elem();
}


// trace = Re(Trace(source1))
template<class T>
struct UnaryReturn<RScalar<T>, FnRealTrace > {
  typedef RScalar<typename UnaryReturn<T, FnRealTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnRealTrace>::Type_t
realTrace(const RScalar<T1>& s1)
{
//  return trace_real(s1.elem());

  /*! NOTE: removed trace here !!!!!  */
  return s1.elem();
}


// trace = Im(Trace(source1))
template<class T>
struct UnaryReturn<RScalar<T>, FnImagTrace > {
  typedef RScalar<typename UnaryReturn<T, FnImagTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnImagTrace>::Type_t
imagTrace(const RScalar<T1>& s1)
{
//  return trace_imag(s1.elem());

  /*! NOTE: removed trace here !!!!!  */
  return s1.elem();
}

//! RScalar = trace(RScalar * RScalar)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnTraceMultiply>::Type_t
traceMultiply(const RScalar<T1>& l, const RScalar<T2>& r)
{
//  return traceMultiply(l.elem(), r.elem());

  /*! NOTE: removed trace here !!!!!  */
  return l.elem() * r.elem();
}


// RScalar = Re(RScalar)  [identity]
template<class T>
inline typename UnaryReturn<RScalar<T>, FnReal>::Type_t
real(const RScalar<T>& s1)
{
  return s1.elem();
}


// RScalar = Im(RScalar) [this is zero]
template<class T>
inline typename UnaryReturn<RScalar<T>, FnImag>::Type_t
imag(const RScalar<T>& s1)
{
  typedef typename InternalScalar<T>::Type_t  S;
  return S(0);
}


// ArcCos
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnArcCos>::Type_t
acos(const RScalar<T1>& s1)
{
  return acos(s1.elem());
}

// ArcSin
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnArcSin>::Type_t
asin(const RScalar<T1>& s1)
{
  return asin(s1.elem());
}

// ArcTan
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnArcTan>::Type_t
atan(const RScalar<T1>& s1)
{
  return atan(s1.elem());
}

// Ceil(ing)
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnCeil>::Type_t
ceil(const RScalar<T1>& s1)
{
  return ceil(s1.elem());
}

// Cos
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnCos>::Type_t
cos(const RScalar<T1>& s1)
{
  return cos(s1.elem());
}

// Cosh
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnHypCos>::Type_t
cosh(const RScalar<T1>& s1)
{
  return cosh(s1.elem());
}

// Exp
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnExp>::Type_t
exp(const RScalar<T1>& s1)
{
  return exp(s1.elem());
}

// Fabs
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnFabs>::Type_t
fabs(const RScalar<T1>& s1)
{
  return fabs(s1.elem());
}

// Floor
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnFloor>::Type_t
floor(const RScalar<T1>& s1)
{
  return floor(s1.elem());
}

// Log
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnLog>::Type_t
log(const RScalar<T1>& s1)
{
  return log(s1.elem());
}

// Log10
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnLog10>::Type_t
log10(const RScalar<T1>& s1)
{
  return log10(s1.elem());
}

// Sin
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnSin>::Type_t
sin(const RScalar<T1>& s1)
{
  return sin(s1.elem());
}

// Sinh
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnHypSin>::Type_t
sinh(const RScalar<T1>& s1)
{
  return sinh(s1.elem());
}

// Sqrt
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnSqrt>::Type_t
sqrt(const RScalar<T1>& s1)
{
  return sqrt(s1.elem());
}

// Tan
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnTan>::Type_t
tan(const RScalar<T1>& s1)
{
  return tan(s1.elem());
}

// Tanh
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnHypTan>::Type_t
tanh(const RScalar<T1>& s1)
{
  return tanh(s1.elem());
}


//-----------------------------------------------------------------------------
// These functions always return bool
//! isnan
template<class T1>
inline bool
isnan(const RScalar<T1>& s1)
{
  return isnan(s1.elem());
}

//! isinf
template<class T1>
inline bool
isinf(const RScalar<T1>& s1)
{
  return isinf(s1.elem());
}

//! isnormal
template<class T1>
inline bool
isnormal(const RScalar<T1>& s1)
{
  return isnormal(s1.elem());
}

//! isfinite
template<class T1>
inline bool
isfinite(const RScalar<T1>& s1)
{
  return isfinite(s1.elem());
}


//-----------------------------------------------------------------------------
//! RScalar<T> = pow(RScalar<T> , RScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnPow>::Type_t
pow(const RScalar<T1>& s1, const RScalar<T2>& s2)
{
  return pow(s1.elem(), s2.elem());
}

//! RScalar<T> = atan2(RScalar<T> , RScalar<T>)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnArcTan2>::Type_t
atan2(const RScalar<T1>& s1, const RScalar<T2>& s2)
{
  return atan2(s1.elem(), s2.elem());
}


//! RScalar = outerProduct(RScalar, RScalar)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const RScalar<T1>& l, const RScalar<T2>& r)
{
  return l.elem() * r.elem();
}


//! dest [float type] = source [seed type]
template<class T1>
inline typename UnaryReturn<RScalar<T1>, FnSeedToFloat>::Type_t
seedToFloat(const RScalar<T1>& s1)
{
  return seedToFloat(s1.elem());
}

//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T>
inline typename UnaryReturn<RScalar<T>, FnGetSite>::Type_t
getSite(const RScalar<T>& s1, int innersite)
{
  return getSite(s1.elem(), innersite);
}

//! Extract color vector components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T>
inline typename UnaryReturn<RScalar<T>, FnPeekColorVector>::Type_t
peekColor(const RScalar<T>& l, int row)
{
  return peekColor(l.elem(),row);
}

//! Extract color matrix components 
/*! Generically, this is an identity operation. Defined differently under color */
template<class T>
inline typename UnaryReturn<RScalar<T>, FnPeekColorMatrix>::Type_t
peekColor(const RScalar<T>& l, int row, int col)
{
  return peekColor(l.elem(),row,col);
}

//! Extract spin vector components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<RScalar<T>, FnPeekSpinVector>::Type_t
peekSpin(const RScalar<T>& l, int row)
{
  return peekSpin(l.elem(),row);
}

//! Extract spin matrix components 
/*! Generically, this is an identity operation. Defined differently under spin */
template<class T>
inline typename UnaryReturn<RScalar<T>, FnPeekSpinMatrix>::Type_t
peekSpin(const RScalar<T>& l, int row, int col)
{
  return peekSpin(l.elem(),row,col);
}

//-----------------------------------------------------------------------------
//! QDP Int to int primitive in conversion routine
template<class T> 
inline int 
toInt(const RScalar<T>& s) 
{
  return toInt(s.elem());
}

//! QDP Real to float primitive in conversion routine
template<class T> 
inline float
toFloat(const RScalar<T>& s) 
{
  return toFloat(s.elem());
}

//! QDP Double to double primitive in conversion routine
template<class T> 
inline double
toDouble(const RScalar<T>& s) 
{
  return toDouble(s.elem());
}

//! QDP Boolean to bool primitive in conversion routine
template<class T> 
inline bool
toBool(const RScalar<T>& s) 
{
  return toBool(s.elem());
}

//! QDP Wordtype to primitive wordtype
template<class T> 
inline typename WordType< RScalar<T> >::Type_t
toWordType(const RScalar<T>& s) 
{
  return toWordType(s.elem());
}



//------------------------------------------
//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline
void copymask(RScalar<T>& d, const RScalar<T1>& mask, const RScalar<T>& s1) 
{
  copymask(d.elem(),mask.elem(),s1.elem());
}

//! dest [float type] = source [int type]
template<class T, class T1>
inline
void cast_rep(T& d, const RScalar<T1>& s1)
{
  cast_rep(d, s1.elem());
}


//! dest [float type] = source [int type]
template<class T, class T1>
inline
void recast_rep(RScalar<T>& d, const RScalar<T1>& s1)
{
  cast_rep(d.elem(), s1.elem());
}


//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
copy_site(RScalar<T>& d, int isite, const RScalar<T1>& s1)
{
  copy_site(d.elem(), isite, s1.elem());
}


//! gather several inner sites together
template<class T, class T1>
inline void 
gather_sites(RScalar<T>& d, 
	     const RScalar<T1>& s0, int i0, 
	     const RScalar<T1>& s1, int i1,
	     const RScalar<T1>& s2, int i2,
	     const RScalar<T1>& s3, int i3)
{
  gather_sites(d.elem(), 
	       s0.elem(), i0, 
	       s1.elem(), i1, 
	       s2.elem(), i2, 
	       s3.elem(), i3);
}


#if 1
// Global sum over site indices only
template<class T>
struct UnaryReturn<RScalar<T>, FnSum > {
  typedef RScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RScalar<T>, FnSum>::Type_t
sum(const RScalar<T>& s1)
{
  return sum(s1.elem());
}
#endif


// Global max
template<class T>
struct UnaryReturn<RScalar<T>, FnGlobalMax> {
  typedef RScalar<typename UnaryReturn<T, FnGlobalMax>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RScalar<T>, FnGlobalMax>::Type_t
globalMax(const RScalar<T>& s1)
{
  return globalMax(s1.elem());
}


// Global min
template<class T>
struct UnaryReturn<RScalar<T>, FnGlobalMin> {
  typedef RScalar<typename UnaryReturn<T, FnGlobalMin>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RScalar<T>, FnGlobalMin>::Type_t
globalMin(const RScalar<T>& s1)
{
  return globalMin(s1.elem());
}



//------------------------------------------
// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<RScalar<T>, FnNorm2 > {
  typedef RScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<RScalar<T>, FnLocalNorm2 > {
  typedef RScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RScalar<T>, FnLocalNorm2>::Type_t
localNorm2(const RScalar<T>& s1)
{
  return localNorm2(s1.elem());
}



//! RScalar<T> = InnerProduct(adj(RScalar<T1>)*RScalar<T2>)
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, FnInnerProduct > {
  typedef RScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, FnLocalInnerProduct > {
  typedef RScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const RScalar<T1>& s1, const RScalar<T2>& s2)
{
  return localInnerProduct(s1.elem(), s2.elem());
}


//! RScalar<T> = InnerProductReal(adj(PMatrix<T1>)*PMatrix<T1>)
// Real-ness is eaten at this level
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, FnInnerProductReal > {
  typedef RScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, FnLocalInnerProductReal > {
  typedef RScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const RScalar<T1>& s1, const RScalar<T2>& s2)
{
  return localInnerProduct(s1.elem(), s2.elem());
}


//! RScalar<T> = where(RScalar, RScalar, RScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<RScalar<T1>, RScalar<T2>, RScalar<T3>, FnWhere> {
  typedef RScalar<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<RScalar<T1>, RScalar<T2>, RScalar<T3>, FnWhere>::Type_t
where(const RScalar<T1>& a, const RScalar<T2>& b, const RScalar<T3>& c)
{
  return where(a.elem(), b.elem(), c.elem());
}



//-----------------------------------------------------------------------------
// Broadcast operations
//! dest = 0
template<class T> 
inline
void zero_rep(RScalar<T>& dest) 
{
  zero_rep(dest.elem());
}


//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
copy_site(RComplex<T>& d, int isite, const RComplex<T1>& s1)
{
  copy_site(d.real(), isite, s1.real());
  copy_site(d.imag(), isite, s1.imag());
}

#if 0
//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
copy_site(RComplex<T>& d, int isite, const RScalar<T1>& s1)
{
  copy_site(d.real(), isite, s1.elem());
  zero_rep(d.imag());   // this is wrong - want zero only at a site. Fix when needed.
}
#endif


//! gather several inner sites together
template<class T, class T1>
inline void 
gather_sites(RComplex<T>& d, 
	     const RComplex<T1>& s0, int i0, 
	     const RComplex<T1>& s1, int i1,
	     const RComplex<T1>& s2, int i2,
	     const RComplex<T1>& s3, int i3)
{
  gather_sites(d.real(), 
	       s0.real(), i0, 
	       s1.real(), i1, 
	       s2.real(), i2, 
	       s3.real(), i3);

  gather_sites(d.imag(), 
	       s0.imag(), i0, 
	       s1.imag(), i1, 
	       s2.imag(), i2, 
	       s3.imag(), i3);
}


//! dest  = random  
template<class T, class T1, class T2>
inline void
fill_random(RScalar<T>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  fill_random(d.elem(), seed, skewed_seed, seed_mult);
}



//! dest  = gaussian  
/*! Real form of complex polar method */
template<class T>
inline void
fill_gaussian(RScalar<T>& d, RScalar<T>& r1, RScalar<T>& r2)
{
  typedef typename InternalScalar<T>::Type_t  S;

  // r1 and r2 are the input random numbers needed

  /* Stage 2: get the cos of the second number  */
  T  g_r;

  r2.elem() *= S(6.283185307);
  g_r = cos(r2.elem());
    
  /* Stage 4: get  sqrt(-2.0 * log(u1)) */
  r1.elem() = sqrt(-S(2.0) * log(r1.elem()));

  /* Stage 5:   g_r = sqrt(-2*log(u1))*cos(2*pi*u2) */
  /* Stage 5:   g_i = sqrt(-2*log(u1))*sin(2*pi*u2) */
  d.elem() = r1.elem() * g_r;
}

/*! @} */   // end of group rscalar



//-----------------------------------------------------------------------------
// Complex Reality
//-----------------------------------------------------------------------------

/*! \addtogroup rcomplex 
 * @{ 
 */

//! RComplex = +RComplex
template<class T1>
inline typename UnaryReturn<RComplex<T1>, OpUnaryPlus>::Type_t
operator+(const RComplex<T1>& l)
{
  typedef typename UnaryReturn<RComplex<T1>, OpUnaryPlus>::Type_t  Ret_t;

  return Ret_t(+l.real(),
	       +l.imag());
}


//! RComplex = -RComplex
template<class T1>
inline typename UnaryReturn<RComplex<T1>, OpUnaryMinus>::Type_t
operator-(const RComplex<T1>& l)
{
  typedef typename UnaryReturn<RComplex<T1>, OpUnaryMinus>::Type_t  Ret_t;

  return Ret_t(-l.real(),
	       -l.imag());
}


//! RComplex = RComplex + RComplex
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpAdd>::Type_t
operator+(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpAdd>::Type_t  Ret_t;

  return Ret_t(l.real()+r.real(),
	       l.imag()+r.imag());
}

//! RComplex = RComplex + RScalar
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpAdd>::Type_t
operator+(const RComplex<T1>& l, const RScalar<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpAdd>::Type_t  Ret_t;

  return Ret_t(l.real()+r.elem(),
	       l.imag());
}

//! RComplex = RScalar + RComplex
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpAdd>::Type_t
operator+(const RScalar<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpAdd>::Type_t  Ret_t;

  return Ret_t(l.elem()+r.real(),
	       r.imag());
}


//! RComplex = RComplex - RComplex
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpSubtract>::Type_t
operator-(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpSubtract>::Type_t  Ret_t;

  return Ret_t(l.real() - r.real(),
	       l.imag() - r.imag());
}

//! RComplex = RComplex - RScalar
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpSubtract>::Type_t
operator-(const RComplex<T1>& l, const RScalar<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpSubtract>::Type_t  Ret_t;

  return Ret_t(l.real() - r.elem(),
	       l.imag());
}

//! RComplex = RScalar - RComplex
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpSubtract>::Type_t
operator-(const RScalar<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpSubtract>::Type_t  Ret_t;

  return Ret_t(l.elem() - r.real(),
	       - r.imag());
}


//! RComplex = RComplex * RComplex
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiply>::Type_t
operator*(const RComplex<T1>& __restrict__ l, const RComplex<T2>& __restrict__ r) 
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiply>::Type_t  Ret_t;

  return Ret_t(l.real()*r.real() - l.imag()*r.imag(),
	       l.real()*r.imag() + l.imag()*r.real());
}

//! RComplex = RScalar * RComplex
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpMultiply>::Type_t
operator*(const RScalar<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpMultiply>::Type_t  Ret_t;

  return Ret_t(l.elem()*r.real(), 
	       l.elem()*r.imag());
}

//! RComplex = RComplex * RScalar
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpMultiply>::Type_t
operator*(const RComplex<T1>& l, const RScalar<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpMultiply>::Type_t  Ret_t;

  return Ret_t(l.real()*r.elem(), 
	       l.imag()*r.elem());
}


// Optimized  adj(RComplex)*RComplex
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpAdjMultiply>::Type_t
adjMultiply(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpAdjMultiply>::Type_t  Ret_t;

  // The complex conjugate nature has been eaten here leaving simple multiples
  // involving transposes - which are probably null
  
//  d.real() = transpose(l.real())*r.real() + transpose(l.imag())*r.imag();
//  d.imag() = transpose(l.real())*r.imag() - transpose(l.imag())*r.real();
//  return d;

  /*! NOTE: removed transpose here !!!!!  */
  return Ret_t(l.real()*r.real() + l.imag()*r.imag(),
	       l.real()*r.imag() - l.imag()*r.real());
}

// Optimized  RComplex*adj(RComplex)
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiplyAdj>::Type_t  Ret_t;

  // The complex conjugate nature has been eaten here leaving simple multiples
  // involving transposes - which are probably null
//  d.real() = l.real()*transpose(r.real()) + l.imag()*transpose(r.imag());
//  d.imag() = l.imag()*transpose(r.real()) - l.real()*transpose(r.imag());
//  return d;

  /*! NOTE: removed transpose here !!!!!  */
  return Ret_t(l.real()*r.real() + l.imag()*r.imag(),
	       l.imag()*r.real() - l.real()*r.imag());
}

// Optimized  adj(RComplex)*adj(RComplex)
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpAdjMultiplyAdj>::Type_t  Ret_t;

  // The complex conjugate nature has been eaten here leaving simple multiples
  // involving transposes - which are probably null
//  d.real() = transpose(l.real())*transpose(r.real()) - transpose(l.imag())*transpose(r.imag());
//  d.imag() = -(transpose(l.real())*transpose(r.imag()) + transpose(l.imag())*transpose(r.real()));
//  return d;

  /*! NOTE: removed transpose here !!!!!  */
  return Ret_t(l.real()*r.real() - l.imag()*r.imag(),
	       -(l.real()*r.imag() + l.imag()*r.real()));
}


//! RComplex = RComplex / RComplex
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpDivide>::Type_t
operator/(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpDivide>::Type_t  Ret_t;

  T2 tmp = T2(1.0) / (r.real()*r.real() + r.imag()*r.imag());

  return Ret_t((l.real()*r.real() + l.imag()*r.imag()) * tmp,
	       (l.imag()*r.real() - l.real()*r.imag()) * tmp);
}

//! RComplex = RComplex / RScalar
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpDivide>::Type_t
operator/(const RComplex<T1>& l, const RScalar<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RScalar<T2>, OpDivide>::Type_t  Ret_t;

  T2 tmp = T2(1.0) / r.elem();

  return Ret_t(l.real() * tmp, 
	       l.imag() * tmp);
}

//! RComplex = RScalar / RComplex
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpDivide>::Type_t
operator/(const RScalar<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RScalar<T1>, RComplex<T2>, OpDivide>::Type_t  Ret_t;

  T2 tmp = T2(1.0) / (r.real()*r.real() + r.imag()*r.imag());

  return Ret_t(l.elem() * r.real() * tmp,
	       -l.elem() * r.imag() * tmp);
}



//-----------------------------------------------------------------------------
// These functions always return bool
//! isnan
template<class T1>
inline bool
isnan(const RComplex<T1>& s1)
{
  return isnan(s1.real()) | isnan(s1.imag());
}

//! isinf
template<class T1>
inline bool
isinf(const RComplex<T1>& s1)
{
  return isinf(s1.real()) | isinf(s1.imag());
}

//! isnormal
template<class T1>
inline bool
isnormal(const RComplex<T1>& s1)
{
  return isnormal(s1.real()) & isnormal(s1.imag());
}

//! isfinite
template<class T1>
inline bool
isfinite(const RComplex<T1>& s1)
{
  return isfinite(s1.real()) & isfinite(s1.imag());
}


//-----------------------------------------------------------------------------
// Functions

// Adjoint
template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnAdjoint>::Type_t
adj(const RComplex<T1>& l)
{
  typedef typename UnaryReturn<RComplex<T1>, FnAdjoint>::Type_t  Ret_t;

  // The complex conjugate nature has been eaten here leaving transpose
//  d.real() = transpose(l.real());
//  d.imag() = -transpose(l.imag());
//  return d;

  /*! NOTE: removed transpose here !!!!!  */
  return Ret_t(l.real(),
	       -l.imag());
}

// Conjugate
template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnConjugate>::Type_t
conj(const RComplex<T1>& l)
{
  typedef typename UnaryReturn<RComplex<T1>, FnConjugate>::Type_t  Ret_t;

  return Ret_t(l.real(),
	       -l.imag());
}

// Transpose
template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnTranspose>::Type_t
transpose(const RComplex<T1>& l)
{
  typedef typename UnaryReturn<RComplex<T1>, FnTranspose>::Type_t  Ret_t;

//  d.real() = transpose(l.real());
//  d.imag() = transpose(l.imag());
//  return d;

  /*! NOTE: removed transpose here !!!!!  */
  return Ret_t(l.real(), 
	       l.imag());
}

// TRACE
// trace = Trace(source1)
template<class T>
struct UnaryReturn<RComplex<T>, FnTrace > {
  typedef RComplex<typename UnaryReturn<T, FnTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnTrace>::Type_t
trace(const RComplex<T1>& s1)
{
  typedef typename UnaryReturn<RComplex<T1>, FnTrace>::Type_t  Ret_t;

  /*! NOTE: removed trace here !!!!!  */
  return Ret_t(s1.real(),
	       s1.imag());
}


// trace = Re(Trace(source1))
template<class T>
struct UnaryReturn<RComplex<T>, FnRealTrace > {
  typedef RScalar<typename UnaryReturn<T, FnRealTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnRealTrace>::Type_t
realTrace(const RComplex<T1>& s1)
{
  /*! NOTE: removed trace here !!!!!  */
  return s1.real();
}


// trace = Im(Trace(source1))
template<class T>
struct UnaryReturn<RComplex<T>, FnImagTrace > {
  typedef RScalar<typename UnaryReturn<T, FnImagTrace>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnImagTrace>::Type_t
imagTrace(const RComplex<T1>& s1)
{
  /*! NOTE: removed trace here !!!!!  */
  return s1.imag();
}

//! RComplex = trace(RComplex * RComplex)
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiply>::Type_t
traceMultiply(const RComplex<T1>& l, const RComplex<T2>& r)
{
//  return traceMultiply(l.elem(), r.elem());

  /*! NOTE: removed trace here !!!!!  */
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, OpMultiply>::Type_t  Ret_t;

  return Ret_t(l.real()*r.real() - l.imag()*r.imag(),
	       l.real()*r.imag() + l.imag()*r.real());
}


// RScalar = Re(RComplex)
template<class T>
struct UnaryReturn<RComplex<T>, FnReal > {
  typedef RScalar<typename UnaryReturn<T, FnReal>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnReal>::Type_t
real(const RComplex<T1>& s1)
{
  return s1.real();
}

// RScalar = Im(RComplex)
template<class T>
struct UnaryReturn<RComplex<T>, FnImag > {
  typedef RScalar<typename UnaryReturn<T, FnImag>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<RComplex<T1>, FnImag>::Type_t
imag(const RComplex<T1>& s1)
{
  return s1.imag();
}


//! RComplex<T> = (RScalar<T> , RScalar<T>)
template<class T1, class T2>
struct BinaryReturn<RScalar<T1>, RScalar<T2>, FnCmplx > {
  typedef RComplex<typename BinaryReturn<T1, T2, FnCmplx>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnCmplx>::Type_t
cmplx(const RScalar<T1>& s1, const RScalar<T2>& s2)
{
  typedef typename BinaryReturn<RScalar<T1>, RScalar<T2>, FnCmplx>::Type_t  Ret_t;

  return Ret_t(s1.elem(),
	       s2.elem());
}



// RComplex = i * RScalar
template<class T>
struct UnaryReturn<RScalar<T>, FnTimesI > {
  typedef RComplex<typename UnaryReturn<T, FnTimesI>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RScalar<T>, FnTimesI>::Type_t
timesI(const RScalar<T>& s1)
{
  typename UnaryReturn<RScalar<T>, FnTimesI>::Type_t  d;

  zero_rep(d.real());
  d.imag() = s1.elem();
  return d;
}

// RComplex = i * RComplex
template<class T>
inline typename UnaryReturn<RComplex<T>, FnTimesI>::Type_t
timesI(const RComplex<T>& s1)
{
  typedef typename UnaryReturn<RComplex<T>, FnTimesI>::Type_t  Ret_t;

  return Ret_t(-s1.imag(),
	       s1.real());
}


// RComplex = -i * RScalar
template<class T>
struct UnaryReturn<RScalar<T>, FnTimesMinusI > {
  typedef RComplex<typename UnaryReturn<T, FnTimesMinusI>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RScalar<T>, FnTimesMinusI>::Type_t
timesMinusI(const RScalar<T>& s1)
{
  typename UnaryReturn<RScalar<T>, FnTimesMinusI>::Type_t  d;

  zero_rep(d.real());
  d.imag() = -s1.elem();
  return d;
}


// RComplex = -i * RComplex
template<class T>
inline typename UnaryReturn<RComplex<T>, FnTimesMinusI>::Type_t
timesMinusI(const RComplex<T>& s1)
{
  typedef typename UnaryReturn<RComplex<T>, FnTimesMinusI>::Type_t  Ret_t;

  return Ret_t(s1.imag(),
	       -s1.real());
}


//! RComplex = outerProduct(RComplex, RComplex)
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, FnOuterProduct>::Type_t
outerProduct(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, FnOuterProduct>::Type_t  Ret_t;

  // Return   l*conj(r)
  return Ret_t(l.real()*r.real() + l.imag()*r.imag(),
	       l.imag()*r.real() - l.real()*r.imag());
}

//! RComplex = outerProduct(RComplex, RScalar)
template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const RComplex<T1>& l, const RScalar<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RScalar<T2>, FnOuterProduct>::Type_t  Ret_t;

  // Return   l*conj(r)
  return Ret_t(l.real()*r.elem(),
	       l.imag()*r.elem());
}

//! RComplex = outerProduct(RScalar, RComplex)
template<class T1, class T2>
inline typename BinaryReturn<RScalar<T1>, RComplex<T2>, FnOuterProduct>::Type_t
outerProduct(const RScalar<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RScalar<T1>, RComplex<T2>, FnOuterProduct>::Type_t  Ret_t;

  // Return   l*conj(r)
  return Ret_t( l.elem()*r.real(),
	       -l.elem()*r.imag());
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T>
inline typename UnaryReturn<RComplex<T>, FnGetSite>::Type_t
getSite(const RComplex<T>& s1, int innersite)
{
  typedef typename UnaryReturn<RComplex<T>, FnGetSite>::Type_t  Ret_t;

  return Ret_t(getSite(s1.real(), innersite), 
	       getSite(s1.imag(), innersite));
}


//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline
void copymask(RComplex<T>& d, const RScalar<T1>& mask, const RComplex<T>& s1) 
{
  copymask(d.real(),mask.elem(),s1.real());
  copymask(d.imag(),mask.elem(),s1.imag());
}


#if 1
// Global sum over site indices only
template<class T>
struct UnaryReturn<RComplex<T>, FnSum> {
  typedef RComplex<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RComplex<T>, FnSum>::Type_t
sum(const RComplex<T>& s1)
{
  typedef typename UnaryReturn<RComplex<T>, FnSum>::Type_t  Ret_t;

  return Ret_t(sum(s1.real()),
	       sum(s1.imag()));
}
#endif


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<RComplex<T>, FnNorm2 > {
  typedef RScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<RComplex<T>, FnLocalNorm2 > {
  typedef RScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<RComplex<T>, FnLocalNorm2>::Type_t
localNorm2(const RComplex<T>& s1)
{
  return localNorm2(s1.real()) + localNorm2(s1.imag());
}



//! RComplex<T> = InnerProduct(adj(RComplex<T1>)*RComplex<T2>)
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, FnInnerProduct > {
  typedef RComplex<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, FnLocalInnerProduct > {
  typedef RComplex<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const RComplex<T1>& l, const RComplex<T2>& r)
{
  typedef typename BinaryReturn<RComplex<T1>, RComplex<T2>, FnLocalInnerProduct>::Type_t  Ret_t;

  return Ret_t(localInnerProduct(l.real(),r.real()) + localInnerProduct(l.imag(),r.imag()),
	       localInnerProduct(l.real(),r.imag()) - localInnerProduct(l.imag(),r.real()));
}


//! RScalar<T> = InnerProductReal(adj(RComplex<T1>)*RComplex<T1>)
// Real-ness is eaten at this level
template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, FnInnerProductReal > {
  typedef RScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<RComplex<T1>, RComplex<T2>, FnLocalInnerProductReal > {
  typedef RScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<RComplex<T1>, RComplex<T2>, FnLocalInnerProductReal>::Type_t
localInnerProductReal(const RComplex<T1>& l, const RComplex<T2>& r)
{
  return localInnerProduct(l.real(),r.real()) + localInnerProduct(l.imag(),r.imag());
}


//! RComplex<T> = where(RScalar, RComplex, RComplex)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<RScalar<T1>, RComplex<T2>, RComplex<T3>, FnWhere> {
  typedef RComplex<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<RScalar<T1>, RComplex<T2>, RComplex<T3>, FnWhere>::Type_t
where(const RScalar<T1>& a, const RComplex<T2>& b, const RComplex<T3>& c)
{
  typedef typename TrinaryReturn<RScalar<T1>, RComplex<T2>, RComplex<T3>, FnWhere>::Type_t  Ret_t;

  // Not optimal - want to have where outside assignment
  return Ret_t(where(a.elem(), b.real(), c.real()),
	       where(a.elem(), b.imag(), c.imag()));
}

//! RComplex<T> = where(RScalar, RComplex, RScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<RScalar<T1>, RComplex<T2>, RScalar<T3>, FnWhere> {
  typedef RComplex<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<RScalar<T1>, RComplex<T2>, RComplex<T3>, FnWhere>::Type_t
where(const RScalar<T1>& a, const RComplex<T2>& b, const RScalar<T3>& c)
{
  typedef typename TrinaryReturn<RScalar<T1>, RComplex<T2>, RScalar<T3>, FnWhere>::Type_t  Ret_t;
  typedef typename InternalScalar<T3>::Type_t  S;

  // Not optimal - want to have where outside assignment
  return Ret_t(where(a.elem(), b.real(), c.real()),
	       where(a.elem(), b.imag(), S(0)));
}

//! RComplex<T> = where(RScalar, RScalar, RComplex)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<RScalar<T1>, RScalar<T2>, RComplex<T3>, FnWhere> {
  typedef RComplex<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<RScalar<T1>, RScalar<T2>, RComplex<T3>, FnWhere>::Type_t
where(const RScalar<T1>& a, const RScalar<T2>& b, const RComplex<T3>& c)
{
  typedef typename TrinaryReturn<RScalar<T1>, RScalar<T2>, RComplex<T3>, FnWhere>::Type_t  Ret_t;
  typedef typename InternalScalar<T2>::Type_t  S;

  // Not optimal - want to have where outside assignment
  return Ret_t(where(a.elem(), b.real(), c.real()),
	       where(a.elem(), S(0), c.imag()));
}


//-----------------------------------------------------------------------------
// Broadcast operations
//! dest = 0
template<class T> 
inline
void zero_rep(RComplex<T>& dest) 
{
  zero_rep(dest.real());
  zero_rep(dest.imag());
}


//! dest  = random  
template<class T, class T1, class T2>
inline void
fill_random(RComplex<T>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  fill_random(d.real(), seed, skewed_seed, seed_mult);
  fill_random(d.imag(), seed, skewed_seed, seed_mult);
}


//! dest  = gaussian
/*! RComplex polar method */
template<class T>
inline void
fill_gaussian(RComplex<T>& d, RComplex<T>& r1, RComplex<T>& r2)
{
  typedef typename InternalScalar<T>::Type_t  S;

  // r1 and r2 are the input random numbers needed

  /* Stage 2: get the cos of the second number  */
  T  g_r, g_i;

  r2.real() *= S(6.283185307);
  g_r = cos(r2.real());
  g_i = sin(r2.real());
    
  /* Stage 4: get  sqrt(-2.0 * log(u1)) */
  r1.real() = sqrt(-S(2.0) * log(r1.real()));

  /* Stage 5:   g_r = sqrt(-2*log(u1))*cos(2*pi*u2) */
  /* Stage 5:   g_i = sqrt(-2*log(u1))*sin(2*pi*u2) */
  d.real() = r1.real() * g_r;
  d.imag() = r1.real() * g_i;
}

/*! @} */  // end of group rcomplex

} // namespace QDP

#endif
