// -*- C++ -*-

/*! \file
 * \brief Inner grid
 */

#ifndef QDP_INNER_H
#define QDP_INNER_H

namespace QDP {

//-------------------------------------------------------------------------------------
/*! \addtogroup iscalar Inner grid scalar
 * \ingroup fiberbundle
 *
 * Inner grid scalar means sites are the slowest varying index. There can
 * still be an Outer grid.
 *
 * @{
 */

//! Scalar inner lattice
/*! All inner lattices are of IScalar or ILattice type */
template<class T> class IScalar
{
public:
  IScalar() {}
  ~IScalar() {}

  //---------------------------------------------------------
  //! construct dest = const
  IScalar(const typename WordType<T>::Type_t& rhs) : F(rhs) {}

  //! construct dest = rhs
  template<class T1>
  IScalar(const IScalar<T1>& rhs) : F(rhs.elem()) {}

  //! construct dest = rhs
  template<class T1>
  IScalar(const T1& rhs) : F(rhs) {}

  //---------------------------------------------------------
#if 0
  //! dest = const
  /*! Fill with an integer constant. Will be promoted to underlying word type */
  inline
  IScalar& operator=(const typename WordType<T>::Type_t& rhs)
    {
      elem() = rhs;
      return *this;
    }
#endif

  //! IScalar = IScalar
  /*! Set equal to another IScalar */
  template<class T1>
  inline
  IScalar& operator=(const IScalar<T1>& rhs) 
    {
      elem() = rhs.elem();
      return *this;
    }

  //! IScalar += IScalar
  template<class T1>
  inline
  IScalar& operator+=(const IScalar<T1>& rhs) 
    {
      elem() += rhs.elem();
      return *this;
    }

  //! IScalar -= IScalar
  template<class T1>
  inline
  IScalar& operator-=(const IScalar<T1>& rhs) 
    {
      elem() -= rhs.elem();
      return *this;
    }

  //! IScalar *= IScalar
  template<class T1>
  inline
  IScalar& operator*=(const IScalar<T1>& rhs) 
    {
      elem() *= rhs.elem();
      return *this;
    }

  //! IScalar /= IScalar
  template<class T1>
  inline
  IScalar& operator/=(const IScalar<T1>& rhs) 
    {
      elem() /= rhs.elem();
      return *this;
    }

  //! IScalar %= IScalar
  template<class T1>
  inline
  IScalar& operator%=(const IScalar<T1>& rhs) 
    {
      elem() %= rhs.elem();
      return *this;
    }

  //! IScalar |= IScalar
  template<class T1>
  inline
  IScalar& operator|=(const IScalar<T1>& rhs) 
    {
      elem() |= rhs.elem();
      return *this;
    }

  //! IScalar &= IScalar
  template<class T1>
  inline
  IScalar& operator&=(const IScalar<T1>& rhs) 
    {
      elem() &= rhs.elem();
      return *this;
    }

  //! IScalar ^= IScalar
  template<class T1>
  inline
  IScalar& operator^=(const IScalar<T1>& rhs) 
    {
      elem() ^= rhs.elem();
      return *this;
    }

  //! IScalar <<= IScalar
  template<class T1>
  inline
  IScalar& operator<<=(const IScalar<T1>& rhs) 
    {
      elem() <<= rhs.elem();
      return *this;
    }

  //! IScalar >>= IScalar
  template<class T1>
  inline
  IScalar& operator>>=(const IScalar<T1>& rhs) 
    {
      elem() >>= rhs.elem();
      return *this;
    }

  //! Deep copies here
  IScalar(const IScalar& a): F(a.F) {/* fprintf(stderr,"copy IScalar\n"); */}

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
std::istream& operator>>(std::istream& s, IScalar<T>& d)
{
  s >> d.elem();
  return s;
}

//! Ascii input
template<class T>
StandardInputStream& operator>>(StandardInputStream& is, IScalar<T>& d)
{
  return is >> d.elem();
}

//! Ascii output
template<class T> 
inline  
std::ostream& operator<<(std::ostream& s, const IScalar<T>& d)
{
  return s << d.elem();
}

//! Ascii output
template<class T>
StandardOutputStream& operator<<(StandardOutputStream& os, const IScalar<T>& d)
{
  return os << d.elem();
}


//! Text input
template<class T>
inline
TextReader& operator>>(TextReader& s, IScalar<T>& d)
{
  return s >> d.elem();
}

//! Text output
template<class T> 
inline  
TextWriter& operator<<(TextWriter& s, const IScalar<T>& d)
{
  return s << d.elem();
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T>
inline
XMLWriter& operator<<(XMLWriter& xml, const IScalar<T>& d)
{
  return xml << d.elem();
}

//! XML input
template<class T>
inline
void read(XMLReader& xml, const std::string& path, IScalar<T>& d)
{
  read(xml, path, d.elem());
}
#endif

/*! @} */  // end of group iscalar




//-------------------------------------------------------------------------------------
//! Inner lattice class
/*!
 * Lattice class for vector like architectures. 
 * Also mixed mode super-scalar architectures with vector extensions

 * \addtogroup ilattice Lattice inner grid
 * \ingroup fiberbundle 
 * @{
 */
template<class T, int N> class ILattice
{
public:
  ILattice() {}
  ~ILattice() {}

  //---------------------------------------------------------
  //! construct dest = const
  ILattice(const typename WordType<T>::Type_t& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs;
    }

  //! construct dest = rhs
  template<class T1>
  ILattice(const IScalar<T1>& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem();
    }

  //! construct dest = rhs
  template<class T1>
  ILattice(const ILattice<T1,N>& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);
    }

  //! construct dest = rhs
  template<class T1>
  ILattice(const T1& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs;
    }

  //---------------------------------------------------------
#if 0
  //! dest = const
  /*! Fill with an integer constant. Will be promoted to underlying word type */
  inline
  ILattice& operator=(const typename WordType<T>::Type_t& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs;

      return *this;
    }
#endif

  //---------------------------------------------------------
  //! ILattice = IScalar
  /*! Set equal to an IScalar */
  template<class T1>
  inline
  ILattice& operator=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem();

      return *this;
    }

  //! ILattice += IScalar
  template<class T1>
  inline
  ILattice& operator+=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem();

      return *this;
    }

  //! ILattice -= IScalar
  template<class T1>
  inline
  ILattice& operator-=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem();

      return *this;
    }

  //! ILattice *= IScalar
  template<class T1>
  inline
  ILattice& operator*=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) *= rhs.elem();

      return *this;
    }

  //! ILattice /= IScalar
  template<class T1>
  inline
  ILattice& operator/=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) /= rhs.elem();

      return *this;
    }

  //! ILattice %= IScalar
  template<class T1>
  inline
  ILattice& operator%=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) %= rhs.elem();

      return *this;
    }

  //! ILattice |= IScalar
  template<class T1>
  inline
  ILattice& operator|=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) |= rhs.elem();

      return *this;
    }

  //! ILattice &= IScalar
  template<class T1>
  inline
  ILattice& operator&=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) &= rhs.elem();

      return *this;
    }

  //! ILattice ^= IScalar
  template<class T1>
  inline
  ILattice& operator^=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) ^= rhs.elem();

      return *this;
    }

  //! ILattice <<= IScalar
  template<class T1>
  inline
  ILattice& operator<<=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) <<= rhs.elem();

      return *this;
    }

  //! ILattice >>= IScalar
  template<class T1>
  inline
  ILattice& operator>>=(const IScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) >>= rhs.elem();

      return *this;
    }

  //---------------------------------------------------------
  //! ILattice = ILattice
  /*! Set equal to another ILattice */
  template<class T1>
  inline
  ILattice& operator=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! ILattice += ILattice
  template<class T1>
  inline
  ILattice& operator+=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem(i);

      return *this;
    }

  //! ILattice -= ILattice
  template<class T1>
  inline
  ILattice& operator-=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem(i);

      return *this;
    }

  //! ILattice *= ILattice
  template<class T1>
  inline
  ILattice& operator*=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) *= rhs.elem(i);

      return *this;
    }

  //! ILattice /= ILattice
  template<class T1>
  inline
  ILattice& operator/=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) /= rhs.elem(i);

      return *this;
    }

  //! ILattice %= ILattice
  template<class T1>
  inline
  ILattice& operator%=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) %= rhs.elem(i);

      return *this;
    }

  //! ILattice |= ILattice
  template<class T1>
  inline
  ILattice& operator|=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) |= rhs.elem(i);

      return *this;
    }

  //! ILattice &= ILattice
  template<class T1>
  inline
  ILattice& operator&=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) &= rhs.elem(i);

      return *this;
    }

  //! ILattice ^= ILattice
  template<class T1>
  inline
  ILattice& operator^=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) ^= rhs.elem(i);

      return *this;
    }

  //! ILattice <<= ILattice
  template<class T1>
  inline
  ILattice& operator<<=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) <<= rhs.elem(i);

      return *this;
    }

  //! ILattice >>= ILattice
  template<class T1>
  inline
  ILattice& operator>>=(const ILattice<T1,N>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) >>= rhs.elem(i);

      return *this;
    }

#if 0
  // NOTE: intentially avoid defining a copy constructor - let the compiler
  // generate one via the bit copy mechanism. This effectively achieves
  // the first form of the if below (QDP_USE_ARRAY_INITIALIZER) without having
  // to use that syntax which is not strictly legal in C++.
#endif

  //! Deep copy constructor
#if defined(QDP_USE_ARRAY_INITIALIZER)
  ILattice(const ILattice& a) : F(a.F) {}
#else
  /*! This is a copy form - legal but not necessarily efficient */
  ILattice(const ILattice& a)
    {
      // fprintf(stderr,"copy ILattice\n");
      for(int i=0; i < N; ++i)
	F[i] = a.F[i];
    }
#endif

public:
  //! The backdoor
  /*! 
   * Used by optimization routines (e.g., SSE) that need the memory address of data.
   * BTW: to make this a friend would be a real pain since functions are templatized.
   */
  inline T* data() {return F;}


public:
  T& elem(int i) {return F[i];}
  const T& elem(int i) const {return F[i];}

private:
  /*! For now a fixed representation */
  T F[N];

};



//! Stream input
template<class T, int N>
inline
std::istream& operator>>(std::istream& s, ILattice<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s >> d.elem(i);

  return s;
}

//! Stream output
template<class T, int N>
inline
std::ostream& operator<<(std::ostream& s, const ILattice<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s << d.elem(i) << "\n";
  return s;
}


//! Text input
template<class T, int N>
inline
TextReader& operator>>(TextReader& s, ILattice<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s >> d.elem(i);
  return s;
}

//! Text output
template<class T, int N> 
inline  
TextWriter& operator<<(TextWriter& s, const ILattice<T,N>& d)
{
  for(int i=0; i < N; ++i)
    s << d.elem(i) << "\n";
  return s;
}


/*! @} */   // end of group ilattice


//-----------------------------------------------------------------------------
// Traits classes to support operations of simple scalars (floating constants, 
// etc.) on QDPTypes
//-----------------------------------------------------------------------------

template<class T>
struct WordType<IScalar<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};


template<class T, int N>
struct WordType<ILattice<T,N> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};

// Get the fixed single precision version of this Type (recursively)
template<class T>
struct SinglePrecType< IScalar<T> >
{
  typedef IScalar< typename SinglePrecType<T>::Type_t > Type_t;
};

// Get the fixed single precision version of this Type (recursively)
template<class T, int N>
struct SinglePrecType<ILattice<T,N> >
{
  typedef ILattice< typename SinglePrecType<T>::Type_t, N> Type_t;
};

// Get the fixed double precision version of this Type (recursively)
template<class T>
struct DoublePrecType< IScalar<T> >
{
  typedef IScalar< typename DoublePrecType<T>::Type_t> Type_t;
};

// Get the fixed double precision version of this Type (recursively)
template<class T, int N>
struct DoublePrecType<ILattice<T,N> >
{
  typedef ILattice< typename DoublePrecType<T>::Type_t, N> Type_t;
};


// Internally used scalars
template<class T>
struct InternalScalar<IScalar<T> > {
  typedef IScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

template<class T, int N>
struct InternalScalar<ILattice<T,N> > {
  typedef IScalar<typename InternalScalar<T>::Type_t>  Type_t;
};


// Trait to make a primitive scalar leaving grid along
template<class T>
struct PrimitiveScalar<IScalar<T> > {
  typedef IScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

template<class T, int N>
struct PrimitiveScalar<ILattice<T,N> > {
  typedef ILattice<typename PrimitiveScalar<T>::Type_t, N>  Type_t;
};


// Trait to make a lattice scalar leaving other indices alone
template<class T>
struct LatticeScalar<IScalar<T> > {
  typedef IScalar<typename LatticeScalar<T>::Type_t>  Type_t;
};

template<class T, int N>
struct LatticeScalar<ILattice<T,N> > {
  typedef IScalar<typename LatticeScalar<T>::Type_t>  Type_t;
};


// Internally used real scalars
template<class T>
struct RealScalar<IScalar<T> > {
  typedef IScalar<typename RealScalar<T>::Type_t>  Type_t;
};

template<class T, int N>
struct RealScalar<ILattice<T,N> > {
  typedef ILattice<typename RealScalar<T>::Type_t, N>  Type_t;
};



//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Default unary(IScalar) -> IScalar
template<class T1, class Op>
struct UnaryReturn<IScalar<T1>, Op> {
  typedef IScalar<typename UnaryReturn<T1, Op>::Type_t>  Type_t;
};

// Default unary(ILattice) -> ILattice
template<class T1, int N, class Op>
struct UnaryReturn<ILattice<T1,N>, Op> {
  typedef ILattice<typename UnaryReturn<T1, Op>::Type_t, N>  Type_t;
};

// Default binary(IScalar,IScalar) -> IScalar
template<class T1, class T2, class Op>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, Op> {
  typedef IScalar<typename BinaryReturn<T1, T2, Op>::Type_t>  Type_t;
};

// Currently, the only trinary operator is ``where'', so return 
// based on T2 and T3
// Default trinary(IScalar,IScalar,IScalar) -> IScalar
template<class T1, class T2, class T3, class Op>
struct TrinaryReturn<IScalar<T1>, IScalar<T2>, IScalar<T3>, Op> {
  typedef IScalar<typename BinaryReturn<T2, T3, Op>::Type_t>  Type_t;
};

// Default binary(ILattice,ILattice) -> ILattice
template<class T1, class T2, int N, class Op>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, Op> {
  typedef ILattice<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(IScalar,ILattice) -> ILattice
template<class T1, class T2, int N, class Op>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, Op> {
  typedef ILattice<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};

// Default binary(ILattice,IScalar) -> ILattice
template<class T1, class T2, int N, class Op>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, Op> {
  typedef ILattice<typename BinaryReturn<T1, T2, Op>::Type_t, N>  Type_t;
};


// Currently, the only trinary operator is ``where'', so return 
// based on T2 and T3

// Default trinary(ILattice,ILattice,ILattice) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, ILattice<T3,N>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};


// Default trinary(ILattice,IScalar,ILattice) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<ILattice<T1,N>, IScalar<T2>, ILattice<T3,N>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};

// Default trinary(ILattice,ILattice,IScalar) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, IScalar<T3>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};

// Default trinary(IScalar,ILattice,ILattice) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<IScalar<T1>, ILattice<T2,N>, ILattice<T3,N>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};


// Default trinary(IScalar,IScalar,ILattice) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<IScalar<T1>, IScalar<T2>, ILattice<T3,N>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};

// Default trinary(OSscalar,ILattice,IScalar) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<IScalar<T1>, ILattice<T2,N>, IScalar<T3>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};

// Default trinary(ILattice,IScalar,IScalar) -> ILattice
template<class T1, class T2, class T3, int N, class Op>
struct TrinaryReturn<ILattice<T1,N>, IScalar<T2>, IScalar<T3>, Op> {
  typedef ILattice<typename BinaryReturn<T2, T3, Op>::Type_t, N>  Type_t;
};




// Specific IScalar cases
// Global operations
template<class T>
struct UnaryReturn<IScalar<T>, FnPeekSite> {
  typedef IScalar<typename UnaryReturn<T, FnPeekSite>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<IScalar<T>, FnSumMulti> {
  typedef IScalar<typename UnaryReturn<T, FnSumMulti>::Type_t>  Type_t;
};


// Gamma algebra
template<int N, int m, class T2, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, IScalar<T2>, OpGammaConstMultiply> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<IScalar<T2>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, IScalar<T2>, OpGammaTypeMultiply> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaType>
struct BinaryReturn<IScalar<T2>, GammaType<N>, OpMultiplyGammaType> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};


// Gamma algebra
template<int N, int m, class T2, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, IScalar<T2>, OpGammaConstDPMultiply> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<IScalar<T2>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, IScalar<T2>, OpGammaTypeDPMultiply> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<IScalar<T2>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef IScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};



// Local operations

#if 0
template<class T1, class T2>
struct UnaryReturn<IScalar<T2>, OpCast<T1> > {
  typedef IScalar<typename UnaryReturn<T, OpCast>::Type_t>  Type_t;
//  typedef T1 Type_t;
};
#endif
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpAddAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpSubtractAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpMultiplyAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpDivideAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpModAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpModAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpBitwiseOrAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpBitwiseAndAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpBitwiseXorAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpLeftShiftAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpRightShiftAssign> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t>  &Type_t;
};
 
template<class T1, class T2, class T3>
struct TrinaryReturn<IScalar<T1>, IScalar<T2>, IScalar<T3>, FnColorContract> {
  typedef IScalar<typename TrinaryReturn<T1, T2, T3, FnColorContract>::Type_t>  Type_t;
};

// Specific ILattice cases
// Global operations
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnSumMulti> {
  typedef IScalar<typename UnaryReturn<T, FnSumMulti>::Type_t>  Type_t;
};


// Gamma algebra
template<int M, int m, class T2, int N, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<M,m>, ILattice<T2,N>, OpGammaConstMultiply> {
  typedef ILattice<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, int M, int m, class OpMultiplyGammaConst>
struct BinaryReturn<ILattice<T2,N>, GammaConst<M,m>, OpMultiplyGammaConst> {
  typedef ILattice<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<int M, class T2, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<M>, ILattice<T2,N>, OpGammaTypeMultiply> {
  typedef ILattice<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};

template<class T2, int N, int M, class OpMultiplyGammaType>
struct BinaryReturn<ILattice<T2,N>, GammaType<M>, OpMultiplyGammaType> {
  typedef ILattice<typename UnaryReturn<T2, OpUnaryPlus>::Type_t, N>  Type_t;
};



// Local operations
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, OpNot> {
  typedef ILattice<typename UnaryReturn<T, OpNot>::Type_t, N>  Type_t;
};


#if 0
template<class T1, class T2, int N>
struct UnaryReturn<ILattice<T2,N>, OpCast<T1>> {
  typedef ILattice<typename UnaryReturn<T, OpCast>::Type_t, N>  Type_t;
//  typedef T1 Type_t;
};
#endif

template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAddAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpSubtractAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMultiplyAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpDivideAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpModAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpModAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseOrAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseAndAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseXorAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLeftShiftAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpRightShiftAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, ILattice<T3,N>, FnColorContract> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnColorContract>::Type_t, N>  Type_t;
};


// Mixed ILattice & IScalar cases
// Global operations
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnPeekSite> {
  typedef IScalar<typename UnaryReturn<T, FnPeekSite>::Type_t>  Type_t;
};

// Local operations
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAddAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpAddAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpSubtractAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpSubtractAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMultiplyAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpMultiplyAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpDivideAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpDivideAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpModAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpModAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseOrAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseOrAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseAndAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseAndAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseXorAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseXorAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLeftShiftAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLeftShiftAssign>::Type_t, N>  &Type_t;
};
 
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpRightShiftAssign> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpRightShiftAssign>::Type_t, N>  &Type_t;
};




//-----------------------------------------------------------------------------
// Inner-grid scalar operations
//-----------------------------------------------------------------------------

/*! \addtogroup iscalar */
/*! @{ */

// Inner grid scalar

// ! IScalar
template<class T>
struct UnaryReturn<IScalar<T>, OpNot> {
  typedef IScalar<typename UnaryReturn<T, OpNot>::Type_t>  Type_t;
};

template<class T1>
inline typename UnaryReturn<IScalar<T1>, OpNot>::Type_t
operator!(const IScalar<T1>& l)
{
  return ! l.elem();
}


// +IScalar
template<class T1>
inline typename UnaryReturn<IScalar<T1>, OpUnaryPlus>::Type_t
operator+(const IScalar<T1>& l)
{
  return +l.elem();
}


// -IScalar
template<class T1>
inline typename UnaryReturn<IScalar<T1>, OpUnaryMinus>::Type_t
operator-(const IScalar<T1>& l)
{
  return -l.elem();
}


// IScalar + IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpAdd>::Type_t
operator+(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem()+r.elem();
}


// IScalar - IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpSubtract>::Type_t
operator-(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() - r.elem();
}


// IScalar * IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpMultiply>::Type_t
operator*(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() * r.elem();
}

// Optimized  adj(IScalar)*IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const IScalar<T1>& l, const IScalar<T2>& r)
{
  // Do not pass on the transpose or the conj
  return l.elem() * r.elem();
}

// Optimized  IScalar*adj(IScalar)
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const IScalar<T1>& l, const IScalar<T2>& r)
{
  // Do not pass on the transpose or the conj
  return l.elem() * r.elem();
}

// Optimized  adj(IScalar)*adj(IScalar)
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const IScalar<T1>& l, const IScalar<T2>& r)
{
  // Do not pass on the transpose or the conj
  return l.elem() * r.elem();
}


// IScalar / IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpDivide>::Type_t
operator/(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() / r.elem();
}


// IScalar << IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpLeftShift> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpLeftShift>::Type_t
operator<<(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() << r.elem();
}


// IScalar >> IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpRightShift> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpRightShift>::Type_t>  Type_t;
};
 

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpRightShift>::Type_t
operator>>(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() >> r.elem();
}


// IScalar % IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpMod>::Type_t
operator%(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() % r.elem();
}


// IScalar ^ IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpBitwiseXor>::Type_t
operator^(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() ^ r.elem();
}


// IScalar & IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() & r.elem();
}


// IScalar | IScalar
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpBitwiseOr>::Type_t
operator|(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() | r.elem();
}



// Comparisons

// IScalar < IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpLT> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpLT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpLT>::Type_t
operator<(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() < r.elem();
}


// IScalar <= IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpLE> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpLE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpLE>::Type_t
operator<=(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() <= r.elem();
}


// IScalar > IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpGT> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpGT>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpGT>::Type_t
operator>(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() > r.elem();
}


// IScalar >= IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpGE> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpGE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpGE>::Type_t
operator>=(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() >= r.elem();
}


// IScalar == IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpEQ> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpEQ>::Type_t
operator==(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() == r.elem();
}


// IScalar != IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpNE> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpNE>::Type_t
operator!=(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() != r.elem();
}


// IScalar && IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpAnd> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpAnd>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpAnd>::Type_t
operator&&(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() && r.elem();
}


// IScalar || IScalar
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, OpOr> {
  typedef IScalar<typename BinaryReturn<T1, T2, OpOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, OpOr>::Type_t
operator||(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return l.elem() || r.elem();
}



//-----------------------------------------------------------------------------
// Functions

// IScalar = adj(IScalar)
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnAdjoint>::Type_t
adj(const IScalar<T1>& s1)
{
  return adj(s1.elem());
}


// IScalar = conj(IScalar)
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnConjugate>::Type_t
conj(const IScalar<T1>& s1)
{
  return conj(s1.elem());
}


// IScalar = transpose(IScalar)
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnTranspose>::Type_t
transpose(const IScalar<T1>& s1)
{
  return transpose(s1.elem());
}


// IScalar = Trace(IScalar)
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnTrace>::Type_t
trace(const IScalar<T1>& s1)
{
  return trace(s1.elem());
}


// IScalar = Re(Trace(IScalar))
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnRealTrace>::Type_t
trace_real(const IScalar<T1>& s1)
{
  return trace_real(s1.elem());
}


// IScalar = Im(Trace(IScalar))
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnImagTrace>::Type_t
trace_imag(const IScalar<T1>& s1)
{
  return trace_imag(s1.elem());
}


// IScalar = Re(IScalar)
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnReal>::Type_t
real(const IScalar<T1>& s1)
{
  return real(s1.elem());
}


// IScalar = Im(IScalar)
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnImag>::Type_t
imag(const IScalar<T1>& s1)
{
  return imag(s1.elem());
}


// ArcCos
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnArcCos>::Type_t
acos(const IScalar<T1>& s1)
{
  return acos(s1.elem());
}


// ArcSin
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnArcSin>::Type_t
asin(const IScalar<T1>& s1)
{
  return asin(s1.elem());
}


// ArcTan
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnArcTan>::Type_t
atan(const IScalar<T1>& s1)
{
  return atan(s1.elem());
}


// Cos
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnCos>::Type_t
cos(const IScalar<T1>& s1)
{
  return cos(s1.elem());
}


// Exp
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnExp>::Type_t
exp(const IScalar<T1>& s1)
{
  return exp(s1.elem());
}


// Fabs
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnFabs>::Type_t
fabs(const IScalar<T1>& s1)
{
  return fabs(s1.elem());
}


// Log
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnLog>::Type_t
log(const IScalar<T1>& s1)
{
  return log(s1.elem());
}


// Sin
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnSin>::Type_t
sin(const IScalar<T1>& s1)
{
  return sin(s1.elem());
}


// Sqrt
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnSqrt>::Type_t
sqrt(const IScalar<T1>& s1)
{
  return sqrt(s1.elem());
}


// Tan
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnTan>::Type_t
tan(const IScalar<T1>& s1)
{
  return tan(s1.elem());
}


//! IScalar = pow(IScalar, IScalar)
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, FnPow>::Type_t
pow(const IScalar<T1>& s1, const IScalar<T2>& s2)
{
  return pow(s1.elem(), s2.elem());
}


//! IScalar = atan2(IScalar, IScalar)
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, FnArcTan2>::Type_t
atan2(const IScalar<T1>& s1, const IScalar<T2>& s2)
{
  return atan2(s1.elem(), s2.elem());
}


//! dest [float type] = source [seed type]
template<class T1>
inline typename UnaryReturn<IScalar<T1>, FnSeedToFloat>::Type_t
seedToFloat(const IScalar<T1>& s1)
{
  return seedToFloat(s1.elem());
}


//! IScalar = outerProduct(IScalar, IScalar)
template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const IScalar<T1>& l, const IScalar<T2>& r)
{
  return outerProduct(l.elem(),r.elem());
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T>
struct UnaryReturn<IScalar<T>, FnGetSite> {
  typedef IScalar<typename UnaryReturn<T, FnGetSite>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<IScalar<T>, FnGetSite>::Type_t
getSite(const IScalar<T>& s1, int innersite)
{
  return s1.elem();
}


//-----------------------------------------------------------------------------
//! QDP Int to int primitive in conversion routine
template<class T> 
inline int 
toInt(const IScalar<T>& s) 
{
  return toInt(s.elem());
}

//! QDP Real to float primitive in conversion routine
template<class T> 
inline float
toFloat(const IScalar<T>& s) 
{
  return toFloat(s.elem());
}

//! QDP Double to double primitive in conversion routine
template<class T> 
inline double
toDouble(const IScalar<T>& s) 
{
  return toDouble(s.elem());
}

//! QDP Boolean to bool primitive in conversion routine
template<class T> 
inline bool
toBool(const IScalar<T>& s) 
{
  return toBool(s.elem());
}

//! QDP Wordtype to primitive wordtype
template<class T> 
inline WordType<IScalar<T> > 
toWordType(const IScalar<T>& s) 
{
  return toWordType(s.elem());
}

//------------------------------------------
//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline
void copymask(IScalar<T>& d, const IScalar<T1>& mask, const IScalar<T>& s1) 
{
  copymask(d.elem(),mask.elem(),s1.elem());
}

//! dest [float type] = source [int type]
template<class T, class T1>
inline
void cast_rep(T& d, const IScalar<T1>& s1)
{
  cast_rep(d, s1.elem());
}


//! dest [float type] = source [int type]
template<class T, class T1>
inline
void recast_rep(IScalar<T>& d, const IScalar<T1>& s1)
{
  cast_rep(d.elem(), s1.elem());
}


//! dest [some type] = source [some type]
template<class T, class T1>
inline void 
copy_site(IScalar<T>& d, int isite, const IScalar<T1>& s1)
{
  d.elem() = s1.elem();
}


#if 0
// This should never be used and is probably an error if needed

//! gather several inner sites together
template<class T, class T1>
inline void 
gather_sites(IScalar<T>& d, 
	     const IScalar<T1>& s0, int i0, 
	     const IScalar<T1>& s1, int i1,
	     const IScalar<T1>& s2, int i2,
	     const IScalar<T1>& s3, int i3)
{
  gather_sites(d.elem(), 
	       s0.elem(), i0, 
	       s1.elem(), i1, 
	       s2.elem(), i2, 
	       s3.elem(), i3);
}
#endif


//------------------------------------------
// Global sum over site indices only
template<class T>
struct UnaryReturn<IScalar<T>, FnSum > {
  typedef IScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<IScalar<T>, FnSum>::Type_t
sum(const IScalar<T>& s1)
{
//  return sum(s1.elem());
  return s1.elem();
}


//------------------------------------------
// Global max
template<class T>
struct UnaryReturn<IScalar<T>, FnGlobalMax> {
  typedef IScalar<typename UnaryReturn<T, FnGlobalMax>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<IScalar<T>, FnGlobalMax>::Type_t
globalMax(const IScalar<T>& s1)
{
//  return globalMax(s1.elem());
  return s1.elem();
}


//------------------------------------------
// Global min
template<class T>
struct UnaryReturn<IScalar<T>, FnGlobalMin> {
  typedef IScalar<typename UnaryReturn<T, FnGlobalMin>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<IScalar<T>, FnGlobalMin>::Type_t
globalMin(const IScalar<T>& s1)
{
//  return globalMin(s1.elem());
  return s1.elem();
}


// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T>
struct UnaryReturn<IScalar<T>, FnNorm2> {
  typedef IScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T>
struct UnaryReturn<IScalar<T>, FnLocalNorm2> {
  typedef IScalar<typename UnaryReturn<T, FnLocalNorm2>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<IScalar<T>, FnLocalNorm2>::Type_t
localNorm2(const IScalar<T>& s1)
{
  return localNorm2(s1.elem());
}



//! IScalar = InnerProduct(adj(IScalar)*IScalar)
template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, FnInnerProduct> {
  typedef IScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
struct BinaryReturn<IScalar<T1>, IScalar<T2>, FnLocalInnerProduct> {
  typedef IScalar<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<IScalar<T1>, IScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const IScalar<T1>& s1, const IScalar<T2>& s2)
{
  return localInnerProduct(s1.elem(), s2.elem());
}


//! IScalar = where(IScalar, IScalar, IScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3>
struct TrinaryReturn<IScalar<T1>, IScalar<T2>, IScalar<T3>, FnWhere> {
  typedef IScalar<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t>  Type_t;
};

template<class T1, class T2, class T3>
inline typename TrinaryReturn<IScalar<T1>, IScalar<T2>, IScalar<T3>, FnWhere>::Type_t
where(const IScalar<T1>& a, const IScalar<T2>& b, const IScalar<T3>& c)
{
  return where(a.elem(), b.elem(), c.elem());
}



//-----------------------------------------------------------------------------
// Broadcast operations
//! dest = 0
template<class T> 
inline
void zero_rep(IScalar<T>& dest) 
{
  zero_rep(dest.elem());
}


//-----------------------------------------------------------------------------
// Random numbers
//! dest  = random  
template<class T, class T1, class T2>
inline void
fill_random(IScalar<T>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  fill_random(d.elem(), seed, skewed_seed, seed_mult);
}


//! dest  = gaussian  
template<class T>
inline void
fill_gaussian(IScalar<T>& d, IScalar<T>& r1, IScalar<T>& r2)
{
  fill_gaussian(d.elem(), r1.elem(), r2.elem());
}


/*! @} */  // end of group iscalar

//-----------------------------------------------------------------------------
// Inner-grid lattice operations
//-----------------------------------------------------------------------------

/*! \addtogroup ilattice */
/*! @{ */

// Ilattice operations

//! ILattice = ! ILattice
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, OpNot>::Type_t
operator!(const ILattice<T1,N>& l)
{
  typename UnaryReturn<ILattice<T1,N>, OpNot>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = ! l.elem(i);
  return d;
}


//! ILattice = +ILattice
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, OpUnaryPlus>::Type_t
operator+(const ILattice<T1,N>& l)
{
  typename UnaryReturn<ILattice<T1,N>, OpUnaryPlus>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = +l.elem(i);
  return d;
}


//! ILattice = -ILattice
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, OpUnaryMinus>::Type_t
operator-(const ILattice<T1,N>& l)
{
  typename UnaryReturn<ILattice<T1,N>, OpUnaryMinus>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = -l.elem(i);
  return d;
}


// ILattice + ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAdd>::Type_t
operator+(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) + r.elem(i);
  return d;
}

// ILattice + IScalar
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAdd>::Type_t
operator+(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) + r.elem();
  return d;
}

// IScalar + ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAdd>::Type_t
operator+(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAdd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() + r.elem(i);
  return d;
}


// ILattice - ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpSubtract>::Type_t
operator-(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) - r.elem(i);
  return d;
}

// ILattice - IScalar
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpSubtract>::Type_t
operator-(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) - r.elem();
  return d;
}

// IScalar - ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpSubtract>::Type_t
operator-(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpSubtract>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() - r.elem(i);
  return d;
}


// ILattice * ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMultiply>::Type_t
operator*(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}

// ILattice * IScalar
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMultiply>::Type_t
operator*(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

// IScalar * ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMultiply>::Type_t
operator*(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMultiply>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}


// Optimized  adj(ILattice)*ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAdjMultiply>::Type_t
adjMultiply(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAdjMultiply>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}

// Optimized  adj(ILattice)*IScalar
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAdjMultiply>::Type_t
adjMultiply(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAdjMultiply>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

// Optimized  adj(IScalar)*ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAdjMultiply>::Type_t
adjMultiply(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAdjMultiply>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}


// Optimized  ILattice*adj(ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMultiplyAdj>::Type_t
multiplyAdj(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMultiplyAdj>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}

// Optimized  ILattice*adj(IScalar)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMultiplyAdj>::Type_t
multiplyAdj(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMultiplyAdj>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem();
  return d;
}

// Optimized  IScalar*adj(ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMultiplyAdj>::Type_t
multiplyAdj(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMultiplyAdj>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() * r.elem(i);
  return d;
}


// Optimized  adj(ILattice)*adj(ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAdjMultiplyAdj>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}

// Optimized  adj(ILattice)*adj(IScalar)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAdjMultiplyAdj>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}

// Optimized  adj(IScalar)*adj(ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAdjMultiplyAdj>::Type_t  d;

  // Do not pass on the transpose or the conj
  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) * r.elem(i);
  return d;
}


// ILattice / ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpDivide>::Type_t
operator/(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpDivide>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) / r.elem(i);
  return d;
}

// ILattice / IScalar
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpDivide>::Type_t
operator/(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpDivide>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) / r.elem();
  return d;
}

// IScalar / ILattice
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpDivide>::Type_t
operator/(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpDivide>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() / r.elem(i);
  return d;
}


//! ILattice << ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLeftShift> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLeftShift>::Type_t
operator<<(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLeftShift>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) << r.elem(i);
  return d;
}

//! ILattice << IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLeftShift> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLeftShift>::Type_t
operator<<(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLeftShift>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) << r.elem();
  return d;
}

//! IScalar << ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLeftShift> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLeftShift>::Type_t
operator<<(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLeftShift>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() << r.elem(i);
  return d;
}


//! ILattice >> ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpRightShift> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpRightShift>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpRightShift>::Type_t
operator>>(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpRightShift>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) >> r.elem(i);
  return d;
}

// ILattice >> IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpRightShift> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpRightShift>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpRightShift>::Type_t
operator>>(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpRightShift>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) >> r.elem();
  return d;
}

// IScalar >> ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpRightShift> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpRightShift>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpRightShift>::Type_t
operator>>(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpRightShift>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() >> r.elem(i);
  return d;
}


// ILattice % ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMod> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpMod>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMod>::Type_t
operator%(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpMod>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) % r.elem(i);
  return d;
}

// ILattice % IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMod> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpMod>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMod>::Type_t
operator%(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpMod>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) % r.elem();
  return d;
}

// IScalar % ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMod> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpMod>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMod>::Type_t
operator%(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpMod>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() % r.elem(i);
  return d;
}


// ILattice ^ ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseXor> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseXor>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseXor>::Type_t
operator^(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseXor>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) ^ r.elem(i);
  return d;
}

// ILattice ^ IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseXor> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseXor>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseXor>::Type_t
operator^(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseXor>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) ^ r.elem();
  return d;
}

// IScalar ^ ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseXor> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseXor>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseXor>::Type_t
operator^(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseXor>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() ^ r.elem(i);
  return d;
}


// ILattice & IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseAnd> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseAnd>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseAnd>::Type_t
operator&(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseAnd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) & r.elem(i);
  return d;
}

// ILattice & IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseAnd> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseAnd>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseAnd>::Type_t
operator&(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseAnd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) & r.elem();
  return d;
}

// IScalar & ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseAnd> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseAnd>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseAnd>::Type_t
operator&(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseAnd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() & r.elem(i);
  return d;
}


// ILattice | IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseOr> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseOr>::Type_t
operator|(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpBitwiseOr>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) | r.elem(i);
  return d;
}

// ILattice | IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseOr> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseOr>::Type_t
operator|(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpBitwiseOr>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) | r.elem();
  return d;
}

// IScalar | ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseOr> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t, N>  Type_t;
};
 
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseOr>::Type_t
operator|(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpBitwiseOr>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() | r.elem(i);
  return d;
}


// Comparisons

// ILattice < ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLT> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLT>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLT>::Type_t
operator<(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) < r.elem(i);
  return d;
}

// ILattice < IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLT> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLT>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLT>::Type_t
operator<(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) < r.elem();
  return d;
}

// IScalar < ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLT> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLT>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLT>::Type_t
operator<(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() < r.elem(i);
  return d;
}


// ILattice <= ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLE>::Type_t
operator<=(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpLE>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) <= r.elem(i);
  return d;
}

// ILattice <= IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLE>::Type_t
operator<=(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpLE>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) <= r.elem();
  return d;
}

// IScalar <= ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpLE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLE>::Type_t
operator<=(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpLE>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() <= r.elem(i);
  return d;
}


// ILattice > ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGT> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpGT>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGT>::Type_t
operator>(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) > r.elem(i);
  return d;
}

// ILattice > IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGT> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpGT>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGT>::Type_t
operator>(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) > r.elem();
  return d;
}

// ILattice > ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGT> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpGT>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGT>::Type_t
operator>(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() > r.elem(i);
  return d;
}



// ILattice >= ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpGE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGE>::Type_t
operator>=(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) >= r.elem(i);
  return d;
}

// ILattice >= IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpGE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGE>::Type_t
operator>=(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) >= r.elem();
  return d;
}

// IScalar >= ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpGE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGE>::Type_t
operator>=(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() >= r.elem(i);
  return d;
}


// ILattice == ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpEQ> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpEQ>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpEQ>::Type_t
operator==(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) == r.elem(i);
  return d;
}

// ILattice == IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpEQ> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpEQ>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpEQ>::Type_t
operator==(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) == r.elem();
  return d;
}

// IScalar == ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpEQ> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpEQ>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpEQ>::Type_t
operator==(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpGT>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() == r.elem(i);
  return d;
}


// ILattice != ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpNE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpNE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpNE>::Type_t
operator!=(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpNE>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) != r.elem(i);
  return d;
}

// ILattice != IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpNE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpNE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpNE>::Type_t
operator!=(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpNE>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) != r.elem();
  return d;
}

// ILattice != ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpNE> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpNE>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpNE>::Type_t
operator!=(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpNE>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() != r.elem(i);
  return d;
}


// ILattice && ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAnd> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpAnd>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAnd>::Type_t
operator&&(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpAnd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) && r.elem(i);
  return d;
}

// ILattice && IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAnd> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpAnd>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAnd>::Type_t
operator&&(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpAnd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) && r.elem();
  return d;
}

// IScalar && ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAnd> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpAnd>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAnd>::Type_t
operator&&(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpAnd>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() && r.elem(i);
  return d;
}


// ILattice || ILattice
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpOr> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpOr>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpOr>::Type_t
operator||(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, OpOr>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) || r.elem(i);
  return d;
}

// ILattice || IScalar
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpOr> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpOr>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpOr>::Type_t
operator||(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, OpOr>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem(i) || r.elem();
  return d;
}

// IScalar || ILattice
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpOr> {
  typedef ILattice<typename BinaryReturn<T1, T2, OpOr>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpOr>::Type_t
operator||(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, OpOr>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = l.elem() || r.elem(i);
  return d;
}


//-----------------------------------------------------------------------------
// Functions

// ILattice = adj(ILattice)
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnAdjoint>::Type_t
adj(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnAdjoint>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = adj(s1.elem(i));
  return d;
}


// ILattice = conj(ILattice)
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnConjugate>::Type_t
conj(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnConjugate>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = conj(s1.elem(i));
  return d;
}


// ILattice = transpose(ILattice)
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnTranspose>::Type_t
transpose(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnTranspose>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = transpose(s1.elem(i));
  return d;
}


// ILattice = Trace(ILattice)
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnTrace>::Type_t
trace(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnTrace>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = trace(s1.elem(i));
  return d;
}


// ILattice = Re(Trace(ILattice))
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnRealTrace>::Type_t
trace_real(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnRealTrace>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = trace_real(s1.elem(i));
  return d;
}


// ILattice = Im(Trace(ILattice))
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnImagTrace>::Type_t
trace_imag(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnImagTrace>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = trace_imag(s1.elem(i));
  return d;
}


// ILattice = Re(ILattice)
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnReal>::Type_t
real(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnReal>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = real(s1.elem(i));
  return d;
}


// ILattice = Im(ILattice)
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnImag>::Type_t
imag(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnImag>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = imag(s1.elem(i));
  return d;
}


// ArcCos
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnArcCos>::Type_t
acos(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnArcCos>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = acos(s1.elem(i));
  return d;
}

// ArcSin
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnArcSin>::Type_t
asin(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnArcSin>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = asin(s1.elem(i));
  return d;
}

// ArcTan
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnArcTan>::Type_t
atan(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnArcTan>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = atan(s1.elem(i));
  return d;
}

// Cos
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnCos>::Type_t
cos(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnCos>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = cos(s1.elem(i));
  return d;
}

// Exp
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnExp>::Type_t
exp(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnExp>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = exp(s1.elem(i));
  return d;
}

// Fabs
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnFabs>::Type_t
fabs(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnFabs>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = fabs(s1.elem(i));
  return d;
}

// Log
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnLog>::Type_t
log(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnLog>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = log(s1.elem(i));
  return d;
}

// Sin
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnSin>::Type_t
sin(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnSin>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = sin(s1.elem(i));
  return d;
}

// Sqrt
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnSqrt>::Type_t
sqrt(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnSqrt>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = sqrt(s1.elem(i));
  return d;
}

// Tan
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnTan>::Type_t
tan(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnTan>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = tan(s1.elem(i));
  return d;
}


//! ILattice = pow(ILattice, ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnPow>::Type_t
pow(const ILattice<T1,N>& s1, const ILattice<T2,N>& s2)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnPow>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem(i));
  return d;
}

//! ILattice = pow(ILattice, IScalar)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnPow>::Type_t
pow(const ILattice<T1,N>& s1, const IScalar<T2>& s2)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnPow>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem());
  return d;
}

//! ILattice = pow(IScalar, ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnPow>::Type_t
pow(const IScalar<T1>& s1, const ILattice<T2,N>& s2)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnPow>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = pow(s1.elem(i), s2.elem(i));
  return d;
}


//! ILattice = atan2(ILattice, ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnArcTan2>::Type_t
atan2(const ILattice<T1,N>& s1, const ILattice<T2,N>& s2)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnArcTan2>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = atan2(s1.elem(i), s2.elem(i));
  return d;
}

//! ILattice = atan2(ILattice, IScalar)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnArcTan2>::Type_t
atan2(const ILattice<T1,N>& s1, const IScalar<T2>& s2)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnArcTan2>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = atan2(s1.elem(i), s2.elem());
  return d;
}

//! ILattice = atan2(IScalar, ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnArcTan2>::Type_t
atan2(const IScalar<T1>& s1, const ILattice<T2,N>& s2)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnArcTan2>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = atan2(s1.elem(), s2.elem(i));
  return d;
}


//! dest [float type] = source [seed type]
template<class T1, int N>
inline typename UnaryReturn<ILattice<T1,N>, FnSeedToFloat>::Type_t
seedToFloat(const ILattice<T1,N>& s1)
{
  typename UnaryReturn<ILattice<T1,N>, FnSeedToFloat>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = seedToFloat(s1.elem(i));
  return d;
}


//! ILattice = outerProduct(ILattice, ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnOuterProduct>::Type_t
outerProduct(const ILattice<T1,N>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnOuterProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = outerProduct(l.elem(i), r.elem(i));
  return d;
}

//! ILattice = outerProduct(ILattice, IScalar)
template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnOuterProduct>::Type_t
outerProduct(const ILattice<T1,N>& l, const IScalar<T2>& r)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnOuterProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = outerProduct(l.elem(i), r.elem());
  return d;
}

//! ILattice = outerProduct(IScalar, ILattice)
template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnOuterProduct>::Type_t
outerProduct(const IScalar<T1>& l, const ILattice<T2,N>& r)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnOuterProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = outerProduct(l.elem(), r.elem(i));
  return d;
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
// Global operations
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnGetSite> {
  typedef IScalar<typename UnaryReturn<T, FnGetSite>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<ILattice<T,N>, FnGetSite>::Type_t
getSite(const ILattice<T,N>& s1, int innersite)
{
  return s1.elem(innersite);
}


//! dest = 0
template<class T, int N> 
inline void 
zero_rep(ILattice<T,N>& dest) 
{
  for(int i=0; i < N; ++i)
    zero_rep(dest.elem(i));
}


//! dest [some type] = source [some type]
template<class T, class T1, int N>
inline void 
copy_site(ILattice<T,N>& d, int isite, const IScalar<T1>& s1)
{
  d.elem(isite) = s1.elem();
}


//! gather several inner sites together
/*! This version is built for inner length of 2 */
template<class T, class T1>
inline void 
gather_sites(ILattice<T,2>& d, 
	     const ILattice<T1,2>& s0, int i0, 
	     const ILattice<T1,2>& s1, int i1)
{
  d.elem(0) = s0.elem(i0);
  d.elem(1) = s0.elem(i1);
}


//! gather several inner sites together
/*! This version is built for inner length of 4 */
template<class T, class T1>
inline void 
gather_sites(ILattice<T,4>& d, 
	     const ILattice<T1,4>& s0, int i0, 
	     const ILattice<T1,4>& s1, int i1,
	     const ILattice<T1,4>& s2, int i2,
	     const ILattice<T1,4>& s3, int i3)
{
  d.elem(0) = s0.elem(i0);
  d.elem(1) = s1.elem(i1);
  d.elem(2) = s2.elem(i2);
  d.elem(3) = s3.elem(i3);
}


//! dest = (mask) ? s1 : dest
template<class T, class T1, int N> 
inline void 
copymask(ILattice<T,N>& d, const ILattice<T1,N>& mask, const ILattice<T,N>& s1) 
{
  for(int i=0; i < N; ++i)
    copymask(d.elem(i),mask.elem(i),s1.elem(i));
}


//! dest  = random  
template<class T, int N, class T1, class T2>
inline void
fill_random(ILattice<T,N>& d, T1& seed, T2& skewed_seed, const T1& seed_mult)
{
  fill_random<T1,T2,N>(d.data(), seed, skewed_seed, seed_mult);
}


#if 1
// Global sum over site indices only
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnSum> {
  typedef IScalar<typename UnaryReturn<T, FnSum>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<ILattice<T,N>, FnSum>::Type_t
sum(const ILattice<T,N>& s1)
{
  typename UnaryReturn<ILattice<T,N>, FnSum>::Type_t  d;

  d.elem() = s1.elem(0);
  for(int i=1; i < N; ++i)
    d.elem() += s1.elem(i);

  return d;
}
#endif


//--------------------------------------------
// Global max
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnGlobalMax> {
  typedef IScalar<typename UnaryReturn<T, FnGlobalMax>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<ILattice<T,N>, FnGlobalMax>::Type_t
sum(const ILattice<T,N>& s1)
{
  typename UnaryReturn<ILattice<T,N>, FnGlobalMax>::Type_t  d;

  d.elem() = s1.elem(0);
  for(int i=1; i < N; ++i)
  {
    if (toBool(s1.elem() > d.elem()))
      d.elem() = s1.elem(i);
  }

  return d;
}


//--------------------------------------------
// Global min
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnGlobalMin> {
  typedef IScalar<typename UnaryReturn<T, FnGlobalMin>::Type_t>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<ILattice<T,N>, FnGlobalMin>::Type_t
sum(const ILattice<T,N>& s1)
{
  typename UnaryReturn<ILattice<T,N>, FnGlobalMin>::Type_t  d;

  d.elem() = s1.elem(0);
  for(int i=1; i < N; ++i)
  {
    if (toBool(s1.elem() > d.elem()))
      d.elem() = s1.elem(i);
  }

  return d;
}


//--------------------------------------------
// InnerProduct (norm-seq) global sum = sum(tr(adj(s1)*s1))
template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnNorm2> {
  typedef IScalar<typename UnaryReturn<T, FnNorm2>::Type_t>  Type_t;
};

template<class T, int N>
struct UnaryReturn<ILattice<T,N>, FnLocalNorm2> {
  typedef ILattice<typename UnaryReturn<T, FnLocalNorm2>::Type_t, N>  Type_t;
};

template<class T, int N>
inline typename UnaryReturn<ILattice<T,N>, FnLocalNorm2>::Type_t
localNorm2(const ILattice<T,N>& s1)
{
  typename UnaryReturn<ILattice<T,N>, FnLocalNorm2>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = localNorm2(s1.elem(i));

  return d;
}


//! IScalar = InnerProduct(adj(ILattice)*ILattice)
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnInnerProduct> {
  typedef IScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnLocalInnerProduct> {
  typedef ILattice<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnLocalInnerProduct>::Type_t
localInnerProduct(const ILattice<T1,N>& s1, const ILattice<T2,N>& s2)
{
  typename BinaryReturn<ILattice<T1,N>, ILattice<T2,N>, FnLocalInnerProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = localInnerProduct(s1.elem(i), s2.elem(i));

  return d;
}

//! IScalar = InnerProduct(adj(ILattice)*IScalar)
template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnInnerProduct> {
  typedef IScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnLocalInnerProduct> {
  typedef ILattice<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnLocalInnerProduct>::Type_t
localInnerProduct(const ILattice<T1,N>& s1, const IScalar<T2>& s2)
{
  typename BinaryReturn<ILattice<T1,N>, IScalar<T2>, FnLocalInnerProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = localInnerProduct(s1.elem(i), s2.elem());

  return d;
}

//! IScalar = InnerProduct(adj(IScalar)*ILattice)
template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnInnerProduct> {
  typedef IScalar<typename BinaryReturn<T1, T2, FnInnerProduct>::Type_t>  Type_t;
};

template<class T1, class T2, int N>
struct BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnLocalInnerProduct> {
  typedef ILattice<typename BinaryReturn<T1, T2, FnLocalInnerProduct>::Type_t, N>  Type_t;
};

template<class T1, class T2, int N>
inline typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnLocalInnerProduct>::Type_t
localInnerProduct(const IScalar<T1>& s1, const ILattice<T2,N>& s2)
{
  typename BinaryReturn<IScalar<T1>, ILattice<T2,N>, FnLocalInnerProduct>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = localInnerProduct(s1.elem(), s2.elem(i));

  return d;
}


//! ILattice = where(ILattice, ILattice, ILattice)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, ILattice<T3,N>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, ILattice<T3,N>, FnWhere>::Type_t
where(const ILattice<T1,N>& a, const ILattice<T2,N>& b, const ILattice<T3,N>& c)
{
  typename TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, ILattice<T3,N>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(i), b.elem(i), c.elem(i));

  return d;
}

//! ILattice = where(ILattice, ILattice, IScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, IScalar<T3>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, IScalar<T3>, FnWhere>::Type_t
where(const ILattice<T1,N>& a, const ILattice<T2,N>& b, const IScalar<T3>& c)
{
  typename TrinaryReturn<ILattice<T1,N>, ILattice<T2,N>, IScalar<T3>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(i), b.elem(i), c.elem());

  return d;
}

//! ILattice = where(ILattice, IScalar, ILattice)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<ILattice<T1,N>, IScalar<T2>, ILattice<T3,N>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<ILattice<T1,N>, IScalar<T2>, ILattice<T3,N>, FnWhere>::Type_t
where(const ILattice<T1,N>& a, const IScalar<T2>& b, const ILattice<T3,N>& c)
{
  typename TrinaryReturn<ILattice<T1,N>, IScalar<T2>, ILattice<T3,N>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(i), b.elem(), c.elem(i));

  return d;
}

//! ILattice = where(ILattice, IScalar, IScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<ILattice<T1,N>, IScalar<T2>, IScalar<T3>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<ILattice<T1,N>, IScalar<T2>, IScalar<T3>, FnWhere>::Type_t
where(const ILattice<T1,N>& a, const IScalar<T2>& b, const IScalar<T3>& c)
{
  typename TrinaryReturn<ILattice<T1,N>, IScalar<T2>, IScalar<T3>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(i), b.elem(), c.elem());

  return d;
}

//! ILattice = where(IScalar, ILattice, ILattice)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<IScalar<T1>, ILattice<T2,N>, ILattice<T3,N>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<IScalar<T1>, ILattice<T2,N>, ILattice<T3,N>, FnWhere>::Type_t
where(const IScalar<T1>& a, const ILattice<T2,N>& b, const ILattice<T3,N>& c)
{
  typename TrinaryReturn<IScalar<T1>, ILattice<T2,N>, ILattice<T3,N>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem(i));

  return d;
}

//! ILattice = where(IScalar, ILattice, IScalar)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<IScalar<T1>, ILattice<T2,N>, IScalar<T3>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<IScalar<T1>, ILattice<T2,N>, IScalar<T3>, FnWhere>::Type_t
where(const IScalar<T1>& a, const ILattice<T2,N>& b, const IScalar<T3>& c)
{
  typename TrinaryReturn<IScalar<T1>, ILattice<T2,N>, IScalar<T3>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(), b.elem(i), c.elem());

  return d;
}

//! ILattice = where(IScalar, IScalar, ILattice)
/*!
 * Where is the ? operation
 * returns  (a) ? b : c;
 */
template<class T1, class T2, class T3, int N>
struct TrinaryReturn<IScalar<T1>, IScalar<T2>, ILattice<T3,N>, FnWhere> {
  typedef ILattice<typename TrinaryReturn<T1, T2, T3, FnWhere>::Type_t, N>  Type_t;
};

template<class T1, class T2, class T3, int N>
inline typename TrinaryReturn<IScalar<T1>, IScalar<T2>, ILattice<T3,N>, FnWhere>::Type_t
where(const IScalar<T1>& a, const IScalar<T2>& b, const ILattice<T3,N>& c)
{
  typename TrinaryReturn<IScalar<T1>, IScalar<T2>, ILattice<T3,N>, FnWhere>::Type_t  d;

  for(int i=0; i < N; ++i)
    d.elem(i) = where(a.elem(), b.elem(), c.elem(i));

  return d;
}


/*! @} */  // end of group ilattice

} // namespace QDP

#endif
