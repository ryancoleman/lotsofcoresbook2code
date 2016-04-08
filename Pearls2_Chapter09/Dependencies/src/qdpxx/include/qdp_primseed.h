// -*- C++ -*-

/*! \file
 * \brief Primitive Seed
 */


#ifndef QDP_PRIMSEED_H
#define QDP_PRIMSEED_H

namespace QDP {

//-------------------------------------------------------------------------------------
/*! \addtogroup primseed Seed primitive
 * \ingroup fiber
 *
 * Primitive type for supporting random numbers. This is really a 
 * big integer class.
 *
 * @{
 */


//! Primitive Seed class
/*!
   * Seed primitives exist to facilitate seed multiplication - namely
   * multiplication of a large integer represented as 4 smaller integers
   *
   * NOTE: the class is still treated as a template since there may be
   * an inner lattice so the type could be represented as an length 4
   * array of lattice integers
   */
template <class T> class PSeed
{
public:
  PSeed() {}
  ~PSeed() {}

  //! construct dest = const
  template<class T1>
  PSeed(const PScalar<T1>& rhs) 
    {
      assign(rhs);
    }


  //! PSeed = PScalar
  /*! Set equal to input scalar (an integer) */
  template<class T1>
  inline
  PSeed& assign(const PScalar<T1>& rhs) 
    {
      typedef typename InternalScalar<T1>::Type_t  S;

      elem(0) = rhs.elem() & S(4095);
      elem(1) = (rhs.elem() >> S(12)) & S(4095);
      elem(2) = (rhs.elem() >> S(24)) & S(4095);
//      elem(3) = (rhs.elem() >> S(36)) & S(2047);  // This probably will never be nonzero
      zero_rep(elem(3));    // assumes 32 bit integers

      return *this;
    }

  //! PSeed = PScalar
  /*! Set equal to input scalar (an integer) */
  template<class T1>
  inline
  PSeed& operator=(const PScalar<T1>& rhs) 
    {
      return assign(rhs);
    }

  //! PSeed = PSeed
  /*! Set equal to another PSeed */
  template<class T1>
  inline
  PSeed& operator=(const PSeed<T1>& rhs) 
    {
      for(int i=0; i < 4; ++i)
	elem(i) = rhs.elem(i);

      return *this;
    }

  //! Deep copy constructor
#if defined(QDP_USE_ARRAY_INITIALIZER)
  /*! This is an array initializer form - may not be strictly legal */
  PSeed(const PSeed& a) : F(a.F) {}
#else
  /*! This is a copy form - legal but not necessarily efficient */
  PSeed(const PSeed& a)
    {
      for(int i=0; i < 4; ++i)
	F[i] = a.F[i];
    }
#endif

public:
  T& elem(int i) {return F[i];}
  const T& elem(int i) const {return F[i];}

private:
  T F[4];
};


//! Text input
template<class T>
inline
std::istream& operator>>(std::istream& s, PSeed<T>& d)
{
  for(int i=0; i < 4; ++i)
    s >> d.elem(i);

  return s;
}

//! Text input
template<class T>
inline
StandardInputStream& operator>>(StandardInputStream& s, PSeed<T>& d)
{
  for(int i=0; i < 4; ++i)
    s >> d.elem(i);

  return s;
}

//! Text output
template<class T>
inline
std::ostream& operator<<(std::ostream& s, const PSeed<T>& d)
{
  s << d.elem(0) << " " << d.elem(1) << " " << d.elem(2) << " " << d.elem(3) << "\n";
  return s;
}

//! Text output
template<class T>
inline
StandardOutputStream& operator<<(StandardOutputStream& s, const PSeed<T>& d)
{
  s << d.elem(0) << " " << d.elem(1) << " " << d.elem(2) << " " << d.elem(3) << "\n";
  return s;
}

//! Text input
template<class T>
inline
TextReader& operator>>(TextReader& txt, PSeed<T>& d)
{
  for(int i=0; i < 4; ++i)
    txt >> d.elem(i);

  return txt;
}

//! Text output
template<class T>
inline
TextWriter& operator<<(TextWriter& txt, const PSeed<T>& d)
{
  for(int i=0; i < 4; ++i)
    txt << d.elem(i) << "\n";

  return txt;
}

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T>
inline
XMLWriter& operator<<(XMLWriter& xml, const PSeed<T>& d)
{
  xml.openTag("Seed");

  // Copy into another array first
  for(int i=0; i < 4; ++i)
  {
    xml.openTag("elem");
    xml << d.elem(i);
    xml.closeTag();
  }

  xml.closeTag();  // Seed
  return xml;
}


//! XML input
template<class T>
inline
void read(XMLReader& xml, const std::string& path, PSeed<T>& d)
{
  typedef typename PrimitiveScalar<T>::Type_t  S;
  multi1d<S> ff(4);

  read(xml, path + "/Seed", ff);
  
  for(int i=0; i < 4; ++i)
  {
    d.elem(i) = S(ff[i]);
  }
}

#endif
/*! @} */   // end of group primseed

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Underlying word type
template<class T1>
struct WordType<PSeed<T1> > 
{
  typedef typename WordType<T1>::Type_t  Type_t;
};

// Fixed Precision versions (do these even make sense? )

template<class T1>
struct SinglePrecType<PSeed<T1> >
{
  typedef PSeed< typename SinglePrecType<T1>::Type_t > Type_t;
};

template<class T1>
struct DoublePrecType<PSeed<T1> >
{
  typedef PSeed< typename DoublePrecType<T1>::Type_t > Type_t;
};


// Internally used scalars
template<class T>
struct InternalScalar<PSeed<T> > {
  typedef PScalar<typename InternalScalar<T>::Type_t>  Type_t;
};

// Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar<PSeed<T> > {
  typedef PScalar<typename PrimitiveScalar<T>::Type_t>  Type_t;
};

// Makes a lattice scalar leaving primitive indices alone
template<class T>
struct LatticeScalar<PSeed<T> > {
  typedef PSeed<typename LatticeScalar<T>::Type_t>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// Assignment is different
template<class T1, class T2 >
struct BinaryReturn<PSeed<T1>, PSeed<T2>, OpAssign > {
  typedef PSeed<T1> &Type_t;
};

 

//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

// PScalar = (PSeed == PSeed)
template<class T1, class T2>
struct BinaryReturn<PSeed<T1>, PSeed<T2>, OpEQ> {
  typedef PScalar<typename BinaryReturn<T1, T2, OpEQ>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PSeed<T1>, PSeed<T2>, OpEQ>::Type_t
operator==(const PSeed<T1>& l, const PSeed<T2>& r)
{
  return 
    (l.elem(0) == r.elem(0)) && 
    (l.elem(1) == r.elem(1)) && 
    (l.elem(2) == r.elem(2)) && 
    (l.elem(3) == r.elem(3));
}


// PScalar = (Seed != Seed)
template<class T1, class T2>
struct BinaryReturn<PSeed<T1>, PSeed<T2>, OpNE> {
  typedef PScalar<typename BinaryReturn<T1, T2, OpNE>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PSeed<T1>, PSeed<T2>, OpNE>::Type_t
operator!=(const PSeed<T1>& l, const PSeed<T2>& r)
{
  return 
    (l.elem(0) != r.elem(0)) ||
    (l.elem(1) != r.elem(1)) || 
    (l.elem(2) != r.elem(2)) || 
    (l.elem(3) != r.elem(3));
}


/*! \addtogroup primseed
 * @{ 
 */

// Primitive Seeds

//! PSeed<T> = PSeed<T> * PSeed<T>
/*!
 * A 47 bit seed multiplication is represented as the multiplication
 * of three 12bit and one 11bit integer
 *
 * i3 = s1(3)*s2(0) + s1(2)*s2(1)
 *    + s1(1)*s2(2) + s1(0)*s2(3);
 * i2 = s1(2)*s2(0) + s1(1)*s2(1)
 *    + s1(0)*s2(2);
 * i1 = s1(1)*s2(0) + s1(0)*s2(1);
 * i0 = s1(0)*s2(0);
 *
 * dest(0) = mod(i0, 4096);
 * i1      = i1 + i0/4096;
 * dest(1) = mod(i1, 4096);
 * i2      = i2 + i1/4096;
 * dest(2) = mod(i2, 4096);
 * i3      = i3 + i2/4096
 * dest(3) = mod(i3, 2048);
 */
template<class T1, class T2>
struct BinaryReturn<PSeed<T1>, PSeed<T2>, OpMultiply> {
  typedef PSeed<typename BinaryReturn<T1, T2, OpMultiply>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PSeed<T1>, PSeed<T2>, OpMultiply>::Type_t
operator*(const PSeed<T1>& s1, const PSeed<T2>& s2)
{
  typename BinaryReturn<PSeed<T1>, PSeed<T2>, OpMultiply>::Type_t  d;
  typedef typename BinaryReturn<T1, T2, OpMultiply>::Type_t  T;
  typedef typename InternalScalar<T>::Type_t  S;
  T  i0, i1, i2, i3;

  /* i3 = s1(3)*s2(0) + s1(2)*s2(1) + s1(1)*s2(2) + s1(0)*s2(3) */
  i3  = s1.elem(3) * s2.elem(0);
  i3 += s1.elem(2) * s2.elem(1);
  i3 += s1.elem(1) * s2.elem(2);
  i3 += s1.elem(0) * s2.elem(3);

  /* i2 = s1(2)*s2(0) + s1(1)*s2(1) + s1(0)*s2(2) */
  i2  = s1.elem(2) * s2.elem(0);
  i2 += s1.elem(1) * s2.elem(1);
  i2 += s1.elem(0) * s2.elem(2);

  /* i1 = s1(1)*s2(0) + s1(0)*s2(1) */
  i1  = s1.elem(1) * s2.elem(0);
  i1 += s1.elem(0) * s2.elem(1);

  /* i0 = s1(0)*s2(0) */
  i0 = s1.elem(0) * s2.elem(0);

  /* dest(0) = mod(i0, 4096) */
  d.elem(0) = i0 & S(4095);

  /* i1 = i1 + i0/4096 */
  i1 += i0 >> S(12);

  /* dest(1) = mod(i1, 4096) */
  d.elem(1) = i1 & S(4095);

  /* i2 = i2 + i1/4096 */
  i2 += i1 >> S(12);

  /* dest(2) = mod(i2, 4096) */
  d.elem(2) = i2 & S(4095);
  /* i3 = i3 + i2/4096 */
  i3 += i2 >> S(12);

  /* dest(3) = mod(i3, 2048) */
  d.elem(3) = i3 & S(2047);

  return d;
}


template<class T1, class T2>
struct BinaryReturn<PSeed<T1>, PSeed<T2>, OpBitwiseOr> {
  typedef PSeed<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PSeed<T1>, PSeed<T2>, OpBitwiseOr>::Type_t
operator|(const PSeed<T1>& l, const PSeed<T2>& r)
{
  typename BinaryReturn<PSeed<T1>, PSeed<T2>, OpBitwiseOr>::Type_t  d;

  d.elem(0) = l.elem(0) | r.elem(0);
  d.elem(1) = l.elem(1) | r.elem(1);
  d.elem(2) = l.elem(2) | r.elem(2);
  d.elem(3) = l.elem(3) | r.elem(3);

  return d;
}



// Mixed versions
template<class T1, class T2>
struct BinaryReturn<PSeed<T1>, PScalar<T2>, OpBitwiseOr> {
  typedef PSeed<typename BinaryReturn<T1, T2, OpBitwiseOr>::Type_t>  Type_t;
};
 
template<class T1, class T2>
inline typename BinaryReturn<PSeed<T1>, PScalar<T2>, OpBitwiseOr>::Type_t
operator|(const PSeed<T1>& l, const PScalar<T2>& r)
{
  // Lazy implementation

  PSeed<T2>  d;
  d = r;

  return (l | d);
}



/*! 
 * This left shift implementation will not work properly for shifts
 * greater than 12
 */
template<class T1, class T2>
struct BinaryReturn<PSeed<T1>, PScalar<T2>, OpLeftShift> {
  typedef PSeed<typename BinaryReturn<T1, T2, OpLeftShift>::Type_t>  Type_t;
};

template<class T1, class T2>
inline typename BinaryReturn<PSeed<T1>, PScalar<T2>, OpLeftShift>::Type_t
operator<<(const PSeed<T1>& s1, const PScalar<T2>& s2)
{
  typename BinaryReturn<PSeed<T1>, PScalar<T2>, OpLeftShift>::Type_t  d;
  typedef typename BinaryReturn<T1, T2, OpLeftShift>::Type_t  T;
  typedef typename InternalScalar<T>::Type_t  S;
  T  i0, i1, i2, i3;

  i0 = s1.elem(0) << s2.elem();
  i1 = s1.elem(1) << s2.elem();
  i2 = s1.elem(2) << s2.elem();
  i3 = s1.elem(3) << s2.elem();

  d.elem(0) = i0 & S(4095);
  i0 >>= S(12);
  i1 |= i0 & S(4095);
  d.elem(1) = i1 & S(4095);
  i1 >>= S(12);
  i2 |= i1 & S(4095);
  d.elem(2) = i2 & S(4095);
  i2 >>= S(12);
  i3 |= i2 & S(4095);
  d.elem(3) = i3 & S(2047);

  return d;
}


//! dest [float type] = source [seed type]
template<class T>
struct UnaryReturn<PSeed<T>, FnSeedToFloat> {
  typedef PScalar<typename UnaryReturn<T, FnSeedToFloat>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<PSeed<T>, FnSeedToFloat>::Type_t
seedToFloat(const PSeed<T>& s1)
{
  typename UnaryReturn<PSeed<T>, FnSeedToFloat>::Type_t  d;
  typedef typename RealScalar<T>::Type_t  S;

  S  twom11(1.0 / 2048.0);
  S  twom12(1.0 / 4096.0);
  S  fs1, fs2;

//  recast_rep(fs1, s1.elem(0));
  fs1 = S(s1.elem(0));
  d.elem() = twom12 * S(s1.elem(0));

//  recast_rep(fs1, s1.elem(1));
  fs1 = S(s1.elem(1));
  fs2 = fs1 + d.elem();
  d.elem() = twom12 * fs2;

//  recast_rep(fs1, s1.elem(2));
  fs1 = S(s1.elem(2));
  fs2 = fs1 + d.elem();
  d.elem() = twom12 * fs2;

//  recast_rep(fs1, s1.elem(3));
  fs1 = S(s1.elem(3));
  fs2 = fs1 + d.elem();
  d.elem() = twom11 * fs2;

  return d;
}


//! dest [some type] = source [some type]
/*! Portable (internal) way of returning a single site */
template<class T>
struct UnaryReturn<PSeed<T>, FnGetSite> {
  typedef PSeed<typename UnaryReturn<T, FnGetSite>::Type_t>  Type_t;
};

template<class T>
inline typename UnaryReturn<PSeed<T>, FnGetSite>::Type_t
getSite(const PSeed<T>& s1, int innersite)
{ 
  typename UnaryReturn<PSeed<T>, FnGetSite>::Type_t  d;

  for(int i=0; i < 4; ++i)
    d.elem(i) = getSite(s1.elem(i), innersite);

  return d;
}


// Functions
//! dest = 0
template<class T> 
inline void 
zero_rep(PSeed<T>& dest) 
{
  for(int i=0; i < 4; ++i)
    zero_rep(dest.elem(i));
}


//! dest = (mask) ? s1 : dest
template<class T, class T1> 
inline void 
copymask(PSeed<T>& d, const PScalar<T1>& mask, const PSeed<T>& s1) 
{
  for(int i=0; i < 4; ++i)
    copymask(d.elem(i),mask.elem(),s1.elem(i));
}

/*! @} */

} // namespace QDP

#endif
