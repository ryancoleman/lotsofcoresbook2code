// -*- C++ -*-
// $Id: qdp_scalarsite_sse_vector.h,v 1.5 2007-08-20 17:08:14 uid4709 Exp $

/*! @file
 * @brief Vector optimizations
 *
 * Vector extension optimizations of basic operations
 *
 * NOT CURRENTLY USED
 */

#ifndef QDP_SCALARSITE_SSE_VECTOR_H
#define QDP_SCALARSITE_SSE_VECTOR_H

#error "THIS IS UNUSED AT THE MOMENT"

// These SSE asm instructions are only supported under GCC/G++
#if defined(__GNUC__)

#include "qdp_sse_intrin.h"

namespace QDP {

/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */

// Use this def just to safe some typing later on in the file
typedef RComplex<PScalar<float> >  RComplexFloat;
typedef PDWVector<float,4>         PDWVectorFloat4;
typedef RComplex<PDWVectorFloat4>  RComplexFloat4;

//-------------------------------------------------------------------------
// Start of PDWVector optimizations
#if 0


// Use this def just to safe some typing later on in the file
//#define PVectorFloat  PDWVector<float,4>
//#define RComplexFloat  RComplex<ILattice<float,4> >


#if 0
// NOTE: the   operator+(v4sf,v4sf) first exists in gcc 3.3.X, not 3.2.X

// v4sf + v4sf
inline v4sf
operator+(v4sf l, v4sf r)
{
  v4sf tmp = _mm_add_ps(l, r);
  return tmp;
}


// v4sf - v4sf
inline v4sf
operator-(v4sf l, v4sf r)
{
  return _mm_sub_ps(l, r);
}


// v4sf * v4sf
inline v4sf
operator*(v4sf l, v4sf r)
{
  return _mm_mul_ps(l, r);
}


// v4sf / v4sf
inline v4sf
operator/(v4sf l, v4sf r)
{
  return _mm_div_ps(l, r);
}
#endif





#if 1
//! Primitive domain-wall vector class
/*! 
 * Supports domain-wall manipulations
 */
template<> class PDWVector<float,4> : public PVector<float, 4, PDWVector>
{
public:
  typedef float  T;
  static const int N = 4;

  PDWVector() {}
  ~PDWVector() {}

  //---------------------------------------------------------
  //! construct dest = const
  PDWVector(const WordType<float>::Type_t& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs;
    }

  //! construct dest = rhs
  template<class T1>
  PDWVector(const PDWVector<T1,N>& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem(i);
    }

  //! construct dest = rhs
  template<class T1>
  PDWVector(const T1& rhs)
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs;
    }


  //! construct dest = rhs
  PDWVector(const v4sf& rhs)
    {
      F.v = rhs;
    }


  //---------------------------------------------------------
  //! PDWVector = PScalar
  /*! Set equal to an PScalar */
  template<class T1>
  inline
  PDWVector& operator=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) = rhs.elem();

      return *this;
    }

  //! PDWVector += PScalar
  template<class T1>
  inline
  PDWVector& operator+=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) += rhs.elem();

      return *this;
    }

  //! PDWVector -= PScalar
  template<class T1>
  inline
  PDWVector& operator-=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) -= rhs.elem();

      return *this;
    }

  //! PDWVector *= PScalar
  template<class T1>
  inline
  PDWVector& operator*=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) *= rhs.elem();

      return *this;
    }

  //! PDWVector /= PScalar
  template<class T1>
  inline
  PDWVector& operator/=(const PScalar<T1>& rhs) 
    {
      for(int i=0; i < N; ++i)
	elem(i) /= rhs.elem();

      return *this;
    }


  //---------------------------------------------------------
  //! PDWVector = PDWVector
  /*! Set equal to another PDWVector */
  inline
  PDWVector& operator=(const PDWVector& rhs) 
    {
      F.v = rhs.F.v;
      return *this;
    }

  //! PDWVector += PDWVector
  inline
  PDWVector& operator+=(const PDWVector& rhs) 
    {
      F.v = _mm_add_ps(F.v, rhs.F.v);
      return *this;
    }

  //! PDWVector -= PDWVector
  inline
  PDWVector& operator-=(const PDWVector& rhs) 
    {
      F.v = _mm_sub_ps(F.v, rhs.F.v);
      return *this;
    }

  //! PDWVector *= PDWVector
  inline
  PDWVector& operator*=(const PDWVector& rhs) 
    {
      F.v = _mm_mul_ps(F.v, rhs.F.v);
      return *this;
    }

  //! PDWVector /= PDWVector
  inline
  PDWVector& operator/=(const PDWVector& rhs) 
    {
      F.v = _mm_div_ps(F.v, rhs.F.v);
      return *this;
    }


  //! Deep copy constructor
  PDWVector(const PDWVector& a)
    {
      // fprintf(stderr,"copy PDWVector\n");
      F.v = a.F.v;
    }


public:
  //! The backdoor
  /*! 
   * Used by optimization routines (e.g., SSE) that need the memory address of data.
   * BTW: to make this a friend would be a real pain since functions are templatized.
   */
  inline T* data() {return F.a;}


public:
  T& elem(int i) {return F.a[i];}
  const T& elem(int i) const {return F.a[i];}

  v4sf& elem_v() {return F.v;}
  const v4sf elem_v() const {return F.v;}

private:
  // SSE attributes
  union {
    v4sf v;
    T    a[4];
  } F  QDP_ALIGN16;

};
#endif




//--------------------------------------------------------------------------------------
// Optimized version of  
//    PDWVectorFloat4 <- PDWVectorFloat4 + PDWVectorFloat4
//template<>
inline PDWVectorFloat4
operator+(const PDWVectorFloat4& l, const PDWVectorFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "DWV+DWV" << endl;
#endif

  return _mm_add_ps(l.elem_v(), r.elem_v());
}


// Optimized version of  
//    PDWVectorFloat4 <- PDWVectorFloat4 - PDWVectorFloat4
//template<>
inline PDWVectorFloat4
operator-(const PDWVectorFloat4& l, const PDWVectorFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "DWV-DWV" << endl;
#endif

  return _mm_sub_ps(l.elem_v(), r.elem_v());
}


// Optimized version of  
//    PDWVectorFloat4 <- PDWVectorFloat4 * PDWVectorFloat4
//template<>
inline PDWVectorFloat4
operator*(const PDWVectorFloat4& l, const PDWVectorFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "DWV * DWV" << endl;
#endif

  return _mm_mul_ps(l.elem_v(), r.elem_v());
}

// Optimized version of  
//    PDWVectorFloat4 <- PScalar * PDWVectorFloat4
inline PDWVectorFloat4
operator*(const PScalar<float>& l, const PDWVectorFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "P * DWV" << endl;
#endif

   v4sf x = _mm_load_ss((float*)&(l.elem()));
  x = _mm_shuffle_ps(x,x,0);

  return _mm_mul_ps(x, r.elem_v());
}

// Optimized version of  
//    PDWVectorFloat4 <- PDWVectorFloat4 * PScalar
inline PDWVectorFloat4
operator*(const PDWVectorFloat4& l, const PScalar<float>& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "DWV * P" << endl;
#endif

   v4sf x = _mm_load_ss((float*)&(r.elem()));

  x = _mm_shuffle_ps(x,x,0);
  return _mm_mul_ps(l.elem_v(), x);
}



// Optimized version of  
//    PDWVectorFloat4 <- PDWVectorFloat4 / PDWVectorFloat4
//template<>
inline PDWVectorFloat4
operator/(const PDWVectorFloat4& l, const PDWVectorFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "DWV / DWV" << endl;
#endif

  return _mm_mul_ps(l.elem_v(), r.elem_v());
}


//--------------------------------------------------------------------------------------
// Optimized version of  
//    RComplexFloat4 <- RComplexFloat4 + RComplexFloat4
template<>
inline RComplexFloat4
operator+(const RComplexFloat4& l, const RComplexFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "C<DWV> + C<DWV>" << endl;
#endif

  return RComplexFloat4(_mm_add_ps(l.real().elem_v(), r.real().elem_v()),
			_mm_add_ps(l.imag().elem_v(), r.imag().elem_v()));
}


// Optimized version of  
//    RComplexFloat4 <- RComplexFloat4 - RComplexFloat4
template<>
inline RComplexFloat4
operator-(const RComplexFloat4& l, const RComplexFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "C<DWV> - C<DWV" << endl;
#endif

  return RComplexFloat4(_mm_sub_ps(l.real().elem_v(), r.real().elem_v()),
			_mm_sub_ps(l.imag().elem_v(), r.imag().elem_v()));
}


// Optimized version of  
//    RComplexFloat4 <- RComplexFloat4 * RComplexFloat4
template<>
inline RComplexFloat4
operator*(const RComplexFloat4& l, const RComplexFloat4& r)
{
  RComplexFloat4 d;

#if defined(QDP_SCALARSITE_DEBUG)
  cout << "C<DWV> * C<DWV>" << endl;
#endif

  v4sf tmp1 = _mm_mul_ps(l.real().elem_v(), r.real().elem_v());
  v4sf tmp2 = _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v());
  d.real().elem_v() = _mm_sub_ps(tmp1, tmp2);

  v4sf tmp3 = _mm_mul_ps(l.real().elem_v(), r.imag().elem_v());
  v4sf tmp4 = _mm_mul_ps(l.imag().elem_v(), r.real().elem_v());
  d.imag().elem_v() = _mm_add_ps(tmp3, tmp4);

  return d;
}

// Optimized version of  
//    RComplexFloat4 <- adj(RComplexFloat4) * RComplexFloat4
template<>
inline BinaryReturn<RComplexFloat4, RComplexFloat4, OpAdjMultiply>::Type_t
adjMultiply(const RComplexFloat4& l, const RComplexFloat4& r)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "adj(C<DWV>) * C<DWV>" << endl;
#endif

  BinaryReturn<RComplexFloat4, RComplexFloat4, OpAdjMultiply>::Type_t  d;

  v4sf tmp1 = _mm_mul_ps(l.real().elem_v(), r.real().elem_v());
  v4sf tmp2 = _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v());
  d.real().elem_v() = _mm_add_ps(tmp1, tmp2);

  v4sf tmp3 = _mm_mul_ps(l.real().elem_v(), r.imag().elem_v());
  v4sf tmp4 = _mm_mul_ps(l.imag().elem_v(), r.real().elem_v());
  d.imag().elem_v() = _mm_sub_ps(tmp3, tmp4);

  return d;
}

// Optimized  RComplex*adj(RComplex)
template<>
inline BinaryReturn<RComplexFloat4, RComplexFloat4, OpMultiplyAdj>::Type_t
multiplyAdj(const RComplexFloat4& l, const RComplexFloat4& r)
{
  BinaryReturn<RComplexFloat4, RComplexFloat4, OpMultiplyAdj>::Type_t  d;

  v4sf tmp1 = _mm_mul_ps(l.real().elem_v(), r.real().elem_v());
  v4sf tmp2 = _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v());
  d.real().elem_v() = _mm_add_ps(tmp1, tmp2);

  v4sf tmp3 = _mm_mul_ps(l.imag().elem_v(), r.real().elem_v());
  v4sf tmp4 = _mm_mul_ps(l.real().elem_v(), r.imag().elem_v());
  d.imag().elem_v() = _mm_sub_ps(tmp3, tmp4);

  return d;
}

// Optimized  adj(RComplex)*adj(RComplex)
template<>
inline BinaryReturn<RComplexFloat4, RComplexFloat4, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const RComplexFloat4& l, const RComplexFloat4& r)
{
  BinaryReturn<RComplexFloat4, RComplexFloat4, OpAdjMultiplyAdj>::Type_t  d;

  typedef struct
  {
    unsigned int c[4];
  } sse_mask __attribute__ ((aligned (16)));
  
  static sse_mask _sse_sgn __attribute__ ((unused)) ={0x80000000, 0x80000000, 0x80000000, 0x80000000};

  v4sf tmp1 = _mm_mul_ps(l.real().elem_v(), r.real().elem_v());
  v4sf tmp2 = _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v());
  d.real().elem_v() = _mm_sub_ps(tmp1, tmp2);

  v4sf tmp3 = _mm_mul_ps(l.real().elem_v(), r.imag().elem_v());
  v4sf tmp4 = _mm_mul_ps(l.imag().elem_v(), r.real().elem_v());
  v4sf tmp5 = _mm_add_ps(tmp3, tmp4);
//  d.imag().elem_v() = _mm_xor_ps(tmp5, _sse_sgn.v);
  v4sf tmp6 = _mm_load_ps((float*)&_sse_sgn);
  d.imag().elem_v() = _mm_xor_ps(tmp5, tmp6);

  return d;
}

#endif


/*! @} */   // end of group optimizations

#if defined(QDP_SCALARSITE_DEBUG)
#undef QDP_SCALARSITE_DEBUG
#endif

} // namespace QDP;

#endif  // defined(__GNUC__)

#endif
