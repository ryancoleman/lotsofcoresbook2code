// -*- C++ -*-

/*! @file
 * @brief Intel SSE optimizations
 *
 * SSE optimizations of basic operations
 */

#ifndef QDP_SCALARVECSITE_SSE_H
#define QDP_SCALARVECSITE_SSE_H


// These SSE asm instructions are only supported under GCC/G++
// Only supported on gcc >= 3.2
#if defined(__GNUC__) && __GNUC_MINOR__ >= 2

namespace QDP {

typedef REAL32 vReal32 __attribute__ ((aligned (16),mode(V4SF)));


static inline vReal32 vmk1(REAL32 a) 
{
  vReal32 v = _mm_load_ss((float *)&a);
  // asm("shufps\t$0,%0,%0" : "+x" (v));
  v = _mm_shuffle_ps(v,v,0);
  return v;
}

static inline vReal32 vmk4(REAL32 a0, REAL32 a1, REAL32 a2, REAL32 a3) 
{
  vReal32 v;
  REAL32 *r = (REAL32 *)&v;
  r[0] = a0;
  r[1] = a1;
  r[2] = a2;
  r[3] = a3;
  return v;
}


//! Specialized Inner lattice class
/*! Uses sse  */
template<> class ILattice<REAL32,4>
{
public:
  typedef REAL32  T;
  static const int N = 4;

  ILattice() {}
  ~ILattice() {}

  //---------------------------------------------------------
  //! construct dest = 4 consts
  ILattice(REAL32 a0, REAL32 a1, REAL32 a2, REAL32 a3) : v(vmk4(a0,a1,a2,a3)) {}
  
  //! construct dest = const
  ILattice(REAL32 rhs) : v(vmk1(rhs)) {}

  //! construct dest = rhs
  template<class T1>
  ILattice(const ILattice<T1,N>& rhs)
    {
      vput_0(rhs.elem(0)); 
      vput_1(rhs.elem(1)); 
      vput_2(rhs.elem(2));
      vput_3(rhs.elem(3));
    }

  //! construct dest = rhs
  ILattice(vReal32 rhs) : v(rhs) {}

  //! construct dest = rhs
  template<class T1>
  ILattice(const IScalar<T1>& rhs) : v(vmk1(REAL32(rhs.elem()))) {}

  //! construct dest = rhs
  template<class T1>
  ILattice(const T1& rhs) : v(vmk1(REAL32(rhs))) {}


  //---------------------------------------------------------
  //! ILattice = IScalar
  /*! Set equal to an IScalar */
  template<class T1>
  inline
  ILattice& operator=(const IScalar<T1>& rhs) 
    {
      v = vmk1(REAL32(rhs.elem()));
      return *this;
    }

  //! ILattice += IScalar
  template<class T1>
  inline
  ILattice& operator+=(const IScalar<T1>& rhs) 
    {
      v = _mm_add_ps(v, vmk1(REAL32(rhs.elem())));
      return *this;
    }

  //! ILattice -= IScalar
  template<class T1>
  inline
  ILattice& operator-=(const IScalar<T1>& rhs) 
    {
      v = _mm_sub_ps(v, vmk1(REAL32(rhs.elem())));
      return *this;
    }

  //! ILattice *= IScalar
  template<class T1>
  inline
  ILattice& operator*=(const IScalar<T1>& rhs) 
    {
      v = _mm_mul_ps(v, vmk1(REAL32(rhs.elem())));
      return *this;
    }

  //! ILattice /= IScalar
  template<class T1>
  inline
  ILattice& operator/=(const IScalar<T1>& rhs) 
    {
      v = _mm_div_ps(v, vmk1(REAL32(rhs.elem())));
      return *this;
    }


  //---------------------------------------------------------
  //! ILattice = ILattice
  /*! Set equal to another ILattice */
  inline
  ILattice& operator=(const ILattice& rhs) 
    {
      v = rhs.v;
      return *this;
    }

  //! ILattice += ILattice
  inline
  ILattice& operator+=(const ILattice& rhs) 
    {
      v = _mm_add_ps(v, rhs.v);
      return *this;
    }

  //! ILattice -= ILattice
  inline
  ILattice& operator-=(const ILattice& rhs) 
    {
      v = _mm_sub_ps(v, rhs.v);
      return *this;
    }

  //! ILattice *= ILattice
  inline
  ILattice& operator*=(const ILattice& rhs) 
    {
      v = _mm_mul_ps(v, rhs.v);
      return *this;
    }

  //! ILattice /= ILattice
  inline
  ILattice& operator/=(const ILattice& rhs) 
    {
      v = _mm_div_ps(v, rhs.v);
      return *this;
    }


  //! Deep copy constructor
  ILattice(const ILattice& a) : v(a.v) {}


public:
  //! The backdoor
  /*! 
   * Used by optimization routines (e.g., SSE) that need the memory address of data.
   * BTW: to make this a friend would be a real pain since functions are templatized.
   */
  inline T* data() {return (REAL32*)&v;}


public:
  void vput_0(REAL32 a) { ((REAL32 *)&v)[0] = a; }
  void vput_1(REAL32 a) { ((REAL32 *)&v)[1] = a; }
  void vput_2(REAL32 a) { ((REAL32 *)&v)[2] = a; }
  void vput_3(REAL32 a) { ((REAL32 *)&v)[3] = a; }  
  
  T& elem(int i) {return ((REAL32*)&v)[i];}
  const T& elem(int i) const {return ((REAL32*)&v)[i];}

  vReal32& elem_v() {return v;}
  const vReal32 elem_v() const {return v;}

private:
  // SSE attributes
  vReal32 v;

} QDP_ALIGN16;


} // namespace QDP

// Use SSE specific Linalg stuff (inline assembler etc)
#include "scalarvecsite_sse/qdp_scalarvecsite_sse_linalg.h"

// Use SSE specific blas stuff (inline assembler etc)
#include "scalarvecsite_sse/qdp_scalarvecsite_sse_blas.h"

#else
#error "This is not a GNUC 3.2 or greater compiler, and therefore does not support the GNU specific asm directives."
#endif  // gnuc

#endif  // guard

