// -*- C++ -*-
// $Id: qdp_scalarvecsite_sse_linalg.h,v 1.8 2007-06-10 14:32:11 edwards Exp $

/*! @file
 * @brief Intel SSE optimizations
 *
 * SSE optimizations of basic operations
 */

#ifndef QDP_SCALARVECSITE_SSE_LINALG_H
#define QDP_SCALARVECSITE_SSE_LINALG_H

// These SSE asm instructions are only supported under GCC/G++ 3.2 or greater
#if defined(__GNUC__) && __GNUC_MINOR__ >= 2

#include "qdp_sse_intrin.h"

namespace QDP {

// #define QDP_SCALARVECSITE_DEBUG

#define QDP_SCALARVECSITE_USE_EVALUATE



/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */

// Use this def just to safe some typing later on in the file
typedef IScalar<REAL32>                IScalarFloat;
typedef ILattice<REAL32,4>             ILatticeFloat;
typedef RComplex<ILattice<REAL32,4> >  RComplexFloat; 

#include "scalarvecsite_sse/ssevec_mult_nn.h"

//--------------------------------------------------------------------------------------
// Optimized version of  
// ILatticeFloat <- ILatticeFloat + ILatticeFloat
template<>
inline BinaryReturn<ILatticeFloat, ILatticeFloat, OpAdd>::Type_t
operator+(const ILatticeFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, ILatticeFloat, OpAdd>::Type_t  Ret_t;

//  cout << "I+I" << endl; 
  return Ret_t(_mm_add_ps(l.elem_v(), r.elem_v()));
}

// ILatticeFloat <- ILatticeFloat + IScalarFloat
template<>
inline BinaryReturn<ILatticeFloat, IScalarFloat, OpAdd>::Type_t
operator+(const ILatticeFloat& l, const IScalarFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, IScalarFloat, OpAdd>::Type_t  Ret_t;

//  cout << "I+I" << endl; 
  return Ret_t(_mm_add_ps(l.elem_v(), vmk1(r.elem())));
}

// ILatticeFloat <- IScalarFloat + ILatticeFloat
template<>
inline BinaryReturn<IScalarFloat, ILatticeFloat, OpAdd>::Type_t
operator+(const IScalarFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<IScalarFloat, ILatticeFloat, OpAdd>::Type_t  Ret_t;

//  cout << "I+I" << endl; 
  return Ret_t(_mm_add_ps(vmk1(l.elem()), r.elem_v()));
}


// Optimized version of  
// ILatticeFloat <- ILatticeFloat - ILatticeFloat
template<>
inline BinaryReturn<ILatticeFloat, ILatticeFloat, OpSubtract>::Type_t
operator-(const ILatticeFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, ILatticeFloat, OpSubtract>::Type_t  Ret_t;

//  cout << "I-I" << endl;
  return Ret_t(_mm_sub_ps(l.elem_v(), r.elem_v()));
}

// ILatticeFloat <- ILatticeFloat - IScalarFloat
template<>
inline BinaryReturn<ILatticeFloat, IScalarFloat, OpSubtract>::Type_t
operator-(const ILatticeFloat& l, const IScalarFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, IScalarFloat, OpSubtract>::Type_t  Ret_t;

//  cout << "I-I" << endl;
  return Ret_t(_mm_sub_ps(l.elem_v(), vmk1(r.elem())));
}

// ILatticeFloat <- IScalarFloat - ILatticeFloat
template<>
inline BinaryReturn<IScalarFloat, ILatticeFloat, OpSubtract>::Type_t
operator-(const IScalarFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<IScalarFloat, ILatticeFloat, OpSubtract>::Type_t  Ret_t;

//  cout << "I-I" << endl;
  return Ret_t(_mm_sub_ps(vmk1(l.elem()), r.elem_v()));
}


// Optimized version of  
//    ILatticeFloat <- ILatticeFloat * ILatticeFloat
template<>
inline BinaryReturn<ILatticeFloat, ILatticeFloat, OpMultiply>::Type_t
operator*(const ILatticeFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, ILatticeFloat, OpMultiply>::Type_t  Ret_t;

//  cout << "I*I" << endl;
  return Ret_t(_mm_mul_ps(l.elem_v(), r.elem_v()));
}

// ILatticeFloat <- ILatticeFloat * IScalarFloat
template<>
inline BinaryReturn<ILatticeFloat, IScalarFloat, OpMultiply>::Type_t
operator*(const ILatticeFloat& l, const IScalarFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, IScalarFloat, OpMultiply>::Type_t  Ret_t;

//  cout << "I*I" << endl;
  return Ret_t(_mm_mul_ps(l.elem_v(), vmk1(r.elem())));
}

// ILatticeFloat <- IScalarFloat * ILatticeFloat
template<>
inline BinaryReturn<IScalarFloat, ILatticeFloat, OpMultiply>::Type_t
operator*(const IScalarFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<IScalarFloat, ILatticeFloat, OpMultiply>::Type_t  Ret_t;

//  cout << "I*I" << endl;
  return Ret_t(_mm_mul_ps(vmk1(l.elem()), r.elem_v()));
}


// Optimized version of  
// ILatticeFloat <- ILatticeFloat / ILatticeFloat
template<>
inline BinaryReturn<ILatticeFloat, ILatticeFloat, OpDivide>::Type_t
operator/(const ILatticeFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, ILatticeFloat, OpDivide>::Type_t  Ret_t;

//  cout << "I/I" << endl;
  return Ret_t(_mm_div_ps(l.elem_v(), r.elem_v()));
}

// ILatticeFloat <- ILatticeFloat / IScalarFloat
template<>
inline BinaryReturn<ILatticeFloat, IScalarFloat, OpDivide>::Type_t
operator/(const ILatticeFloat& l, const IScalarFloat& r)
{
  typedef BinaryReturn<ILatticeFloat, IScalarFloat, OpDivide>::Type_t  Ret_t;

//  cout << "I/I" << endl;
  return Ret_t(_mm_div_ps(l.elem_v(), vmk1(r.elem())));
}

// ILatticeFloat <- IScalarFloat / ILatticeFloat
template<>
inline BinaryReturn<IScalarFloat, ILatticeFloat, OpDivide>::Type_t
operator/(const IScalarFloat& l, const ILatticeFloat& r)
{
  typedef BinaryReturn<IScalarFloat, ILatticeFloat, OpDivide>::Type_t  Ret_t;

//  cout << "I/I" << endl;
  return Ret_t(_mm_div_ps(vmk1(l.elem()), r.elem_v()));
}



//--------------------------------------------------------------------------------------
// Optimized version of  
//   RComplexFloat <- RComplexFloat + RComplexFloat
template<>
inline BinaryReturn<RComplexFloat, RComplexFloat, OpAdd>::Type_t
operator+(const RComplexFloat& l, const RComplexFloat& r)
{
  typedef BinaryReturn<RComplexFloat, RComplexFloat, OpAdd>::Type_t  Ret_t;

//  cout << "C+C" << endl;
  return Ret_t(_mm_add_ps(l.real().elem_v(), r.real().elem_v()),
	       _mm_add_ps(l.imag().elem_v(), r.imag().elem_v()));
}


// Optimized version of  
//    RComplexFloat <- RComplexFloat - RComplexFloat
template<>
inline BinaryReturn<RComplexFloat, RComplexFloat, OpSubtract>::Type_t
operator-(const RComplexFloat& l, const RComplexFloat& r)
{
  typedef BinaryReturn<RComplexFloat, RComplexFloat, OpSubtract>::Type_t  Ret_t;

//  cout << "C-C" << endl;
  return Ret_t(_mm_sub_ps(l.real().elem_v(), r.real().elem_v()),
	       _mm_sub_ps(l.imag().elem_v(), r.imag().elem_v()));
}


// Optimized version of  
//    RComplexFloat <- RComplexFloat * RComplexFloat
template<>
inline BinaryReturn<RComplexFloat, RComplexFloat, OpMultiply>::Type_t
operator*(const RComplexFloat& l, const RComplexFloat& r)
{
  typedef BinaryReturn<RComplexFloat, RComplexFloat, OpMultiply>::Type_t  Ret_t;

//  cout << "C*C" << endl;
  return Ret_t(_mm_sub_ps(_mm_mul_ps(l.real().elem_v(), r.real().elem_v()),
				    _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v())),
	       _mm_add_ps(_mm_mul_ps(l.real().elem_v(), r.imag().elem_v()),
				    _mm_mul_ps(l.imag().elem_v(), r.real().elem_v())));
}

// Optimized version of  
//    RComplexFloat <- adj(RComplexFloat) * RComplexFloat
template<>
inline BinaryReturn<RComplexFloat, RComplexFloat, OpAdjMultiply>::Type_t
adjMultiply(const RComplexFloat& l, const RComplexFloat& r)
{
  typedef BinaryReturn<RComplexFloat, RComplexFloat, OpAdjMultiply>::Type_t  Ret_t;

//  cout << "adj(C)*C" << endl;
  return Ret_t(_mm_add_ps(_mm_mul_ps(l.real().elem_v(), r.real().elem_v()),
				    _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v())),
	       _mm_sub_ps(_mm_mul_ps(l.real().elem_v(), r.imag().elem_v()),
				    _mm_mul_ps(l.imag().elem_v(), r.real().elem_v())));
}

// Optimized  RComplex*adj(RComplex)
template<>
inline BinaryReturn<RComplexFloat, RComplexFloat, OpMultiplyAdj>::Type_t
multiplyAdj(const RComplexFloat& l, const RComplexFloat& r)
{
  typedef BinaryReturn<RComplexFloat, RComplexFloat, OpMultiplyAdj>::Type_t  Ret_t;

//  cout << "C*adj(C)" << endl;
  return Ret_t(_mm_add_ps(_mm_mul_ps(l.real().elem_v(), r.real().elem_v()),
				    _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v())),
	       _mm_sub_ps(_mm_mul_ps(l.imag().elem_v(), r.real().elem_v()),
				    _mm_mul_ps(l.real().elem_v(), r.imag().elem_v())));
}

// Optimized  adj(RComplex)*adj(RComplex)
template<>
inline BinaryReturn<RComplexFloat, RComplexFloat, OpAdjMultiplyAdj>::Type_t
adjMultiplyAdj(const RComplexFloat& l, const RComplexFloat& r)
{
  typedef BinaryReturn<RComplexFloat, RComplexFloat, OpAdjMultiplyAdj>::Type_t  Ret_t;
  REAL32 zero = 0.0;

//  cout << "adj(C)*adj(C)" << endl;
  return Ret_t(_mm_sub_ps(_mm_mul_ps(l.real().elem_v(), r.real().elem_v()),
				    _mm_mul_ps(l.imag().elem_v(), r.imag().elem_v())),
	       _mm_sub_ps(vmk1(zero),
				    _mm_add_ps(_mm_mul_ps(l.real().elem_v(), r.imag().elem_v()),
							 _mm_mul_ps(l.imag().elem_v(), r.real().elem_v()))));
}


//--------------------------------------------------------------------------------------


// Optimized version of  
//    PColorMatrix<RComplexFloat,3> <- PColorMatrix<RComplexFloat,3> * PColorMatrix<RComplexFloat,3>
template<>
inline BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
  PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiply>::Type_t
operator*(const PMatrix<RComplexFloat,3,PColorMatrix>& l, 
	  const PMatrix<RComplexFloat,3,PColorMatrix>& r)
{
//  cout << "M*M" << endl;

  BinaryReturn<PMatrix<RComplexFloat,3,PColorMatrix>, 
    PMatrix<RComplexFloat,3,PColorMatrix>, OpMultiply>::Type_t  d;

  REAL32 *dd = (REAL32*)&d;
  REAL32 *ll = (REAL32*)&l;
  REAL32 *rr = (REAL32*)&r;

  _inline_ssevec_mult_su3_nn(dd,ll,rr,0);
  _inline_ssevec_mult_su3_nn(dd,ll,rr,1);
  _inline_ssevec_mult_su3_nn(dd,ll,rr,2);

  return d;
}



#if defined(QDP_SCALARVECSITE_USE_EVALUATE)
// Specialization to optimize the case   
//    LatticeColorMatrix = LatticeColorMatrix * LatticeColorMatrix
// NOTE: let this be a subroutine to save space
template<>
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s);
#endif



/*! @} */   // end of group optimizations

#if defined(QDP_SCALARVECSITE_DEBUG)
#undef QDP_SCALARVECSITE_DEBUG
#endif

#if defined(QDP_SCALARVECSITE_USE_EVALUATE)
#undef QDP_SCALARVECSITE_USE_EVALUATE
#endif


} // namespace QDP;

#endif  // defined(__GNUC__)

#endif
