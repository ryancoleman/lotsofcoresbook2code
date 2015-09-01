// -*- C++ -*-
//
//
// QDP data parallel interface
//

#ifndef QDP_PRIMGAMMA_H
#define QDP_PRIMGAMMA_H

namespace QDP {

//-------------------------------------------------------------------------------------
//! Gamma matrices
//! Simple interface for gamma matrices. These are run-time constructable
template<int N> class GammaType
{
public:
  //! Main constructor 
  GammaType() {}
  //! Destructor
  ~GammaType() {}

  //! Index  from an integer
  GammaType(int mm) : m(mm) {}

public:
  //! The integer representation for which product of gamma matrices
  int elem() const {return m;}  

private:
  //! Representation
  /*! 
   * The integer in the range 0 to Ns*Ns-1 that represents which product
   * gamma matrices
   */
  const int m;
};


//-------------------------------------------------------------------------------------
//! Gamma matrices
//! Simple interface for gamma matrices. These are run-time constructable
template<int N> class GammaTypeDP
{
public:
  //! Main constructor 
  GammaTypeDP() {}
  //! Destructor
  ~GammaTypeDP() {}

  //! Index  from an integer
  GammaTypeDP(int mm) : m(mm) {}

public:
  //! The integer representation for which product of gamma matrices
  int elem() const {return m;}  

private:
  //! Representation
  /*! 
   * The integer in the range 0 to Ns*Ns-1 that represents which product
   * gamma matrices
   */
  const int m;
};


//-------------------------------------------------------------------------------------
//! Gamma matrices
//! Simple interface for gamma matrices. These are compile-time fixed
template<int N, int m> class GammaConst
{
};


//-------------------------------------------------------------------------------------
//! Gamma matrices
//! Simple interface for gamma matrices. These are compile-time fixed
template<int N, int m> class GammaConstDP
{
};


//-----------------------------------------------------------------------------
// Traits classes to support return types
//-----------------------------------------------------------------------------

// PScalar
template<int N, int m, class T2, class OpGammaConstMultiply>
struct BinaryReturn<GammaConst<N,m>, PScalar<T2>, OpGammaConstMultiply> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConst>
struct BinaryReturn<PScalar<T2>, GammaConst<N,m>, OpMultiplyGammaConst> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeMultiply>
struct BinaryReturn<GammaType<N>, PScalar<T2>, OpGammaTypeMultiply> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaType>
struct BinaryReturn<PScalar<T2>, GammaType<N>, OpMultiplyGammaType> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

 
// PScalar
template<int N, int m, class T2, class OpGammaConstDPMultiply>
struct BinaryReturn<GammaConstDP<N,m>, PScalar<T2>, OpGammaConstDPMultiply> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, int m, class OpMultiplyGammaConstDP>
struct BinaryReturn<PScalar<T2>, GammaConstDP<N,m>, OpMultiplyGammaConstDP> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpGammaTypeDPMultiply>
struct BinaryReturn<GammaTypeDP<N>, PScalar<T2>, OpGammaTypeDPMultiply> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

template<class T2, int N, class OpMultiplyGammaTypeDP>
struct BinaryReturn<PScalar<T2>, GammaTypeDP<N>, OpMultiplyGammaTypeDP> {
  typedef PScalar<typename UnaryReturn<T2, OpUnaryPlus>::Type_t>  Type_t;
};

 
//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

// Wrappers for standard binary multiplication with gamma matrix constants
template<class T, int N> inline T Multiply_gamma0_arg(const T& s1) {return (GammaConst<N,0>() * s1);}
template<class T, int N> inline T Multiply_gamma1_arg(const T& s1) {return (GammaConst<N,1>() * s1);}
template<class T, int N> inline T Multiply_gamma2_arg(const T& s1) {return (GammaConst<N,2>() * s1);}
template<class T, int N> inline T Multiply_gamma3_arg(const T& s1) {return (GammaConst<N,3>() * s1);}
template<class T, int N> inline T Multiply_gamma4_arg(const T& s1) {return (GammaConst<N,4>() * s1);}
template<class T, int N> inline T Multiply_gamma5_arg(const T& s1) {return (GammaConst<N,5>() * s1);}
template<class T, int N> inline T Multiply_gamma6_arg(const T& s1) {return (GammaConst<N,6>() * s1);}
template<class T, int N> inline T Multiply_gamma7_arg(const T& s1) {return (GammaConst<N,7>() * s1);}
template<class T, int N> inline T Multiply_gamma8_arg(const T& s1) {return (GammaConst<N,8>() * s1);}
template<class T, int N> inline T Multiply_gamma9_arg(const T& s1) {return (GammaConst<N,9>() * s1);}
template<class T, int N> inline T Multiply_gamma10_arg(const T& s1) {return (GammaConst<N,10>() * s1);}
template<class T, int N> inline T Multiply_gamma11_arg(const T& s1) {return (GammaConst<N,11>() * s1);}
template<class T, int N> inline T Multiply_gamma12_arg(const T& s1) {return (GammaConst<N,12>() * s1);}
template<class T, int N> inline T Multiply_gamma13_arg(const T& s1) {return (GammaConst<N,13>() * s1);}
template<class T, int N> inline T Multiply_gamma14_arg(const T& s1) {return (GammaConst<N,14>() * s1);}
template<class T, int N> inline T Multiply_gamma15_arg(const T& s1) {return (GammaConst<N,15>() * s1);}

// Wrappers for standard binary multiplication with gamma matrix constants
template<class T, int N> inline T Multiply_arg_gamma0(const T& s1) {return (s1 * GammaConst<N,0>());}
template<class T, int N> inline T Multiply_arg_gamma1(const T& s1) {return (s1 * GammaConst<N,1>());}
template<class T, int N> inline T Multiply_arg_gamma2(const T& s1) {return (s1 * GammaConst<N,2>());}
template<class T, int N> inline T Multiply_arg_gamma3(const T& s1) {return (s1 * GammaConst<N,3>());}
template<class T, int N> inline T Multiply_arg_gamma4(const T& s1) {return (s1 * GammaConst<N,4>());}
template<class T, int N> inline T Multiply_arg_gamma5(const T& s1) {return (s1 * GammaConst<N,5>());}
template<class T, int N> inline T Multiply_arg_gamma6(const T& s1) {return (s1 * GammaConst<N,6>());}
template<class T, int N> inline T Multiply_arg_gamma7(const T& s1) {return (s1 * GammaConst<N,7>());}
template<class T, int N> inline T Multiply_arg_gamma8(const T& s1) {return (s1 * GammaConst<N,8>());}
template<class T, int N> inline T Multiply_arg_gamma9(const T& s1) {return (s1 * GammaConst<N,9>());}
template<class T, int N> inline T Multiply_arg_gamma10(const T& s1) {return (s1 * GammaConst<N,10>());}
template<class T, int N> inline T Multiply_arg_gamma11(const T& s1) {return (s1 * GammaConst<N,11>());}
template<class T, int N> inline T Multiply_arg_gamma12(const T& s1) {return (s1 * GammaConst<N,12>());}
template<class T, int N> inline T Multiply_arg_gamma13(const T& s1) {return (s1 * GammaConst<N,13>());}
template<class T, int N> inline T Multiply_arg_gamma14(const T& s1) {return (s1 * GammaConst<N,14>());}
template<class T, int N> inline T Multiply_arg_gamma15(const T& s1) {return (s1 * GammaConst<N,15>());}



template<int N, class T>
inline T OpGammaTypeMultiply:: operator()(const GammaType<N>& a, const T &b) const
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_gamma0_arg<T,N>, 
		 &Multiply_gamma1_arg<T,N>,
		 &Multiply_gamma2_arg<T,N>,
		 &Multiply_gamma3_arg<T,N>,
		 &Multiply_gamma4_arg<T,N>,
		 &Multiply_gamma5_arg<T,N>,
		 &Multiply_gamma6_arg<T,N>,
		 &Multiply_gamma7_arg<T,N>,
		 &Multiply_gamma8_arg<T,N>,
		 &Multiply_gamma9_arg<T,N>,
		 &Multiply_gamma10_arg<T,N>,
		 &Multiply_gamma11_arg<T,N>,
		 &Multiply_gamma12_arg<T,N>,
		 &Multiply_gamma13_arg<T,N>,
		 &Multiply_gamma14_arg<T,N>,
		 &Multiply_gamma15_arg<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[a.elem()](b);
}


template<class T, int N>
inline T OpMultiplyGammaType::operator()(const T &a, const GammaType<N>& b) const
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_arg_gamma0<T,N>, 
		 &Multiply_arg_gamma1<T,N>,
		 &Multiply_arg_gamma2<T,N>,
		 &Multiply_arg_gamma3<T,N>,
		 &Multiply_arg_gamma4<T,N>,
		 &Multiply_arg_gamma5<T,N>,
		 &Multiply_arg_gamma6<T,N>,
		 &Multiply_arg_gamma7<T,N>,
		 &Multiply_arg_gamma8<T,N>,
		 &Multiply_arg_gamma9<T,N>,
		 &Multiply_arg_gamma10<T,N>,
		 &Multiply_arg_gamma11<T,N>,
		 &Multiply_arg_gamma12<T,N>,
		 &Multiply_arg_gamma13<T,N>,
		 &Multiply_arg_gamma14<T,N>,
		 &Multiply_arg_gamma15<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[b.elem()](a);
}

//-----------------------------------------------------------------------------
// Operators
//-----------------------------------------------------------------------------

// Wrappers for standard binary multiplication with gamma matrix constants
template<class T, int N> inline T Multiply_gamma0_DP_arg(const T& s1) {return (GammaConstDP<N,0>() * s1);}
template<class T, int N> inline T Multiply_gamma1_DP_arg(const T& s1) {return (GammaConstDP<N,1>() * s1);}
template<class T, int N> inline T Multiply_gamma2_DP_arg(const T& s1) {return (GammaConstDP<N,2>() * s1);}
template<class T, int N> inline T Multiply_gamma3_DP_arg(const T& s1) {return (GammaConstDP<N,3>() * s1);}
template<class T, int N> inline T Multiply_gamma4_DP_arg(const T& s1) {return (GammaConstDP<N,4>() * s1);}
template<class T, int N> inline T Multiply_gamma5_DP_arg(const T& s1) {return (GammaConstDP<N,5>() * s1);}
template<class T, int N> inline T Multiply_gamma6_DP_arg(const T& s1) {return (GammaConstDP<N,6>() * s1);}
template<class T, int N> inline T Multiply_gamma7_DP_arg(const T& s1) {return (GammaConstDP<N,7>() * s1);}
template<class T, int N> inline T Multiply_gamma8_DP_arg(const T& s1) {return (GammaConstDP<N,8>() * s1);}
template<class T, int N> inline T Multiply_gamma9_DP_arg(const T& s1) {return (GammaConstDP<N,9>() * s1);}
template<class T, int N> inline T Multiply_gamma10_DP_arg(const T& s1) {return (GammaConstDP<N,10>() * s1);}
template<class T, int N> inline T Multiply_gamma11_DP_arg(const T& s1) {return (GammaConstDP<N,11>() * s1);}
template<class T, int N> inline T Multiply_gamma12_DP_arg(const T& s1) {return (GammaConstDP<N,12>() * s1);}
template<class T, int N> inline T Multiply_gamma13_DP_arg(const T& s1) {return (GammaConstDP<N,13>() * s1);}
template<class T, int N> inline T Multiply_gamma14_DP_arg(const T& s1) {return (GammaConstDP<N,14>() * s1);}
template<class T, int N> inline T Multiply_gamma15_DP_arg(const T& s1) {return (GammaConstDP<N,15>() * s1);}

// Wrappers for standard binary multiplication with gamma matrix constants
template<class T, int N> inline T Multiply_DP_arg_gamma0(const T& s1) {return (s1 * GammaConstDP<N,0>());}
template<class T, int N> inline T Multiply_DP_arg_gamma1(const T& s1) {return (s1 * GammaConstDP<N,1>());}
template<class T, int N> inline T Multiply_DP_arg_gamma2(const T& s1) {return (s1 * GammaConstDP<N,2>());}
template<class T, int N> inline T Multiply_DP_arg_gamma3(const T& s1) {return (s1 * GammaConstDP<N,3>());}
template<class T, int N> inline T Multiply_DP_arg_gamma4(const T& s1) {return (s1 * GammaConstDP<N,4>());}
template<class T, int N> inline T Multiply_DP_arg_gamma5(const T& s1) {return (s1 * GammaConstDP<N,5>());}
template<class T, int N> inline T Multiply_DP_arg_gamma6(const T& s1) {return (s1 * GammaConstDP<N,6>());}
template<class T, int N> inline T Multiply_DP_arg_gamma7(const T& s1) {return (s1 * GammaConstDP<N,7>());}
template<class T, int N> inline T Multiply_DP_arg_gamma8(const T& s1) {return (s1 * GammaConstDP<N,8>());}
template<class T, int N> inline T Multiply_DP_arg_gamma9(const T& s1) {return (s1 * GammaConstDP<N,9>());}
template<class T, int N> inline T Multiply_DP_arg_gamma10(const T& s1) {return (s1 * GammaConstDP<N,10>());}
template<class T, int N> inline T Multiply_DP_arg_gamma11(const T& s1) {return (s1 * GammaConstDP<N,11>());}
template<class T, int N> inline T Multiply_DP_arg_gamma12(const T& s1) {return (s1 * GammaConstDP<N,12>());}
template<class T, int N> inline T Multiply_DP_arg_gamma13(const T& s1) {return (s1 * GammaConstDP<N,13>());}
template<class T, int N> inline T Multiply_DP_arg_gamma14(const T& s1) {return (s1 * GammaConstDP<N,14>());}
template<class T, int N> inline T Multiply_DP_arg_gamma15(const T& s1) {return (s1 * GammaConstDP<N,15>());}



template<int N, class T>
inline T OpGammaTypeDPMultiply:: operator()(const GammaTypeDP<N>& a, const T &b) const
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_gamma0_DP_arg<T,N>, 
		 &Multiply_gamma1_DP_arg<T,N>,
		 &Multiply_gamma2_DP_arg<T,N>,
		 &Multiply_gamma3_DP_arg<T,N>,
		 &Multiply_gamma4_DP_arg<T,N>,
		 &Multiply_gamma5_DP_arg<T,N>,
		 &Multiply_gamma6_DP_arg<T,N>,
		 &Multiply_gamma7_DP_arg<T,N>,
		 &Multiply_gamma8_DP_arg<T,N>,
		 &Multiply_gamma9_DP_arg<T,N>,
		 &Multiply_gamma10_DP_arg<T,N>,
		 &Multiply_gamma11_DP_arg<T,N>,
		 &Multiply_gamma12_DP_arg<T,N>,
		 &Multiply_gamma13_DP_arg<T,N>,
		 &Multiply_gamma14_DP_arg<T,N>,
		 &Multiply_gamma15_DP_arg<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[a.elem()](b);
}


template<class T, int N>
inline T OpMultiplyGammaTypeDP::operator()(const T &a, const GammaTypeDP<N>& b) const
{
  typedef T (*Ptrfunc)(const T&);
  Ptrfunc s[] = {&Multiply_DP_arg_gamma0<T,N>, 
		 &Multiply_DP_arg_gamma1<T,N>,
		 &Multiply_DP_arg_gamma2<T,N>,
		 &Multiply_DP_arg_gamma3<T,N>,
		 &Multiply_DP_arg_gamma4<T,N>,
		 &Multiply_DP_arg_gamma5<T,N>,
		 &Multiply_DP_arg_gamma6<T,N>,
		 &Multiply_DP_arg_gamma7<T,N>,
		 &Multiply_DP_arg_gamma8<T,N>,
		 &Multiply_DP_arg_gamma9<T,N>,
		 &Multiply_DP_arg_gamma10<T,N>,
		 &Multiply_DP_arg_gamma11<T,N>,
		 &Multiply_DP_arg_gamma12<T,N>,
		 &Multiply_DP_arg_gamma13<T,N>,
		 &Multiply_DP_arg_gamma14<T,N>,
		 &Multiply_DP_arg_gamma15<T,N>};

  // Using the function pointer table, call via indirection the desired function
  // based on the argument to Gamma. It is assumed the argument must be a valid 
  // integer since Gamma was constructed (and should check it).
  return s[b.elem()](a);
}

} // namespace QDP

#endif
