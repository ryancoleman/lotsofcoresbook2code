// -*- C++ -*-

/*! \file
 * \brief QDP Operations on built-in types
 */

#ifndef QDP_SIMPLEWORD_H
#define QDP_SIMPLEWORD_H

#include <cmath>

namespace QDP {

using std::abs;
using std::div;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::ceil;
using std::cos;
using std::cosh;
using std::exp;
using std::fabs;
using std::floor;
using std::fmod;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;

using std::isnan;
using std::isinf;
using std::isfinite;
using std::isnormal;


/*! \addtogroup simpleword Builtin word operations
 * \ingroup fiber
 *
 * Fiber operations on builtin types. These operations are in support of QDP.
 * Namely, the builtin types are the lowest level (most deeply nested) types
 * within the QDP type composition
 *
 * @{
 */

// All these are explicit to avoid any general template clashes

//! dest = 0
inline
void zero_rep(int& dest) 
{
  dest = 0;
}

//! dest = 0
inline
void zero_rep(float& dest) 
{
  dest = 0;
}

//! dest = 0
inline
void zero_rep(double& dest) 
{
  dest = 0;
}

//! No bool(dest) = 0


//! d = (mask) ? s1 : d;
inline
void copymask(int& d, bool mask, int s1) 
{
  if (mask)
    d = s1;
}

//! d = (mask) ? s1 : d;
inline
void copymask(float& d, bool mask, float s1) 
{
  if (mask)
    d = s1;
}

//! d = (mask) ? s1 : d;
inline
void copymask(double& d, bool mask, double s1) 
{
  if (mask)
    d = s1;
}


//---------------------------
//! dest [float type] = source [int type]
inline
void cast_rep(float& d, int s1)
{
  d = float(s1);
}

//! dest [float type] = source [float type]
inline
void cast_rep(float& d, float s1)
{
  d = float(s1);
}

//! dest [float type] = source [double type]
inline
void cast_rep(float& d, double s1)
{
  d = float(s1);
}


//! dest [float type] = source [int type]
inline
void cast_rep(double& d, int s1)
{
  d = double(s1);
}

//! dest [double type] = source [float type]
inline
void cast_rep(double& d, float s1)
{
  d = double(s1);
}

//! dest [double type] = source [double type]
inline
void cast_rep(double& d, double s1)
{
  d = double(s1);
}


//---------------------------
//! dest [float type] = source [int type]
inline
void recast_rep(float& d, int s1)
{
  cast_rep(d,s1);
}

//! dest [float type] = source [float type]
inline
void recast_rep(float& d, float s1)
{
  cast_rep(d,s1);
}

//! dest [float type] = source [double type]
inline
void recast_rep(float& d, double s1)
{
  cast_rep(d,s1);
}


//! dest [float type] = source [int type]
inline
void recast_rep(double& d, int s1)
{
  cast_rep(d,s1);
}

//! dest [double type] = source [float type]
inline
void recast_rep(double& d, float s1)
{
  cast_rep(d,s1);
}

//! dest [double type] = source [double type]
inline
void recast_rep(double& d, double s1)
{
  cast_rep(d,s1);
}





//-------------------------------------------------------
// Functions

// Conjugate
inline 
float conj(float l)
{
  return l;
}

// Conjugate
inline 
double conj(double l)
{
  return l;
}

// Conjugate
inline 
int conj(int l)
{
  return l;
}


// Transpose
inline 
float transpose(float l)
{
  return l;
}

// Transpose
inline 
double transpose(double l)
{
  return l;
}

// Transpose
inline 
int transpose(int l)
{
  return l;
}



// TRACE
// trace = Trace(source1)
inline 
float trace(float s1)
{
  return s1;
}


// trace = Trace(source1)
inline 
double trace(double s1)
{
  return s1;
}


// trace = Trace(source1)
inline 
int trace(int s1)
{
  return s1;
}


// GetSite is only non-trivial at an inner grid sit
inline float  getSite(float s1, int innersite) {return s1;}
inline double getSite(double s1, int innersite) {return s1;}
inline int    getSite(int s1, int innersite) {return s1;}
inline bool   getSite(bool s1, int innersite) {return s1;}



// int = toInt(source1)
inline 
int toInt(int s1)
{
  return s1;
}

// float = toFloat(source1)
inline 
float toFloat(float s1)
{
  return s1;
}

// double = toDouble(source1)
inline 
double toDouble(double s1)
{
  return s1;
}

// bool = toBool(source1)
inline 
bool toBool(bool s1)
{
  return s1;
}


// int = toWordType(<int>)
inline 
int toWordType(int s1)
{
  return s1;
}

// float = toWordType(<float>)
inline 
float toWordType(float s1)
{
  return s1;
}

// double = toWordType(<double>)
inline 
double toWordType(double s1)
{
  return s1;
}

// bool = toWordType(<bool>)
inline 
bool toWordType(bool s1)
{
  return s1;
}



// Where is the ? operator
inline 
int where(bool a, int b, int c)
{
  if (a) return b; else return c;
}

inline 
float where(bool a, float b, float c)
{
  if (a) return b; else return c;
}

inline 
double where(bool a, double b, double c)
{
  if (a) return b; else return c;
}



// Global sum over site indices only
inline
int sum(int s1)
{
  return s1;
}

inline
int localNorm2(int s1)
{
  return s1*s1;
}

inline
int localInnerProduct(int s1, int s2)
{
  return s1*s2;
}

inline
unsigned int sum(unsigned int s1)
{
  return s1;
}

inline
unsigned int localNorm2(unsigned int s1)
{
  return s1*s1;
}

inline
unsigned int localInnerProduct(unsigned int s1, unsigned int s2)
{
  return s1*s2;
}

inline
double sum(float s1)
{
  return double(s1);
}

inline
double localNorm2(float s1)
{
  return double(s1*s1);
}

inline
double localInnerProduct(float s1, float s2)
{
  return double(s1*s2);
}

inline
double sum(double s1)
{
  return s1;
}

inline
double localNorm2(double s1)
{
  return s1*s1;
}

inline
double localInnerProduct(float s1, double s2)
{
  return double(s1)*s2;
}

inline
double localInnerProduct(double s1, float s2)
{
  return s1*double(s2);
}

inline
double localInnerProduct(double s1, double s2)
{
  return s1*s2;
}



/*! @} */  // end of group simpleword

//-----------------------------------------------------------------------------
// Traits classes 
//-----------------------------------------------------------------------------

// Return types
template<>
struct UnaryReturn<int, FnSeedToFloat> {
  typedef float  Type_t;
};

template<>
struct UnaryReturn<int, FnSum> {
  typedef int  Type_t;
};

template<>
struct UnaryReturn<int, FnGlobalMax> {
  typedef int  Type_t;
};

template<>
struct UnaryReturn<int, FnGlobalMin> {
  typedef int  Type_t;
};

template<>
struct UnaryReturn<int, FnSumMulti> {
  typedef int  Type_t;
};

template<>
struct UnaryReturn<int, FnNorm2> {
  typedef int  Type_t;
};

template<>
struct BinaryReturn<int, int, FnInnerProduct> {
  typedef int  Type_t;
};

template<>
struct UnaryReturn<int, FnLocalNorm2> {
  typedef int  Type_t;
};

template<>
struct BinaryReturn<int, int, FnLocalInnerProduct> {
  typedef int  Type_t;
};

template<>
struct TrinaryReturn<bool, int, int, FnWhere> {
  typedef int  Type_t;
};


template<>
struct UnaryReturn<float, FnSum> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<float, FnGlobalMax> {
  typedef float  Type_t;
};

template<>
struct UnaryReturn<float, FnGlobalMin> {
  typedef float  Type_t;
};

template<>
struct UnaryReturn<float, FnSumMulti> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<float, FnNorm2> {
  typedef double  Type_t;
};

template<>
struct BinaryReturn<float, float, FnInnerProduct> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<float, FnLocalNorm2> {
  typedef double  Type_t;
};

template<>
struct BinaryReturn<float, float, FnLocalInnerProduct> {
  typedef double  Type_t;
};

template<>
struct TrinaryReturn<bool, float, float, FnWhere> {
  typedef float  Type_t;
};


template<>
struct UnaryReturn<double, FnSum> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<double, FnGlobalMax> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<double, FnGlobalMin> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<double, FnSumMulti> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<double, FnNorm2> {
  typedef double  Type_t;
};

template<>
struct BinaryReturn<double, double, FnInnerProduct> {
  typedef double  Type_t;
};

template<>
struct UnaryReturn<double, FnLocalNorm2> {
  typedef double  Type_t;
};

template<>
struct BinaryReturn<double, double, FnLocalInnerProduct> {
  typedef double  Type_t;
};

template<>
struct TrinaryReturn<bool, double, double, FnWhere> {
  typedef double  Type_t;
};



template<>
struct BinaryReturn<int, int, OpAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpAddAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpSubtractAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpMultiplyAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpDivideAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpModAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpBitwiseOrAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpBitwiseAndAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpBitwiseXorAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpLeftShiftAssign> {
  typedef int  Type_t;
};
 
template<>
struct BinaryReturn<int, int, OpRightShiftAssign> {
  typedef int  Type_t;
};
 

template<>
struct BinaryReturn<float, float, OpAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpAddAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpSubtractAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpMultiplyAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpDivideAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpModAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpBitwiseOrAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpBitwiseAndAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpBitwiseXorAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpLeftShiftAssign> {
  typedef float  Type_t;
};
 
template<>
struct BinaryReturn<float, float, OpRightShiftAssign> {
  typedef float  Type_t;
};
 

template<>
struct BinaryReturn<double, double, OpAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpAddAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpSubtractAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpMultiplyAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpDivideAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpModAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpBitwiseOrAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpBitwiseAndAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpBitwiseXorAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpLeftShiftAssign> {
  typedef double  Type_t;
};
 
template<>
struct BinaryReturn<double, double, OpRightShiftAssign> {
  typedef double  Type_t;
};
 

template<>
struct BinaryReturn<bool, bool, OpAssign> {
  typedef bool  Type_t;
};
 


template<>
struct TrinaryReturn<float, float, float, FnColorContract> {
  typedef float Type_t;
};

} // namespace QDP

#endif
