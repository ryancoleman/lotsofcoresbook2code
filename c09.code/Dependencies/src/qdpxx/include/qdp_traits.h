// -*- C++ -*-

/*! @file
 * @brief Traits classes
 * 
 * Traits classes needed internally
 */

#ifndef QDP_TRAITS_H
#define QDP_TRAITS_H

namespace QDP {

//-----------------------------------------------------------------------------
// Traits class for returning the subset-ted class name of a outer grid class
//-----------------------------------------------------------------------------

template<class T>
struct QDPSubTypeTrait {};

//-----------------------------------------------------------------------------
// Traits classes to support operations of simple scalars (floating constants, 
// etc.) on QDPTypes
//-----------------------------------------------------------------------------

//! Find the underlying word type of a field
template<class T>
struct WordType
{
  typedef T  Type_t;
};


//-----------------------------------------------------------------------------
// Traits Classes to support getting fixed precision versions of floating 
// precision classes
// ----------------------------------------------------------------------------
template<class T>
struct SinglePrecType
{
  typedef T Type_t; // This should never be instantiated as such
};

template<class T>
struct DoublePrecType
{
  typedef T Type_t; // This should never be instantiated as such
};

// Now we need to specialise to the bit whose precisions float
// The single prec types for both REAL32 and REAL64 are REAL32
template<>
struct SinglePrecType<REAL32>
{
  typedef REAL32 Type_t;
};

template<>
struct SinglePrecType<REAL64>
{
  typedef REAL32 Type_t;
};

// The Double prec types for both REAL32 and REAL64 are REAL64
template<>
struct DoublePrecType<REAL32>
{
  typedef REAL64 Type_t;
};

template<>
struct DoublePrecType<REAL64>
{
  typedef REAL64 Type_t;
};

  template<typename T> 
  struct DoublePrecType< multi1d< T > >
  {
    typedef multi1d< typename DoublePrecType< T >::Type_t > Type_t;
  };


//-----------------------------------------------------------------------------
// Constructors for simple word types
//-----------------------------------------------------------------------------

//! Construct simple word type. Default behavior is empty
template<class T>
struct SimpleScalar {};


//! Construct simple word type used at some level within primitives
template<class T>
struct InternalScalar {};


//! Makes a primitive scalar leaving grid alone
template<class T>
struct PrimitiveScalar {};


//! Makes a lattice scalar leaving primitive indices alone
template<class T>
struct LatticeScalar {};


//! Construct simple word type used at some level within primitives
template<class T>
struct RealScalar {};


//! Construct primitive type of input but always RScalar complex type
template<class T>
struct NoComplex {};


//! Simple zero tag
struct Zero {};

//! Put zero in some unnamed space
namespace {
 Zero zero;
}

} // namespace QDP

#endif

