// -*- C++ -*-
//
// QDP data parallel interface
//

#ifndef QDP_SPECIALIZATIONS_H
#define QDP_SPECIALIZATIONS_H

namespace QDP {


//
// Conversion routines. These cannot be implicit conversion functions
// since they foul up the PETE defs in QDPOperators.h using primitive
// types
//

//! Make an int from an Integer
inline int 
toInt(const Integer& s) 
{
  return toInt(s.elem());
}

//! Make a float from a Real32
inline float
toFloat(const Real32& s) 
{
  return toFloat(s.elem());
}

//! Make a double from a Real64
inline double
toDouble(const Real64& s) 
{
  return toDouble(s.elem());
}

//! Make a bool from a Boolean
inline bool
toBool(const Boolean& s) 
{
  return toBool(s.elem());
}

#ifdef QDP_USE_LIBXML2
// XML readers
template<>
void read(XMLReader& xml, const std::string& s, multi1d<Integer>& d);
template<>
void read(XMLReader& xml, const std::string& s, multi1d<Real32>& d);
template<>
void read(XMLReader& xml, const std::string& s, multi1d<Real64>& d);
template<>
void read(XMLReader& xml, const std::string& s, multi1d<Boolean>& d);

template<>
void read(XMLReader& xml, const std::string& s, std::vector<Integer>& d);
template<>
void read(XMLReader& xml, const std::string& s, std::vector<Real32>& d);
template<>
void read(XMLReader& xml, const std::string& s, std::vector<Real64>& d);
template<>
void read(XMLReader& xml, const std::string& s, std::vector<Boolean>& d);

template<>
void read(XMLReader& xml, const std::string& s, std::list<Integer>& d);
template<>
void read(XMLReader& xml, const std::string& s, std::list<Real32>& d);
template<>
void read(XMLReader& xml, const std::string& s, std::list<Real64>& d);
template<>
void read(XMLReader& xml, const std::string& s, std::list<Boolean>& d);

// XML writers
template<>
void write(XMLWriter& xml, const std::string& s, const multi1d<Integer>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const multi1d<Real32>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const multi1d<Real64>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const multi1d<Boolean>& d);

template<>
void write(XMLWriter& xml, const std::string& s, const std::vector<Integer>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const std::vector<Real32>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const std::vector<Real64>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const std::vector<Boolean>& d);

template<>
void write(XMLWriter& xml, const std::string& s, const std::list<Integer>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const std::list<Real32>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const std::list<Real64>& d);
template<>
void write(XMLWriter& xml, const std::string& s, const std::list<Boolean>& d);

#endif

//
// Return an equivalent QDP type given some simple machine type
//
template<>
struct SimpleScalar<float>
{
  typedef Real32   Type_t;
};

// Construct simple float word
template<>
struct SimpleScalar<int>
{
  typedef Integer   Type_t;
};

// Construct simple double word
template<>
struct SimpleScalar<double>
{
  typedef Real64   Type_t;
};

// Construct simple boolean word
template<>
struct SimpleScalar<bool>
{
  typedef Boolean   Type_t;
};


//
// Type constructors for QDP types within the type system. Namely,
// at some level like a primitive, sometimes scalar temporaries are needed
// These are the bottom most constructors given a machine type
//
// Construct simple word
template<>
struct InternalScalar<float>
{
  typedef float  Type_t;
};

template<>
struct InternalScalar<int>
{
  typedef int   Type_t;
};

template<>
struct InternalScalar<double>
{
  typedef double   Type_t;
};

template<>
struct InternalScalar<bool>
{
  typedef bool  Type_t;
};


// Makes a primitive scalar leaving grid alone
template<>
struct PrimitiveScalar<float>
{
  typedef float  Type_t;
};

template<>
struct PrimitiveScalar<int>
{
  typedef int   Type_t;
};

template<>
struct PrimitiveScalar<double>
{
  typedef double   Type_t;
};

template<>
struct PrimitiveScalar<bool>
{
  typedef bool  Type_t;
};



// Makes a lattice scalar leaving primitive indices alone
template<>
struct LatticeScalar<float>
{
  typedef float  Type_t;
};

template<>
struct LatticeScalar<int>
{
  typedef int   Type_t;
};

template<>
struct LatticeScalar<double>
{
  typedef double   Type_t;
};

template<>
struct LatticeScalar<bool>
{
  typedef bool  Type_t;
};



// Internally used real scalars
template<>
struct RealScalar<int> {
  typedef REAL32  Type_t;
};

template<>
struct RealScalar<float> {
  typedef REAL32  Type_t;
};

template<>
struct RealScalar<double> {
  typedef REAL64  Type_t;
};



//
// Leaf constructors for simple machine types. These are a specialization over the
// default constructors. The point is to avoid wrapping the simple types in
// references which do not help much (primitive objects are word size!), 
// but get in the way of type computations.
//
template<>
struct CreateLeaf<int>
{
  typedef int Inp_t;
  typedef Integer  Leaf_t;
  inline static
  Leaf_t make(const Inp_t &a) { return Leaf_t(a); }
};

template<>
struct CreateLeaf<float>
{
  typedef float Inp_t;
  typedef Real32  Leaf_t;
  inline static
  Leaf_t make(const Inp_t &a) { return Leaf_t(a); }
};

template<>
struct CreateLeaf<double>
{
  typedef double Inp_t;
  typedef Real64  Leaf_t;
  inline static
  Leaf_t make(const Inp_t &a) { return Leaf_t(a); }
};

template<>
struct CreateLeaf<bool>
{
  typedef bool Inp_t;
  typedef Boolean  Leaf_t;
  inline static
  Leaf_t make(const Inp_t &a) { return Leaf_t(a); }
};


template<>
struct CreateLeaf<OScalar<IntReal32> >
{
  typedef OScalar<IntReal32> Leaf_t;
  inline static
  Leaf_t make(const OScalar<IntReal32> &a) { return Leaf_t(a); }
};

template<>
struct CreateLeaf<OScalar<IntReal64> >
{
  typedef OScalar<IntReal64> Leaf_t;
  inline static
  Leaf_t make(const OScalar<IntReal64> &a) { return Leaf_t(a); }
};

template<>
struct CreateLeaf<OScalar<IntInteger> >
{
  typedef OScalar<IntInteger> Leaf_t;
  inline static
  Leaf_t make(const OScalar<IntInteger> &a) { return Leaf_t(a); }
};


template<>
struct CreateLeaf<OScalar<IntBoolean> >
{
  typedef OScalar<IntBoolean> Leaf_t;
  inline static
  Leaf_t make(const OScalar<IntBoolean> &a) { return Leaf_t(a); }
};


} // namespace QDP

#endif
