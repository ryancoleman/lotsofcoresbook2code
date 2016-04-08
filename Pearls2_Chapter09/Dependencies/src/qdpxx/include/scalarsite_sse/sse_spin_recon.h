#ifndef SSE_SPIN_RECON_H
#define SSE_SPIN_RECON_H

namespace QDP {

typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 4> Spin4_32;
typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 2> Spin2_32;

template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir0Minus> {
  typedef Spin4_32  Type_t;
};

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir0Minus>::Type_t
spinReconstructDir0Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir0Minus>::Type_t  d;

  inlineSpinReconDir0Minus( &(s1.elem(0).elem(0).real()),
			    &(d.elem(0).elem(0).real()),
			    1);
 

  return d;
}

template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir0Plus> {
  typedef Spin4_32  Type_t;
};

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir0Plus>::Type_t
spinReconstructDir0Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir0Plus>::Type_t  d;

  inlineSpinReconDir0Plus( &(s1.elem(0).elem(0).real()),
			   &(d.elem(0).elem(0).real()),
			   1);


  return d;
}


template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir1Minus> {
  typedef Spin4_32  Type_t;
};


template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir1Minus>::Type_t
spinReconstructDir1Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir1Minus>::Type_t  d;

  inlineSpinReconDir1Minus( &(s1.elem(0).elem(0).real()),
			    &(d.elem(0).elem(0).real()),
			    1);

 

  return d;
}

template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir1Plus> {
  typedef Spin4_32  Type_t;
};

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir1Plus>::Type_t
spinReconstructDir1Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir1Plus>::Type_t  d;

  inlineSpinReconDir1Plus( &(s1.elem(0).elem(0).real()),
			   &(d.elem(0).elem(0).real()),
			   1);


  return d;
}


template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir2Minus> {
  typedef Spin4_32  Type_t;
};


template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir2Minus>::Type_t
spinReconstructDir2Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir2Minus>::Type_t  d;

  inlineSpinReconDir2Minus( &(s1.elem(0).elem(0).real()),
			    &(d.elem(0).elem(0).real()),
			    1);

  return d;
}


template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir2Plus> {
  typedef Spin4_32  Type_t;
};

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir2Plus>::Type_t
spinReconstructDir2Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir2Plus>::Type_t  d;
 
  inlineSpinReconDir2Plus( &(s1.elem(0).elem(0).real()),
			   &(d.elem(0).elem(0).real()),
			   1);


  return d;
}


template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir3Minus> {
  typedef Spin4_32  Type_t;
};

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir3Minus>::Type_t
spinReconstructDir3Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir3Minus>::Type_t  d;
  
  inlineSpinReconDir3Minus( &(s1.elem(0).elem(0).real()),
			    &(d.elem(0).elem(0).real()),
			    1);

  return d;
}

template<>
struct UnaryReturn<Spin2_32, FnSpinReconstructDir3Plus> {
  typedef Spin4_32  Type_t;
};

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir3Plus>::Type_t
spinReconstructDir3Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir3Plus>::Type_t  d;
  
  inlineSpinReconDir3Plus( &(s1.elem(0).elem(0).real()),
			   &(d.elem(0).elem(0).real()),
			   1);
  
  return d;
}


} // namespace QDP


#endif
