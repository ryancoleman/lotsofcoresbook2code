#ifndef SSE_SPIN_PROJ_H
#define SSE_SPIN_PROJ_H

#include <stdio.h>

namespace QDP {
typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 4> Spin4_32;
typedef PSpinVector<PColorVector<RComplex<REAL32>, 3>, 2> Spin2_32;

//-------------------------
// Proj Dir=0 Minus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir0Minus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir0Minus>::Type_t
spinProjectDir0Minus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir0Minus>::Type_t  d;

  inlineSpinProjDir0Minus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);
			   
  
  return d;
}

//-------------------------
// Proj Dir=0 Plus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir0Plus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir0Plus>::Type_t
spinProjectDir0Plus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir0Plus>::Type_t  d;

 
  inlineSpinProjDir0Plus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);
			   
 
  
  return d;
}

//-------------------------
// Proj Dir=1 Minus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir1Minus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir1Minus>::Type_t
spinProjectDir1Minus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir1Minus>::Type_t  d;

  inlineSpinProjDir1Minus(&(s1.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);

  
  return d;
}


//-------------------------
// Proj Dir=1 Plus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir1Plus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir1Plus>::Type_t
spinProjectDir1Plus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir1Plus>::Type_t  d;

  inlineSpinProjDir1Plus(&(s1.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);
  return d;
}


//-------------------------
// Proj Dir=2 Minus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir2Minus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir2Minus>::Type_t
spinProjectDir2Minus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir2Minus>::Type_t  d;

  /* 1 - \gamma_2 =  1  0  -i  0 
                     0  1  0  +i
                    +i  0  1   0
                     0 -i  0   1 


   *      ( b0r + i b0i )  =  ( {a0r + a2i} + i{a0i - a2r} )
   *      ( b1r + i b1i )     ( {a1r - a3i} + i{a1i + a3r} )
   */
  inlineSpinProjDir2Minus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);
			  
  
  return d;
}


//-------------------------
// Proj Dir=2 Plus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir2Plus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir2Plus>::Type_t
spinProjectDir2Plus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir2Plus>::Type_t  d;


  inlineSpinProjDir2Plus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);

  
  return d;
}

//-------------------------
// Proj Dir=3 Minus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir3Minus> {
 typedef Spin2_32  Type_t;
};


template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir3Minus>::Type_t
spinProjectDir3Minus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir3Minus>::Type_t  d;

  inlineSpinProjDir3Minus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);
  
  return d;
}


//-------------------------
// Proj Dir=3 Plus
//-------------------------
template<>
struct  UnaryReturn<Spin4_32, FnSpinProjectDir3Plus> {
 typedef Spin2_32  Type_t;
};

template<>
inline UnaryReturn<Spin4_32, FnSpinProjectDir3Plus>::Type_t
spinProjectDir3Plus(const Spin4_32& s1)
{
  UnaryReturn<Spin4_32, FnSpinProjectDir3Plus>::Type_t  d;

  inlineSpinProjDir3Plus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);
  
  return d;
}

#if 0
template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir0Minus>::Type_t
spinReconstructDir0Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir0Minus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }

  /*                              ( 1  0  0 -i)  ( a0 )    ( a0 - i a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -i  0)  ( a1 )  = ( a1 - i a2 )
   *                    0         ( 0  i  1  0)  ( a2 )    ( a2 + i a1 )
   *                              ( i  0  0  1)  ( a3 )    ( a3 + i a0 )

   * The bottom components of be may be reconstructed using the formula

   *   ( b2r + i b2i )  =  ( {a2r - a1i} + i{a2i + a1r} )  =  ( - b1i + i b1r )
   *   ( b3r + i b3i )     ( {a3r - a0i} + i{a3i + a0r} )     ( - b0i + i b0r ) 
   */
  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() = -s1.elem(1).elem(col).imag();
    d.elem(2).elem(col).imag() =  s1.elem(1).elem(col).real();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() = -s1.elem(0).elem(col).imag();
    d.elem(3).elem(col).imag() =  s1.elem(0).elem(col).real();
  }
 

  return d;
}

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir0Plus>::Type_t
spinReconstructDir0Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir0Plus>::Type_t  d;

  inlineSpinReconDir0Plus(&(s1.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);


  return d;
}


template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir1Minus>::Type_t
spinReconstructDir1Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir1Minus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }

  /*                              ( 1  0  0  1)  ( a0 )    ( a0 + a3 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1 -1  0)  ( a1 )  = ( a1 - a2 )
   *                    1         ( 0 -1  1  0)  ( a2 )    ( a2 - a1 )
   *                              ( 1  0  0  1)  ( a3 )    ( a3 + a0 )
	 
   * The bottom components of be may be reconstructed using the formula

   *  ( b2r + i b2i )  =  ( {a2r - a1r} + i{a2i - a1i} )  =  ( - b1r - i b1i )
   *  ( b3r + i b3i )     ( {a3r + a0r} + i{a3i + a0i} )     (   b0r + i b0i ) 
   */
  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() = -s1.elem(1).elem(col).real();
    d.elem(2).elem(col).imag() = -s1.elem(1).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() =  s1.elem(0).elem(col).real();
    d.elem(3).elem(col).imag() =  s1.elem(0).elem(col).imag();
  }
 

  return d;
}

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir1Plus>::Type_t
spinReconstructDir1Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir1Plus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }
  /* 1 + \gamma_1 =  1  0  0 -1 
                     0  1  1  0
                     0  1  1  0
                    -1  0  0  1 
 

   *   ( b2r + i b2i )  =  ( {a2r + a1r} + i{a2i + a1i} )  =  (   b1r + i b1i )
   *   ( b3r + i b3i )     ( {a3r - a0r} + i{a3i - a0i} )     ( - b0r - i b0i ) 
  
  */
  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() =  s1.elem(1).elem(col).real();
    d.elem(2).elem(col).imag() =  s1.elem(1).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() = -s1.elem(0).elem(col).real();
    d.elem(3).elem(col).imag() = -s1.elem(0).elem(col).imag();
  }
 

  return d;
}


template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir2Minus>::Type_t
spinReconstructDir2Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir2Minus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }

  /*                              ( 1  0 -i  0)  ( a0 )    ( a0 - i a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0  i)  ( a1 )  = ( a1 + i a3 )
   *                    2         ( i  0  1  0)  ( a2 )    ( a2 + i a0 )
   *                              ( 0 -i  0  1)  ( a3 )    ( a3 - i a1 )

   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0i} + i{a2i + a0r} )  =  ( - b0i + i b0r )
   *  ( b3r + i b3i )     ( {a3r + a1i} + i{a3i - a1r} )     (   b1i - i b1r )
   */

  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() = -s1.elem(0).elem(col).imag();
    d.elem(2).elem(col).imag() =  s1.elem(0).elem(col).real();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() =  s1.elem(1).elem(col).imag();
    d.elem(3).elem(col).imag() = -s1.elem(1).elem(col).real();
  }
 

  return d;
}

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir2Plus>::Type_t
spinReconstructDir2Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir2Plus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }


  /* 1 + \gamma_2 =  1  0  i  0 
                     0  1  0 -i
                    -i  0  1  0
                     0  i  0  1 
		     
   *  ( b2r + i b2i )  =  ( {a2r + a0i} + i{a2i - a0r} )  =  (   b0i - i b0r )
   *  ( b3r + i b3i )     ( {a3r - a1i} + i{a3i + a1r} )     ( - b1i + i b1r ) 
  */
  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() =  s1.elem(0).elem(col).imag();
    d.elem(2).elem(col).imag() = -s1.elem(0).elem(col).real();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() = -s1.elem(1).elem(col).imag();
    d.elem(3).elem(col).imag() =  s1.elem(1).elem(col).real();
  }
 

  return d;
}


template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir3Minus>::Type_t
spinReconstructDir3Minus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir3Minus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }

  /*                              ( 1  0 -1  0)  ( a0 )    ( a0 - a2 )
   *  B  :=  ( 1 - Gamma  ) A  =  ( 0  1  0 -1)  ( a1 )  = ( a1 - a3 )
   *                    3         (-1  0  1  0)  ( a2 )    ( a2 - a0 )
   *                              ( 0 -1  0  1)  ( a3 )    ( a3 - a1 )
      
   * The bottom components of be may be reconstructed using the formula
   *  ( b2r + i b2i )  =  ( {a2r - a0r} + i{a2i - a0i} )  =  ( - b0r - i b0i )
   *  ( b3r + i b3i )     ( {a3r - a1r} + i{a3i - a1i} )     ( - b1r - i b1i ) 
   */
  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() = -s1.elem(0).elem(col).real();
    d.elem(2).elem(col).imag() = -s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() = -s1.elem(1).elem(col).real();
    d.elem(3).elem(col).imag() = -s1.elem(1).elem(col).imag();
  }
 

  return d;
}

template<>
inline UnaryReturn<Spin2_32, FnSpinReconstructDir3Plus>::Type_t
spinReconstructDir3Plus(const Spin2_32& s1)
{
  UnaryReturn<Spin2_32, FnSpinReconstructDir3Plus>::Type_t  d;

  for(int col=0; col < 3; col++) { 
    d.elem(0).elem(col).real() = s1.elem(0).elem(col).real();
    d.elem(0).elem(col).imag() = s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(1).elem(col).real() = s1.elem(1).elem(col).real();
    d.elem(1).elem(col).imag() = s1.elem(1).elem(col).imag();
  }


  /*                              ( 1  0  1  0)  ( a0 )    ( a0 + a2 )
   *  B  :=  ( 1 + Gamma  ) A  =  ( 0  1  0  1)  ( a1 )  = ( a1 + a3 )
   *                    3         ( 1  0  1  0)  ( a2 )    ( a2 + a0 )
   *                              ( 0  1  0  1)  ( a3 )    ( a3 + a1 )
   
   * The bottom components of be may be reconstructed using the formula
   
   *   ( b2r + i b2i )  =  ( {a2r + a0r} + i{a2i + a0i} )  =  ( b0r + i b0i )
   *   ( b3r + i b3i )     ( {a3r + a1r} + i{a3i + a1i} )     ( b1r + i b1i ) 
   */
  
  for(int col=0; col < 3; col++) { 
    d.elem(2).elem(col).real() =  s1.elem(0).elem(col).real();
    d.elem(2).elem(col).imag() =  s1.elem(0).elem(col).imag();
  }

  for(int col=0; col < 3; col++) { 
    d.elem(3).elem(col).real() =  s1.elem(1).elem(col).real();
    d.elem(3).elem(col).imag() =  s1.elem(1).elem(col).imag();
  }
 

  return d;
}

#endif

} // namespace QDP;

#endif
