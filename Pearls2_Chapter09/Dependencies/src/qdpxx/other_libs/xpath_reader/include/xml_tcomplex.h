/*! $Id: xml_tcomplex.h,v 1.1 2005-04-27 19:48:33 edwards Exp $
 *
 * File: tcomplex.h
 * 
 * Description: This file defines a basic complex container.
 *  IT DOES NOT DEFINE COMPLEX ALGEBRA. It is intended so that
 *  you can use it to read complex things, and then copy them 
 *  to your own classes.
 *
 */

#ifndef COMPLEX_H
#define COMPLEX_H

namespace XMLTComplex {

// This is a templated class for complex things (not necessarily numbers
template < typename T > class TComplex {
 public:
  // Destructor
  ~TComplex() {};
  
  // Empty constructor
  TComplex() {};
  
  // Other onstructor
  TComplex(const T& _real, const T& _imag) {
    // COPY the "contents" of the reference into r
    r = T(_real);
    
    // COPY the "contents" of the reference into r
    i = T(_imag);
  }

  // Copy Constructor 
  TComplex(const TComplex& c)
    {
      r = c.r;
      i = c.i;
    }
  
  // Get a reference to real part
  // not constant so that x.real() = T 
  // assignment can work.
  T& real() { return r; }
  
  // Get 'imag' part
  // not constant sot that x.imag() = T 
  // can work.
  T& imag() { return i; }
  
 private:
  T r;
  T i;
};

} // namespace XMLTComplex

#endif
