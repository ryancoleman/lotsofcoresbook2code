// -*- C++ -*-
/* $Id: xml_array.h,v 1.11 2009-05-02 03:59:33 edwards Exp $ 
 *
 * File: xml_array.h
 *
 * Description: 
 *
 * This file has been blatantly stolen from QDP++ (where it was known 
 * as multi.h.) It provides templated array functionality. 
 * Unlike multi.h which defined 1, 2, 3 and 4D arrays, here I am 
 * restricting myself to 1D arrays. These can hold other 1D arrays 
 *
 */

#ifndef XML_ARRAY_H
#define XML_ARRAY_H

#include <iostream>
#include <cmath>
#include <cstdlib>

//#include "xml_array2d.h"

namespace XMLArray {

  template<typename T> class Array2d ;

  //! Container for a multi-dimensional 1D array
  template<typename T> class Array
  {
  public:
    Array() {F=0;n1=0;}
    explicit Array(int ns1) {F=0;n1=0;resize(ns1);}
    ~Array() {delete[] F;}
  
    //! Copy constructor
    Array(const Array& s): n1(s.n1), F(0)
      {
	resize(n1);
      
	for(int i=0; i < n1; ++i)
	  F[i] = s.F[i];
      }
  
    //! Allocate mem for the array
    void resize(int ns1) 
      {
	if (ns1 < 0) 
	{
	  std::cerr<<"Invalid resize in 1d: size= " << ns1 << std::endl;
	  std::exit(1);
	}
	delete[] F; 
	n1 = ns1;
     
	F = new(std::nothrow) T[n1];
	if( F == 0x0 ) { 
	  std::cerr << "Failed to resize array: n1= " << n1 << std::endl;
	  std::exit(-1);
	}
      
      }
  
    //! Size of array
    int size() const {return n1;}
    int size1() const {return n1;}
  
    //! Equal operator uses underlying = of T
    Array<T>& operator=(const Array<T>& s1)
      {
	if (size() != s1.size())   // a simple check avoids resizing always
	  resize(s1.size());
    
	for(int i=0; i < n1; ++i)
	  F[i] = s1.F[i];
	return *this;
      }

    //! Equal operator: converts a 2d array to an 1d
    // can't make this work... although I  the 1d to 2d conversion works
    //Array<T>& operator=(const Array2d<T>& s1) ;
    
    //! Equal operator uses underlying = of T
    template<typename T1>
    Array<T>& operator=(const T1& s1)
      {
	if (F == 0)
	{
	  std::cerr << __func__ << ": Array - left hand side not initialized in T1=s1\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] = s1;
	return *this;
      }
  
    //! Set equal to a old-style C 1-D array
    Array<T>& operator=(const T s1[])
      {
	if (F == 0)
	{
	  std::cerr << __func__ << ": Array - left hand side not initialized = s1[]\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] = s1[i];
	return *this;
      }
  
    //! Add-replace on each element
    /*! Uses underlying += */
    Array<T>& operator+=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in +=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] += s1.F[i];
	return *this;
      }
  
    //! Subtract-replace on each element
    /*! Uses underlying -= */
    Array<T>& operator-=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in -=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] -= s1.F[i];
	return *this;
      }
  
    //! Mult-replace on each element
    /*! Uses underlying *= */
    Array<T>& operator*=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in *=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] *= s1.F[i];
	return *this;
      }
  
    //! Divide-replace on each element
    /*! Uses underlying /= */
    Array<T>& operator/=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in /=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] /= s1.F[i];
	return *this;
      }

    //! And-replace on each element
    /*! Uses underlying &= */
    Array<T>& operator&=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in &=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] &= s1.F[i];
	return *this;
      }

    //! Or-replace on each element
    /*! Uses underlying |= */
    Array<T>& operator|=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in |=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] |= s1.F[i];
	return *this;
      }

    //! ^-replace on each element
    /*! Uses underlying ^= */
    Array<T>& operator^=(const Array<T>& s1)
      {
	if (size() != s1.size())
	{
	  std::cerr << "Sizes incompatible in ^=\n";
	  std::exit(1);
	}
    
	for(int i=0; i < n1; ++i)
	  F[i] ^= s1.F[i];
	return *this;
      }

    //! add to all elements
    /*! Uses the underlying += */
    Array<T>& operator+=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] += s1 ;
	return *this;
      }

    //! subtruct from all elements
    /*! Uses the underlying -= */
    Array<T>& operator-=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] -= s1 ;
	return *this;
      }
  
    //! multiply  all elements
    /*! Uses the underlying *= */
    Array<T>& operator*=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] *= s1 ;
	return *this;
      }

    //! divide all elements
    /*! Uses the underlying /= */
    Array<T>& operator/=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] /= s1 ;
	return *this;
      }

    //! and all elements
    /*! Uses the underlying &= */
    Array<T>& operator&=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] &= s1 ;
	return *this;
      }

    //! or all elements
    /*! Uses the underlying |= */
    Array<T>& operator|=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] |= s1 ;
	return *this;
      }

    //! ^ all elements
    /*! Uses the underlying ^= */
    Array<T>& operator^=(const T& s1)
      {
	for(int i=0; i < n1; ++i)
	  F[i] ^= s1 ;
	return *this;
      }

    //! Return ref to a column slice
    const T* slice() const {return F;}
  
    //! Return ref to an element
    T& operator()(int i) {return F[i];}
  
    //! Return const ref to an element
    const T& operator()(int i) const {return F[i];}
  
    //! Return ref to an element
    T& operator[](int i) {return F[i];}
  
    //! Return const ref to an element
    const T& operator[](int i) const {return F[i];}
  
  private:
    int n1;
    T *F;
  };


  //---------------------------------------------------------------
  // Comparisons/recombinations

  //! Concatenate two Array's
  template<typename T>
  inline Array<T> concat(const Array<T>& l, const Array<T>& r)
  {
    Array<T> nz(l.size() + r.size());
    int j = 0;
    for(int i=0; i < l.size(); ++i)
      nz[j++] = l[i];
  
    for(int i=0; i < r.size(); ++i)
      nz[j++] = r[i];

    return nz;
  }

  //! Check if two Array's are the same
  template<typename T>
  inline bool operator==(const Array<T>& n1, const Array<T>& n2)
  {
    if (n1.size() == 0 || n1.size() != n2.size())
      return false;
    
    for(int i=0; i < n1.size(); ++i)
      if (n2[i] != n1[i])
	return false;

    return true;
  }
  
  //! Check if two Array's are no the same
  template<typename T>
  inline bool operator!=(const Array<T>& n1, const Array<T>& n2)
  {
    return ! (n1 == n2);
  }

  //! a < b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator<(const Array<T>& a, const Array<T>& b)
  {
    bool ret = false;
    int  len = (a.size() < b.size()) ? a.size() : b.size();

    for(int i=0; i < len; ++i)
    {
      if (a[i] != b[i])
	return (a[i] < b[i]) ? true : false;
    }
    
    return (a.size() == b.size()) ? false : (a.size() < b.size()) ? true : false;
  }

  //! a > b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator>(const Array<T>& a, const Array<T>& b)
  {
    bool ret = false;
    int  len = (a.size() < b.size()) ? a.size() : b.size();

    for(int i=0; i < len; ++i)
    {
      if (a[i] != b[i])
	return (a[i] > b[i]) ? true : false;
    }
    
    return (a.size() == b.size()) ? false : (a.size() > b.size()) ? true : false;
  }

  //! a <= b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator<=(const Array<T>& a, const Array<T>& b)
  {
    return (a < b) || (a == b);
  }

  //! a >= b
  /*! This definition follows that of string comparison */
  template<typename T>
  inline bool operator>=(const Array<T>& a, const Array<T>& b)
  {
    return (a > b) || (a == b);
  }


  //---------------------------------------------------------------
  // Basic math support
  //
  //!unary -
  template< typename T> 
  inline
  Array<T> operator-(const Array<T>& a)
  {
    Array<T> d(a.size());
    for(int i=0; i < d.size(); ++i)
      d[i] = -a[i];
    return d;
  }

  //! add Arrays
  template< typename T> 
  inline
  Array<T> operator+(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c+=b;
    return c;
  }
  
  //! subtract Arrays
  template< typename T> 
  inline
  Array<T> operator-(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c-=b;
    return c;
  }
  
  //! multiply Arrays
  template< typename T> 
  inline
  Array<T> operator*(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c*=b;
    return c;
  }
  
  //! divide Arrays
  template< typename T> 
  inline
  Array<T> operator/(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c/=b;
    return c;
  }

  //! and Arrays
  template< typename T> 
  inline
  Array<T> operator&(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c&=b;
    return c;
  }

  //! or Arrays
  template< typename T> 
  inline
  Array<T> operator|(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c|=b;
    return c;
  }

  //! ^ Arrays
  template< typename T> 
  inline
  Array<T> operator^(const Array<T>& a, const Array<T>& b)
  {
    Array<T> c(a); 
    c^=b;
    return c;
  }

  //! scalar + Array
  template< typename T> 
  inline
  Array<T> operator+(const T& s, const Array<T>& a)
  {
    Array<T> c(a); 
    c+=s;
    return c;
  }

  //! Array + scalar
  template< typename T> 
  inline
  Array<T> operator+(const Array<T>& a, const T& s)
  {
    Array<T> c(a); 
    c+=s;
    return c;
  }
  
  //! scalar - Array
  template< typename T> 
  inline
  Array<T> operator-(const T& s, const Array<T>& a)
  {
    Array<T> c(-a); 
    c+=s;
    return c;
  }
  //! Array - scalar
  template< typename T> 
  inline
  Array<T> operator-(const Array<T>& a, const T& s)
  {
    Array<T> c(a); 
    c-=s;
    return c;
  }

  //! scalar * Array
  template< typename T> 
  inline
  Array<T> operator*(const T& s, const Array<T>& a)
  {
    Array<T> c(a); 
    c*=s;
    return c;
  }

  //! Array * scalar
  template< typename T> 
  inline
  Array<T> operator*(const Array<T>& a, const T& s)
  {
    Array<T> c(a); 
    c*=s;
    return c;
  }

  //! scalar / Array
  template< typename T> 
  inline
  Array<T> operator/(const T& s, const Array<T>& a)
  {
    Array<T> c(a.size());
    c = s;
    c/= a;
    return c;
  }

  //! Array / scalar
  template< typename T> 
  inline
  Array<T> operator/(const Array<T>& a, const T& s)
  {
    Array<T> c(a); 
    c/=s;
    return c;
  }

  //! sqrt
  template< typename T> 
  inline
  Array<T> sqrt(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::sqrt(a[i]);
    }
    return c;
  }

 //! log
  template< typename T> 
  inline
  Array<T> log(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::log(a[i]);
    }
    return c;
  }

 //! sin
  template< typename T> 
  inline
  Array<T> sin(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::sin(a[i]);
    }
    return c;
  }


 //! cos
  template< typename T> 
  inline
  Array<T> cos(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::cos(a[i]);
    }
    return c;
  }

 //! tan
  template< typename T> 
  inline
  Array<T> tan(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::tan(a[i]);
    }
    return c;
  }

//! asin
  template< typename T> 
  inline
  Array<T> asin(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::asin(a[i]);
    }
    return c;
  }


 //! acos
  template< typename T> 
  inline
  Array<T> acos(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::acos(a[i]);
    }
    return c;
  }

 //! atan
  template< typename T> 
  inline
  Array<T> atan(const Array<T>& a)
  {
    Array<T> c(a.size()); 
    for(int i(0);i<a.size();i++)
    {
      T tt;
      tt = a[i];
      c[i] = std::atan(a[i]);
    }
    return c;
  }



  //! norm2 of an array
  template< typename T> 
  inline
  T norm2(const Array<T>& a)
  {
    T nn = a[0]*a[0];  // assumes at least 1 element
    for(int i=1; i < a.size(); ++i)
      nn += a[i]*a[i];

    return nn;
  }

} // namespace XMLArray

#endif
