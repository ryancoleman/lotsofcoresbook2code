// -*- C++ -*-
/*
 * File: xml_array3d.h
 *
 * Description: 
 *
 * A 4d array class
 */

#ifndef XML_ARRAY3D_H
#define XML_ARRAY3D_H

#include <iostream>
#include <cmath>
#include "xml_array.h"

namespace XMLArray 
{
  
  //! Container for a Array-dimensional 3D array
  template<typename T> class Array3d
  {
  public:
    Array3d() {F=0;n1=n2=n3=sz=0;}
    explicit Array3d(int ns3, int ns2, int ns1) {F=0;resize(ns3,ns2,ns1);}
    ~Array3d() {delete[] F;}

    //! Copy constructor
    Array3d(const Array3d& s): n1(s.n1), n2(s.n2), n3(s.n3), sz(s.sz), F(0)
      {
	resize(n3,n2,n1);

	for(int i=0; i < sz; ++i)
	  F[i] = s.F[i];
      }

    //! Allocate mem for the array
    void resize(int ns3, int ns2, int ns1) 
      {
	delete[] F; 
	n1 = ns1; 
	n2 = ns2;  
	n3 = ns3;  
	sz = n1*n2*n3; 
	F = new(std::nothrow) T[sz];
	if( F == 0x0 ) { 
          std::cerr << "Unable to new memory in Array3d::resize(): n_left= " << ns3 << " n_middle= " << ns2 << " n_right= " << ns1 << "\n";
	  exit(1);
	}
      }

    //! Size of array
    int size1() const {return n1;}
    int size2() const {return n2;}
    int size3() const {return n3;}

    //! Another variant on the size of the 3d array
    int leftSize()   const {return n3;}
    int middleSize() const {return n2;}
    int rightSize()  const {return n1;}

    //! Equal operator uses underlying = of T
    Array3d<T>& operator=(const Array3d<T>& s1)
      {
	resize(s1.size3(), s1.size2(), s1.size1());   // always resize

	for(int i=0; i < sz; ++i)
	  F[i] = s1.F[i];
	return *this;
      }

    //! Equal operator uses underlying = of T
    template<typename T1>
    Array3d<T>& operator=(const T1& s1)
      {
	if (F == 0)
	{
	  std::cerr << __func__ << ": Array3d - left hand side not initialized\n";
	  exit(1);
	}
    
	for(int i=0; i < sz; ++i)
	  F[i] = s1;
	return *this;
      }
  
    //! Add-replace on each element
    /*! Uses underlying += */
    Array3d<T>& operator+=(const Array3d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3())
	{
	  std::cerr << "Sizes incompatible in +=\n";
	  exit(1);
	}
    
	for(int i=0; i < sz; ++i)
	  F[i] += s1.F[i];
	return *this;
      }
  
    //! Subtract-replace on each element
    /*! Uses underlying -= */
    Array3d<T>& operator-=(const Array3d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3())
	{
	  std::cerr << "Sizes incompatible in -=\n";
	  exit(1);
	}
    
	for(int i=0; i < sz; ++i)
	  F[i] -= s1.F[i];
	return *this;
      }
  
    //! Mult-replace on each element
    /*! Uses underlying *= */
    Array3d<T>& operator*=(const Array3d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3())
	{
	  std::cerr << "Sizes incompatible in *=\n";
	  exit(1);
	}
    
	for(int i=0; i < sz; ++i)
	  F[i] *= s1.F[i];
	return *this;
      }
  
    //! Divide-replace on each element
    /*! Uses underlying /= */
    Array3d<T>& operator/=(const Array3d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3())
	{
	  std::cerr << "Sizes incompatible in /=\n";
	  exit(1);
	}
    
	for(int i=0; i < sz; ++i)
	  F[i] /= s1.F[i];
	return *this;
      }

    //! multiply   all ellements
    /*! Uses the underlying *= */
    Array3d<T>& operator*=(const T& s1)
      {
	for(int i=0; i < sz; ++i)
	  F[i] *= s1 ;
	return *this;
      }

    //! devide  all ellements
    /*! Uses the underlying *= */
    Array3d<T>& operator/=(const T& s1)
      {
	for(int i=0; i < sz; ++i)
	  F[i] /= s1 ;
	return *this;
      }

    //!unary -
    /*! Uses the underlying unary - */
    Array3d<T> operator-() const
      {
	Array3d<T> d(n3,n2,n1);
	for(int i=0; i < sz; ++i)
	  d.F[i] = -F[i] ;
	return d;
      }

    //! Return ref to an element
    T& operator()(int k, int j, int i) {return F[i+n1*(j+n2*(k))];}

    //! Return const ref to an element
    const T& operator()(int k, int j, int i) const {return F[i+n1*(j+n2*(k))];}

    //! Raw pointer
    inline operator T*()
    {
      return F;
    }

    //! Raw pointer
    inline operator const T*() const
    {
      return F;
    }

  private:
    int n1;
    int n2;
    int n3;
    int sz;
    T *F;
  };


  //---------------------------------------------------------------
  // Basic math support
  //
  //!unary -
  template< typename T> 
  inline
  Array3d<T> operator-(const Array3d<T>& a)
  {
    Array3d<T> d(a);
    d *= T(-1);
    return d;
  }

  //! add Arrays
  template< typename T> 
  inline
  Array3d<T> operator+(const Array3d<T>& a, const Array3d<T>& b)
  {
    Array3d<T> c(a); 
    c+=b;
    return c;
  }
  
  //! subtract Arrays
  template< typename T> 
  inline
  Array3d<T> operator-(const Array3d<T>& a, const Array3d<T>& b)
  {
    Array3d<T> c(a); 
    c-=b;
    return c;
  }
  
  //! scalar + Array
  template< typename T> 
  inline
  Array3d<T> operator+(const T& s, const Array3d<T>& a)
  {
    Array3d<T> c(a); 
    c+=s;
    return c;
  }

  //! Array + scalar
  template< typename T> 
  inline
  Array3d<T> operator+(const Array3d<T>& a, const T& s)
  {
    Array3d<T> c(a); 
    c+=s;
    return c;
  }
  
  //! scalar - Array
  template< typename T> 
  inline
  Array3d<T> operator-(const T& s, const Array3d<T>& a)
  {
    Array3d<T> c(-a); 
    c+=s;
    return c;
  }
  //! Array - scalar
  template< typename T> 
  inline
  Array3d<T> operator-(const Array3d<T>& a, const T& s)
  {
    Array3d<T> c(a); 
    c-=s;
    return c;
  }

  //! scalar * Array
  template< typename T> 
  inline
  Array3d<T> operator*(const T& s, const Array3d<T>& a)
  {
    Array3d<T> c(a);
    c*=s;
    return c;
  }

  //! Array * scalar
  template< typename T> 
  inline
  Array3d<T> operator*(const Array3d<T>& a, const T& s)
  {
    Array3d<T> c(a); 
    c*=s;
    return c;
  }

  //! Array / scalar
  template< typename T> 
  inline
  Array3d<T> operator/(const Array3d<T>& a, const T& s)
  {
    Array3d<T> c(a); 
    c/=s;
    return c;
  }


} // namespace XMLArray

#endif
