#ifndef XML_ARRAY4D_H
#define XML_ARRAY4D_H

#include <iostream>
#include <cmath>
#include "xml_array.h"

namespace XMLArray 
{
  
  //! Container for a Array-dimensional 4D array
  template<typename T> class Array4d
  {
  public:
    Array4d() {F=0;n1=n2=n3=n4=sz=0;}
    explicit Array4d(int ns4, int ns3, int ns2, int ns1) {F=0;resize(ns4,ns3,ns2,ns1);}
    ~Array4d() {delete[] F;}

    //! Copy constructor
    Array4d(const Array4d& s): n1(s.n1), n2(s.n2), n3(s.n3), n4(s.n4), sz(s.sz), F(0)
      {
	resize(n4,n3,n2,n1);

	for(int i=0; i < sz; ++i)
	  F[i] = s.F[i];
      }

    //! Allocate mem for the array
    void resize(int ns4, int ns3, int ns2, int ns1) 
      {
	delete[] F; 
	n1 = ns1; 
	n2 = ns2;  
	n3 = ns3;  
	n4 = ns4;
	sz = n1*n2*n3*n4; 
	F = new(std::nothrow) T[sz];
	if( F == 0x0 ) { 
	  std::cerr << "Unable to new memory in Array4d::resize()\n";
	  exit(1);
	}
      }

    //! Size of array
    int size1() const {return n1;}
    int size2() const {return n2;}
    int size3() const {return n3;}
    int size4() const {return n4;}

    //! Equal operator uses underlying = of T
    Array4d<T>& operator=(const Array4d<T>& s1)
      {
	resize(s1.size4(), s1.size3(), s1.size2(), s1.size1());   // always resize

	for(int i=0; i < sz; ++i)
	  F[i] = s1.F[i];
	return *this;
      }

    //! Equal operator uses underlying = of T
    template<typename T1>
    Array4d<T>& operator=(const T1& s1)
      {
	if (F == 0)
	{
	  std::cerr << "left hand side not initialized\n";
	  exit(1);
	}
    
	for(int i=0; i < sz; ++i)
	  F[i] = s1;
	return *this;
      }
  
    //! Add-replace on each element
    /*! Uses underlying += */
    Array4d<T>& operator+=(const Array4d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3() || size4() != s1.size4() )
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
    Array4d<T>& operator-=(const Array4d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3()  || size4() != s1.size4() )
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
    Array4d<T>& operator*=(const Array4d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3()  || size4() != s1.size4() )
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
    Array4d<T>& operator/=(const Array4d<T>& s1)
      {
	if (size1() != s1.size1() || size2() != s1.size2() || size3() != s1.size3()  || size4() != s1.size4()  )
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
    Array4d<T>& operator*=(const T& s1)
      {
	for(int i=0; i < sz; ++i)
	  F[i] *= s1 ;
	return *this;
      }

    //! devide  all ellements
    /*! Uses the underlying *= */
    Array4d<T>& operator/=(const T& s1)
      {
	for(int i=0; i < sz; ++i)
	  F[i] /= s1 ;
	return *this;
      }

    //!unary -
    /*! Uses the underlying unary - */
    Array4d<T> operator-()
      {
	Array4d<T> d(n3,n2,n1);
	for(int i=0; i < sz; ++i)
	  d.F[i] = -F[i] ;
	return d;
      }

    //! Return ref to an element
    //T& operator()(int k, int j, int i) {return F[i+n1*(j+n2*(k))];}
    T& operator()(int l, int k, int j, int i) {return F[ i + n1*( j + n2*(k + n3*l )   )];}

    //! Return const ref to an element
    const T& operator()(int l, int k, int j, int i) const {return F[i + n1*( j + n2*( k + n3*l )    )];}

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
    int n4;
    int sz;
    T *F;
  };


} // namespace XMLArray

#endif
