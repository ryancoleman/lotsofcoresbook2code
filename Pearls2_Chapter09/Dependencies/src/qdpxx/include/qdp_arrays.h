// -*- C++ -*-
/*! \file
 *  \brief One based arrays
 */

#ifndef __QDP_ARRAYS_H__
#define __QDP_ARRAYS_H__

#include <vector>

namespace QDP 
{

  //----------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------
  //! 1d array, with 1-based indices
  template<typename T> class Array1dO
  {
  public:
    Array1dO() {}
    explicit Array1dO(int n) {d.resize(n);}   // initialize
    ~Array1dO() {}

    //! Copy constructor
    Array1dO(const Array1dO& s): d(s.d) {}

    //! Copy constructor from a zero-based array
    Array1dO(const multi1d<T>& s) {
      d.resize(s.size());
      for(int i=0; i < s.size(); ++i)
	d[i] = s[i];
    }

    //! Copy constructor from a zero-based array
    Array1dO(const std::vector<T>& s) : d(s) {}

    void resize(int n) {d.resize(n);}
    int size() const {return d.size();}

    Array1dO& operator=(const Array1dO& s) {d = s.d; return *this;}

    //! Return ref to an element
    T& operator[](int i) {return d[i-1];}
    
    //! Return const ref to an element
    const T& operator[](int i) const {return d[i-1];}

    //! Append to a vector - it will grow
    void push_back(const T& s1) {d.push_back(s1);}

    //! Return ref to underlying array
    const std::vector<T>& ref() const {return d;}

    //! Return ref to underlying array
    std::vector<T>& ref() {return d;}

  private:
    std::vector<T> d;
  };


  //----------------------------------------------------------------------------------
  //! 2d array, with 1-based indices
  template<typename T> class Array2dO
  {
  public:
    Array2dO() {}
    explicit Array2dO(int n) {d.resize(n,n);}   // initialize
    explicit Array2dO(int n2, int n1) {d.resize(n2,n1);}   // initialize
    ~Array2dO() {}

    //! Copy constructor
    Array2dO(const Array2dO& s): d(s.d) {}

    void resize(int n2, int n1) {d.resize(n2,n1);}
    void resize(int n) {d.resize(n,n);}
    int size() const {return d.size1();}
    int size1() const {return d.size1();}
    int size2() const {return d.size2();}

    //! Another variant on the size of the 2d array
    int nrows() const {return d.size2();}
    int ncols() const {return d.size1();}

    Array2dO& operator=(const Array2dO& s) {d = s.d; return *this;}

    //! Equal operator uses underlying = of T
    template<typename T1>
    Array2dO<T>& operator=(const T1& s1) {d = s1; return *this;}

    //! Return ref to an element
    T& operator()(int i, int j) {return d(i-1,j-1);}
    
    //! Return const ref to an element
    const T& operator()(int i, int j) const {return d(i-1,j-1);}

    //! Return ref to underlying array
    const multi2d<T>& ref() const {return d;}

    //! Return ref to underlying array
    multi2d<T>& ref() {return d;}

  private:
    multi2d<T> d;
  };


  //----------------------------------------------------------------------------------
  //! 3d array, with 1-based indices
  template<typename T> class Array3dO
  {
  public:
    Array3dO() {}
    explicit Array3dO(int n) {d.resize(n,n,n);}   // initialize
    explicit Array3dO(int n3, int n2, int n1) {d.resize(n3,n2,n1);}   // initialize
    ~Array3dO() {}

    //! Copy constructor
    Array3dO(const Array3dO& s): d(s.d) {}

    void resize(int n3, int n2, int n1) {d.resize(n3,n2,n1);}
    void resize(int n) {d.resize(n,n,n);}
    int size() const {return d.size1();}
    int size1() const {return d.size1();}
    int size2() const {return d.size2();}
    int size3() const {return d.size3();}

    //! Another variant on the size of the 3d array
    int leftSize()   const {return d.size3();}
    int middleSize() const {return d.size2();}
    int rightSize()  const {return d.size1();}

    Array3dO& operator=(const Array3dO& s) {d = s.d; return *this;}

    //! Equal operator uses underlying = of T
    template<typename T1>
    Array3dO<T>& operator=(const T1& s1) {d = s1; return *this;}

    //! Return ref to an element
    T& operator()(int i, int j, int k) {return d(i-1,j-1,k-1);}
    
    //! Return const ref to an element
    const T& operator()(int i, int j, int k) const {return d(i-1,j-1,k-1);}

    //! Return ref to underlying array
    const multi3d<T>& ref() const {return d;}

    //! Return ref to underlying array
    multi3d<T>& ref() {return d;}


  private:
    multi3d<T> d;
  };


  //----------------------------------------------------------------------------------
  //! 4d array, with 1-based indices
  template<typename T> class Array4dO
  {
  public:
    Array4dO() {}
    explicit Array4dO(int n) {d.resize(n,n,n,n);}   // initialize
    explicit Array4dO(int n4, int n3, int n2, int n1) {d.resize(n4,n3,n2,n1);}   // initialize
    ~Array4dO() {}

    //! Copy constructor
    Array4dO(const Array4dO& s): d(s.d) {}

    void resize(int n4, int n3, int n2, int n1) {d.resize(n4,n3,n2,n1);}
    void resize(int n) {d.resize(n,n,n,n);}
    int size() const {return d.size1();}
    int size1() const {return d.size1();}
    int size2() const {return d.size2();}
    int size3() const {return d.size3();}
    int size4() const {return d.size4();}

    Array4dO& operator=(const Array4dO& s) {d = s.d; return *this;}

    //! Equal operator uses underlying = of T
    template<typename T1>
    Array4dO<T>& operator=(const T1& s1) {d = s1; return *this;}

    //! Return ref to an element
    T& operator()(int i, int j, int k, int l) {return d(i-1,j-1,k-1,l-1);}
    
    //! Return const ref to an element
    const T& operator()(int i, int j, int k, int l) const {return d(i-1,j-1,k-1,l-1);}

    //! Return ref to underlying array
    const multi4d<T>& ref() const {return d;}

    //! Return ref to underlying array
    multi4d<T>& ref() {return d;}


  private:
    multi4d<T> d;
  };


  //----------------------------------------------------------------------------------
  //! 5d array, with 1-based indices
  template<typename T> class Array5dO
  {
  public:
    Array5dO() {}
    explicit Array5dO(int n) {d.resize(n,n,n,n,n);}   // initialize
    explicit Array5dO(int n5,int n4, int n3, int n2, int n1) {d.resize(n5,n4,n3,n2,n1);}   // initialize
    ~Array5dO() {}

    //! Copy constructor
    Array5dO(const Array5dO& s): d(s.d) {}

    void resize(int n5, int n4, int n3, int n2, int n1) {d.resize(n5,n4,n3,n2,n1);}
    void resize(int n) {d.resize(n,n,n,n,n);}
    int size() const {return d.size1();}
    int size1() const {return d.size1();}
    int size2() const {return d.size2();}
    int size3() const {return d.size3();}
    int size4() const {return d.size4();}
    int size5() const {return d.size5();}

    Array5dO& operator=(const Array5dO& s) {d = s.d; return *this;}

    //! Equal operator uses underlying = of T
    template<typename T1>
    Array5dO<T>& operator=(const T1& s1) {d = s1; return *this;}

    //! Return ref to an element
    T& operator()(int i, int j, int k, int l, int m) {return d(i-1,j-1,k-1,l-1,m-1);}
    
    //! Return const ref to an element
    const T& operator()(int i, int j, int k, int l, int m) const {return d(i-1,j-1,k-1,l-1,m-1);}

    //! Return ref to underlying array
    const multi5d<T>& ref() const {return d;}

    //! Return ref to underlying array
    multi5d<T>& ref() {return d;}


  private:
    multi5d<T> d;
  };


  //---------------------------------------------------------------
  // Basic math support
  //
  //!unary -
  template< typename T> 
  inline
  Array1dO<T> operator-(const Array1dO<T>& a)
  {
    Array1dO<T> d(a.size());
    for(int i=0; i < d.size(); ++i)
      d.ref()[i] = -a.ref()[i];
    return d;
  }

  //! scalar * Array
  template< typename T> 
  inline
  Array1dO<T> operator*(const T& s, const Array1dO<T>& a)
  {
    Array1dO<T> c(a.size());
    for(int i=0; i < c.size(); ++i)
      c.ref()[i] = s * a.ref()[i];
    return c;
  }

  //! Array / scalar
  template< typename T> 
  inline
  Array1dO<T> operator/(const Array1dO<T>& a, const T& s)
  {
    Array1dO<T> c(a.size());
    for(int i=0; i < c.size(); ++i)
      c.ref()[i] = a.ref()[i] / s;
    return c;
  }


#if 0
  //------------------------------------------------
  //!unary -
  template< typename T> 
  inline
  Array3dO<T> operator-(const Array3dO<T>& a)
  {
    Array3dO<T> d;
    d.ref() = -a.ref();
    return d;
  }

  //! scalar * Array
  template< typename T> 
  inline
  Array3dO<T> operator*(const T& s, const Array3dO<T>& a)
  {
    Array3dO<T> c;
    c = s * a.ref();
    return c;
  }

  //! Array / scalar
  template< typename T> 
  inline
  Array3dO<T> operator/(const Array3dO<T>& a, const T& s)
  {
    Array3dO<T> c;
    c = a.ref() / s;
    return c;
  }
#endif


}  // end namespace QDP

#endif
