/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

    The source code contained or described herein and all documents related
    to the source code ("Material") are owned by Intel Corporation or its
    suppliers or licensors.  Title to the Material remains with Intel
    Corporation or its suppliers and licensors.  The Material is protected
    by worldwide copyright laws and treaty provisions.  No part of the
    Material may be used, copied, reproduced, modified, published, uploaded,
    posted, transmitted, distributed, or disclosed in any way without
    Intel's prior express written permission.

    No license under any patent, copyright, trade secret or other
    intellectual property right is granted to or conferred upon you by
    disclosure or delivery of the Materials, either expressly, by
    implication, inducement, estoppel or otherwise.  Any license under such
    intellectual property rights must be express and approved by Intel in
    writing.
*/

#ifndef MATRIX_GENERAL_HPP
#define MATRIX_GENERAL_HPP
#include <limits>
#include <vector>
template<typename T>
inline bool check_diff(int n, int m, const T* restrict a, const T* restrict b)
{
  for(int i=0; i<16; ++i)
    std::cout << a[i*n+i] << " " << b[i*n+i] << std::endl;

  for(int i=0; i<n*n; ++i)
    if(abs(a[i]-b[i])>std::numeric_limits<T>::epsilon()) 
    {
      std::cout << "Something is wrong " << i << " " << a[i] << " " << b[i] << std::endl;
      return false;
    }
  std::cout << "OK " << std::endl;
  return true;
}

template<typename T, unsigned ALIGN>
struct Matrix
{

  typedef T value_type;
  int Mrows;
  int Mcols;
  int Mcols_max;
  int Nalloc;
  T* Base;

  explicit Matrix(int m=0, int n=-1): Nalloc(0),Base(nullptr)
  {
    resize(m,n);
  }

  ~Matrix()
  {
    if(Base!=nullptr) _mm_free(Base);
  }

  void resize(int m, int n)
  {
    if(m==0) return;
    Mrows=m;
    Mcols=(n<0)?Mrows:n;
    Mcols_max=((Mcols*sizeof(T)+ALIGN-1)/ALIGN)*(ALIGN/sizeof(T));
    if(Base!=nullptr) _mm_free(Base);
    Base = static_cast<T*>(_mm_malloc((sizeof(T))*Mrows*Mcols_max+16,ALIGN));
  }
  
  inline T* data() { return Base;}
  inline const T* data() const { return Base;}

  inline int rows() const { return Mrows;}
  inline int cols() const { return Mcols;}
  inline int cols_max() const { return Mcols_max;}

  inline T* operator[](int i) 
  {
    return Base+i*Mcols_max;
  }

  inline const T* operator[](int i) const
  {
    return Base+i*Mcols_max;
  }

  inline void init(T row, T col, T off)
  {
//#pragma omp parallel for 
    const T x=1.0/static_cast<T>(Mrows);
    for(int i=0; i<Mrows; i++)
    {
      T* b=Base+i*Mcols_max;
      for(int j=0; j<Mcols; j++) b[j]=(row*i+col*j+off)*x;
      for(int j=Mcols; j<Mcols_max; ++j) b[j]=T();
    }
  }


  inline void set(T val)
  {
    for(int i=0; i<Mrows; i++)
    {
      T* b=Base+i*Mcols_max;
      std::fill(b,b+Mcols,val);
      std::fill(b+Mcols,b+Mcols_max,T());
    }
  }
};

/** specialization of unaligned */
template<typename T>
struct Matrix<T,1>
{

  typedef T value_type;
  int Mrows;
  int Mcols;
  int Mcols_max;
  T* Base;

  explicit Matrix(int m, int n=-1): Mrows(m)
  {
    Mcols_max=Mcols=(n<0)?Mrows:n;
    Base=new T[Mrows*Mcols];
  }

  ~Matrix()
  {
    delete [] Base;
  }
  
  inline T* data() { return Base;}
  inline const T* data() const { return Base;}

  inline int rows() const { return Mrows;}
  inline int cols() const { return Mcols;}
  inline int cols_max() const { return Mcols;}

  inline T* operator[](int i) 
  {
    return Base+i*Mcols;
  }

  inline const T* operator[](int i) const
  {
    return Base+i*Mcols;
  }

  inline void init(T row, T col, T off)
  {
//#pragma omp parallel for 
    const T x=1.0/static_cast<T>(Mrows);
    for(int i=0; i<Mrows; i++)
    {
      T* b=Base+i*Mcols;
      for(int j=0; j<Mcols; j++) b[j]=(row*i+col*j+off)*x;
    }
  }

  inline void set(T val)
  {
    std::fill(Base,Base+Mrows*Mcols,val);
  }
};

#include <mkl.h>

template<typename MT>
MUST_INLINE void multiply_mkl(const MT& a, const MT& b, MT& c)
{
  double alpha = 1.0, beta = 0.;
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,c.rows(),c.cols(),a.cols(),
      alpha,a.data(),a.Mcols_max,b.data(),b.Mcols_max,beta,c.data(),c.Mcols_max);
}

template<typename T, unsigned NB>
MUST_INLINE void docopy(const T* restrict a, T* restrict b)
{
//#pragma unroll(8)
//#pragma ivdep
#pragma omp simd
  for(int k=0; k<NB;++k) b[k]=a[k];
}


#endif
