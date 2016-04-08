//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_MKL_RANDOM_H
#define OHMMS_MKL_RANDOM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <ctime>
#include <sstream>
#include <string>
#include <limits>
#include <vector>
#include <mkl_vsl.h>

/** MKL generator declaration to be specialized for double and float 
 * 
 * MKLRandom uses the static functions to handle datatype
 * - uniform(stream) return a single random value [0,1)
 * - uniform(stream,d,n)  return n random values [0,1)
 * - gaussian(stream,d,n,center,sigma)  return n Gaussian random values with center & sigma
 */
template<typename T>
struct mkl_generator { };

template<>
struct mkl_generator<double>
{
  static inline double uniform(VSLStreamStatePtr stream) 
  {
    double r;
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,&r,0.0,1.0);
    return r;
  }

  static inline void uniform(VSLStreamStatePtr stream, double* restrict d, int n)
  {
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,n,d,0.0,1.0);
  }

  static inline void gaussian(VSLStreamStatePtr stream, double* restrict d,
      int n, double center, double sigma) 
  {
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2,stream,n,d,center,sigma);
  }
};

template<>
struct mkl_generator<float>
{
  static inline double uniform(VSLStreamStatePtr stream) 
  {
    float r;
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,&r,0.0f,1.0f);
    return r;
  }

  static inline void uniform(VSLStreamStatePtr stream, float* restrict d, int n)
  {
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,n,d,0.0,1.0);
  }

  static inline void gaussian(VSLStreamStatePtr stream, float* restrict d,
      int n, double center, double sigma) 
  {
    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2,stream,n,d,center,sigma);
  }
};

/** random number generator using mkl vectorize library
 * @tparam T real data type
 * @tparam RNG type of generator, default is MT
 */
template<typename T, int RNG=VSL_BRNG_SFMT19937>
class MKLRandom
{

public:
  /// real result type
  typedef T result_type;
  /// unsigned integer type
  typedef MKL_UINT uint_type;

  ///context number
  int myContext;
  ///number of contexts
  int nContexts;
  ///offset of the random seed
  int baseOffset;
  ///stream of this generator
  VSLStreamStatePtr myStream;

  std::string ClassName;
  std::string EngineName;

  ///default constructor
  explicit MKLRandom(uint_type iseed=911, const std::string& aname="mt19937")
    : ClassName("mkl"), EngineName(aname),
    myContext(0), nContexts(1), baseOffset(0)
  {
    vslNewStream( &myStream, RNG, iseed);
  }

  ///copy constructor
  MKLRandom(const MKLRandom& rng): ClassName(rng.ClassName), EngineName(rng.EngineName),
    myContext(rng.myContext), nContexts(rng.nContexts), baseOffset(rng.baseOffset)
  {
    vslCopyStream(&myStream,rng.myStream);
  }

  ///copy operator (unnecessary but why not)
  MKLRandom<T,RNG>& operator=(const MKLRandom& r)
  {
    ClassName=r.ClassName;
    EngineName=r.EngineName;
    myContext=r.myContext;
    nContexts=r.nContexts;
    baseOffset=r.baseOffset;
    vslCopyStream(&myStream,rng.myStream);
    return *this;
  }

  ~MKLRandom() {vslDeleteStream(&myStream); }

  /** initialize the generator
   * @param i thread index
   * @param nstr number of threads
   * @param iseed_in input seed
   *
   * Initialize generator with the seed.
   */
  void init(int i, int nstr, int iseed_in, uint_type offset=1)
  {
    uint_type baseSeed=iseed_in;
    myContext=i;
    nContexts=nstr;
    if(iseed_in<=0) baseSeed=make_seed(0,nstr);
    baseOffset=offset;
    this->seed(baseSeed);
  }

  ///get baseOffset
  inline int offset() const
  {
    return baseOffset;
  }
  ///assign baseOffset
  inline int& offset()
  {
    return baseOffset;
  }

  ///assign seed
  inline void seed(uint_type aseed)
  {
    VSLStreamStatePtr astream;
    vslNewStream( &astream, RNG+myContext, aseed);
    vslCopyStreamState(myStream,astream);
    vslDeleteStream(&astream);
  }

  /** return a random number [0,1)
   */
  inline result_type rand()
  {
    return mkl_generator<T>::uniform(myStream);
  }

  /** return a random number [0,1)
   */
  inline result_type operator()()
  {
    return mkl_generator<T>::uniform(myStream);
  }

  /** return a random integer
   */
  inline uint_type irand()
  {
    uint_type i;
    virnguniformbits32(VSL_RNG_METHOD_UNIFORM_STD, myStream, 1, &i )
    return i;
  }

  inline void generate_uniform(T* restrict d, int n)
  {
    mkl_generator<T>::uniform(myStream,d,n);
  }

  inline void generate_gaussian(T* restrict d, int n, T center=0.0, T sigma=1.0)
  {
    mkl_generator<T>::gaussian(myStream,d,n,center,sigma);
  }

  inline int state_size() const
  {
    return 0;
  }

  inline void read(std::istream& rin)
  {
    int memsize=vslGetStreamSize( myStream );
    char *memptr=new char[memsize];
    rin.read(memptr,memsize);
    vslLoadStreamM(&myStream, memptr);
    delete [] memptr;
  }

  inline void write(std::ostream& rout) const
  {
    int memsize=vslGetStreamSize( myStream );
    char *memptr=new char[memsize];
    vslSaveStreamM( myStream, memptr);
    rout.write(memptr,memsize);
    delete [] memptr;
  }

  inline void save(std::vector<uint_type>& curstate) const
  {
  }

  inline void load(const std::vector<uint_type>& newstate)
  {
  }

private:
};
#endif

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: MKLRandom.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
