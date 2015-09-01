// -*- C++ -*-
//
// QDP data parallel interface
//
// Outer lattice routines specific to a scalar platform 

#ifndef QDP_SCALAR_SPECIFIC_H
#define QDP_SCALAR_SPECIFIC_H

namespace QDP {

// Use separate defs here. This will cause subroutine calls under g++

//-----------------------------------------------------------------------------
// Layout stuff specific to a scalar architecture
namespace Layout
{
  //! coord[mu]  <- mu  : fill with lattice coord in mu direction
  LatticeInteger latticeCoordinate(int mu);
}


//-----------------------------------------------------------------------------
// Internal ops designed to look like those in parscalar
// These dummy routines exist just to make code more portable
namespace QDPInternal
{
  //! Dummy array sum accross all nodes
  template<class T>
  inline void globalSumArray(T* dest, int n) {}

  //! Dummy global sum on a multi1d
  template<class T>
  inline void globalSumArray(multi1d<T>& dest) {}

  //! Dummy global sum on a multi2d
  template<class T>
  inline void globalSumArray(multi2d<T>& dest) {}

  //! Dummy sum across all nodes
  template<class T>
  inline void globalSum(T& dest) {}

  //! Dummy global And
  inline void globalAnd(bool& in) {}

  //! Dummy global Or
  inline void globalOr(bool& in) {}

  //! Dummy broadcast from primary node to all other nodes
  template<class T>
  inline void broadcast(T& dest) {}

  //! Dummy broadcast from primary node to all other nodes
  template<>
  inline void broadcast(std::string& dest) {}

  //! Dummy broadcast a string from primary node to all other nodes
  inline void broadcast_str(std::string& dest) {}

  //! Dummy broadcast from primary node to all other nodes
  inline void broadcast(void* dest, size_t nbytes) {}
}

/////////////////////////////////////////////////////////
// Threading evaluate with openmp and qmt implementation
//
// by Xu Guo, EPCC, 16 June 2008
/////////////////////////////////////////////////////////

//! user argument for the evaluate function:
// "OLattice Op Scalar(Expression(source)) under an Subset"
//
template<class T, class T1, class Op, class RHS>
struct u_arg{
        OLattice<T>& d;
        const QDPExpr<RHS,OScalar<T1> >& r;
        const Op& op;
        const int *tab;
  u_arg( OLattice<T>& d_,
	 const QDPExpr<RHS, OScalar<T1> >& r_,
	 const Op& op_,
	 const int *tab_ ) : d(d_), r(r_), op(op_), tab(tab_) {}
   };

//! user function for the evaluate function:
// "OLattice Op Scalar(Expression(source)) under an Subset"
//
template<class T, class T1, class Op, class RHS>
void ev_userfunc(int lo, int hi, int myId, u_arg<T,T1,Op,RHS> *a)
{
   OLattice<T>& dest = a->d;
   const QDPExpr<RHS,OScalar<T1> >&rhs = a->r;
   const int* tab = a->tab;
   const Op& op= a->op;

      
   for(int j=lo; j < hi; ++j)
   {
     int i = tab[j];
     op(dest.elem(i), forEach(rhs, EvalLeaf1(0), OpCombine()));
   }
}


//! user argument for the evaluate function:
// "OLattice Op OLattice(Expression(source)) under an Subset"
//
template<class T, class T1, class Op, class RHS>
struct user_arg{
        OLattice<T>& d;
        const QDPExpr<RHS,OLattice<T1> >& r;
        const Op& op;
        const int *tab;
  user_arg(OLattice<T>& d_,
	   const QDPExpr<RHS,OLattice<T1> >& r_,
	   const Op& op_,
	   const int *tab_) : d(d_), r(r_), op(op_), tab(tab_) {}

   };

//! user function for the evaluate function:
// "OLattice Op OLattice(Expression(source)) under an Subset"
//
template<class T, class T1, class Op, class RHS>
void evaluate_userfunc(int lo, int hi, int myId, user_arg<T,T1,Op,RHS> *a)
{

   OLattice<T>& dest = a->d;
   const QDPExpr<RHS,OLattice<T1> >&rhs = a->r;
   const int* tab = a->tab;
   const Op& op= a->op;

      
   for(int j=lo; j < hi; ++j)
   {
     int i = tab[j];
     op(dest.elem(i), forEach(rhs, EvalLeaf1(i), OpCombine()));
   }
}

//! include the header file for dispatch
#include "qdp_dispatch.h"

//-----------------------------------------------------------------------------
//! OLattice Op Scalar(Expression(source)) under an Subset
/*! 
 * OLattice Op Expression, where Op is some kind of binary operation 
 * involving the destination 
 */
template<class T, class T1, class Op, class RHS>
//inline
void evaluate(OLattice<T>& dest, const Op& op, const QDPExpr<RHS,OScalar<T1> >& rhs,
	      const Subset& s)
{
//  cerr << "In evaluateUnorderedSubet(olattice,oscalar)\n";

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest, op, rhs);
  prof.time -= getClockTime();
#endif

  int numSiteTable = s.numSiteTable();
  
  u_arg<T,T1,Op,RHS> a(dest, rhs, op, s.siteTable().slice());

  dispatch_to_threads< u_arg<T,T1,Op,RHS> >(numSiteTable, a, ev_userfunc);
 
  ///////////////////
  // Original code
  //////////////////
  //const int *tab = s.siteTable().slice();
  //for(int j=0; j < s.numSiteTable(); ++j) 
  //{
  //int i = tab[j];
//    fprintf(stderr,"eval(olattice,oscalar): site %d\n",i);
//    op(dest.elem(i), forEach(rhs, ElemLeaf(), OpCombine()));
  //op(dest.elem(i), forEach(rhs, EvalLeaf1(0), OpCombine()));
  //}

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif
}



//! OLattice Op OLattice(Expression(source)) under an Subset
/*! 
 * OLattice Op Expression, where Op is some kind of binary operation 
 * involving the destination 
 */
template<class T, class T1, class Op, class RHS>
//inline
void evaluate(OLattice<T>& dest, const Op& op, const QDPExpr<RHS,OLattice<T1> >& rhs,
	      const Subset& s)
{
//  cerr << "In evaluateSubset(olattice,olattice)" << endl;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest, op, rhs);
  prof.time -= getClockTime();
#endif

  int numSiteTable = s.numSiteTable();

  user_arg<T,T1,Op,RHS> a(dest, rhs, op, s.siteTable().slice());

  dispatch_to_threads<user_arg<T,T1,Op,RHS> >(numSiteTable, a, evaluate_userfunc);

  ////////////////////
  // Original code
  ///////////////////

  // General form of loop structure
  //const int *tab = s.siteTable().slice();
  //for(int j=0; j < s.numSiteTable(); ++j) 
  //{
  //int i = tab[j];
//    fprintf(stderr,"eval(olattice,olattice): site %d\n",i);
  //op(dest.elem(i), forEach(rhs, EvalLeaf1(i), OpCombine()));
  //}

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif
}



//-----------------------------------------------------------------------------
template<class T, class T1, class Op, class RHS>
//inline
void evaluate_F(T* dest, const Op& op, const QDPExpr<RHS,OLattice<T1> >& rhs,
	      const Subset& s)
{
  //cerr << "In evaluate_F(olattice,olattice)" << endl;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest, op, rhs);
  prof.time -= getClockTime();
#endif

  // int numSiteTable = s.numSiteTable();
  // user_arg<T,T1,Op,RHS> a(dest, rhs, op, s.siteTable().slice());
  // dispatch_to_threads< user_arg<T,T1,Op,RHS> >(numSiteTable, a, evaluate_userfunc);

  ////////////////////
  // Original code
  ///////////////////

  //QDP_info("eval_F %d sites",s.numSiteTable());

  // General form of loop structure
  const int *tab = s.siteTable().slice();

#pragma omp parallel for
  for(int j=0; j < s.numSiteTable(); ++j) 
  {
    int i = tab[j];
    //fprintf(stderr,"eval(olattice,olattice): site %d\n",i);
    op( dest[j], forEach(rhs, EvalLeaf1(i), OpCombine()));
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif
}


//-----------------------------------------------------------------------------
//! dest = (mask) ? s1 : dest
template<class T1, class T2> 
void 
copymask(OSubLattice<T2> d, const OLattice<T1>& mask, const OLattice<T2>& s1) 
{
  OLattice<T2>& dest = d.field();
  const Subset& s = d.subset();

  const int *tab = s.siteTable().slice();
#pragma omp parallel for
  for(int j=0; j < s.numSiteTable(); ++j) 
  {
    int i = tab[j];
    copymask(dest.elem(i), mask.elem(i), s1.elem(i));
  }
}


//! dest = (mask) ? s1 : dest
template<class T1, class T2> 
void 
copymask(OLattice<T2>& dest, const OLattice<T1>& mask, const OLattice<T2>& s1) 
{
  const int vvol = Layout::vol();

#pragma omp parallel for
  for(int i=0; i < vvol; ++i) 
    copymask(dest.elem(i), mask.elem(i), s1.elem(i));
}



//-----------------------------------------------------------------------------
// Random numbers
namespace RNG
{
  extern Seed ran_seed;
  extern Seed ran_mult;
  extern Seed ran_mult_n;
  extern LatticeSeed *lattice_ran_mult;
}


//! dest  = random  
/*! This implementation is correct for no inner grid */
template<class T>
void 
random(OScalar<T>& d)
{
  Seed seed = RNG::ran_seed;
  Seed skewed_seed = RNG::ran_seed * RNG::ran_mult;

  fill_random(d.elem(), seed, skewed_seed, RNG::ran_mult);

  RNG::ran_seed = seed;  // The seed from any site is the same as the new global seed
}


//! dest  = random    under a subset
template<class T>
void 
random(OLattice<T>& d, const Subset& s)
{
  // Grab table array
  const int *tab = s.siteTable().slice();

#pragma omp parallel 
  {
    Seed seed;
    Seed skewed_seed;

#pragma omp for // need the barrier to avoid that RNG::ran_seed is changed too early
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      seed = RNG::ran_seed;
      skewed_seed.elem() = RNG::ran_seed.elem() * RNG::lattice_ran_mult->elem(i);
      fill_random(d.elem(i), seed, skewed_seed, RNG::ran_mult_n);
    }

#pragma omp critical (random)
    {
      RNG::ran_seed = seed;  // The seed from any site is the same as the new global seed
    }
  }

}


//! dest  = random   under a subset
template<class T>
void random(OSubLattice<T> dd)
{
  OLattice<T>& d = dd.field();
  const Subset& s = dd.subset();

  random(d,s);
}


//! dest  = random  
template<class T>
void random(OLattice<T>& d)
{
  random(d,all);
}


//! dest  = gaussian   under a subset
template<class T>
void gaussian(OLattice<T>& d, const Subset& s)
{
  OLattice<T>  r1, r2;

  random(r1,s);
  random(r2,s);

  const int *tab = s.siteTable().slice();

#pragma omp parallel for
  for(int j=0; j < s.numSiteTable(); ++j) 
  {
    int i = tab[j];
    fill_gaussian(d.elem(i), r1.elem(i), r2.elem(i));
  }
}



//! dest  = gaussian   under a subset
template<class T>
void gaussian(OSubLattice<T> dd)
{
  OLattice<T>& d = dd.field();
  const Subset& s = dd.subset();

  gaussian(d,s);
}


//! dest  = gaussian
template<class T>
void gaussian(OLattice<T>& d)
{
  gaussian(d,all);
}



//-----------------------------------------------------------------------------
// Broadcast operations
//! dest  = 0 
template<class T> 
void zero_rep(OLattice<T>& dest, const Subset& s) 
{
  const int *tab = s.siteTable().slice();

#pragma omp parallel for
  for(int j=0; j < s.numSiteTable(); ++j) 
  {
    int i = tab[j];
    zero_rep(dest.elem(i));
  }
}


#if 0
//! dest  = 0 
template<class T, class S>
void zero_rep(OSubLattice<T> dd) 
{
  OLattice<T>& d = dd.field();
  const Subset& s = dd.subset();
  
  zero_rep(d,s);
}
#endif

//! dest  = 0 
template<class T> 
void zero_rep(OLattice<T>& dest) 
{
  const int vvol = Layout::vol();
#pragma omp parallel for
  for(int i=0; i < vvol; ++i) 
    zero_rep(dest.elem(i));
}



//-----------------------------------------------
// Global sums
//! OScalar = sum(OScalar) under an explicit subset
/*!
 * Allow a global sum that sums over the lattice, but returns an object
 * of the same primitive type. E.g., contract only over lattice indices
 */
template<class RHS, class T>
typename UnaryReturn<OScalar<T>, FnSum>::Type_t
sum(const QDPExpr<RHS,OScalar<T> >& s1, const Subset& s)
{
  typename UnaryReturn<OScalar<T>, FnSum>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  evaluate(d,OpAssign(),s1,all);   // since OScalar, no global sum needed

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! OScalar = sum(OScalar)
/*!
 * Allow a global sum that sums over the lattice, but returns an object
 * of the same primitive type. E.g., contract only over lattice indices
 */
template<class RHS, class T>
typename UnaryReturn<OScalar<T>, FnSum>::Type_t
sum(const QDPExpr<RHS,OScalar<T> >& s1)
{
  typename UnaryReturn<OScalar<T>, FnSum>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  evaluate(d,OpAssign(),s1,all);

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}



//! OScalar = sum(OLattice)  under an explicit subset
/*!
 * Allow a global sum that sums over the lattice, but returns an object
 * of the same primitive type. E.g., contract only over lattice indices
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSum>::Type_t
sum(const QDPExpr<RHS,OLattice<T> >& s1, const Subset& s)
{
  typename UnaryReturn<OLattice<T>, FnSum>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  // Must initialize to zero since we do not know if the loop will be entered
  zero_rep(d.elem());

  const int *tab = s.siteTable().slice();

#pragma omp parallel 
  {
    typename UnaryReturn<OLattice<T>, FnSum>::Type_t dthread;
    zero_rep(dthread.elem());

#pragma omp for nowait
    for(int j=0; j < s.numSiteTable(); ++j) {
      
      int i = tab[j];
      dthread.elem() += forEach(s1, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW
    }

#pragma omp critical
    {
      d.elem() += dthread.elem();
    }

  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! OScalar = sum(OLattice)
/*!
 * Allow a global sum that sums over the lattice, but returns an object
 * of the same primitive type. E.g., contract only over lattice indices
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSum>::Type_t
sum(const QDPExpr<RHS,OLattice<T> >& s1)
{
  typename UnaryReturn<OLattice<T>, FnSum>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  // Loop always entered - could unroll
  zero_rep(d.elem());
  
  const int vvol = Layout::vol();

#pragma omp parallel
  {
    typename UnaryReturn<OLattice<T>, FnSum>::Type_t	dthread;
    zero_rep(dthread.elem());
    
#pragma omp for nowait
    for(int i=0; i < vvol; ++i)  {
      dthread.elem() += forEach(s1, EvalLeaf1(i), OpCombine());
    }

#pragma omp critical
    {
      d.elem() += dthread.elem();
    }
  }

#if defined(QDP_USE_PROFILING)	 
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif
  
  return d;
}


#if 0
//! OScalar = sum(OLattice)
/*!
 * Allow a global sum that sums over the lattice, but returns an object
 * of the same primitive type. E.g., contract only over lattice indices
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSum>::Type_t
sum(const QDPExpr<RHS,OLattice<T> >& s1)
{
  typename UnaryReturn<OLattice<T>, FnSum>::Type_t	d;
  
#if defined(QDP_USE_PROFILING)	 
  static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif
  
  // Loop always entered - could unroll
  zero_rep(d.elem());
  const int nodeSites = Layout::sitesOnNode();

#pragma omp parallel
  {
    typename UnaryReturn<OLattice<T>, FnSum>::Type_t	dthread;
    zero_rep(dthread.elem());
    
#pragma omp for nowait
    for(int i=0; i < vvol; ++i) {
      d.elem() += forEach(s1, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW
    }

#pragma omp critical
    {
      d.elem() += dthread.elem();
    }
  }
	
#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}
#endif

//-----------------------------------------------------------------------------
// Multiple global sums 
//! multi1d<OScalar> dest  = sumMulti(OScalar,Set) 
/*!
 * Compute the global sum on multiple subsets specified by Set 
 *
 * This implementation is specific to a purely olattice like
 * types. The scalar input value is replicated to all the
 * slices
 */
template<class RHS, class T>
typename UnaryReturn<OScalar<T>, FnSum>::Type_t
sumMulti(const QDPExpr<RHS,OScalar<T> >& s1, const Set& ss)
{
  typename UnaryReturn<OScalar<T>, FnSumMulti>::Type_t  dest(ss.numSubsets());

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest[0], OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  // lazy - evaluate repeatedly
  for(int i=0; i < ss.numSubsets(); ++i)
    evaluate(dest[i],OpAssign(),s1,all);

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return dest;
}


//! multi1d<OScalar> dest  = sumMulti(OLattice,Set) 
/*!
 * Compute the global sum on multiple subsets specified by Set 
 *
 * This is a very simple implementation. There is no need for
 * anything fancier unless global sums are just so extraordinarily
 * slow. Otherwise, generalized sums happen so infrequently the slow
 * version is fine.
 */
  template<class RHS, class T>
  struct SumMultiOLatticeThreadArgs {
    const multi1d<int>& lat_color;
    const QDPExpr<RHS,OLattice<T> >& s;
    multi1d<typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t>& dest;
    SumMultiOLatticeThreadArgs(const multi1d<int>& lat_color_,
			       const QDPExpr<RHS,OLattice<T> >& s_,
			       multi1d<typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t>& dest_) : lat_color(lat_color_), s(s_), dest(dest_) {}

  };

  template<class RHS, class T>
  void sumMultiKernel(int lo, int hi, int my_id, SumMultiOLatticeThreadArgs<RHS,T>* a)
  {
    const multi1d<int>& lat_color = a->lat_color;
    const  QDPExpr<RHS,OLattice<T> >& s=a->s;
    multi1d<typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t>& dest=a->dest;
    for(int i=lo; i < hi; ++i) { 
      int j = lat_color[i];
      (dest[my_id])[j].elem() += forEach(s, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW
    }
  }

template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t
sumMulti(const QDPExpr<RHS,OLattice<T> >& s1, const Set& ss)
{
  typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t  dest(ss.numSubsets());

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest[0], OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  multi1d< typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t > pdest(qdpNumThreads());

  // Initialize result with zero
  for(int thread=0; thread < qdpNumThreads(); ++thread) {
    pdest[thread].resize(ss.numSubsets());

    for(int k=0; k < ss.numSubsets(); ++k) {
      zero_rep(pdest[thread][k]);
    }
  }

  // Loop over all sites and accumulate based on the coloring 
  const multi1d<int>& lat_color =  ss.latticeColoring();
  SumMultiOLatticeThreadArgs<RHS,T> args(lat_color,s1,pdest);

  const int vvol = Layout::vol();
  dispatch_to_threads(vvol, args, sumMultiKernel<RHS,T>);

  for(int k=0; k< ss.numSubsets(); ++k) { 
    dest[k] = pdest[0][k];
  }

  for(int thread=1; thread < qdpNumThreads(); thread++) { 
    for(int k=0; k< ss.numSubsets(); ++k) { 
      dest[k] += pdest[thread][k];
    }
  }
  
#if 0
  for(int i=0; i < vvol; ++i) 
  {
    int j = lat_color[i];
    dest[j].elem() += forEach(s1, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW
  }
#endif

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return dest;
}

#if 0
  // Original code
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t
sumMulti(const QDPExpr<RHS,OLattice<T> >& s1, const Set& ss)
{
  typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t  dest(ss.numSubsets());

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest[0], OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  // Initialize result with zero
  for(int k=0; k < ss.numSubsets(); ++k)
    zero_rep(dest[k]);

  // Loop over all sites and accumulate based on the coloring 
  const multi1d<int>& lat_color =  ss.latticeColoring();

  const int vvol = Layout::vol();
  for(int i=0; i < vvol; ++i) 
  {
    int j = lat_color[i];
    dest[j].elem() += forEach(s1, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return dest;
}
#endif

//-----------------------------------------------------------------------------
// Multiple global sums 
//! multi2d<OScalar> dest  = sumMulti(multi1d<OScalar>,Set) 
/*!
 * Compute the global sum on multiple subsets specified by Set 
 *
 * This implementation is specific to a purely olattice like
 * types. The scalar input value is replicated to all the
 * slices
 */
template<class T>
multi2d<typename UnaryReturn<OScalar<T>, FnSum>::Type_t>
sumMulti(const multi1d< OScalar<T> >& s1, const Set& ss)
{
  multi2d<typename UnaryReturn<OScalar<T>, FnSum>::Type_t>  dest(s1.size(),ss.numSubsets());

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest(0,0), OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  // lazy - evaluate repeatedly
  for(int i=0; i < dest.size1(); ++i)
    for(int j=0; j < dest.size2(); ++j)
      dest(j,i) = s1[j];

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return dest;
}


//! multi2d<OScalar> dest  = sumMulti(multi1d<OLattice>,Set) 
/*!
 * Compute the global sum on multiple subsets specified by Set 
 *
 * This is a very simple implementation. There is no need for
 * anything fancier unless global sums are just so extraordinarily
 * slow. Otherwise, generalized sums happen so infrequently the slow
 * version is fine.
 */
template<class T>
multi2d<typename UnaryReturn<OLattice<T>, FnSum>::Type_t>
sumMulti(const multi1d< OLattice<T> >& s1, const Set& ss)
{
  multi2d<typename UnaryReturn<OLattice<T>, FnSum>::Type_t>  dest(s1.size(),ss.numSubsets());

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(dest(0,0), OpAssign(), FnSum(), s1);
  prof.time -= getClockTime();
#endif

  // Initialize result with zero
  for(int i=0; i < dest.size1(); ++i)
    for(int j=0; j < dest.size2(); ++j)
      zero_rep(dest(j,i));

  // Loop over all sites and accumulate based on the coloring 
  const multi1d<int>& lat_color =  ss.latticeColoring();

  const int vvol = Layout::vol();
  for(int k=0; k < s1.size(); ++k)
  {
    const OLattice<T>& ss1 = s1[k];

    for(int i=0; i < vvol; ++i) 
    {
      int j = lat_color[i];
      dest(k,j).elem() += ss1.elem(i);
    }
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return dest;
}


//-----------------------------------------------------------------------------
//! OScalar = norm2(trace(adj(multi1d<source>)*multi1d<source>))
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnNorm2>::Type_t
norm2(const multi1d< OScalar<T> >& s1)
{
  typename UnaryReturn<OScalar<T>, FnNorm2>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnNorm2(), s1[0]);
  prof.time -= getClockTime();
#endif

  // Possibly loop entered
  zero_rep(d.elem());

  for(int n=0; n < s1.size(); ++n)
  {
    OScalar<T>& ss1 = s1[n];
    d.elem() += localNorm2(ss1.elem());
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}

//! OScalar = sum(OScalar)  under an explicit subset
/*! Discards subset */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnNorm2>::Type_t
norm2(const multi1d< OScalar<T> >& s1, const Subset& s)
{
  return norm2(s1);
}



//! OScalar = norm2(multi1d<OLattice>) under an explicit subset
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t
norm2(const multi1d< OLattice<T> >& s1, const Subset& s)
{
  typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnNorm2(), s1[0]);
  prof.time -= getClockTime();
#endif

  // Possibly loop entered
  zero_rep(d.elem());

  const int *tab = s.siteTable().slice();

  for(int n=0; n < s1.size(); ++n) {
    
    const OLattice<T>& ss1 = s1[n];
    
#pragma omp parallel
    {
      typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t	dthread;
      zero_rep(dthread.elem());
      
#pragma omp for
      for(int j=0; j < s.numSiteTable(); ++j) {
	
	int i = tab[j];
	dthread.elem() += localNorm2(ss1.elem(i));
      }
      
#pragma omp critical
      {
	d.elem() += dthread.elem();
      }
    }
  }


#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}



//! OScalar = norm2(multi1d<OLattice>)
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t
norm2(const multi1d< OLattice<T> >& s1)
{
  return norm2(s1,all);
}



//-----------------------------------------------------------------------------
//! OScalar = innerProduct(multi1d<source1>,multi1d<source2>))
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OScalar<T1> >& s1, const multi1d< OScalar<T2> >& s2)
{
  typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProduct>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnInnerProduct(), s1[0]);
  prof.time -= getClockTime();
#endif

  // Possibly loop entered
  zero_rep(d.elem());

  for(int n=0; n < s1.size(); ++n)
  {
    OScalar<T1>& ss1 = s1[n];
    OScalar<T2>& ss2 = s2[n];
    d.elem() += localInnerProduct(ss1.elem(),ss2.elem());
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}

//! OScalar = sum(OScalar)  under an explicit subset
/*! Discards subset */
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OScalar<T1> >& s1, const multi1d< OScalar<T2> >& s2,
	     const Subset& s)
{
  return innerProduct(s1,s2);
}


//! OScalar = innerProduct(multi1d<OLattice>,multi1d<OLattice>) under an explicit subset
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OLattice<T1> >& s1, const multi1d< OLattice<T2> >& s2,
	     const Subset& s)
{
  typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProduct>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnInnerProduct(), s1[0]);
  prof.time -= getClockTime();
#endif

  // Possibly loop entered
  zero_rep(d.elem());

  const int *tab = s.siteTable().slice();

  for(int n=0; n < s1.size(); ++n) {
    const OLattice<T1>& ss1 = s1[n];
    const OLattice<T2>& ss2 = s2[n];
		
#pragma omp parallel
    {
      typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProduct>::Type_t	 dthread;
      zero_rep(dthread.elem());
      
#pragma omp for
      for(int j=0; j < s.numSiteTable(); ++j) 
	{
	  int i = tab[j];
	  dthread.elem() += localInnerProduct(ss1.elem(i),ss2.elem(i));
	}
#pragma omp critical
      {
	d.elem() += dthread.elem();
      }
    }
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! OScalar = innerProduct(multi1d<OLattice>,multi1d<OLattice>)
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OLattice<T1> >& s1, const multi1d< OLattice<T2> >& s2)
{
  return innerProduct(s1,s2,all);
}



//-----------------------------------------------------------------------------
//! OScalar = innerProductReal(multi1d<source1>,multi1d<source2>))
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OScalar<T1> >& s1, const multi1d< OScalar<T2> >& s2)
{
  typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProductReal>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnInnerProductReal(), s1[0]);
  prof.time -= getClockTime();
#endif

  // Possibly loop entered
  zero_rep(d.elem());

  for(int n=0; n < s1.size(); ++n)
  {
    OScalar<T1>& ss1 = s1[n];
    OScalar<T2>& ss2 = s2[n];
    d.elem() += localInnerProductReal(ss1.elem(),ss2.elem());
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}

//! OScalar = sum(OScalar)  under an explicit subset
/*! Discards subset */
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OScalar<T1> >& s1, const multi1d< OScalar<T2> >& s2,
		 const Subset& s)
{
  return innerProductReal(s1,s2);
}



//! OScalar = innerProductReal(multi1d<OLattice>,multi1d<OLattice>) under an explicit subset
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<T1> >& s1, const multi1d< OLattice<T2> >& s2,
		 const Subset& s)
{
  typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProductReal>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnInnerProductReal(), s1[0]);
  prof.time -= getClockTime();
#endif

  // Possibly loop entered
  zero_rep(d.elem());

  const int *tab = s.siteTable().slice();

  for(int n=0; n < s1.size(); ++n) {

    const OLattice<T1>& ss1 = s1[n];
    const OLattice<T2>& ss2 = s2[n];
		
#pragma omp parallel
    {
      typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProductReal>::Type_t	 dthread;
      zero_rep(dthread.elem());
      
#pragma omp for
      for(int j=0; j < s.numSiteTable(); ++j) 
	{
	  int i = tab[j];
	  dthread.elem() += localInnerProductReal(ss1.elem(i),ss2.elem(i));
	}
#pragma omp critical
      {
	d.elem() += dthread.elem();
      }
    }
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! OScalar = innerProductReal(multi1d<OLattice>,multi1d<OLattice>)
/*!
 * return  \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<T1> >& s1, const multi1d< OLattice<T2> >& s2)
{
  return innerProductReal(s1,s2,all);
}


//-----------------------------------------------
// Global max and min
// NOTE: there are no subset version of these operations. It is very problematic
// and QMP does not support them.
//! OScalar = globalMax(OScalar)
/*!
 * Find the maximum an object has across the lattice
 */
template<class RHS, class T>
typename UnaryReturn<OScalar<T>, FnGlobalMax>::Type_t
globalMax(const QDPExpr<RHS,OScalar<T> >& s1)
{
  typename UnaryReturn<OScalar<T>, FnGlobalMax>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnGlobalMax(), s1);
  prof.time -= getClockTime();
#endif

  evaluate(d,OpAssign(),s1,all);   // since OScalar, no global max needed

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}



//! OScalar = globalMax(OLattice)
/*!
 * Find the maximum an object has across the lattice
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnGlobalMax>::Type_t
globalMax(const QDPExpr<RHS,OLattice<T> >& s1)
{
  typename UnaryReturn<OLattice<T>, FnGlobalMax>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnGlobalMax(), s1);
  prof.time -= getClockTime();
#endif

  // Loop always entered so unroll
  d.elem() = forEach(s1, EvalLeaf1(0), OpCombine());   // SINGLE NODE VERSION FOR NOW

  const int vvol = Layout::vol();
  for(int i=1; i < vvol; ++i) 
  {
    typename UnaryReturn<T, FnGlobalMax>::Type_t  dd = 
      forEach(s1, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW

    if (toBool(dd > d.elem()))
      d.elem() = dd;
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! OScalar = globalMin(OScalar)
/*!
 * Find the minimum an object has across the lattice
 */
template<class RHS, class T>
typename UnaryReturn<OScalar<T>, FnGlobalMin>::Type_t
globalMin(const QDPExpr<RHS,OScalar<T> >& s1)
{
  typename UnaryReturn<OScalar<T>, FnGlobalMin>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnGlobalMin(), s1);
  prof.time -= getClockTime();
#endif

  evaluate(d,OpAssign(),s1,all);   // since OScalar, no global min needed

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}



//! OScalar = globalMin(OLattice)
/*!
 * Find the minimum of an object under a subset of the lattice
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnGlobalMin>::Type_t
globalMin(const QDPExpr<RHS,OLattice<T> >& s1)
{
  typename UnaryReturn<OLattice<T>, FnGlobalMin>::Type_t  d;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnGlobalMin(), s1);
  prof.time -= getClockTime();
#endif

  // Loop always entered so unroll
  d.elem() = forEach(s1, EvalLeaf1(0), OpCombine());   // SINGLE NODE VERSION FOR NOW

  const int vvol = Layout::vol();
  for(int i=1; i < vvol; ++i) 
  {
    typename UnaryReturn<T, FnGlobalMin>::Type_t  dd = 
      forEach(s1, EvalLeaf1(i), OpCombine());   // SINGLE NODE VERSION FOR NOW

    if (toBool(dd < d.elem()))
      d.elem() = dd;
  }

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//-----------------------------------------------
// Test badness/goodness of floating point numbers.
// These functions always return bool
//! bool = isnan(OScalar)
/*!
 * Return true if there is a NaN anywhere in the source
 */
template<class T>
inline bool
isnan(const QDPType<T,OScalar<T> >& s1)
{
  return isnan(s1.elem());
}


//! bool = isnan(OLattice)
/*!
 * Return true if there is a NaN anywhere in the source
 */
template<class T>
inline bool
isnan(const OLattice<T>& s1)
{
  bool d = false;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnIsNan(), s1);
  prof.time -= getClockTime();
#endif

  const int vvol = Layout::vol();
  for(int i=0; i < vvol; ++i) 
  {
    d |= isnan(s1.elem(i));
  }

  QDPInternal::globalOr(d);

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! bool = isinf(OScalar)
/*!
 * Return true if there is a NaN anywhere in the source
 */
template<class T>
inline bool
isinf(const QDPType<T,OScalar<T> >& s1)
{
  return isinf(s1.elem());
}


//! bool = isinf(OLattice)
/*!
 * Return true if there is an Inf anywhere in the source
 */
template<class T>
inline bool
isinf(const OLattice<T>& s1)
{
  bool d = false;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnIsInf(), s1);
  prof.time -= getClockTime();
#endif

  const int vvol = Layout::vol();
  for(int i=0; i < vvol; ++i) 
  {
    d |= isinf(s1.elem(i));
  }

  QDPInternal::globalOr(d);

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! bool = isfinite(OScalar)
/*!
 * Return true if all the values in source are finite floating point numbers
 */
template<class T>
inline bool
isfinite(const QDPType<T,OScalar<T> >& s1)
{
  return isfinite(s1.elem());
}


//! bool = isfinite(OLattice)
/*!
 * Return true if all the values in source are finite floating point numbers
 */
template<class T>
inline bool
isfinite(const OLattice<T>& s1)
{
  bool d = true;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnIsFinite(), s1);
  prof.time -= getClockTime();
#endif

  const int vvol = Layout::vol();
  for(int i=0; i < vvol; ++i) 
  {
    d &= isfinite(s1.elem(i));
  }

  QDPInternal::globalAnd(d);

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//! bool = isnormal(OScalar)
/*!
 * Return true if all the values in source are normal floating point numbers
 */
template<class T>
inline bool
isnormal(const QDPType<T,OScalar<T> >& s1)
{
  return isnormal(s1.elem());
}


//! bool = isnormal(OLattice)
/*!
 * Return true if all the values in source are normal floating point numbers
 */
template<class T>
inline bool
isnormal(const OLattice<T>& s1)
{
  bool d = true;

#if defined(QDP_USE_PROFILING)   
  static QDPProfile_t prof(d, OpAssign(), FnIsNormal(), s1);
  prof.time -= getClockTime();
#endif

  const int vvol = Layout::vol();
  for(int i=0; i < vvol; ++i) 
  {
    d &= isnormal(s1.elem(i));
  }

  QDPInternal::globalAnd(d);

#if defined(QDP_USE_PROFILING)   
  prof.time += getClockTime();
  prof.count++;
  prof.print();
#endif

  return d;
}


//-----------------------------------------------------------------------------
// Peek and poke at individual sites. This is very architecture specific
// NOTE: these two routines assume there is no underlying inner grid

//! Extract site element
/*! @ingroup group1
  @param l  source to examine
  @param coord Nd lattice coordinates to examine
  @return single site object of the same primitive type
  @ingroup group1
  @relates QDPType */
template<class T1>
inline OScalar<T1>
peekSite(const OScalar<T1>& l, const multi1d<int>& coord)
{
  return l;
}

//! Extract site element
/*! @ingroup group1
  @param l  source to examine
  @param coord Nd lattice coordinates to examine
  @return single site object of the same primitive type
  @ingroup group1
  @relates QDPType */
template<class RHS, class T1>
inline OScalar<T1>
peekSite(const QDPExpr<RHS,OScalar<T1> > & l, const multi1d<int>& coord)
{
  // For now, simply evaluate the expression and then call the function
  typedef OScalar<T1> C1;
  
  return peekSite(C1(l), coord);
}


//! Extract site element
/*! @ingroup group1
  @param l  source to examine
  @param coord Nd lattice coordinates to examine
  @return single site object of the same primitive type
  @ingroup group1
  @relates QDPType */
template<class T1>
inline OScalar<T1>
peekSite(const OLattice<T1>& l, const multi1d<int>& coord)
{
  OScalar<T1> dest;

  dest.elem() = l.elem(Layout::linearSiteIndex(coord));
  return dest;
}

//! Extract site element
/*! @ingroup group1
  @param l  source to examine
  @param coord Nd lattice coordinates to examine
  @return single site object of the same primitive type
  @ingroup group1
  @relates QDPType */
template<class RHS, class T1>
inline OScalar<T1>
peekSite(const QDPExpr<RHS,OLattice<T1> > & l, const multi1d<int>& coord)
{
  // For now, simply evaluate the expression and then call the function
  typedef OLattice<T1> C1;
  
  return peekSite(C1(l), coord);
}


//! Insert site element
/*! @ingroup group1
  @param l  target to update
  @param r  source to insert
  @param coord Nd lattice coordinates where to insert
  @return object of the same primitive type but of promoted lattice type
  @ingroup group1
  @relates QDPType */
template<class T1>
inline OLattice<T1>&
pokeSite(OLattice<T1>& l, const OScalar<T1>& r, const multi1d<int>& coord)
{
  l.elem(Layout::linearSiteIndex(coord)) = r.elem();
  return l;
}


//! Copy data values from field src to array dest
/*! @ingroup group1
  @param dest  target to update
  @param src   QDP source to insert
  @param s     subset
  @ingroup group1
  @relates QDPType */
template<class T>
inline void 
QDP_extract(multi1d<OScalar<T> >& dest, const OLattice<T>& src, const Subset& s)
{
  const int *tab = s.siteTable().slice();

#pragma omp parallel for
  for(int j=0; j < s.numSiteTable(); ++j) 
  {
    int i = tab[j];
    dest[i].elem() = src.elem(i);
  }
}

//! Inserts data values from site array src.
/*! @ingroup group1
  @param dest  QDP target to update
  @param src   source to insert
  @param s     subset
  @ingroup group1
  @relates QDPType */
template<class T>
inline void 
QDP_insert(OLattice<T>& dest, const multi1d<OScalar<T> >& src, const Subset& s)
{
  const int *tab = s.siteTable().slice();

#pragma omp parallel for 
  for(int j=0; j < s.numSiteTable(); ++j) 
  {
    int i = tab[j];
    dest.elem(i) = src[i].elem();
  }
}


//-----------------------------------------------------------------------------
// This is the PETE version of a map, namely return an expression
struct FnMap
{
  PETE_EMPTY_CONSTRUCTORS(FnMap)

  const int *goff;
  FnMap(const int *goffsets): goff(goffsets)
    {
//    fprintf(stderr,"FnMap(): goff=0x%x\n",goff);
    }
  
  template<class T>
  inline typename UnaryReturn<T, FnMap>::Type_t
  operator()(const T &a) const
  {
    return (a);
  }
};


#if defined(QDP_USE_PROFILING)   
template <>
struct TagVisitor<FnMap, PrintTag> : public ParenPrinter<FnMap>
{ 
  static void visit(FnMap op, PrintTag t) 
    { t.os_m << "shift"; }
};
#endif


// Specialization of ForEach deals with maps. 
template<class A, class CTag>
struct ForEach<UnaryNode<FnMap, A>, EvalLeaf1, CTag>
{
  typedef typename ForEach<A, EvalLeaf1, CTag>::Type_t TypeA_t;
  typedef typename Combine1<TypeA_t, FnMap, CTag>::Type_t Type_t;
  inline static
  Type_t apply(const UnaryNode<FnMap, A> &expr, const EvalLeaf1 &f, 
    const CTag &c) 
  {
    EvalLeaf1 ff(expr.operation().goff[f.val1()]);
//  fprintf(stderr,"ForEach<Unary<FnMap>>: site = %d, new = %d\n",f.val1(),ff.val1());

    return Combine1<TypeA_t, FnMap, CTag>::
      combine(ForEach<A, EvalLeaf1, CTag>::apply(expr.child(), ff, c),
              expr.operation(), c);
  }
};



//-----------------------------------------------------------------------------
//! General permutation map class for communications
class Map
{
public:
  //! Constructor - does nothing really
  Map() {}

  //! Destructor
  ~Map() {}

  //! Constructor from a function object
  Map(const MapFunc& fn) {make(fn);}

  //! Actual constructor from a function object
  /*! The semantics are   source_site = func(dest_site) */
  void make(const MapFunc& func);

  //! Function call operator for a shift
  /*! 
   * map(source)
   *
   * Implements:  dest(x) = s1(x+offsets)
   *
   * Shifts on a OLattice are non-trivial.
   * Notice, there may be an ILattice underneath which requires shift args.
   * This routine is very architecture dependent.
   */
  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPType<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPType<T1,C1> & l)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPType<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(goffsets.slice()),
	CreateLeaf<QDPType<T1,C1> >::make(l)));
    }


  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPExpr<T1,C1> & l)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(goffsets.slice()),
	CreateLeaf<QDPExpr<T1,C1> >::make(l)));
    }


public:
  //! Accessor to offsets
  const multi1d<int>& Offsets() const {return goffsets;}

private:
  //! Hide copy constructor
  Map(const Map&) {}

  //! Hide operator=
  void operator=(const Map&) {}

private:
  //! Offset table used for communications. 
  multi1d<int> goffsets;
};


//-----------------------------------------------------------------------------
//! Array of general permutation map class for communications
class ArrayMap
{
public:
  //! Constructor - does nothing really
  ArrayMap() {}

  //! Destructor
  ~ArrayMap() {}

  //! Constructor from a function object
  ArrayMap(const ArrayMapFunc& fn) {make(fn);}

  //! Actual constructor from a function object
  /*! The semantics are   source_site = func(dest_site,sign) */
  void make(const ArrayMapFunc& func);

  //! Function call operator for a shift
  /*! 
   * map(source,dir)
   *
   * Implements:  dest(x) = source(map(x,dir))
   *
   * Shifts on a OLattice are non-trivial.
   * Notice, there may be an ILattice underneath which requires shift args.
   * This routine is very architecture dependent.
   */
  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPType<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPType<T1,C1> & l, int dir)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPType<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(mapsa[dir].Offsets().slice()),
	CreateLeaf<QDPType<T1,C1> >::make(l)));
    }


  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPExpr<T1,C1> & l, int dir)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(mapsa[dir].Offsets().slice()),
	CreateLeaf<QDPExpr<T1,C1> >::make(l)));
    }


private:
  //! Hide copy constructor
  ArrayMap(const ArrayMap&) {}

  //! Hide operator=
  void operator=(const ArrayMap&) {}

private:
  multi1d<Map> mapsa;
  
};


//-----------------------------------------------------------------------------
//! Bi-directional version of general permutation map class for communications
class BiDirectionalMap
{
public:
  //! Constructor - does nothing really
  BiDirectionalMap() {}

  //! Destructor
  ~BiDirectionalMap() {}

  //! Constructor from a function object
  BiDirectionalMap(const MapFunc& fn) {make(fn);}

  //! Actual constructor from a function object
  /*! The semantics are   source_site = func(dest_site,sign) */
  void make(const MapFunc& func);

  //! Function call operator for a shift
  /*! 
   * map(source,isign)
   *
   * Implements:  dest(x) = source(map(x,isign))
   *
   * Shifts on a OLattice are non-trivial.
   * Notice, there may be an ILattice underneath which requires shift args.
   * This routine is very architecture dependent.
   */
  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPType<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPType<T1,C1> & l, int isign)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPType<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(bimaps[(isign+1)>>1].Offsets().slice()),
	CreateLeaf<QDPType<T1,C1> >::make(l)));
    }


  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPExpr<T1,C1> & l, int isign)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(bimaps[(isign+1)>>1].Offsets().slice()),
	CreateLeaf<QDPExpr<T1,C1> >::make(l)));
    }


private:
  //! Hide copy constructor
  BiDirectionalMap(const BiDirectionalMap&) {}

  //! Hide operator=
  void operator=(const BiDirectionalMap&) {}

private:
  multi1d<Map> bimaps;
  
};


//-----------------------------------------------------------------------------
//! Bi-directional version of general permutation map class for communications
class ArrayBiDirectionalMap
{
public:
  //! Constructor - does nothing really
  ArrayBiDirectionalMap() {}

  //! Destructor
  ~ArrayBiDirectionalMap() {}

  //! Constructor from a function object
  ArrayBiDirectionalMap(const ArrayMapFunc& fn) {make(fn);}

  //! Actual constructor from a function object
  /*! The semantics are   source_site = func(dest_site,sign) */
  void make(const ArrayMapFunc& func);

  //! Function call operator for a shift
  /*! 
   * map(source,isign,dir)
   *
   * Implements:  dest(x) = source(map(x,isign,dir))
   *
   * Shifts on a OLattice are non-trivial.
   * Notice, there may be an ILattice underneath which requires shift args.
   * This routine is very architecture dependent.
   */
  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPType<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPType<T1,C1> & l, int isign, int dir)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPType<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(bimapsa((isign+1)>>1,dir).Offsets().slice()),
	CreateLeaf<QDPType<T1,C1> >::make(l)));
    }


  template<class T1,class C1>
  inline typename MakeReturn<UnaryNode<FnMap,
    typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t>, C1>::Expression_t
  operator()(const QDPExpr<T1,C1> & l, int isign, int dir)
    {
      typedef UnaryNode<FnMap,
	typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t> Tree_t;
      return MakeReturn<Tree_t,C1>::make(Tree_t(FnMap(bimapsa((isign+1)>>1,dir).Offsets().slice()),
	CreateLeaf<QDPExpr<T1,C1> >::make(l)));
    }


private:
  //! Hide copy constructor
  ArrayBiDirectionalMap(const ArrayBiDirectionalMap&) {}

  //! Hide operator=
  void operator=(const ArrayBiDirectionalMap&) {}

private:
  multi2d<Map> bimapsa;
  
};


//-----------------------------------------------------------------------------
// Input and output of various flavors that are architecture specific

//! Decompose a lexicographic site into coordinates
multi1d<int> crtesn(int ipos, const multi1d<int>& latt_size);

#ifdef QDP_USE_LIBXML2
//! XML output
template<class T>  
XMLWriter& operator<<(XMLWriter& xml, const OLattice<T>& d)
{
  xml.openTag("OLattice");

  XMLWriterAPI::AttributeList alist;

  const int vvol = Layout::vol();
  for(int site=0; site < vvol; ++site) 
  { 
    multi1d<int> coord = crtesn(site, Layout::lattSize());
    std::ostringstream os;
    os << coord[0];
    for(int i=1; i < coord.size(); ++i)
      os << " " << coord[i];

    alist.clear();
    alist.push_back(XMLWriterAPI::Attribute("site", site));
    alist.push_back(XMLWriterAPI::Attribute("coord", os.str()));

    xml.openTag("elem", alist);
    xml << d.elem(Layout::linearSiteIndex(site));
    xml.closeTag();
  }

  xml.closeTag(); // OLattice

  return xml;
}
#endif


//! Binary output
/*! Assumes no inner grid */
template<class T>
inline
void write(BinaryWriter& bin, const OScalar<T>& d)
{
  if (Layout::primaryNode()) 
    bin.writeArray((const char *)&(d.elem()), 
		   sizeof(typename WordType<T>::Type_t), 
		   sizeof(T) / sizeof(typename WordType<T>::Type_t));
}

//! Binary output
/*! Assumes no inner grid */
template<class T>  
void write(BinaryWriter& bin, const OLattice<T>& d)
{
  const int vvol = Layout::vol();
  for(int site=0; site < vvol; ++site) 
  {
    int i = Layout::linearSiteIndex(site);
    bin.writeArray((const char*)&(d.elem(i)), 
		   sizeof(typename WordType<T>::Type_t), 
		   sizeof(T) / sizeof(typename WordType<T>::Type_t));
  }
}

//! Binary output
/*! Assumes no inner grid */
template<class T>  
void write(BinaryWriter& bin, OSubLattice<T> dd)
{
  // Single node code
  const Subset& sub = dd.subset();
  const Set& set    = sub.getSet();

  const OLattice<T>& d = dd.field();

  const multi1d<int>& lat_color = set.latticeColoring();
  const int color = sub.color();

  // Choose only this color within a lexicographic loop
  const int vvol = Layout::vol();
  for(int site=0; site < vvol; ++site) 
  {
    int i = Layout::linearSiteIndex(site);
    if (lat_color[i] == color)
    {
      bin.writeArray((const char*)&(d.elem(i)), 
		     sizeof(typename WordType<T>::Type_t), 
		     sizeof(T) / sizeof(typename WordType<T>::Type_t));
    }
  }
}

//! Write a single site of a lattice quantity at coord
/*! Assumes no inner grid */
template<class T>  
void write(BinaryWriter& bin, const OLattice<T>& d, const multi1d<int>& coord)
{
  int i = Layout::linearSiteIndex(coord);
  bin.writeArray((const char*)&(d.elem(i)), 
		 sizeof(typename WordType<T>::Type_t), 
		 sizeof(T) / sizeof(typename WordType<T>::Type_t));
}

//! Binary input
/*! Assumes no inner grid */
template<class T>
void read(BinaryReader& bin, OScalar<T>& d)
{
  bin.readArray((char*)&(d.elem()), 
		sizeof(typename WordType<T>::Type_t), 
		sizeof(T) / sizeof(typename WordType<T>::Type_t)); 
}

//! Binary input
/*! Assumes no inner grid */
template<class T>  
void read(BinaryReader& bin, OLattice<T>& d)
{
  const int vvol = Layout::vol();
  for(int site=0; site < vvol; ++site) 
  {
    int i = Layout::linearSiteIndex(site);
    bin.readArray((char*)&(d.elem(i)), 
		  sizeof(typename WordType<T>::Type_t), 
		  sizeof(T) / sizeof(typename WordType<T>::Type_t));
  }
}

//! Read a single site and place it at coord
/*! Assumes no inner grid */
template<class T>  
void read(BinaryReader& bin, OLattice<T>& d, const multi1d<int>& coord)
{
  int i = Layout::linearSiteIndex(coord);
  bin.readArray((char*)&(d.elem(i)), 
		sizeof(typename WordType<T>::Type_t), 
		sizeof(T) / sizeof(typename WordType<T>::Type_t));
}

//! Binary input
/*! Assumes no inner grid */
template<class T>  
void read(BinaryReader& bin, OSubLattice<T> dd)
{
  // Single node code
  const Subset& sub = dd.subset();
  const Set& set    = sub.getSet();

  OLattice<T>& d = dd.field();

  const multi1d<int>& lat_color = set.latticeColoring();
  const int color = sub.color();

  // Choose only this color within a lexicographic loop
  const int vvol = Layout::vol();
  for(int site=0; site < vvol; ++site) 
  {
    int i = Layout::linearSiteIndex(site);
    if (lat_color[i] == color)
    {
      bin.readArray((char*)&(d.elem(i)), 
		    sizeof(typename WordType<T>::Type_t), 
		    sizeof(T) / sizeof(typename WordType<T>::Type_t));
    }
  }
}

// **************************************************************
// Special support for slices of a lattice
namespace LatticeTimeSliceIO 
{
  template<class T>
  void readSlice(BinaryReader& bin, OLattice<T>& data, 
		 int start_lexico, int stop_lexico)
  {
    for(int site=start_lexico; site < stop_lexico; ++site)
    {
      int i = Layout::linearSiteIndex(site);
      bin.readArray((char*)&(data.elem(i)), 
		    sizeof(typename WordType<T>::Type_t), 
		    sizeof(T) / sizeof(typename WordType<T>::Type_t));
    }
  }

  template<class T>
  void writeSlice(BinaryWriter& bin, OLattice<T>& data, 
		  int start_lexico, int stop_lexico)
  {
    for(int site=start_lexico; site < stop_lexico; ++site)
    {
      int i = Layout::linearSiteIndex(site);
      bin.writeArray((const char*)&(data.elem(i)), 
		     sizeof(typename WordType<T>::Type_t), 
		     sizeof(T) / sizeof(typename WordType<T>::Type_t));
    }
  }

} // namespace LatticeTimeSliceIO
} // namespace QDP

#endif
