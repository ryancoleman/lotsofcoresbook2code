// -*- C++ -*-

/*! @file
 * @brief Outer lattice routines specific to a parallel platform with scalar layout
 */

#ifndef QDP_PARSCALAR_SPECIFIC_H
#define QDP_PARSCALAR_SPECIFIC_H

#include "qmp.h"

namespace QDP {


// Use separate defs here. This will cause subroutine calls under g++

//-----------------------------------------------------------------------------
// Layout stuff specific to a parallel architecture
namespace Layout
{
	//! coord[mu]	 <- mu	: fill with lattice coord in mu direction
	LatticeInteger latticeCoordinate(int mu);
}


//-----------------------------------------------------------------------------
// Internal ops with ties to QMP
namespace QDPInternal
{
	//! Route to another node (blocking)
	void route(void *send_buf, int srce_node, int dest_node, int count);

	//! Wait on send-receive
	void wait(int dir);

	//! Send to another node (wait)
	/*! All nodes participate */
	void sendToWait(void *send_buf, int dest_node, int count);

	//! Receive from another node (wait)
	void recvFromWait(void *recv_buf, int srce_node, int count);

	//! Via some mechanism, get the dest to node 0
	/*! Ultimately, I do not want to use point-to-point */
	template<class T>
	void sendToPrimaryNode(T& dest, int srcnode)
	{
		if (srcnode != 0)
		{
			if (Layout::primaryNode())
	recvFromWait((void *)&dest, srcnode, sizeof(T));

			if (Layout::nodeNumber() == srcnode)
	sendToWait((void *)&dest, 0, sizeof(T));
		}
	}

	//! Unsigned accumulate
	inline void sumAnUnsigned(void* inout, void* in)
	{
		*(unsigned int*)inout += *(unsigned int*)in;
	}

	//! Wrapper to get a functional unsigned global sum
	inline void globalSumArray(unsigned int *dest, int len)
	{
		for(int i=0; i < len; i++, dest++)
			QMP_binary_reduction(dest, sizeof(unsigned int), sumAnUnsigned);
	}

	//! Low level hook to QMP_global_sum
	inline void globalSumArray(int *dest, int len)
	{
		for(int i=0; i < len; i++, dest++)
			QMP_sum_int(dest);
	}

	//! Low level hook to QMP_global_sum
	inline void globalSumArray(float *dest, int len)
	{
		QMP_sum_float_array(dest, len);
	}

	//! Low level hook to QMP_global_sum
	inline void globalSumArray(double *dest, int len)
	{
		QMP_sum_double_array(dest, len);
	}

	//! Global sum on a multi1d
	template<class T>
	inline void globalSumArray(multi1d<T>& dest)
	{
		// The implementation here is relying on the structure being packed
		// tightly in memory - no padding
		typedef typename WordType<T>::Type_t	W;	 // find the machine word type

#if 0
		QDPIO::cout << "sizeof(T) = " << sizeof(T) << endl;
		QDPIO::cout << "sizeof(W) = " << sizeof(W) << endl;
		QDPIO::cout << "Calling multi1d global sum array with length " << dest.size()*sizeof(T)/sizeof(W) << endl;
#endif
		globalSumArray((W *)dest.slice(), dest.size()*sizeof(T)/sizeof(W)); // call appropriate hook
	}

	//! Global sum on a multi2d
	template<class T>
	inline void globalSumArray(multi2d<T>& dest)
	{
		// The implementation here is relying on the structure being packed
		// tightly in memory - no padding
		typedef typename WordType<T>::Type_t	W;	 // find the machine word type

#if 0
		QDPIO::cout << "sizeof(T) = " << sizeof(T) << endl;
		QDPIO::cout << "sizeof(W) = " << sizeof(W) << endl;
		QDPIO::cout << "Calling multi2d global sum array with length " << dest.size1()*dest.size2()*sizeof(T)/sizeof(W) << endl;
#endif
		// call appropriate hook
		globalSumArray((W *)dest.slice(0), dest.size1()*dest.size2()*sizeof(T)/sizeof(W));
	}

	//! Sum across all nodes
	template<class T>
	inline void globalSum(T& dest)
	{
		// The implementation here is relying on the structure being packed
		// tightly in memory - no padding
		typedef typename WordType<T>::Type_t	W;	 // find the machine word type

#if 0 
		QDPIO::cout << "sizeof(T) = " << sizeof(T) << endl;
		QDPIO::cout << "sizeof(W) = " << sizeof(W) << endl;
		QDPIO::cout << "Calling global sum array with length " << sizeof(T)/sizeof(W) << endl;
#endif
		globalSumArray((W *)&dest, int(sizeof(T)/sizeof(W))); // call appropriate hook
	}

  template<>
  inline void globalSum(double& dest)
  {
#if 0 
    QDPIO::cout << "Using simple sum_double" << endl;
#endif
    QMP_sum_double(&dest);
  }


  //! Low level hook to QMP_max_double
  inline void globalMaxValue(float* dest)
  {
    QMP_max_float(dest);
  }

  //! Low level hook to QMP_max_double
  inline void globalMaxValue(double* dest)
  {
    QMP_max_double(dest);
  }


  //! Global max across all nodes
  template<class T>
  inline void globalMax(T& dest)
  {
    typedef typename WordType<T>::Type_t  W;   // find the machine word type

    globalMaxValue((W *)&dest);
  }


  //! Low level hook to QMP_min_float
  inline void globalMinValue(float* dest)
  {
    QMP_min_float(dest);
  }

  //! Low level hook to QMP_min_double
  inline void globalMinValue(double* dest)
  {
    QMP_min_double(dest);
  }


  //! Global min across all nodes
  template<class T>
  inline void globalMin(T& dest)
  {
    typedef typename WordType<T>::Type_t  W;   // find the machine word type

    globalMinValue((W *)&dest);
  }


  //! Global And
  inline void globalCheckAnd(void* inout, void* in)
  {
    *(unsigned int*)inout = *(unsigned int*)inout & *(unsigned int*)in;
  }

  //! Wrapper to get a functional global And
  inline void globalAnd(bool& dest)
  {
    QMP_binary_reduction(&dest, sizeof(bool), globalCheckAnd);
  }



  //! Global Or
  inline void globalCheckOr(void* inout, void* in)
  {
    *(unsigned int*)inout = *(unsigned int*)inout | *(unsigned int*)in;
  }

  //! Wrapper to get a functional global Or
  inline void globalOr(bool& dest)
  {
    QMP_binary_reduction(&dest, sizeof(bool), globalCheckOr);
  }

  //! Broadcast from primary node to all other nodes
  template<class T>
  inline void broadcast(T& dest)
  {
    QMP_broadcast((void *)&dest, sizeof(T));
  }

  //! Broadcast a string from primary node to all other nodes
  void broadcast_str(std::string& dest);

  //! Broadcast from primary node to all other nodes
  inline void broadcast(void* dest, size_t nbytes)
  {
    QMP_broadcast(dest, nbytes);
  }

  //! Broadcast a string from primary node to all other nodes
  template<>
  inline void broadcast(std::string& dest)
  {
    broadcast_str(dest);
  }

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
		u_arg(
				OLattice<T>& d_,
				const QDPExpr<RHS,OScalar<T1> >& r_,
				const Op& op_,
				const int *tab_
		) : d(d_), r(r_), op(op_), tab(tab_) {}
		
		OLattice<T>& d;
		const QDPExpr<RHS,OScalar<T1> >& r;
		const Op& op;
		const int *tab;
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
		user_arg(
				OLattice<T>& d_,
				const QDPExpr<RHS,OLattice<T1> >& r_,
				const Op& op_,
				const int *tab_ ) : d(d_), r(r_), op(op_), tab(tab_) {}

				OLattice<T>& d;
				const QDPExpr<RHS,OLattice<T1> >& r;
				const Op& op;
				const int *tab;
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
// cerr << "In evaluateSubset(olattice,oscalar)\n";

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
//		fprintf(stderr,"eval(olattice,oscalar): site %d\n",i);
//		op(dest.elem(i), forEach(rhs, ElemLeaf(), OpCombine()));
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
//	cerr << "In evaluateSubset(olattice,olattice)" << endl;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(dest, op, rhs);
	prof.time -= getClockTime();
#endif

	int numSiteTable = s.numSiteTable();

	user_arg<T,T1,Op,RHS> a(dest, rhs, op, s.siteTable().slice());

	dispatch_to_threads< user_arg<T,T1,Op,RHS> >(numSiteTable, a, evaluate_userfunc);

	////////////////////
	// Original code
	///////////////////

	// General form of loop structure
	//const int *tab = s.siteTable().slice();
	//for(int j=0; j < s.numSiteTable(); ++j) 
	//{
	//int i = tab[j];
//		fprintf(stderr,"eval(olattice,olattice): site %d\n",i);
	//op(dest.elem(i), forEach(rhs, EvalLeaf1(i), OpCombine()));
	//}

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif
}




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
void copymask(OSubLattice<T2> d, const OLattice<T1>& mask, const OLattice<T2>& s1) 
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
void copymask(OLattice<T2>& dest, const OLattice<T1>& mask, const OLattice<T2>& s1) 
{
	int nodeSites = Layout::sitesOnNode();
	
#pragma omp parallel for
	for(int i=0; i < nodeSites; ++i) 
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


//! dest	= random	
/*! This implementation is correct for no inner grid */
template<class T>
void 
random(OScalar<T>& d)
{
	Seed seed = RNG::ran_seed;
	Seed skewed_seed = RNG::ran_seed * RNG::ran_mult;

	fill_random(d.elem(), seed, skewed_seed, RNG::ran_mult);

	RNG::ran_seed = seed;	 // The seed from any site is the same as the new global seed
}


//! dest	= random		under a subset
template<class T>
void 
random(OLattice<T>& d, const Subset& s)
{
	const int *tab = s.siteTable().slice();

#pragma omp parallel
	{
		Seed seed;
		Seed skewed_seed;
#pragma omp for // need the barrier to avoid that RNG::ran_seed is changed too early!
		for(int j=0; j < s.numSiteTable(); ++j) 
		{
			int i = tab[j];
			seed = RNG::ran_seed;

			skewed_seed.elem() = RNG::ran_seed.elem() * RNG::lattice_ran_mult->elem(i);
			fill_random(d.elem(i), seed, skewed_seed, RNG::ran_mult_n);
		}

#pragma omp critical (random)
		{
			RNG::ran_seed = seed;	 // The seed from any site is the same as the new global seed
		}
	}
}



//! dest	= random	 under a subset
template<class T>
void random(OSubLattice<T> dd)
{
	OLattice<T>& d = dd.field();
	const Subset& s = dd.subset();

	random(d,s);
}


//! dest	= random	
template<class T>
void random(OLattice<T>& d)
{
	random(d,all);
}


//! dest	= gaussian	 under a subset
template<class T>
void gaussian(OLattice<T>& d, const Subset& s)
{
	OLattice<T>	 r1, r2;

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



//! dest	= gaussian	 under a subset
template<class T>
void gaussian(OSubLattice<T> dd)
{
	OLattice<T>& d = dd.field();
	const Subset& s = dd.subset();

	gaussian(d,s);
}


//! dest	= gaussian
template<class T>
void gaussian(OLattice<T>& d)
{
	gaussian(d,all);
}



//-----------------------------------------------------------------------------
// Broadcast operations
//! dest	= 0 
template<class T> 
inline
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
  /* There is remotes/origin/jacques */
//! dest	= 0 
template<class T>
void zero_rep(OSubLattice<T> dd) 
{
	OLattice<T>& d = dd.field();
	const Subset& s = dd.subset();
	
	zero_rep(d,s);
}
#endif



//! dest	= 0 
template<class T> 
void zero_rep(OLattice<T>& dest) 
{
	const int nodeSites = Layout::sitesOnNode();

#pragma omp parallel for
	for(int i=0; i < nodeSites; ++i) 
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
	typename UnaryReturn<OScalar<T>, FnSum>::Type_t	 d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
	prof.time -= getClockTime();
#endif

	evaluate(d,OpAssign(),s1,all);	 // since OScalar, no global sum needed

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
	typename UnaryReturn<OScalar<T>, FnSum>::Type_t	 d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnSum(), s1);
	prof.time -= getClockTime();
#endif

	evaluate(d,OpAssign(),s1,all);	 // since OScalar, no global sum needed

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return d;
}



//! OScalar = sum(OLattice)	 under an explicit subset
/*!
 * Allow a global sum that sums over the lattice, but returns an object
 * of the same primitive type. E.g., contract only over lattice indices
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSum>::Type_t
sum(const QDPExpr<RHS,OLattice<T> >& s1, const Subset& s)
{
	typename UnaryReturn<OLattice<T>, FnSum>::Type_t	d;

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
		for(int j=0; j < s.numSiteTable(); ++j) 
		{
			int i = tab[j];
			dthread.elem() += forEach(s1, EvalLeaf1(i), OpCombine());
		}
#pragma omp critical
		{
			d.elem() += dthread.elem();
		}
	}
	
	// Do a global sum on the result
	QDPInternal::globalSum(d);

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
		for(int i=0; i < nodeSites; ++i) 
			dthread.elem() += forEach(s1, EvalLeaf1(i), OpCombine());
			
#pragma omp critical
		{
			d.elem() += dthread.elem();
		}
	}

	// Do a global sum on the result
	QDPInternal::globalSum(d);

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return d;
}


//-----------------------------------------------------------------------------
// Multiple global sums 
//! multi1d<OScalar> dest	 = sumMulti(OScalar,Set) 
/*!
 * Compute the global sum on multiple subsets specified by Set 
 *
 * This implementation is specific to a purely olattice like
 * types. The scalar input value is replicated to all the
 * slices
 */
template<class RHS, class T>
typename UnaryReturn<OScalar<T>, FnSumMulti>::Type_t
sumMulti(const QDPExpr<RHS,OScalar<T> >& s1, const Set& ss)
{
	typename UnaryReturn<OScalar<T>, FnSumMulti>::Type_t	dest(ss.numSubsets());

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


//! multi1d<OScalar> dest	 = sumMulti(OLattice,Set) 
/*!
 * Compute the global sum on multiple subsets specified by Set 
 *
 * This is a very simple implementation. There is no need for
 * anything fancier unless global sums are just so extraordinarily
 * slow. Otherwise, generalized sums happen so infrequently the slow
 * version is fine.
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t
sumMulti(const QDPExpr<RHS,OLattice<T> >& s1, const Set& ss)
{
	typename UnaryReturn<OLattice<T>, FnSumMulti>::Type_t	 dest(ss.numSubsets());

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(dest[0], OpAssign(), FnSum(), s1);
	prof.time -= getClockTime();
#endif

	// Initialize result with zero
	for(int k=0; k < ss.numSubsets(); ++k)
		zero_rep(dest[k]);

	// Loop over all sites and accumulate based on the coloring 
	const multi1d<int>& lat_color =	 ss.latticeColoring();
	const int nodeSites = Layout::sitesOnNode();

	for(int i=0; i < nodeSites; ++i) 
	{
		int j = lat_color[i];
		dest[j].elem() += forEach(s1, EvalLeaf1(i), OpCombine());
	}

	// Do a global sum on the result
	QDPInternal::globalSumArray(dest);

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return dest;
}


//-----------------------------------------------------------------------------
// Multiple global sums on an array
//! multi2d<OScalar> dest	 = sumMulti(multi1d<OScalar>,Set) 
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
	multi2d<typename UnaryReturn<OScalar<T>, FnSumMulti>::Type_t> dest(s1.size(), ss.numSubsets());

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


//! multi2d<OScalar> dest	 = sumMulti(multi1d<OLattice>,Set) 
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
	multi2d<typename UnaryReturn<OLattice<T>, FnSum>::Type_t> dest(s1.size(), ss.numSubsets());

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(dest[0], OpAssign(), FnSum(), s1);
	prof.time -= getClockTime();
#endif

	// Initialize result with zero
	for(int i=0; i < dest.size1(); ++i)
		for(int j=0; j < dest.size2(); ++j)
			zero_rep(dest(j,i));

	// Loop over all sites and accumulate based on the coloring 
	const multi1d<int>& lat_color =	 ss.latticeColoring();

	for(int k=0; k < s1.size(); ++k)
	{
		const OLattice<T>& ss1 = s1[k];
		const int nodeSites = Layout::sitesOnNode();
		for(int i=0; i < nodeSites; ++i) 
		{
			int j = lat_color[i];
			dest(k,j).elem() += ss1.elem(i);
		}
	}

	// Do a global sum on the result
	QDPInternal::globalSumArray(dest);

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
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnNorm2>::Type_t
norm2(const multi1d< OScalar<T> >& s1)
{
	typename UnaryReturn<OScalar<T>, FnNorm2>::Type_t	 d;

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

//! OScalar = sum(OScalar)	under an explicit subset
/*! Discards subset */
template<class T>
inline typename UnaryReturn<OScalar<T>, FnNorm2>::Type_t
norm2(const multi1d< OScalar<T> >& s1, const Subset& s)
{
	return norm2(s1);
}


//! OScalar = norm2(multi1d<OLattice>) under an explicit subset
/*!
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T>
inline typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t
norm2(const multi1d< OLattice<T> >& s1, const Subset& s)
{
	typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t	d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnNorm2(), s1[0]);
	prof.time -= getClockTime();
#endif

	// Possibly loop entered
	zero_rep(d.elem());

	const int *tab = s.siteTable().slice();

	for(int n=0; n < s1.size(); ++n)
	{
		const OLattice<T>& ss1 = s1[n];

		#pragma omp parallel
		{
			typename UnaryReturn<OLattice<T>, FnNorm2>::Type_t	dthread;
			zero_rep(dthread.elem());

			#pragma omp for
			for(int j=0; j < s.numSiteTable(); ++j)
			{
				int i = tab[j];
				dthread.elem() += localNorm2(ss1.elem(i));
			}
			#pragma omp critical
			{
				d.elem() += dthread.elem();
			}
		}
	}

	// Do a global sum on the result
	QDPInternal::globalSum(d);

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return d;
}


//! OScalar = norm2(multi1d<OLattice>)
/*!
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
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
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OScalar<T1> >& s1, const multi1d< OScalar<T2> >& s2)
{
	typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProduct>::Type_t	 d;

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

//! OScalar = sum(OScalar)	under an explicit subset
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
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProduct>::Type_t
innerProduct(const multi1d< OLattice<T1> >& s1, const multi1d< OLattice<T2> >& s2,
			 const Subset& s)
{
	typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProduct>::Type_t	 d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnInnerProduct(), s1[0]);
	prof.time -= getClockTime();
#endif

	// Possibly loop entered
	zero_rep(d.elem());

	const int *tab = s.siteTable().slice();
	for(int n=0; n < s1.size(); ++n)
	{
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

	// Do a global sum on the result
	QDPInternal::globalSum(d);

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return d;
}


//! OScalar = innerProduct(multi1d<OLattice>,multi1d<OLattice>)
/*!
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
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
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OScalar<T1> >& s1, const multi1d< OScalar<T2> >& s2)
{
	typename BinaryReturn<OScalar<T1>, OScalar<T2>, FnInnerProductReal>::Type_t	 d;

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

//! OScalar = sum(OScalar)	under an explicit subset
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
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
 *
 * Sum over the lattice
 * Allow a global sum that sums over all indices
 */
template<class T1, class T2>
inline typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProductReal>::Type_t
innerProductReal(const multi1d< OLattice<T1> >& s1, const multi1d< OLattice<T2> >& s2,
		 const Subset& s)
{
	typename BinaryReturn<OLattice<T1>, OLattice<T2>, FnInnerProductReal>::Type_t	 d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnInnerProductReal(), s1[0]);
	prof.time -= getClockTime();
#endif

	// Possibly loop entered
	zero_rep(d.elem());

	const int *tab = s.siteTable().slice();
	for(int n=0; n < s1.size(); ++n)
	{
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

	// Do a global sum on the result
	QDPInternal::globalSum(d);

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return d;
}


//! OScalar = innerProductReal(multi1d<OLattice>,multi1d<OLattice>)
/*!
 * return	 \sum_{multi1d} \sum_x(trace(adj(multi1d<source>)*multi1d<source>))
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
	typename UnaryReturn<OScalar<T>, FnGlobalMax>::Type_t	 d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnGlobalMax(), s1);
	prof.time -= getClockTime();
#endif

	evaluate(d,OpAssign(),s1,all);	 // since OScalar, no global max needed

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
	typename UnaryReturn<OLattice<T>, FnGlobalMax>::Type_t	d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnGlobalMax(), s1);
	prof.time -= getClockTime();
#endif

	// Loop always entered so unroll
	d.elem() = forEach(s1, EvalLeaf1(0), OpCombine());	 // SINGLE NODE VERSION FOR NOW

	const int vvol = Layout::sitesOnNode();
	for(int i=1; i < vvol; ++i) 
	{
		typename UnaryReturn<T, FnGlobalMax>::Type_t	dd = 
			forEach(s1, EvalLeaf1(i), OpCombine());		// SINGLE NODE VERSION FOR NOW

		if (toBool(dd > d.elem()))
			d.elem() = dd;
	}

	// Do a global max on the result
	QDPInternal::globalMax(d); 

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
	typename UnaryReturn<OScalar<T>, FnGlobalMin>::Type_t	 d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnGlobalMin(), s1);
	prof.time -= getClockTime();
#endif

	evaluate(d,OpAssign(),s1,all);	 // since OScalar, no global min needed

#if defined(QDP_USE_PROFILING)	 
	prof.time += getClockTime();
	prof.count++;
	prof.print();
#endif

	return d;
}


//! OScalar = globalMin(OLattice)
/*!
 * Find the minimum an object has across the lattice
 */
template<class RHS, class T>
typename UnaryReturn<OLattice<T>, FnGlobalMin>::Type_t
globalMin(const QDPExpr<RHS,OLattice<T> >& s1)
{
	typename UnaryReturn<OLattice<T>, FnGlobalMin>::Type_t	d;

#if defined(QDP_USE_PROFILING)	 
	static QDPProfile_t prof(d, OpAssign(), FnGlobalMin(), s1);
	prof.time -= getClockTime();
#endif

	// Loop always entered so unroll
	d.elem() = forEach(s1, EvalLeaf1(0), OpCombine());	 // SINGLE NODE VERSION FOR NOW

	const int vvol = Layout::sitesOnNode();
	for(int i=1; i < vvol; ++i) 
	{
		typename UnaryReturn<T, FnGlobalMin>::Type_t	dd = 
			forEach(s1, EvalLeaf1(i), OpCombine());		// SINGLE NODE VERSION FOR NOW

		if (toBool(dd < d.elem()))
			d.elem() = dd;
	}

	// Do a global min on the result
	QDPInternal::globalMin(d); 

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

  const int nodeSites = Layout::sitesOnNode();
  for(int i=0; i < nodeSites; ++i) 
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

  const int nodeSites = Layout::sitesOnNode();
  for(int i=0; i < nodeSites; ++i) 
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

  const int nodeSites = Layout::sitesOnNode();
  for(int i=0; i < nodeSites; ++i) 
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

  const int nodeSites = Layout::sitesOnNode();
  for(int i=0; i < nodeSites; ++i) 
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
	@param l	source to examine
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
	@param l	source to examine
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
	@param l	source to examine
	@param coord Nd lattice coordinates to examine
	@return single site object of the same primitive type
	@ingroup group1
	@relates QDPType */
template<class T1>
inline OScalar<T1>
peekSite(const OLattice<T1>& l, const multi1d<int>& coord)
{
	OScalar<T1> dest;
	int nodenum = Layout::nodeNumber(coord);

	// Find the result somewhere within the machine.
	// Then we must get it to node zero so we can broadcast it
	// out to all nodes
	if (Layout::nodeNumber() == nodenum)
		dest.elem() = l.elem(Layout::linearSiteIndex(coord));
	else
		zero_rep(dest.elem());

	// Send result to primary node via some mechanism
	QDPInternal::sendToPrimaryNode(dest, nodenum);

	// Now broadcast back out to all nodes
	QDPInternal::broadcast(dest);

	return dest;
}

//! Extract site element
/*! @ingroup group1
	@param l	source to examine
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
	@param l	target to update
	@param r	source to insert
	@param coord Nd lattice coordinates where to insert
	@return object of the same primitive type but of promoted lattice type
	@ingroup group1
	@relates QDPType */
template<class T1>
inline OLattice<T1>&
pokeSite(OLattice<T1>& l, const OScalar<T1>& r, const multi1d<int>& coord)
{
	if (Layout::nodeNumber() == Layout::nodeNumber(coord))
		l.elem(Layout::linearSiteIndex(coord)) = r.elem();

	return l;
}


//! Copy data values from field src to array dest
/*! @ingroup group1
	@param dest	 target to update
	@param src	 QDP source to insert
	@param s		 subset
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
	@param dest	 QDP target to update
	@param src	 source to insert
	@param s		 subset
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
// Map
//

// Empty map
struct FnMap
{
	PETE_EMPTY_CONSTRUCTORS(FnMap)
};

#if defined(QDP_USE_PROFILING)	 
template <>
struct TagVisitor<FnMap, PrintTag> : public ParenPrinter<FnMap>
{ 
	static void visit(FnMap op, PrintTag t) 
		{ t.os_m << "shift"; }
};
#endif


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
	/*! The semantics are		source_site = func(dest_site,isign) */
	void make(const MapFunc& func);

	//! Function call operator for a shift
	/*! 
	 * map(source)
	 *
	 * Implements:	dest(x) = s1(x+offsets)
	 *
	 * Shifts on a OLattice are non-trivial.
	 *
	 * Notice, this implementation does not allow an Inner grid
	 */
	template<class T1>
	OLattice<T1>
	operator()(const OLattice<T1> & l)
	{
#if QDP_DEBUG >= 3
		QDP_info("Map()");
#endif

		OLattice<T1> d;
		const int nodeSites = Layout::sitesOnNode();

		if (offnodeP)
		{
	// Off-node communications required
#if QDP_DEBUG >= 3
			QDP_info("Map: off-node communications required");
#endif

	// Eventually these declarations should move into d - the return object
			typedef T1 * T1ptr;
			T1 **dest = new(std::nothrow) T1ptr[nodeSites];
			if( dest == 0x0 ) { 
				QDP_error_exit("Unable to new T1ptr in OLattice<T1>::operator()\n");
			}
			QMP_msgmem_t msg[2];
			QMP_msghandle_t mh_a[2], mh;

#if 0
			// This test has been moved into Map::create.
			// For now, will assume there is only 1 destination nod e 
			// and receive from only 1 node
			if (srcenodes.size() != 1)
				QDP_error_exit("Map: for now only allow 1 destination node");
	
			if (destnodes.size() != 1)
				QDP_error_exit("Map: for now only allow receives from 1 node");
#endif

			int dstnum = destnodes_num[0]*sizeof(T1);
			int srcnum = srcenodes_num[0]*sizeof(T1);

	// Try getting fast and communicable memory
			QMP_mem_t *send_buf_mem = QMP_allocate_aligned_memory(dstnum,QDP_ALIGNMENT_SIZE, 
								(QMP_MEM_COMMS|QMP_MEM_FAST) ); // packed data to send
	if( send_buf_mem == 0x0 ) { 
				send_buf_mem = QMP_allocate_aligned_memory(dstnum, QDP_ALIGNMENT_SIZE, 
									QMP_MEM_COMMS);
				if( send_buf_mem == 0x0 ) { 
					QDP_error_exit("Unable to allocate send_buf_mem\n");
				}
			}

			QMP_mem_t *recv_buf_mem = QMP_allocate_aligned_memory(srcnum,QDP_ALIGNMENT_SIZE, 
										(QMP_MEM_COMMS|QMP_MEM_FAST)); // packed receive data
			if( recv_buf_mem == 0x0 ) { 
				recv_buf_mem = QMP_allocate_aligned_memory(srcnum, QDP_ALIGNMENT_SIZE, QMP_MEM_COMMS); 
				if( recv_buf_mem == 0x0 ) { 
					QDP_error_exit("Unable to allocate recv_buf_mem\n");
				}
			}
				
			T1 *send_buf = (T1 *)QMP_get_memory_pointer(send_buf_mem);
			T1 *recv_buf = (T1 *)QMP_get_memory_pointer(recv_buf_mem);

			// Total and utter paranoia
			if ( send_buf == 0x0 ) { 
				QDP_error_exit("QMP_get_memory_pointer returned NULL pointer from non NULL QMP_mem_t (send_buf)\n");
			}

			if ( recv_buf == 0x0 ) { 
				QDP_error_exit("QMP_get_memory_pointer returned NULL pointer from non NULL QMP_mem_t (recv_buf)\n"); 
			}

			const int my_node = Layout::nodeNumber();

			// Gather the face of data to send
			// For now, use the all subset
			for(int si=0; si < soffsets.size(); ++si) 
			{
#if QDP_DEBUG >= 3
				QDP_info("Map_scatter_send(buf[%d],olattice[%d])",si,soffsets[si]);
#endif

				send_buf[si] = l.elem(soffsets[si]);
			}

			// Set the dest gather pointers
			// For now, use the all subset

		// no threading to avoid confusion with order of recv_buf (and fast anyways)
			for(int i=0, ri=0; i < nodeSites; ++i) 
			{
				if (srcnode[i] != my_node)
				{
#if QDP_DEBUG >= 3
					QDP_info("Map_gather_recv(olattice[%d],recv[%d])",i,ri);
#endif

					dest[i] = &(recv_buf[ri++]);
				}
				else
				{
#if QDP_DEBUG >= 3
					QDP_info("Map_gather_onnode(olattice[%d],olattice[%d])",i,goffsets[i]);
#endif

					dest[i] = &(const_cast<T1&>(l.elem(goffsets[i])));
				}
			}

			QMP_status_t err;

#if QDP_DEBUG >= 3
			QDP_info("Map: send = 0x%x	recv = 0x%x",send_buf,recv_buf);
			QDP_info("Map: establish send=%d recv=%d",destnodes[0],srcenodes[0]);
			{
				const multi1d<int>& me = Layout::nodeCoord();
				multi1d<int> scrd = Layout::getLogicalCoordFrom(destnodes[0]);
				multi1d<int> rcrd = Layout::getLogicalCoordFrom(srcenodes[0]);

				QDP_info("Map: establish-info		my_crds=[%d,%d,%d,%d]",me[0],me[1],me[2],me[3]);
				QDP_info("Map: establish-info send_crds=[%d,%d,%d,%d]",scrd[0],scrd[1],scrd[2],scrd[3]);
				QDP_info("Map: establish-info recv_crds=[%d,%d,%d,%d]",rcrd[0],rcrd[1],rcrd[2],rcrd[3]);
			}
#endif

			msg[0]	= QMP_declare_msgmem(recv_buf, srcnum);
			if( msg[0] == (QMP_msgmem_t)NULL ) { 
				QDP_error_exit("QMP_declare_msgmem for msg[0] failed in Map::operator()\n");
			}
			msg[1]	= QMP_declare_msgmem(send_buf, dstnum);
			if( msg[1] == (QMP_msgmem_t)NULL ) {
				QDP_error_exit("QMP_declare_msgmem for msg[1] failed in Map::operator()\n");
			}

			mh_a[0] = QMP_declare_receive_from(msg[0], srcenodes[0], 0);
			if( mh_a[0] == (QMP_msghandle_t)NULL ) { 
				QDP_error_exit("QMP_declare_receive_from for mh_a[0] failed in Map::operator()\n");
			}

			mh_a[1] = QMP_declare_send_to(msg[1], destnodes[0], 0);
			if( mh_a[1] == (QMP_msghandle_t)NULL ) {
				QDP_error_exit("QMP_declare_send_to for mh_a[1] failed in Map::operator()\n");
			}

			mh			= QMP_declare_multiple(mh_a, 2);
			if( mh == (QMP_msghandle_t)NULL ) { 
				QDP_error_exit("QMP_declare_multiple for mh failed in Map::operator()\n");
			}

#if QDP_DEBUG >= 3
			QDP_info("Map: calling start send=%d recv=%d",destnodes[0],srcenodes[0]);
#endif

			// Launch the faces
			if ((err = QMP_start(mh)) != QMP_SUCCESS)
				QDP_error_exit(QMP_error_string(err));

#if QDP_DEBUG >= 3
			QDP_info("Map: calling wait");
#endif

			// Wait on the faces
			if ((err = QMP_wait(mh)) != QMP_SUCCESS)
				QDP_error_exit(QMP_error_string(err));

#if QDP_DEBUG >= 3
			QDP_info("Map: calling free msgs");
#endif

			//	QMP_free_msghandle(mh_a[1]);
			//	QMP_free_msghandle(mh_a[0]);
			QMP_free_msghandle(mh);
			QMP_free_msgmem(msg[1]);
			QMP_free_msgmem(msg[0]);

			// Scatter the data into the destination
			// Some of the data maybe in receive buffers
			// For now, use the all subset
#pragma omp parallel for
			for(int i=0; i < nodeSites; ++i) 
			{
#if QDP_DEBUG >= 3
				QDP_info("Map_scatter(olattice[%d],olattice[0x%x])",i,dest[i]);
#endif
				d.elem(i) = *(dest[i]);
			}

			// Cleanup
			QMP_free_memory(recv_buf_mem);
			QMP_free_memory(send_buf_mem);
			delete[] dest;

#if QDP_DEBUG >= 3
			QDP_info("finished cleanup");
#endif
		}
		else 
		{
			// No off-node communications - copy on node
#if QDP_DEBUG >= 3
			QDP_info("Map: copy on node - no communications, try this");
#endif

			// For now, use the all subset
#pragma omp parallel for
			for(int i=0; i < nodeSites; ++i) 
			{
#if QDP_DEBUG >= 3
				QDP_info("Map(olattice[%d],olattice[%d])",i,goffsets[i]);
#endif
				d.elem(i) = l.elem(goffsets[i]);
			}
		}

#if QDP_DEBUG >= 3
		QDP_info("exiting Map()");
#endif

		return d;
	}


	template<class T1>
	OScalar<T1>
	operator()(const OScalar<T1> & l)
		{
			return l;
		}

	template<class RHS, class T1>
	OScalar<T1>
	operator()(const QDPExpr<RHS,OScalar<T1> > & l)
		{
			// For now, simply evaluate the expression and then do the map
			typedef OScalar<T1> C1;

//		fprintf(stderr,"map(QDPExpr<OScalar>)\n");
			OScalar<T1> d = this->operator()(C1(l));

			return d;
		}

	template<class RHS, class T1>
	OLattice<T1>
	operator()(const QDPExpr<RHS,OLattice<T1> > & l)
		{
			// For now, simply evaluate the expression and then do the map
			typedef OLattice<T1> C1;

//		fprintf(stderr,"map(QDPExpr<OLattice>)\n");
			OLattice<T1> d = this->operator()(C1(l));

			return d;
		}


public:
	//! Accessor to offsets
	const multi1d<int>& goffset() const {return goffsets;}
	const multi1d<int>& soffset() const {return soffsets;}

private:
	//! Hide copy constructor
	Map(const Map&) {}

	//! Hide operator=
	void operator=(const Map&) {}

private:
	//! Offset table used for communications. 
	/*! 
	 * The direction is in the sense of the Map or Shift functions from QDP.
	 * goffsets(position) 
	 */ 
	multi1d<int> goffsets;
	multi1d<int> soffsets;
	multi1d<int> srcnode;
	multi1d<int> dstnode;

	multi1d<int> srcenodes;
	multi1d<int> destnodes;

	multi1d<int> srcenodes_num;
	multi1d<int> destnodes_num;

	// Indicate off-node communications is needed;
	bool offnodeP;
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
	/*! The semantics are		source_site = func(dest_site,isign,dir) */
	void make(const ArrayMapFunc& func);

	//! Function call operator for a shift
	/*! 
	 * map(source,dir)
	 *
	 * Implements:	dest(x) = source(map(x,dir))
	 *
	 * Shifts on a OLattice are non-trivial.
	 * Notice, there may be an ILattice underneath which requires shift args.
	 * This routine is very architecture dependent.
	 */
	template<class T1>
	OLattice<T1>
	operator()(const OLattice<T1> & l, int dir)
		{
#if QDP_DEBUG >= 3
			QDP_info("ArrayMap(OLattice,%d)",dir);
#endif

			return mapsa[dir](l);
		}

	template<class T1>
	OScalar<T1>
	operator()(const OScalar<T1> & l, int dir)
		{
#if QDP_DEBUG >= 3
			QDP_info("ArrayMap(OScalar,%d)",dir);
#endif

			return mapsa[dir](l);
		}


	template<class RHS, class T1>
	OScalar<T1>
	operator()(const QDPExpr<RHS,OScalar<T1> > & l, int dir)
		{
//		fprintf(stderr,"ArrayMap(QDPExpr<OScalar>,%d)\n",dir);

			// For now, simply evaluate the expression and then do the map
			return mapsa[dir](l);
		}

	template<class RHS, class T1>
	OLattice<T1>
	operator()(const QDPExpr<RHS,OLattice<T1> > & l, int dir)
		{
//		fprintf(stderr,"ArrayMap(QDPExpr<OLattice>,%d)\n",dir);

			// For now, simply evaluate the expression and then do the map
			return mapsa[dir](l);
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
//! BiDirectional of general permutation map class for communications
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
	/*! The semantics are		source_site = func(dest_site,isign) */
	void make(const MapFunc& func);

	//! Function call operator for a shift
	/*! 
	 * map(source,isign)
	 *
	 * Implements:	dest(x) = source(map(x,isign))
	 *
	 * Shifts on a OLattice are non-trivial.
	 * Notice, there may be an ILattice underneath which requires shift args.
	 * This routine is very architecture dependent.
	 */
	template<class T1>
	OLattice<T1>
	operator()(const OLattice<T1> & l, int isign)
		{
#if QDP_DEBUG >= 3
			QDP_info("BiDirectionalMap(OLattice,%d)",isign);
#endif

			return bimaps[(isign+1)>>1](l);
		}


	template<class T1>
	OScalar<T1>
	operator()(const OScalar<T1> & l, int isign)
		{
#if QDP_DEBUG >= 3
			QDP_info("BiDirectionalMap(OScalar,%d)",isign);
#endif

			return bimaps[(isign+1)>>1](l);
		}


	template<class RHS, class T1>
	OScalar<T1>
	operator()(const QDPExpr<RHS,OScalar<T1> > & l, int isign)
		{
//		fprintf(stderr,"BiDirectionalMap(QDPExpr<OScalar>,%d)\n",isign);

			// For now, simply evaluate the expression and then do the map
			return bimaps[(isign+1)>>1](l);
		}

	template<class RHS, class T1>
	OLattice<T1>
	operator()(const QDPExpr<RHS,OLattice<T1> > & l, int isign)
		{
//		fprintf(stderr,"BiDirectionalMap(QDPExpr<OLattice>,%d)\n",isign);

			// For now, simply evaluate the expression and then do the map
			return bimaps[(isign+1)>>1](l);
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
//! ArrayBiDirectional of general permutation map class for communications
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
	/*! The semantics are		source_site = func(dest_site,isign,dir) */
	void make(const ArrayMapFunc& func);

	//! Function call operator for a shift
	/*! 
	 * Implements:	dest(x) = source(map(x,isign,dir))
	 *
	 * Syntax:
	 * map(source,isign,dir)
	 *
	 * isign = parity of direction (+1 or -1)
	 * dir	 = array index (could be direction in range [0,...,Nd-1])
	 *
	 * Implements:	dest(x) = s1(x+isign*dir)
	 * There are cpp macros called	FORWARD and BACKWARD that are +1,-1 resp.
	 * that are often used as arguments
	 *
	 * Shifts on a OLattice are non-trivial.
	 * Notice, there may be an ILattice underneath which requires shift args.
	 * This routine is very architecture dependent.
	 */
	template<class T1>
	OLattice<T1>
	operator()(const OLattice<T1> & l, int isign, int dir)
		{
#if QDP_DEBUG >= 3
			QDP_info("ArrayBiDirectionalMap(OLattice,%d,%d)",isign,dir);
#endif

			return bimapsa((isign+1)>>1,dir)(l);
		}

	template<class T1>
	OScalar<T1>
	operator()(const OScalar<T1> & l, int isign, int dir)
		{
#if QDP_DEBUG >= 3
			QDP_info("ArrayBiDirectionalMap(OScalar,%d,%d)",isign,dir);
#endif

			return bimapsa((isign+1)>>1,dir)(l);
		}


	template<class RHS, class T1>
	OScalar<T1>
	operator()(const QDPExpr<RHS,OScalar<T1> > & l, int isign, int dir)
		{
//		fprintf(stderr,"ArrayBiDirectionalMap(QDPExpr<OScalar>,%d,%d)\n",isign,dir);

			// For now, simply evaluate the expression and then do the map
			return bimapsa((isign+1)>>1,dir)(l);
		}

	template<class RHS, class T1>
	OLattice<T1>
	operator()(const QDPExpr<RHS,OLattice<T1> > & l, int isign, int dir)
		{
//		fprintf(stderr,"ArrayBiDirectionalMap(QDPExpr<OLattice>,%d,%d)\n",isign,dir);

			// For now, simply evaluate the expression and then do the map
			return bimapsa((isign+1)>>1,dir)(l);
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

//! Binary output
/*! Assumes no inner grid */
template<class T>
inline
void write(BinaryWriter& bin, const OScalar<T>& d)
{
	bin.writeArray((const char *)&(d.elem()), 
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



// There are 2 main classes of binary/xml reader/writer methods.
// The first is a simple/portable but inefficient method of send/recv
// to/from the destination node.
// The second method (the else) is a more efficient roll-around method.
// However, this method more constrains the data layout - it must be
// close to the original lexicographic order.
// For now, use the direct send method

//! Decompose a lexicographic site into coordinates
multi1d<int> crtesn(int ipos, const multi1d<int>& latt_size);

//! XML output
template<class T>	 
XMLWriter& operator<<(XMLWriter& xml, const OLattice<T>& d)
{
	T recv_buf;

	xml.openTag("OLattice");
	XMLWriterAPI::AttributeList alist;

	// Find the location of each site and send to primary node
	for(int site=0; site < Layout::vol(); ++site)
	{
		multi1d<int> coord = crtesn(site, Layout::lattSize());

		int node	 = Layout::nodeNumber(coord);
		int linear = Layout::linearSiteIndex(coord);

		// Copy to buffer: be really careful since max(linear) could vary among nodes
		if (Layout::nodeNumber() == node)
			recv_buf = d.elem(linear);

		// Send result to primary node. Avoid sending prim-node sending to itself
		if (node != 0)
		{
#if 1
			// All nodes participate
			QDPInternal::route((void *)&recv_buf, node, 0, sizeof(T));
#else
			if (Layout::primaryNode())
	QDPInternal::recvFromWait((void *)&recv_buf, node, sizeof(T));

			if (Layout::nodeNumber() == node)
	QDPInternal::sendToWait((void *)&recv_buf, 0, sizeof(T));
#endif
		}

		if (Layout::primaryNode())
		{
			std::ostringstream os;
			os << coord[0];
			for(int i=1; i < coord.size(); ++i)
	os << " " << coord[i];

			alist.clear();
			alist.push_back(XMLWriterAPI::Attribute("site", site));
			alist.push_back(XMLWriterAPI::Attribute("coord", os.str()));

			xml.openTag("elem", alist);
			xml << recv_buf;
			xml.closeTag();
		}
	}

	xml.closeTag(); // OLattice
	return xml;
}


//! Write a lattice quantity
/*! This code assumes no inner grid */
void writeOLattice(BinaryWriter& bin, 
			 const char* output, size_t size, size_t nmemb);

//! Binary output
/*! Assumes no inner grid */
template<class T>
void write(BinaryWriter& bin, const OLattice<T>& d)
{
	writeOLattice(bin, (const char *)&(d.elem(0)), 
		sizeof(typename WordType<T>::Type_t), 
		sizeof(T) / sizeof(typename WordType<T>::Type_t));
}

//! Write a single site of a lattice quantity
/*! This code assumes no inner grid */
void writeOLattice(BinaryWriter& bin, 
			 const char* output, size_t size, size_t nmemb,
			 const multi1d<int>& coord);

//! Write a single site of a lattice quantity
/*! Assumes no inner grid */
template<class T>
void write(BinaryWriter& bin, const OLattice<T>& d, const multi1d<int>& coord)
{
	writeOLattice(bin, (const char *)&(d.elem(0)), 
		sizeof(typename WordType<T>::Type_t), 
		sizeof(T) / sizeof(typename WordType<T>::Type_t),
		coord);
}

//! Write a single site of a lattice quantity
/*! This code assumes no inner grid */
void writeOLattice(BinaryWriter& bin, 
			 const char* output, size_t size, size_t nmemb,
			 const Subset& sub);

//! Write a single site of a lattice quantity
/*! Assumes no inner grid */
template<class T>
void write(BinaryWriter& bin, OSubLattice<T> dd)
{
	const OLattice<T>& d = dd.field();

	writeOLattice(bin, (const char *)&(d.elem(0)), 
		sizeof(typename WordType<T>::Type_t), 
		sizeof(T) / sizeof(typename WordType<T>::Type_t),
		dd.subset());
}


//! Read a lattice quantity
/*! This code assumes no inner grid */
void readOLattice(BinaryReader& bin, 
			char* input, size_t size, size_t nmemb);

//! Binary input
/*! Assumes no inner grid */
template<class T>
void read(BinaryReader& bin, OLattice<T>& d)
{
	readOLattice(bin, (char *)&(d.elem(0)), 
				 sizeof(typename WordType<T>::Type_t), 
				 sizeof(T) / sizeof(typename WordType<T>::Type_t));
}

//! Read a single site of a lattice quantity
/*! This code assumes no inner grid */
void readOLattice(BinaryReader& bin, 
			char* input, size_t size, size_t nmemb, 
			const multi1d<int>& coord);

//! Read a single site of a lattice quantity
/*! Assumes no inner grid */
template<class T>
void read(BinaryReader& bin, OLattice<T>& d, const multi1d<int>& coord)
{
	readOLattice(bin, (char *)&(d.elem(0)), 
				 sizeof(typename WordType<T>::Type_t), 
				 sizeof(T) / sizeof(typename WordType<T>::Type_t),
				 coord);
}

//! Read a single site of a lattice quantity
/*! This code assumes no inner grid */
void readOLattice(BinaryReader& bin, 
			char* input, size_t size, size_t nmemb, 
			const Subset& sub);

//! Read a single site of a lattice quantity
/*! Assumes no inner grid */
template<class T>
void read(BinaryReader& bin, OSubLattice<T> d)
{
	readOLattice(bin, (char *)(d.field().getF()),
				 sizeof(typename WordType<T>::Type_t), 
				 sizeof(T) / sizeof(typename WordType<T>::Type_t),
				 d.subset());
}



// **************************************************************
// Special support for slices of a lattice
namespace LatticeTimeSliceIO 
{
	//! Lattice time slice reader
	void readOLatticeSlice(BinaryReader& bin, char* data, 
			 size_t size, size_t nmemb,
			 int start_lexico, int stop_lexico);

	void writeOLatticeSlice(BinaryWriter& bin, const char* data, 
				size_t size, size_t nmemb,
				int start_lexico, int stop_lexico);


	// Read a time slice of a lattice quantity (time must be most slowly varying)
	template<class T>
	void readSlice(BinaryReader& bin, OLattice<T>& data, 
		 int start_lexico, int stop_lexico)
	{
		readOLatticeSlice(bin, (char *)&(data.elem(0)), 
					sizeof(typename WordType<T>::Type_t), 
					sizeof(T) / sizeof(typename WordType<T>::Type_t),
					start_lexico, stop_lexico);
	}


	// Write a time slice of a lattice quantity (time must be most slowly varying)
	template<class T>
	void writeSlice(BinaryWriter& bin, const OLattice<T>& data, 
			int start_lexico, int stop_lexico)
	{
		writeOLatticeSlice(bin, (const char *)&(data.elem(0)), 
					 sizeof(typename WordType<T>::Type_t), 
					 sizeof(T) / sizeof(typename WordType<T>::Type_t),
					 start_lexico, stop_lexico);
	}

} // namespace LatticeTimeSliceIO

} // namespace QDP
#endif
