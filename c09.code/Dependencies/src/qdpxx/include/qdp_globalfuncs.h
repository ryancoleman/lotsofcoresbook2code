// -*- C++ -*-

/*! \file
 * \brief Global functions on QDPType
 */

#ifndef QDP_GLOBALFUNCS_H
#define QDP_GLOBALFUNCS_H

namespace QDP
{

  /** \defgroup group3 QDP global reductions
   *  
   *  Global reductions, like sum, norm2, etc.
   *  @{
   */
  //-----------------------------------------------
  // Global sums
  //! OScalar = sum(source)
  /*!
   * Allow a global sum that sums over the lattice, but returns an object
   * of the same primitive type. E.g., contract only over lattice indices
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnSum>::Type_t
  sum(const QDPType<T,C>& s1)
  {
    return sum(PETE_identity(s1));
  }


  //! OScalar = sum(source)  under an explicit subset
  /*!
   * Allow a global sum that sums over the lattice, but returns an object
   * of the same primitive type. E.g., contract only over lattice indices
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnSum>::Type_t
  sum(const QDPType<T,C>& s1, const Subset& s)
  {
    return sum(PETE_identity(s1),s);
  }


  //! OScalar = norm2(trace(adj(source)*source))
  /*!
   * return  num(trace(adj(source)*source))
   *
   * Sum over the lattice
   * Allow a global sum that sums over all indices
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnNorm2>::Type_t
  norm2(const QDPType<T,C>& s1)
  {
    return sum(localNorm2(s1));
  }

  template<class T, class C>
  inline typename UnaryReturn<C, FnNorm2>::Type_t
  norm2(const QDPExpr<T,C>& s1)
  {
    return sum(localNorm2(s1));
  }


  //! OScalar = norm2(trace(adj(source)*source)) under an explicit subset
  /*!
   * return  num(trace(adj(source)*source))
   *
   * Sum over the lattice
   * Allow a global sum that sums over all indices
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnNorm2>::Type_t
  norm2(const QDPType<T,C>& s1, const Subset& s)
  {
    return sum(localNorm2(s1),s);
  }

  template<class T, class C>
  inline typename UnaryReturn<C, FnNorm2>::Type_t
  norm2(const QDPExpr<T,C>& s1, const Subset& s)
  {
    return sum(localNorm2(s1),s);
  }


  //! OScalar = innerProduct(adj(source1)*source2)
  /*!
   * return  sum(trace(adj(source1)*source2))
   *
   * Sum over the lattice
   */
  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPType<T1,C1>& s1, const QDPType<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPType<T1,C1>& s1, const QDPExpr<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPExpr<T1,C1>& s1, const QDPType<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPExpr<T1,C1>& s1, const QDPExpr<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }


  //! OScalar = innerProduct(adj(source1)*source2) under an explicit subset
  /*!
   * return  sum(trace(adj(source1)*source2))
   *
   * Sum over the lattice
   */
  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPType<T1,C1>& s1, const QDPType<T2,C2>& s2,
	       const Subset& s)
  {
    return sum(localInnerProduct(s1,s2),s);
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPType<T1,C1>& s1, const QDPExpr<T2,C2>& s2,
	       const Subset& s)
  {
    return sum(localInnerProduct(s1,s2),s);
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPExpr<T1,C1>& s1, const QDPType<T2,C2>& s2,
	       const Subset& s)
  {
    return sum(localInnerProduct(s1,s2),s);
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPExpr<T1,C1>& s1, const QDPExpr<T2,C2>& s2,
	       const Subset& s)
  {
    return sum(localInnerProduct(s1,s2),s);
  }


  //! OScalar = innerProductReal(adj(source1)*source2)
  /*!
   * return  sum(trace(adj(source1)*source2))
   *
   * Sum over the lattice
   */
  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPType<T1,C1>& s1, const QDPType<T2,C2>& s2)
  {
    return sum(localInnerProductReal(s1,s2));
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPType<T1,C1>& s1, const QDPExpr<T2,C2>& s2)
  {
    return sum(localInnerProductReal(s1,s2));
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPExpr<T1,C1>& s1, const QDPType<T2,C2>& s2)
  {
    return sum(localInnerProductReal(s1,s2));
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPExpr<T1,C1>& s1, const QDPExpr<T2,C2>& s2)
  {
    return sum(localInnerProductReal(s1,s2));
  }


  //! OScalar = innerProductReal(adj(source1)*source2) under an explicit subset
  /*!
   * return  sum(trace(adj(source1)*source2))
   *
   * Sum over the lattice
   */
  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPType<T1,C1>& s1, const QDPType<T2,C2>& s2,
		   const Subset& s)
  {
    return sum(localInnerProductReal(s1,s2),s);
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPType<T1,C1>& s1, const QDPExpr<T2,C2>& s2,
		   const Subset& s)
  {
    return sum(localInnerProductReal(s1,s2),s);
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPExpr<T1,C1>& s1, const QDPType<T2,C2>& s2,
		   const Subset& s)
  {
    return sum(localInnerProductReal(s1,s2),s);
  }

  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProductReal>::Type_t
  innerProductReal(const QDPExpr<T1,C1>& s1, const QDPExpr<T2,C2>& s2,
		   const Subset& s)
  {
    return sum(localInnerProductReal(s1,s2),s);
  }


  //-----------------------------------------------------------------------------
  // Multiple global sums 
  //! dest  = sumMulti(source1,Set) 
  /*!
   * Compute the global sum on multiple subsets specified by Set 
   *
   * This is a very simple implementation. There is no need for
   * anything fancier unless global sums are just so extraordinarily
   * slow. Otherwise, generalized sums happen so infrequently the slow
   * version is fine.
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnSumMulti>::Type_t
  sumMulti(const QDPType<T,C>& s1, const Set& ss)
  {
    return sumMulti(PETE_identity(s1), ss);
  }


  //-----------------------------------------------
  // Global max and min
  //! OScalar = globalMax(source)
  /*!
   * Find the maximum of an object across the lattice
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnGlobalMax>::Type_t
  globalMax(const QDPType<T,C>& s1)
  {
    return globalMax(PETE_identity(s1));
  }


  //! OScalar = globalMin(source)
  /*!
   * Find the minimum of an object across the lattice
   */
  template<class T, class C>
  inline typename UnaryReturn<C, FnGlobalMin>::Type_t
  globalMin(const QDPType<T,C>& s1)
  {
    return globalMin(PETE_identity(s1));
  }


  //-----------------------------------------------------------------------------
  // These functions always return bool
  //! bool = isnan(source)
  /*!
   * Return true if there is a NaN anywhere in the source
   */
  template<class T, class C>
  inline bool
  isnan(const QDPExpr<T,C>& s1)
  {
    return isnan(C(s1));
  }

  //! bool = isinf(source)
  /*!
   * Return true if there is an Inf anywhere in the source
   */
  template<class T, class C>
  inline bool
  isinf(const QDPExpr<T,C>& s1)
  {
    return isinf(C(s1));
  }

  //! bool = isnormal(source)
  /*!
   * Return true if all the values in source are normal floating point numbers
   */
  template<class T, class C>
  inline bool
  isnormal(const QDPExpr<T,C>& s1)
  {
    return isnormal(C(s1));
  }

  //! bool = isfinite(source)
  /*!
   * Return true if all the values in source are finite floating point numbers
   */
  template<class T, class C>
  inline bool
  isfinite(const QDPExpr<T,C>& s1)
  {
    return isfinite(C(s1));
  }


  /** @} */ // end of group3


  /** \defgroup group5 QDP auxilliary functions
   *  
   *  SU(N) operations, spin projections, etc.
   *  @{
   */

  //-----------------------------------------------
  // Spin projection
  //! dest  = spinProject(source1) 
  /*! Boneheaded simple implementation till I get a better one... */
  template<class T, class C>
  inline typename UnaryReturn<C, FnSpinProject>::Type_t
  spinProject(const QDPType<T,C>& s1, int mu, int isign)
  {
    typedef typename UnaryReturn<C, FnSpinProject>::Type_t  Ret_t;
    Ret_t  d;

    switch (isign)
    {
    case -1:
      switch (mu)
      {
      case 0:
	d = spinProjectDir0Minus(s1);
	break;
      case 1:
	d = spinProjectDir1Minus(s1);
	break;
      case 2:
	d = spinProjectDir2Minus(s1);
	break;
      case 3:
	d = spinProjectDir3Minus(s1);
	break;
      default:
	std::cerr << "Spin project: illegal direction\n";
	exit(1);
      }
      break;

    case +1:
      switch (mu)
      {
      case 0:
	d = spinProjectDir0Plus(s1);
	break;
      case 1:
	d = spinProjectDir1Plus(s1);
	break;
      case 2:
	d = spinProjectDir2Plus(s1);
	break;
      case 3:
	d = spinProjectDir3Plus(s1);
	break;
      default:
	std::cerr << "Spin project: illegal direction\n";
	exit(1);
      }
      break;

    default:
      std::cerr << "Spin project: isign must be pos or neg.\n";
      exit(1);
    }

    return d;
  }


  //-----------------------------------------------
  // Spin reconstruction
  //! dest  = spinReconstruct(source1) 
  /*! Boneheaded simple implementation till I get a better one... */
  template<class T, class C>
  inline typename UnaryReturn<C, FnSpinReconstruct>::Type_t
  spinReconstruct(const QDPType<T,C>& s1, int mu, int isign)
  {
    //  typedef typename UnaryReturn<C, FnSpinReconstruct>::Type_t  Ret_t;
    //  Ret_t  d;

    switch (isign)
    {
    case -1:
      switch (mu)
      {
      case 0:
	return spinReconstructDir0Minus(s1);
      case 1:
	return spinReconstructDir1Minus(s1);
      case 2:
	return spinReconstructDir2Minus(s1);
      case 3:
	return spinReconstructDir3Minus(s1);
      default:
	std::cerr << "Spin reconstruct: illegal direction\n";
	exit(1);
      }
      break;

    case +1:
      switch (mu)
      {
      case 0:
	return spinReconstructDir0Plus(s1);
      case 1:
	return spinReconstructDir1Plus(s1);
      case 2:
	return spinReconstructDir2Plus(s1);
      case 3:
	return spinReconstructDir3Plus(s1);
      default:
	std::cerr << "Spin reconstruct: illegal direction\n";
	exit(1);
      }
      break;

    default:
      std::cerr << "Spin reconstruct: isign must be pos or neg.\n";
      exit(1);
    }

    //  return d;
  }


  /** @} */ // end of group5

} // namespace QDP

#endif
