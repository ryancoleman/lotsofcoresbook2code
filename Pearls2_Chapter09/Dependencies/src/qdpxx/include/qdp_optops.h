// -*- C++ -*-

/*! @file
 * @brief PETE optimized operations on QDPTypes
 *
 * These are optimizations that collapse several QDP operations into
 * one mega operation. Assembly/specialization opts do NOT go here.
 */

#ifndef QDP_OPTOPS_H
#define QDP_OPTOPS_H

namespace QDP {

//-----------------------------------------------------------------------------
// Optimization hooks
//-----------------------------------------------------------------------------

struct OpAdjMultiply
{
  PETE_EMPTY_CONSTRUCTORS(OpAdjMultiply)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpAdjMultiply >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//  cerr << "adjMultiply" << endl;
//  return (adj(a)*b);
    return adjMultiply(a,b);
  }
};

struct OpMultiplyAdj
{
  PETE_EMPTY_CONSTRUCTORS(OpMultiplyAdj)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpMultiplyAdj >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//  cerr << "multiplyAdj" << endl;
//  return (a*adj(b));
    return multiplyAdj(a,b);
  }
};

struct OpAdjMultiplyAdj
{
  PETE_EMPTY_CONSTRUCTORS(OpAdjMultiplyAdj)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, OpAdjMultiplyAdj >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//  cerr << "adjMultiplyAdj" << endl;
//  return (adj(a)*adj(b));
    return adjMultiplyAdj(a,b);
  }
};

// adjMultiply(l,r)  <-  adj(l)*r
template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<OpAdjMultiply,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPType<T2,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,OpAdjMultiply>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPType<T2,C2> & r)
{
  //  cerr << "adjMultiply(l,r) <- adj(l)*r" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // The adj does not change container type

  typedef BinaryNode<OpAdjMultiply,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPType<T2,C2> >::Leaf_t> Tree_t;
  typedef typename BinaryReturn<C1,C2,OpAdjMultiply>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPType<T2,C2> >::make(r)));
}

// adjMultiply<l,Expr>  <-  adj(l)*Expr
template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<OpAdjMultiply,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<T2,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,OpAdjMultiply>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<T2,C2> & r)
{
  //cerr << "adjMultiply(l,Expr) <- adj(l)*Expr" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // The adj does not change container type

  typedef BinaryNode<OpAdjMultiply,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<T2,C2> >::Leaf_t> Tree_t;
  typedef typename BinaryReturn<C1,C2,OpAdjMultiply>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<T2,C2> >::make(r)));
}

// multplyAdj(l,r)  <-  l*adj(r)
template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<OpMultiplyAdj,
  typename CreateLeaf<QDPType<T1,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,OpMultiplyAdj>::Type_t >::Expression_t
operator*(const QDPType<T1,C1> & l,
	  const QDPExpr<UnaryNode<FnAdjoint,T2>,C2> & r)
{
// cerr << "multiplyAdj(l,r) <- l*adj(r)" << endl;

  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // The adj does not change container type

  typedef BinaryNode<OpMultiplyAdj,
    typename CreateLeaf<QDPType<T1,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;
  typedef typename BinaryReturn<C1,C2,OpMultiplyAdj>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    CreateLeaf<QDPType<T1,C1> >::make(l),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))));
}

// multiplyAdj(Expr,r)  <-  Expr*adj(r)
template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<OpMultiplyAdj,
  typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,OpMultiplyAdj>::Type_t >::Expression_t
operator*(const QDPExpr<T1,C1> & l,
	  const QDPExpr<UnaryNode<FnAdjoint,T2>,C2> & r)
{
// cerr << "multiplyAdj(Expr,r) <- Expr*adj(r)" << endl;

  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // The adj does not change container type

  typedef BinaryNode<OpMultiplyAdj,
    typename CreateLeaf<QDPExpr<T1,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;
  typedef typename BinaryReturn<C1,C2,OpMultiplyAdj>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    CreateLeaf<QDPExpr<T1,C1> >::make(l),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))));
}

// adjMultiplyAdj(l,r)  <-  adj(l)*adj(r)
template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<OpAdjMultiplyAdj,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,OpAdjMultiplyAdj>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnAdjoint,T2>,C2> & r)
{
// cerr << "adjMultiplyAdj(l,r) <- adj(l)*adj(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // The adj does not change container type
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // The adj does not change container type

  typedef BinaryNode<OpAdjMultiplyAdj,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;
  typedef typename BinaryReturn<C1,C2,OpAdjMultiplyAdj>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))));
}


struct FnTraceMultiply
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceMultiply)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnTraceMultiply >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//    cerr << "FnTraceMultiply()" << endl;
//    return trace(a*b);
    return traceMultiply(a,b);
  }
};

struct FnTraceColorMultiply
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceColorMultiply)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnTraceColorMultiply >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//    cerr << "FnTraceColorMultiply()" << endl;
//    return traceColor(a*b);
    return traceColorMultiply(a,b);
  }
};

struct FnTraceSpinMultiply
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceSpinMultiply)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnTraceSpinMultiply >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//    cerr << "FnTraceSpinMultiply()" << endl;
//    return traceSpin(a*b);
    return traceSpinMultiply(a,b);
  }
};

struct FnTraceSpinOuterProduct
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceSpinOuterProduct)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnTraceSpinOuterProduct >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//    cerr << "FnTraceSpinOuterProduct()" << endl;
//    return traceSpin(outerProduct(a,b));
    return traceSpinOuterProduct(a,b);
  }
};

struct FnTraceSpinQuarkContract13
{
  PETE_EMPTY_CONSTRUCTORS(FnTraceSpinQuarkContract13)
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnTraceSpinQuarkContract13 >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
//    cerr << "FnTraceSpinQuarkContract13()" << endl;
//    return (traceSpin(quarkContract13(a,b)));
    return (traceSpinQuarkContract13(a,b));
  }
};


#if 1

//! traceMultiply(l,r)  <-  trace(l*r)
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnTraceMultiply,T1,T2>,
  typename UnaryReturn<CC,FnTrace>::Type_t>::Expression_t
trace(const QDPExpr<BinaryNode<OpMultiply,T1,T2>,CC> & ll)
{
//  cerr << "traceMultiply(l,r) <- trace(l*r)" << endl;

  typedef BinaryNode<FnTraceMultiply,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnTrace>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}

//! traceColorMultiply(l,r)  <-  traceColor(l*r)
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnTraceColorMultiply,T1,T2>,
  typename UnaryReturn<CC,FnTrace>::Type_t>::Expression_t
traceColor(const QDPExpr<BinaryNode<OpMultiply,T1,T2>,CC> & ll)
{
//  cerr << "traceColorMultiply(l,r) <- traceColor(l*r)" << endl;

  typedef BinaryNode<FnTraceColorMultiply,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnTrace>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}

//! traceSpinMultiply(l,r)  <-  traceSpin(l*r)
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnTraceSpinMultiply,T1,T2>,
  typename UnaryReturn<CC,FnTrace>::Type_t>::Expression_t
traceSpin(const QDPExpr<BinaryNode<OpMultiply,T1,T2>,CC> & ll)
{
//  cerr << "traceSpinMultiply(l,r) <- traceSpin(l*r)" << endl;

  typedef BinaryNode<FnTraceSpinMultiply,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnTrace>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}

//! traceOuterProduct(l,r)  <-  traceSpin(outerProduct(l,r))
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnTraceSpinOuterProduct,T1,T2>,
  typename UnaryReturn<CC,FnTrace>::Type_t>::Expression_t
traceSpin(const QDPExpr<BinaryNode<FnOuterProduct,T1,T2>,CC> & ll)
{
//  cerr << "traceSpinOuterProduct(l,r) <- traceSpin(outerProduct(l,r))" << endl;

  typedef BinaryNode<FnTraceSpinOuterProduct,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnTrace>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}

//! localInnerProduct(l,r)  <-  trace(adj(l)*r)
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnLocalInnerProduct,T1,T2>,
  typename UnaryReturn<CC,FnTrace>::Type_t>::Expression_t
trace(const QDPExpr<BinaryNode<OpAdjMultiply,T1,T2>,CC> & ll)
{
//  cerr << "localInnerProduct(l,r) <- trace(adj(l)*r)" << endl;

  typedef BinaryNode<FnLocalInnerProduct,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnTrace>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}


//! traceSpinQuarkContract13(l,r)  <-  traceSpin(quarkContract13(l,r))
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnTraceSpinQuarkContract13,T1,T2>,
  typename UnaryReturn<CC,FnTrace>::Type_t>::Expression_t
traceSpin(const QDPExpr<BinaryNode<FnQuarkContract13,T1,T2>,CC> & ll)
{
//  cerr << "traceSpinQuarkContract13(l,r) <- traceSpin(quarkContract13(l,r))" << endl;

  typedef BinaryNode<FnTraceSpinQuarkContract13,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnTrace>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}


//! realTrace(l) <- real(trace(l))
template<class T1,class CC>
inline typename MakeReturn<UnaryNode<FnRealTrace,T1>,
  typename UnaryReturn<CC,FnReal>::Type_t>::Expression_t
real(const QDPExpr<UnaryNode<FnTrace,T1>,CC> & ll)
{
//  cerr << "realTrace(ll) <- real(trace(ll))" << endl;

  typedef UnaryNode<FnRealTrace,T1> Tree_t;
  typedef typename UnaryReturn<CC,FnReal>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().child()));
}

//! localInnerProductReal(l,r)  <-  real(trace(adj(l)*r))
template<class T1,class T2,class CC>
inline typename MakeReturn<BinaryNode<FnLocalInnerProductReal,T1,T2>,
  typename UnaryReturn<CC,FnReal>::Type_t>::Expression_t
real(const QDPExpr<BinaryNode<FnLocalInnerProduct,T1,T2>,CC> & ll)
{
//  cerr << "localInnerProductReal(l,r) <- real(localInnerProduct(l,r))" << endl;

  typedef BinaryNode<FnLocalInnerProductReal,T1,T2> Tree_t;
  typedef typename UnaryReturn<CC,FnReal>::Type_t Container_t;
  return MakeReturn<Tree_t,Container_t>::make(Tree_t(
    ll.expression().left(), 
    ll.expression().right()));
}

#endif


} // namespace QDP

#endif
