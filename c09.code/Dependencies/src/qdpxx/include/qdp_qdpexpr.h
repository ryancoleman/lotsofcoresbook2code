// -*- C++ -*-

/*! @file
 * @brief Expression class for QDP
 */

#ifndef QDP_QDPEXPR_H
#define QDP_QDPEXPR_H

namespace QDP {


//! Expression class for QDP
template<class T, class C>
class QDPExpr
{
public:
  //! Type of the first argument.
  typedef T Subtype_t;

  //! Type of the expression.
  typedef T Expression_t;

  //! Type of the container class
  typedef C Container_t;

  //! Construct from an expression.
  QDPExpr(const T& expr) : expr_m(expr)
  { }

  //! Accessor that returns the expression.
  const Expression_t& expression() const
  {
    return expr_m;
  }

#if 0
  //! Conversion (evaluation) to QDPType
  operator C() const
  {
    return C(*this);
  }
#endif

private:
  /*! 
   * Store the expression by value since it is a temporary produced
   * by operator functions.
   */

  T expr_m;
};


//-----------------------------------------------------------------------------
// For Expressions, we strip the QDPExpr<> wrapper since it is intended
// to wrap the whole expression. (QDPExpr<Scalar<>>+QDPExpr<BinaryNode<>>
// becomes QDPExpr<BinaryNode<OpAdd,Scalar<>,BinaryNode<>>>)
//-----------------------------------------------------------------------------

template<class T, class C>
struct MakeReturn
{
  typedef QDPExpr<T, C>  Expression_t;
  inline static
  Expression_t make(const T &a) { return Expression_t(a); }
};

template<class T, class C>
struct CreateLeaf<QDPExpr<T,C> >
{
  typedef QDPExpr<T,C> Input_t;
  typedef typename Input_t::Expression_t Leaf_t;

  inline static
  const Leaf_t &make(const Input_t& a)
  {
    return a.expression();
  }
};


template<class T, class FTag, class CTag, class DTag>
struct ForEach<QDPExpr<T,DTag>, FTag, CTag>
{
  typedef typename ForEach<T, FTag, CTag>::Type_t Type_t;
  inline static
  Type_t apply(const QDPExpr<T,DTag>& expr, const FTag &f, 
	       const CTag &c) 
  {
    return ForEach<T, FTag, CTag>::apply(expr.expression(), f, c);
  }
};


//-----------------------------------------------------------------------------
// Alternative OpCombiner's that deal explicitly with references
//-----------------------------------------------------------------------------

template<class A,class Op>
struct Combine1<Reference<A>, Op, OpCombine>
{
  typedef typename UnaryReturn<A, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, Op op, OpCombine) { return op(a.reference()); }
};

template<class A,class B,class Op>
struct Combine2<Reference<A>, Reference<B>, Op, OpCombine>
{
  typedef typename BinaryReturn<A, B, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, Reference<B> b, Op op, OpCombine)
  {
    return op(a.reference(), b.reference());
  }
};

template<class A,class B,class Op>
struct Combine2<Reference<A>, B, Op, OpCombine>
{
  typedef typename BinaryReturn<A, B, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, B b, Op op, OpCombine)
  {
    return op(a.reference(), b);
  }
};

template<class A,class B,class Op>
struct Combine2<A, Reference<B>, Op, OpCombine>
{
  typedef typename BinaryReturn<A, B, Op>::Type_t Type_t;
  inline static
  Type_t combine(A a, Reference<B> b, Op op, OpCombine)
  {
    return op(a, b.reference());
  }
};


template<class A,class B,class C,class Op>
struct Combine3<Reference<A>, B, C, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, B b, C c, Op op, OpCombine)
  {
    return op(a.reference(), b, c);
  }
};

template<class A,class B,class C,class Op>
struct Combine3<A, Reference<B>, C, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(A a, Reference<B> b, C c, Op op, OpCombine)
  {
    return op(a, b.reference(), c);
  }
};

template<class A,class B,class C,class Op>
struct Combine3<A, B, Reference<C>, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(A a, B b, Reference<C> c, Op op, OpCombine)
  {
    return op(a, b, c.reference());
  }
};


template<class A,class B,class C,class Op>
struct Combine3<Reference<A>, Reference<B>, C, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, Reference<B> b, C c, Op op, OpCombine)
  {
    return op(a.reference(), b.reference(), c);
  }
};

template<class A,class B,class C,class Op>
struct Combine3<A, Reference<B>, Reference<C>, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(A a, Reference<B> b, Reference<C> c, Op op, OpCombine)
  {
    return op(a, b.reference(), c.reference());
  }
};

template<class A,class B,class C,class Op>
struct Combine3<Reference<A>, B, Reference<C>, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, B b, Reference<C> c, Op op, OpCombine)
  {
    return op(a.reference(), b, c.reference());
  }
};

template<class A,class B,class C,class Op>
struct Combine3<Reference<A>, Reference<B>, Reference<C>, Op, OpCombine>
{
  typedef typename TrinaryReturn<A, B, C, Op>::Type_t Type_t;
  inline static
  Type_t combine(Reference<A> a, Reference<B> b, Reference<C> c, Op op, OpCombine)
  {
    return op(a.reference(), b.reference(), c.reference());
  }
};

} // namespace QDP

#endif
