// -*- C++ -*-
// ACL:license
// ----------------------------------------------------------------------
// This software and ancillary information (herein called "SOFTWARE")
// called PETE (Portable Expression Template Engine) is
// made available under the terms described here.  The SOFTWARE has been
// approved for release with associated LA-CC Number LA-CC-99-5.
// 
// Unless otherwise indicated, this SOFTWARE has been authored by an
// employee or employees of the University of California, operator of the
// Los Alamos National Laboratory under Contract No.  W-7405-ENG-36 with
// the U.S. Department of Energy.  The U.S. Government has rights to use,
// reproduce, and distribute this SOFTWARE. The public may copy, distribute,
// prepare derivative works and publicly display this SOFTWARE without 
// charge, provided that this Notice and any statement of authorship are 
// reproduced on all copies.  Neither the Government nor the University 
// makes any warranty, express or implied, or assumes any liability or 
// responsibility for the use of this SOFTWARE.
// 
// If SOFTWARE is modified to produce derivative works, such modified
// SOFTWARE should be clearly marked, so as not to confuse it with the
// version available from LANL.
// 
// For more information about PETE, send e-mail to pete@acl.lanl.gov,
// or visit the PETE web page at http://www.acl.lanl.gov/pete/.
// ----------------------------------------------------------------------
// ACL:license

//-----------------------------------------------------------------------------
// Classes:
// ForEachInOrder
// TagVisitor
//-----------------------------------------------------------------------------

#ifndef POOMA_PETE_FOREACHINORDER_H
#define POOMA_PETE_FOREACHINORDER_H

//-----------------------------------------------------------------------------
// Overview: 
//
//   ForEachInOrder is like ForEach except that it traverses the parse
//   tree "in order", meaning it visits the parts of a TBTree as follows:
//
//         visit left child
//         visit value
//         visit right child
//
//   In addition, it calls a start() function on the visit tag before
//   visit(left) and a finish() function after visit(right). This
//   additional bit of generality allows special actions to be taken,
//   in essence, when the ForEachInOrder::apply moves down and back up
//   the edges of the parse tree (such as printing parentheses).
//
//   An "in order" traversal is not what one does for evaluating
//   expressions, so this may not be useful for much, but I wanted to
//   do it to gain some more experience with PETE.
//
//   This first cut will only do TBTrees.
//
//   The TagFunctor and TagCombine structs from ForEach.h can be reused. 
//
//   TagVisitor is a new class that visits the "value" node prior,
//   between, and after the left and right children are visited.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include "PETE/PETE.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
// The ForEachInOrder struct implements a template meta-program
// traversal of a PETE Expression parse-tree. As explained above, this 
// is done "in order" rather than the "post order" traversal done by
// ForEach. 
//
// The ForEachInOrder struct defines:
//
//   typename ForEachInOrder<Expr,FTag,VTag,CTag>::Type_t
//
//   Type_t::apply(Expr& expr, FTag f, VTag v, CTag c) {...};
//
// where
//
//   Expr is the type of an expression tree.
//   FTag is the type of a functor tag.
//   VTag is the type of a visitor tag.
//   CTag is the type of a combiner tag.
//
// Details:
//
//   Type_t::apply(Expr &e, FTag f, VTag v, CTag c) 
//
// function that traverses the expression tree defined by e, and for
// each binary-tree node it does:
//
//     TagVisitor<...>::start(e.value_m,v);
//
//     left_val = ForEachInOrder<...>::apply(e.left_m,f,v,c),
//
//     TagVisitor<...>::visit(e.value_m,v);
//
//     right_val = ForEachInOrder<...>::apply(e.right_m,f,v,c),
//
//     retval = TagCombineInOrdere<...>::
//               apply(left_val, right_val, e.value_m, c)
//
//     TagVisitor<...>::finish(e.value_m,v);
//
//     return retval;
//
// The TagFunctor is specialized to perform a specific action at the
// leaf nodes. 
//
// The TagVisitor is specialized to perform specific actions both at
// the start and finish of a new TBTree node, and to perform a
// specific operation when it visits the "value" node of the parse
// tree (i.e. this can be specialized to perform specific operations
// for every type of operator). Note that the value returned by
// TagVisitor::apply must be of the type Op. This usually means
// that TagVisitor::apply will just return e.value_m after it does
// its calculation.
//
// The TagCombiner is specialized to combine the results of visiting
// the left, right, and value fields.
//
// The type of object returned is given by: 
//
//    typename ForEachInOrder<Expr,FTag,VTag,CTag>::Type_t 
//
//-----------------------------------------------------------------------------

//
// struct TagVisitor
//
// "Visitor" functor whose apply() method is applied to the value_m
// field of an expression between the left-traversal and the
// right-traversal.
//
// Default is "null" behavior. Just return the op. This should make
// this equivalent to ForEach. This should probably always return the
// unless it is ignored by everything else. But it can take other
// actions as well.
//
// Also includes start() and finish() functions that are called when
// the traversal moves down and back up an edge, respectively.
//

template <class Op, class VTag>
struct TagVisitor 
{
  static void start(Op, VTag) { }
  static void center(Op, VTag) { }
  static void visit(Op, VTag) { }
  static void finish(Op, VTag) { }
};  


// 
// struct ForEachInOrder
//
// Template meta-program for traversing the parse tree.
//
// Default behaviour assumes you're at a leaf, in which case
// it just applies the FTag functor
//

template<class Expr, class FTag, class VTag, class CTag>
struct ForEachInOrder
{
  typedef LeafFunctor<Expr,FTag> Tag_t;
  typedef typename Tag_t::Type_t Type_t;

  static Type_t apply(const Expr &expr, const FTag &f, const VTag &v, 
		      const CTag &c) 
  {
    return Tag_t::apply(expr,f);
  }
};

//
// The Refernce case needs to apply the functor to the wrapped object.
//

template<class T, class FTag, class VTag, class CTag>
struct ForEachInOrder<Reference<T>,FTag,VTag,CTag>
{
  typedef LeafFunctor<T,FTag> Tag_t;
  typedef typename Tag_t::Type_t Type_t;

  static Type_t apply(const Reference<T> &expr, const FTag &f,
		      const VTag &v, const CTag &c) 
  {
    return Tag_t::apply(expr.reference(),f);
  }
};

//
// struct ForEachInOrder
//
// Specialization for a TBTree. This just performs the recursive
// traversal described above.
//

template<class Op, class A, class FTag, class VTag, 
  class CTag>
struct ForEachInOrder<UnaryNode<Op, A>, FTag, VTag, CTag>
{
  typedef ForEachInOrder<A, FTag, VTag, CTag> ForEachA_t;
  typedef TagVisitor<Op, VTag>          Visitor_t;

  typedef typename ForEachA_t::Type_t   TypeA_t;

  typedef Combine1<TypeA_t, Op, CTag>   Combiner_t;

  typedef typename Combiner_t::Type_t   Type_t;

  static Type_t apply(const UnaryNode<Op, A> &expr, const FTag &f, 
		      const VTag &v, const CTag &c)
  {
    Visitor_t::visit(expr.operation(),v);

    Visitor_t::start(expr.operation(),v);

    TypeA_t A_val  = ForEachA_t::apply(expr.child(), f, v, c);
    Type_t val = Combiner_t::combine(A_val, expr.operation(), c);
        
    Visitor_t::finish(expr.operation(),v);

    return val;
  }
};


/*!
 * struct ForEachInOrder for BinaryNode
 */

template<class Op, class A, class B, class FTag, class VTag, 
  class CTag>
struct ForEachInOrder<BinaryNode<Op, A, B>, FTag, VTag, CTag>
{
  typedef ForEachInOrder<A, FTag, VTag, CTag> ForEachA_t;
  typedef ForEachInOrder<B, FTag, VTag, CTag> ForEachB_t;
  typedef TagVisitor<Op, VTag>                Visitor_t;

  typedef typename ForEachA_t::Type_t  TypeA_t;
  typedef typename ForEachB_t::Type_t  TypeB_t;

  typedef Combine2<TypeA_t, TypeB_t, Op, CTag>  Combiner_t;

  typedef typename Combiner_t::Type_t Type_t;

  static Type_t apply(const BinaryNode<Op, A, B> &expr, const FTag &f, 
		      const VTag &v, const CTag &c) 
  {
    Visitor_t::visit(expr.operation(),v);

    Visitor_t::start(expr.operation(),v);

    TypeA_t left_val  = ForEachA_t::apply(expr.left(), f, v, c);

    Visitor_t::center(expr.operation(),v);

    TypeB_t right_val = ForEachB_t::apply(expr.right(), f, v, c);

    Type_t val = Combiner_t::combine(left_val, right_val, expr.operation(), c);
        
    Visitor_t::finish(expr.operation(),v);

    return val;
  }
};


/*!
 * struct ForEachInOrder for BinaryNode
 */

template<class Op, class A, class B, class C, class FTag, class VTag, 
  class CTag>
struct ForEachInOrder<TrinaryNode<Op, A, B, C>, FTag, VTag, CTag>
{
  typedef ForEachInOrder<A, FTag, VTag, CTag> ForEachA_t;
  typedef ForEachInOrder<B, FTag, VTag, CTag> ForEachB_t;
  typedef ForEachInOrder<C, FTag, VTag, CTag> ForEachC_t;
  typedef TagVisitor<Op, VTag>                Visitor_t;

  typedef typename ForEachA_t::Type_t  TypeA_t;
  typedef typename ForEachB_t::Type_t  TypeB_t;
  typedef typename ForEachC_t::Type_t  TypeC_t;

  typedef Combine3<TypeA_t, TypeB_t, TypeC_t, Op, CTag> Combiner_t;

  typedef typename Combiner_t::Type_t Type_t;

  static Type_t apply(const TrinaryNode<Op, A, B, C> &expr, const FTag &f, 
		      const VTag &v, const CTag &c) 
  {
    Visitor_t::visit(expr.operation(),v);

    Visitor_t::start(expr.operation(),v);

    TypeA_t left_val  = ForEachA_t::apply(expr.left(), f, v, c);

    Visitor_t::center(expr.operation(),v);

    TypeB_t middle_val= ForEachB_t::apply(expr.middle(), f, v, c);

    Visitor_t::center(expr.operation(),v);

    TypeC_t right_val = ForEachC_t::apply(expr.right(), f, v, c);

    Type_t val = Combiner_t::combine(left_val, middle_val, right_val, expr.operation(), c);
        
    Visitor_t::finish(expr.operation(),v);

    return val;
  }
};


//
// Specializations of Combine2 for NullTag()
// (seems like something like this should be in ForEach.h)
//

struct NullTag {};

template<class A,class B,class Op>
struct Combine2<A, B, Op, NullTag>
{
  typedef int Type_t;
  static Type_t combine(A, B, Op, NullTag)
    { return 0; }
};

#endif  // PETE_PETE_FOREACHINORDER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ForEachInOrder.h,v $   $Author: edwards $
// $Revision: 1.2 $   $Date: 2004-07-27 05:24:35 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
