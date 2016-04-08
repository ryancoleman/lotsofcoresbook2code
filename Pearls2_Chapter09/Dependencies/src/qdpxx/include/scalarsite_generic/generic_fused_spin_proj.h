#ifndef GENERIC_FUSED_SPIN_PROJ_H
#define GENERIC_FUSED_SPIN_PROJ_H


namespace QDP {

// Convenience Types
typedef PColorVector<RComplex<REAL>, 3> ColVec;
typedef PSpinVector< ColVec, 4> Spin4;
typedef PSpinVector< ColVec, 2> Spin2;
typedef PColorMatrix<RComplex<REAL>, 3> ColMat;

/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir0Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
// This is a struct for adj(x)*spinProjectDir0Plus(y)
struct FnAdjMultSprojDir0Plus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir0Plus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Plus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    // Print Diagnostics
    //    cout << "FnAdjMultSprojDir0Plus" << endl << flush;
    
    // Call the appropriate match
    return (adjMultSprojDir0Plus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir0Plus(r)
// to:    FnAdjMultSprojDir0Plus(l,r)
//
// The spinProjectDir0Plus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir0Plus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir0Plus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir0Plus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir0Plus,T2>,C2> & r)
{
  // Print Diagnostics
  // cout << "FnAdjMultSprojDir0Plus(l,r) <- adj(l)*FnSpinProjectDir0Plus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir0Plus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir0Plus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir0Plus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir0Plus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Plus >::Type_t
adjMultSprojDir0Plus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Plus >::Type_t tmp;

  tmp = spinProjectDir0Plus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir0Plus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus >::Type_t
adjMultSprojDir0Plus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Plus >::Type_t d;

  inlineSpinProjDir0Plus(&(b.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);

  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;

  return ret;
}


/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir0Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

// This is a struct for adj(x)*spinProjectDir0Minus(y)
struct FnAdjMultSprojDir0Minus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir0Minus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Minus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    // Call the appropriate match
    return (adjMultSprojDir0Minus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir0Minus(r)
// to:    FnAdjMultSprojDir0Minus(l,r)
//
// The spinProjectDir0Minus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir0Minus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir0Minus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir0Minus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir0Minus,T2>,C2> & r)
{
  // Print Diagnostics
  //  cout << "FnAdjMultSprojDir0Minus(l,r) <- adj(l)*FnSpinProjectDir0Minus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir0Minus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir0Minus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir0Minus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir0Minus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Minus >::Type_t
adjMultSprojDir0Minus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir0Minus >::Type_t tmp;

  tmp = spinProjectDir0Minus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir0Minus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus >::Type_t
adjMultSprojDir0Minus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir0Minus >::Type_t d;

  inlineSpinProjDir0Minus(&(b.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);

  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;

  return ret;
}


/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir1Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

// This is a struct for adj(x)*spinProjectDir1Plus(y)
struct FnAdjMultSprojDir1Plus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir1Plus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Plus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    // Call the appropriate match
    return (adjMultSprojDir1Plus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir1Plus(r)
// to:    FnAdjMultSprojDir1Plus(l,r)
//
// The spinProjectDir1Plus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir1Plus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir1Plus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir1Plus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir1Plus,T2>,C2> & r)
{
  // Print Diagnostics
  // cout << "FnAdjMultSprojDir1Plus(l,r) <- adj(l)*FnSpinProjectDir1Plus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir1Plus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir1Plus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir1Plus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir1Plus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Plus >::Type_t
adjMultSprojDir1Plus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Plus >::Type_t tmp;

  tmp = spinProjectDir1Plus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir1Plus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus >::Type_t
adjMultSprojDir1Plus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Plus >::Type_t d;

  inlineSpinProjDir1Plus(&(b.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);

  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;

  return ret;
}

/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir1Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */


// This is a struct for adj(x)*spinProjectDir1Minus(y)
struct FnAdjMultSprojDir1Minus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir1Minus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Minus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {

    // Call the appropriate match
    return (adjMultSprojDir1Minus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir1Minus(r)
// to:    FnAdjMultSprojDir1Minus(l,r)
//
// The spinProjectDir1Minus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir1Minus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir1Minus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir1Minus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir1Minus,T2>,C2> & r)
{
  // Print Diagnostics
  //  cout << "FnAdjMultSprojDir1Minus(l,r) <- adj(l)*FnSpinProjectDir1Minus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir1Minus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir1Minus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir1Minus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir1Minus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Minus >::Type_t
adjMultSprojDir1Minus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir1Minus >::Type_t tmp;

  tmp = spinProjectDir1Minus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir1Minus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus >::Type_t
adjMultSprojDir1Minus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir1Minus >::Type_t d;

  inlineSpinProjDir1Minus(&(b.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);

  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;
    
    return ret;
}

/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir2Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */


// This is a struct for adj(x)*spinProjectDir2Plus(y)
struct FnAdjMultSprojDir2Plus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir2Plus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Plus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    // Call the appropriate match
    return (adjMultSprojDir2Plus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir2Plus(r)
// to:    FnAdjMultSprojDir2Plus(l,r)
//
// The spinProjectDir2Plus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir2Plus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir2Plus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir2Plus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir2Plus,T2>,C2> & r)
{
  // Print Diagnostics
  // cout << "FnAdjMultSprojDir2Plus(l,r) <- adj(l)*FnSpinProjectDir2Plus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir2Plus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir2Plus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir2Plus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir2Plus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Plus >::Type_t
adjMultSprojDir2Plus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Plus >::Type_t tmp;

  tmp = spinProjectDir2Plus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir2Plus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus >::Type_t
adjMultSprojDir2Plus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Plus >::Type_t d;

  inlineSpinProjDir2Plus(&(b.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);

  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;

  return ret;
}

/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir2Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
// This is a struct for adj(x)*spinProjectDir2Minus(y)
struct FnAdjMultSprojDir2Minus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir2Minus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Minus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {

    // Call the appropriate match
    return (adjMultSprojDir2Minus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir2Minus(r)
// to:    FnAdjMultSprojDir2Minus(l,r)
//
// The spinProjectDir2Minus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir2Minus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir2Minus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir2Minus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir2Minus,T2>,C2> & r)
{
  // Print Diagnostics
  //  cout << "FnAdjMultSprojDir2Minus(l,r) <- adj(l)*FnSpinProjectDir2Minus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir2Minus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir2Minus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir2Minus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir2Minus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Minus >::Type_t
adjMultSprojDir2Minus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir2Minus >::Type_t tmp;

  tmp = spinProjectDir2Minus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir2Minus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus >::Type_t
adjMultSprojDir2Minus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir2Minus >::Type_t d;

  inlineSpinProjDir2Minus(&(b.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);

  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;
    return ret;
}

/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir3Plus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */
// This is a struct for adj(x)*spinProjectDir3Plus(y)
struct FnAdjMultSprojDir3Plus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir3Plus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Plus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {
    // Call the appropriate match
    return (adjMultSprojDir3Plus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir3Plus(r)
// to:    FnAdjMultSprojDir3Plus(l,r)
//
// The spinProjectDir3Plus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir3Plus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir3Plus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir3Plus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir3Plus,T2>,C2> & r)
{
  // Print Diagnostics
  // cout << "FnAdjMultSprojDir3Plus(l,r) <- adj(l)*FnSpinProjectDir3Plus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir3Plus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir3Plus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir3Plus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir3Plus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Plus >::Type_t
adjMultSprojDir3Plus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Plus >::Type_t tmp;

  tmp = spinProjectDir3Plus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir3Plus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus >::Type_t
adjMultSprojDir3Plus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Plus >::Type_t d;

  inlineSpinProjDir3Plus(&(b.elem(0).elem(0).real()),
			 &(d.elem(0).elem(0).real()),
			 1);

    
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0)) ;
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1)) ;

  return ret;
}

/* **************************************************************************
 * **************************************************************************
 * FUSED: adj(x)*spinProjectDir3Minus(y)
 *
 * All that is below is needed.
 *
 * **************************************************************************
 * ************************************************************************* */

// This is a struct for adj(x)*spinProjectDir3Minus(y)
struct FnAdjMultSprojDir3Minus
{
  // Boilerplate
  PETE_EMPTY_CONSTRUCTORS(FnAdjMultSprojDir3Minus)

  // OK This is an operator() so we can create instances of this 
  // object and treat them as functions()
  template<class T1, class T2>
  inline typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Minus >::Type_t
  operator()(const T1 &a, const T2 &b) const
  {

    // Call the appropriate match
    return (adjMultSprojDir3Minus(a,b));
  }
};


// This is an operator* that rewrites:
//
// from:  adj(l)*spinProjectDir3Minus(r)
// to:    FnAdjMultSprojDir3Minus(l,r)
//
// The spinProjectDir3Minus gets grabbed and turned into an op identity
//     adj(l)              gets grabbed and turned into an op identity
//  the operation is encoded in the FnAdjMultSprojDir3Minus 
//

template<class T1,class C1,class T2,class C2>
inline typename MakeReturn<BinaryNode<FnAdjMultSprojDir3Minus,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
  typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t>,
  typename BinaryReturn<C1,C2,FnAdjMultSprojDir3Minus>::Type_t >::Expression_t
operator*(const QDPExpr<UnaryNode<FnAdjoint,T1>,C1> & l,
	  const QDPExpr<UnaryNode<FnSpinProjectDir3Minus,T2>,C2> & r)
{
  // Print Diagnostics
  //  cout << "FnAdjMultSprojDir3Minus(l,r) <- adj(l)*FnSpinProjectDir3Minus(r)" << endl;

  typedef UnaryNode<OpIdentity,T1> NewExpr1_t; // Shorthand for the new type for adj(l)
  typedef UnaryNode<OpIdentity,T2> NewExpr2_t; // Shorthand for the new type for spinProj

  // A name for the new node: Tree_t
  typedef BinaryNode<FnAdjMultSprojDir3Minus,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T1>,C1> >::Leaf_t,
    typename CreateLeaf<QDPExpr<UnaryNode<OpIdentity,T2>,C2> >::Leaf_t> Tree_t;

  // Create the result:  Tree_t = Binary node
  //                     Return type is the return type defined in structs below
  //                     The CreateLeaf-s do the rewriting - 
  //                     unwrap the previous expression and rewrap it as the new type
  return MakeReturn<Tree_t,typename BinaryReturn<C1,C2,FnAdjMultSprojDir3Minus>::Type_t>::make(Tree_t(
    CreateLeaf<QDPExpr<NewExpr1_t,C1> >::make(NewExpr1_t(l.expression().child())),
    CreateLeaf<QDPExpr<NewExpr2_t,C2> >::make(NewExpr2_t(r.expression().child()))
    )
   );
}

// Return types Fused su3*spinProj(Spin4)->Spin2
template<>
struct BinaryReturn< PScalar< ColMat >, Spin4, FnAdjMultSprojDir3Minus > {
  typedef Spin2 Type_t;
};

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir3Minus struct
template<typename T1, typename T2>
inline 
typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Minus >::Type_t
adjMultSprojDir3Minus(const T1& a, const T2& b)
{
  typename BinaryReturn<T1, T2, FnAdjMultSprojDir3Minus >::Type_t tmp;

  tmp = spinProjectDir3Minus(b);
  return( adj(a)*tmp );
}

// This is what you need to specialise now. It'll get called by 
// the operator() of the FnAdjMultSprojDir3Minus struct
template<>
inline 
BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus >::Type_t
adjMultSprojDir3Minus(const PScalar<ColMat>& a, const Spin4& b)
{
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus >::Type_t ret;
  BinaryReturn<PScalar<ColMat>, Spin4, FnAdjMultSprojDir3Minus >::Type_t d;

  inlineSpinProjDir3Minus(&(b.elem(0).elem(0).real()),
			  &(d.elem(0).elem(0).real()),
			  1);
  
  
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(0),ret.elem(0));
  _inline_mult_adj_su3_mat_vec(a.elem(),d.elem(1),ret.elem(1));
    
    return ret;
}


} // namespace QDP;

#endif
