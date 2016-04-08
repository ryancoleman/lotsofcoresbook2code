// -*- C++ -*-

/*! @file
 * @brief Print profiling info
 *
 * Diagnostics to print profiling.
 */

#ifndef QDP_PROFILE_INCLUDE
#define QDP_PROFILE_INCLUDE

namespace QDP {

typedef unsigned long  QDPTime_t;

//! Get the wallclock time
/*!
  \return The wallclock time (since Epoch) in seconds.
*/
QDPTime_t getClockTime();
void initProfile(const std::string& file, const std::string& caller, int line);
void closeProfile();
void printProfile();
int setProfileLevel(int n);
int setProgramProfileLevel(int n);
int getProfileLevel();
int getProgramProfileLevel();
void pushProfileInfo(int level, const std::string& file, const std::string& caller, int line);
void popProfileInfo();


//--------------------------------------------------------------------------------------
// Selectively turn on profiling
//--------------------------------------------------------------------------------------

#if ! defined(QDP_USE_PROFILING)   
// No profiling
#define QDP_PUSH_PROFILE(a)
#define QDP_POP_PROFILE()

#else   // Profiling enabled

#define QDP_PUSH_PROFILE(a) pushProfileInfo(a, __FILE__, __func__, __LINE__)
#define QDP_POP_PROFILE()  popProfileInfo()

#include <PETE/ForEachInOrder.h>

//-----------------------------------------------------------------------------
// Support of printing
//-----------------------------------------------------------------------------

struct QDPProfile_t;
void registerProfile(QDPProfile_t* qp);


//! Profiling object
/*!
 * Hold profiling state
 */
struct QDPProfile_t
{
  QDPTime_t     time;
  std::string   expr;
  int           count;
  QDPProfile_t* next;

  void print();

  void init();
  QDPProfile_t() {init();}

  //! Profile rhs
  template<class T, class C, class Op, class RHS, class C1>
  QDPProfile_t(const QDPType<T,C>& dest, const Op& op, const QDPExpr<RHS,C1>& rhs)
    {
      init();

      if (getProfileLevel() > 0)
      {
	ostringstream os;
	printExprTree(os, dest, op, rhs);
	expr = os.str();
	registerProfile(this);
      }
    }

  //! Profile  opOuter(rhs)
  template<class T, class C, class Op, class OpOuter, class RHS, class C1>
  QDPProfile_t(const QDPType<T,C>& dest, const Op& op, const OpOuter& opOuter, const QDPExpr<RHS,C1>& rhs)
    {
      init();

      if (getProfileLevel() > 0)
      {
	typedef UnaryNode<OpOuter, typename CreateLeaf<QDPExpr<RHS,C1> >::Leaf_t> Tree_t;
	typedef typename UnaryReturn<C1,OpOuter>::Type_t Container_t;

	ostringstream os;
	printExprTree(os, dest, op, 
		      MakeReturn<Tree_t,Container_t>::make(Tree_t(
			CreateLeaf<QDPExpr<RHS,C1> >::make(rhs))));
	expr = os.str();
	registerProfile(this);
      }
    }

  //! Profile  opOuter(rhs)
  template<class T, class C, class Op, class OpOuter, class T1, class C1>
  QDPProfile_t(const QDPType<T,C>& dest, const Op& op, const OpOuter& opOuter, const QDPType<T1,C1>& rhs)
    {
      init();

      if (getProfileLevel() > 0)
      {
	typedef UnaryNode<OpOuter, typename CreateLeaf<QDPType<T1,C1> >::Leaf_t> Tree_t;
	typedef typename UnaryReturn<C1,OpOuter>::Type_t Container_t;

	ostringstream os;
	printExprTree(os, dest, op, 
		      MakeReturn<Tree_t,Container_t>::make(Tree_t(
			CreateLeaf<QDPType<T1,C1> >::make(rhs))));
	expr = os.str();
	registerProfile(this);
      }
    }
};


struct QDPProfileInfo_t
{
  int           level;
  std::string   file;
  std::string   caller;
  int           line;

  QDPProfileInfo_t() {line=level=0;}

  QDPProfileInfo_t(const QDPProfileInfo_t& a) :
    level(a.level), file(a.file), caller(a.caller), line(a.line) {}

  QDPProfileInfo_t(int _level, const std::string& _file, const std::string& _caller, int _line) :
    level(_level), file(_file), caller(_caller), line(_line) {}
};


struct QDPProfileHead_t
{
  QDPProfileInfo_t  info;
  QDPProfile_t*     start;
  QDPProfile_t*     end;

  QDPProfileHead_t() {start=0; end=0;}
  QDPProfileHead_t(const QDPProfileInfo_t& a) : info(a), start(0), end(0) {}
  QDPProfileHead_t(const QDPProfileHead_t& a) : info(a.info), start(a.start), end(a.end) {}
};



//! Print an expression tree
/*! 
 * PrintAssign traverses the parse-tree and prints the expression 
 * including the assignment to the output stream.
 */
template<class T, class C, class Op, class RHS, class C1>
//inline
void printExprTree(ostream& os, 
		   const QDPType<T,C>& dest, const Op& op, const QDPExpr<RHS,C1>& rhs)
{
  typedef EvalLeaf1    FTag_t;
  typedef OpCombine    CTag_t;
  typedef NullCombine  VTag_t;
  typedef QDPExpr<RHS,C1>  Expr;
    
  typedef typename CreateLeaf<Expr>::Leaf_t Expr_t;
  const Expr_t &e = CreateLeaf<Expr>::make(rhs);

#if 0
  // Compact version
  typedef BinaryNode<Op, QDPType<T,C>, Expr_t> Assign_t;
  typedef ForEachInOrder<Assign_t, PrintTag, PrintTag, NullTag> Print_t;
  Assign_t t(op, dest, e);
  Print_t::apply(t, PrintTag(os), PrintTag(os), NullTag());
#else
  // Makes assignment part special
  typedef ForEachInOrder<RHS, PrintTag, PrintTag, NullTag> Print_t;
  LeafFunctor<C,PrintTag>::apply(static_cast<const C&>(dest),PrintTag(os));
  os << " ";
  TagVisitor<Op,PrintTag>::visit(op, PrintTag(os));
  os << " ";
  Print_t::apply(e, PrintTag(os), PrintTag(os), NullTag());
#endif

  os << ";";
}

//
// struct PrintTag
//
// Tag struct to carry an ostream ref through the parse tree.
//

struct PrintTag
{
  mutable ostream &os_m;
  PrintTag(ostream &os) : os_m(os) {}
};


//
// struct LeafFunctor<QDPType<T,C>, PrintTag>
//
// Specialication of TagFunctor class for applying the PrintTag tag to
// a QDPType. The apply method simply prints the result to the
// ostream carried by PrintTag.
//

template<class T, class C>
struct LeafFunctor<QDPType<T,C>, PrintTag>
{
  typedef int Type_t;
  static int apply(const QDPType<T,C> &s, const PrintTag &f)
    { 
      return LeafFunctor<C,PrintTag>::apply(static_cast<const C&>(s),f);
    }
};

//
// struct LeafFunctor<T, PrintTag>
//
// Specialication of TagFunctor class for applying the PrintTag tag to
// a  "T". The apply method simply prints the result to the
// ostream carried by PrintTag.
//

template<>
struct LeafFunctor<float, PrintTag>
{
  typedef int Type_t;
  static int apply(const float &s, const PrintTag &f)
    { 
      f.os_m << "float"; 
      return 0;
    }
};

template<>
struct LeafFunctor<double, PrintTag>
{
  typedef int Type_t;
  static int apply(const double &s, const PrintTag &f)
    { 
      f.os_m << "double"; 
      return 0;
    }
};

template<>
struct LeafFunctor<int, PrintTag>
{
  typedef int Type_t;
  static int apply(const int &s, const PrintTag &f)
    { 
      f.os_m << "int"; 
      return 0;
    }
};

template<>
struct LeafFunctor<char, PrintTag>
{
  typedef int Type_t;
  static int apply(const char &s, const PrintTag &f)
    { 
      f.os_m << "char"; 
      return 0;
    }
};

template<>
struct LeafFunctor<bool, PrintTag>
{
  typedef int Type_t;
  static int apply(const bool &s, const PrintTag &f)
    { 
      f.os_m << "bool"; 
      return 0;
    }
};

template<class T>
struct LeafFunctor<IScalar<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const IScalar<T> &s, const PrintTag &f)
    { 
      f.os_m << "IScal<";
      LeafFunctor<T,PrintTag>::apply(s.elem(),f);
      f.os_m << ">"; 
      return 0;
    }
};

template<class T, int N>
struct LeafFunctor<ILattice<T,N>, PrintTag>
{
  typedef int Type_t;
  static int apply(const ILattice<T,N> &s, const PrintTag &f)
    { 
      f.os_m << "ILat<";
      LeafFunctor<T,PrintTag>::apply(s.elem(0),f);
      f.os_m << "," << N << ">"; 
      return 0;
    }
};

template<class T>
struct LeafFunctor<RScalar<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const RScalar<T> &s, const PrintTag &f)
    { 
      f.os_m << "RScal<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(),f);
      f.os_m << ">"; 
      return 0;
    }
};

template<class T>
struct LeafFunctor<RComplex<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const RComplex<T> &s, const PrintTag &f)
    { 
      f.os_m << "RCplx<"; 
      LeafFunctor<T,PrintTag>::apply(s.real(),f);
      f.os_m << ">"; 
      return 0;
    }
};

template<class T>
struct LeafFunctor<PScalar<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PScalar<T> &s, const PrintTag &f)
    { 
      f.os_m << "PScal<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(),f);
      f.os_m << ">"; 
      return 0;
    }
};

template<class T>
struct LeafFunctor<PSeed<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PSeed<T> &s, const PrintTag &f)
    { 
      f.os_m << "PSeed<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(0),f);
      f.os_m << ">"; 
      return 0;
    }
};

template <class T, int N>
struct LeafFunctor<PColorMatrix<T,N>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PColorMatrix<T,N> &s, const PrintTag &f)
    { 
      f.os_m << "CMat<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(0,0),f);
      f.os_m << "," << N << ">"; 
      return 0;
    }
};

template <class T, int N>
struct LeafFunctor<PSpinMatrix<T,N>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PSpinMatrix<T,N> &s, const PrintTag &f)
    { 
      f.os_m << "SMat<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(0,0),f);
      f.os_m << "," << N << ">"; 
      return 0;
    }
};

template <class T, int N>
struct LeafFunctor<PColorVector<T,N>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PColorVector<T,N> &s, const PrintTag &f)
    { 
      f.os_m << "ColorVec<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(0),f);
      f.os_m << "," << N << ">"; 
      return 0;
    }
};

template <class T, int N>
struct LeafFunctor<PSpinVector<T,N>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PSpinVector<T,N> &s, const PrintTag &f)
    { 
      f.os_m << "SpinVec<"; 
      LeafFunctor<T,PrintTag>::apply(s.elem(0),f);
      f.os_m << "," << N << ">"; 
      return 0;
    }
};

template <class T, int N, template<class,int> class C>
struct LeafFunctor<PMatrix<T,N,C>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PMatrix<T,N,C> &s, const PrintTag &f)
    { 
      return LeafFunctor<C<T,N>,PrintTag>::apply(static_cast<const C<T,N>&>(s),f);
    }
};

template <class T, int N, template<class,int> class C>
struct LeafFunctor<PVector<T,N,C>, PrintTag>
{
  typedef int Type_t;
  static int apply(const PVector<T,N,C> &s, const PrintTag &f)
    { 
      return LeafFunctor<C<T,N>,PrintTag>::apply(static_cast<const C<T,N>&>(s),f);
    }
};

template<class T>
struct LeafFunctor<OScalar<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const OScalar<T> &s, const PrintTag &f)
    { 
      f.os_m << "OScal<";
      LeafFunctor<T,PrintTag>::apply(s.elem(),f);
      f.os_m << ">"; 
      return 0;
    }
};

template<class T>
struct LeafFunctor<OLattice<T>, PrintTag>
{
  typedef int Type_t;
  static int apply(const OLattice<T> &s, const PrintTag &f)
    { 
      f.os_m << "OLat<";
      LeafFunctor<T,PrintTag>::apply(s.elem(0),f);
      f.os_m << ">"; 
      return 0;
    }
};

template<int N>
struct LeafFunctor<GammaType<N>, PrintTag>
{
  typedef int Type_t;
  static int apply(const GammaType<N> &s, const PrintTag &f)
    { 
      f.os_m << "GammaType";
      return 0;
    }
};

template<int N, int m>
struct LeafFunctor<GammaConst<N,m>, PrintTag>
{
  typedef int Type_t;
  static int apply(const GammaConst<N,m> &s, const PrintTag &f)
    { 
      f.os_m << "GammaConst";
      return 0;
    }
};




//
// struct ParenPrinter
//
// Utility class to provide start and finish functions that print the
// open paren and closed paren for an expression as the ForEach moves
// down and back up an edge of the parse-tree.
//

template<class Op>
struct ParenPrinter
{
  static void start(Op,PrintTag p)
    { p.os_m << "("; }

  static void center(Op,PrintTag p)
    { p.os_m << ","; }

  static void finish(Op,PrintTag p)
    { p.os_m << ")"; }
};


//
// struct TagVisitor<Op, PrintTag>
//
// Specialication of TagVisitor class for applying the PrintTag
// tag to Op nodes. The visit method simply prints a symbol
// to the ostream carried by PrintTag.
//
// FnArcCos
template <>
struct TagVisitor<FnArcCos, PrintTag> : public ParenPrinter<FnArcCos>
{ 
  static void visit(FnArcCos op, PrintTag t) 
    { t.os_m << "acos"; }
};

// FnArcSin
template <>
struct TagVisitor<FnArcSin, PrintTag> : public ParenPrinter<FnArcSin>
{ 
  static void visit(FnArcSin op, PrintTag t) 
    { t.os_m << "asin"; }
};

// FnArcTan
template <>
struct TagVisitor<FnArcTan, PrintTag> : public ParenPrinter<FnArcTan>
{ 
  static void visit(FnArcTan op, PrintTag t) 
    { t.os_m << "atan"; }
};

// FnCeil
template <>
struct TagVisitor<FnCeil, PrintTag> : public ParenPrinter<FnCeil>
{ 
  static void visit(FnCeil op, PrintTag t) 
    { t.os_m << "ceil"; }
};

// FnCos
template <>
struct TagVisitor<FnCos, PrintTag> : public ParenPrinter<FnCos>
{ 
  static void visit(FnCos op, PrintTag t) 
    { t.os_m << "cos"; }
};

// FnHypCos
template <>
struct TagVisitor<FnHypCos, PrintTag> : public ParenPrinter<FnHypCos>
{ 
  static void visit(FnHypCos op, PrintTag t) 
    { t.os_m << "cosh"; }
};

// FnExp
template <>
struct TagVisitor<FnExp, PrintTag> : public ParenPrinter<FnExp>
{ 
  static void visit(FnExp op, PrintTag t) 
    { t.os_m << "exp"; }
};

// FnFabs
template <>
struct TagVisitor<FnFabs, PrintTag> : public ParenPrinter<FnFabs>
{ 
  static void visit(FnFabs op, PrintTag t) 
    { t.os_m << "fabs"; }
};

// FnFloor
template <>
struct TagVisitor<FnFloor, PrintTag> : public ParenPrinter<FnFloor>
{ 
  static void visit(FnFloor op, PrintTag t) 
    { t.os_m << "floor"; }
};

// FnLog
template <>
struct TagVisitor<FnLog, PrintTag> : public ParenPrinter<FnLog>
{ 
  static void visit(FnLog op, PrintTag t) 
    { t.os_m << "log"; }
};

// FnLog10
template <>
struct TagVisitor<FnLog10, PrintTag> : public ParenPrinter<FnLog10>
{ 
  static void visit(FnLog10 op, PrintTag t) 
    { t.os_m << "log10"; }
};

// FnSin
template <>
struct TagVisitor<FnSin, PrintTag> : public ParenPrinter<FnSin>
{ 
  static void visit(FnSin op, PrintTag t) 
    { t.os_m << "sin"; }
};

// FnHypSin
template <>
struct TagVisitor<FnHypSin, PrintTag> : public ParenPrinter<FnHypSin>
{ 
  static void visit(FnHypSin op, PrintTag t) 
    { t.os_m << "sinh"; }
};

// FnSqrt
template <>
struct TagVisitor<FnSqrt, PrintTag> : public ParenPrinter<FnSqrt>
{ 
  static void visit(FnSqrt op, PrintTag t) 
    { t.os_m << "sqrt"; }
};

// FnTan
template <>
struct TagVisitor<FnTan, PrintTag> : public ParenPrinter<FnTan>
{ 
  static void visit(FnTan op, PrintTag t) 
    { t.os_m << "tan"; }
};

// FnHypTan
template <>
struct TagVisitor<FnHypTan, PrintTag> : public ParenPrinter<FnHypTan>
{ 
  static void visit(FnHypTan op, PrintTag t) 
    { t.os_m << "tanh"; }
};

// OpUnaryMinus
template <>
struct TagVisitor<OpUnaryMinus, PrintTag> : public ParenPrinter<OpUnaryMinus>
{ 
  static void visit(OpUnaryMinus op, PrintTag t) 
    { t.os_m << "-"; }
};

// OpUnaryPlus
template <>
struct TagVisitor<OpUnaryPlus, PrintTag> : public ParenPrinter<OpUnaryPlus>
{ 
  static void visit(OpUnaryPlus op, PrintTag t) 
    { t.os_m << "+"; }
};

// OpBitwiseNot
template <>
struct TagVisitor<OpBitwiseNot, PrintTag> : public ParenPrinter<OpBitwiseNot>
{ 
  static void visit(OpBitwiseNot op, PrintTag t) 
    { t.os_m << "~"; }
};

// OpIdentity
template <>
struct TagVisitor<OpIdentity, PrintTag> : public ParenPrinter<OpIdentity>
{ 
  static void visit(OpIdentity op, PrintTag t) 
    { t.os_m << "PETE_identity"; }
};

// OpNot
template <>
struct TagVisitor<OpNot, PrintTag> : public ParenPrinter<OpNot>
{ 
  static void visit(OpNot op, PrintTag t) 
    { t.os_m << "!"; }
};

// OpAdd
template <>
struct TagVisitor<OpAdd, PrintTag> : public ParenPrinter<OpAdd>
{ 
  static void visit(OpAdd op, PrintTag t) 
    { t.os_m << "+"; }
};

// OpSubtract
template <>
struct TagVisitor<OpSubtract, PrintTag> : public ParenPrinter<OpSubtract>
{ 
  static void visit(OpSubtract op, PrintTag t) 
    { t.os_m << "-"; }
};

// OpMultiply
template <>
struct TagVisitor<OpMultiply, PrintTag> : public ParenPrinter<OpMultiply>
{ 
  static void visit(OpMultiply op, PrintTag t) 
    { t.os_m << "*"; }
};

// OpDivide
template <>
struct TagVisitor<OpDivide, PrintTag> : public ParenPrinter<OpDivide>
{ 
  static void visit(OpDivide op, PrintTag t) 
    { t.os_m << "/"; }
};

// OpMod
template <>
struct TagVisitor<OpMod, PrintTag> : public ParenPrinter<OpMod>
{ 
  static void visit(OpMod op, PrintTag t) 
    { t.os_m << "%"; }
};

// OpBitwiseAnd
template <>
struct TagVisitor<OpBitwiseAnd, PrintTag> : public ParenPrinter<OpBitwiseAnd>
{ 
  static void visit(OpBitwiseAnd op, PrintTag t) 
    { t.os_m << "&"; }
};

// OpBitwiseOr
template <>
struct TagVisitor<OpBitwiseOr, PrintTag> : public ParenPrinter<OpBitwiseOr>
{ 
  static void visit(OpBitwiseOr op, PrintTag t) 
    { t.os_m << "|"; }
};

// OpBitwiseXor
template <>
struct TagVisitor<OpBitwiseXor, PrintTag> : public ParenPrinter<OpBitwiseXor>
{ 
  static void visit(OpBitwiseXor op, PrintTag t) 
    { t.os_m << "^"; }
};

// FnLdexp
template <>
struct TagVisitor<FnLdexp, PrintTag> : public ParenPrinter<FnLdexp>
{ 
  static void visit(FnLdexp op, PrintTag t) 
    { t.os_m << "ldexp"; }
};

// FnPow
template <>
struct TagVisitor<FnPow, PrintTag> : public ParenPrinter<FnPow>
{ 
  static void visit(FnPow op, PrintTag t) 
    { t.os_m << "pow"; }
};

// FnFmod
template <>
struct TagVisitor<FnFmod, PrintTag> : public ParenPrinter<FnFmod>
{ 
  static void visit(FnFmod op, PrintTag t) 
    { t.os_m << "fmod"; }
};

// FnArcTan2
template <>
struct TagVisitor<FnArcTan2, PrintTag> : public ParenPrinter<FnArcTan2>
{ 
  static void visit(FnArcTan2 op, PrintTag t) 
    { t.os_m << "atan2"; }
};

// OpLT
template <>
struct TagVisitor<OpLT, PrintTag> : public ParenPrinter<OpLT>
{ 
  static void visit(OpLT op, PrintTag t) 
    { t.os_m << "<"; }
};

// OpLE
template <>
struct TagVisitor<OpLE, PrintTag> : public ParenPrinter<OpLE>
{ 
  static void visit(OpLE op, PrintTag t) 
    { t.os_m << "<="; }
};

// OpGT
template <>
struct TagVisitor<OpGT, PrintTag> : public ParenPrinter<OpGT>
{ 
  static void visit(OpGT op, PrintTag t) 
    { t.os_m << ">"; }
};

// OpGE
template <>
struct TagVisitor<OpGE, PrintTag> : public ParenPrinter<OpGE>
{ 
  static void visit(OpGE op, PrintTag t) 
    { t.os_m << ">="; }
};

// OpEQ
template <>
struct TagVisitor<OpEQ, PrintTag> : public ParenPrinter<OpEQ>
{ 
  static void visit(OpEQ op, PrintTag t) 
    { t.os_m << "=="; }
};

// OpNE
template <>
struct TagVisitor<OpNE, PrintTag> : public ParenPrinter<OpNE>
{ 
  static void visit(OpNE op, PrintTag t) 
    { t.os_m << "!="; }
};

// OpAnd
template <>
struct TagVisitor<OpAnd, PrintTag> : public ParenPrinter<OpAnd>
{ 
  static void visit(OpAnd op, PrintTag t) 
    { t.os_m << "&&"; }
};

// OpOr
template <>
struct TagVisitor<OpOr, PrintTag> : public ParenPrinter<OpOr>
{ 
  static void visit(OpOr op, PrintTag t) 
    { t.os_m << "||"; }
};

// OpLeftShift
template <>
struct TagVisitor<OpLeftShift, PrintTag> : public ParenPrinter<OpLeftShift>
{ 
  static void visit(OpLeftShift op, PrintTag t) 
    { t.os_m << "<<"; }
};

// OpRightShift
template <>
struct TagVisitor<OpRightShift, PrintTag> : public ParenPrinter<OpRightShift>
{ 
  static void visit(OpRightShift op, PrintTag t) 
    { t.os_m << ">>"; }
};

// OpAssign
template <>
struct TagVisitor<OpAssign, PrintTag> : public ParenPrinter<OpAssign>
{ 
  static void visit(OpAssign op, PrintTag t) 
    { t.os_m << "="; }
};

// OpAddAssign
template <>
struct TagVisitor<OpAddAssign, PrintTag> : public ParenPrinter<OpAddAssign>
{ 
  static void visit(OpAddAssign op, PrintTag t) 
    { t.os_m << "+="; }
};

// OpSubtractAssign
template <>
struct TagVisitor<OpSubtractAssign, PrintTag> : public ParenPrinter<OpSubtractAssign>
{ 
  static void visit(OpSubtractAssign op, PrintTag t) 
    { t.os_m << "-="; }
};

// OpMultiplyAssign
template <>
struct TagVisitor<OpMultiplyAssign, PrintTag> : public ParenPrinter<OpMultiplyAssign>
{ 
  static void visit(OpMultiplyAssign op, PrintTag t) 
    { t.os_m << "*="; }
};

// OpDivideAssign
template <>
struct TagVisitor<OpDivideAssign, PrintTag> : public ParenPrinter<OpDivideAssign>
{ 
  static void visit(OpDivideAssign op, PrintTag t) 
    { t.os_m << "/="; }
};

// OpModAssign
template <>
struct TagVisitor<OpModAssign, PrintTag> : public ParenPrinter<OpModAssign>
{ 
  static void visit(OpModAssign op, PrintTag t) 
    { t.os_m << "%="; }
};

// OpBitwiseOrAssign
template <>
struct TagVisitor<OpBitwiseOrAssign, PrintTag> : public ParenPrinter<OpBitwiseOrAssign>
{ 
  static void visit(OpBitwiseOrAssign op, PrintTag t) 
    { t.os_m << "|="; }
};

// OpBitwiseAndAssign
template <>
struct TagVisitor<OpBitwiseAndAssign, PrintTag> : public ParenPrinter<OpBitwiseAndAssign>
{ 
  static void visit(OpBitwiseAndAssign op, PrintTag t) 
    { t.os_m << "&="; }
};

// OpBitwiseXorAssign
template <>
struct TagVisitor<OpBitwiseXorAssign, PrintTag> : public ParenPrinter<OpBitwiseXorAssign>
{ 
  static void visit(OpBitwiseXorAssign op, PrintTag t) 
    { t.os_m << "^="; }
};

// OpLeftShiftAssign
template <>
struct TagVisitor<OpLeftShiftAssign, PrintTag> : public ParenPrinter<OpLeftShiftAssign>
{ 
  static void visit(OpLeftShiftAssign op, PrintTag t) 
    { t.os_m << "<<="; }
};

// OpRightShiftAssign
template <>
struct TagVisitor<OpRightShiftAssign, PrintTag> : public ParenPrinter<OpRightShiftAssign>
{ 
  static void visit(OpRightShiftAssign op, PrintTag t) 
    { t.os_m << ">>="; }
};

// FnWhere
template <>
struct TagVisitor<FnWhere, PrintTag> : public ParenPrinter<FnWhere>
{ 
  static void visit(FnWhere op, PrintTag t) 
    { t.os_m << "where"; }
};

// FnAdjoint
template <>
struct TagVisitor<FnAdjoint, PrintTag> : public ParenPrinter<FnAdjoint>
{ 
  static void visit(FnAdjoint op, PrintTag t) 
    { t.os_m << "adj"; }
};

// FnConjugate
template <>
struct TagVisitor<FnConjugate, PrintTag> : public ParenPrinter<FnConjugate>
{ 
  static void visit(FnConjugate op, PrintTag t) 
    { t.os_m << "conj"; }
};

// FnTranspose
template <>
struct TagVisitor<FnTranspose, PrintTag> : public ParenPrinter<FnTranspose>
{ 
  static void visit(FnTranspose op, PrintTag t) 
    { t.os_m << "transpose"; }
};

// FnTrace
template <>
struct TagVisitor<FnTrace, PrintTag> : public ParenPrinter<FnTrace>
{ 
  static void visit(FnTrace op, PrintTag t) 
    { t.os_m << "trace"; }
};

// FnRealTrace
template <>
struct TagVisitor<FnRealTrace, PrintTag> : public ParenPrinter<FnRealTrace>
{ 
  static void visit(FnRealTrace op, PrintTag t) 
    { t.os_m << "realTrace"; }
};

// FnImagTrace
template <>
struct TagVisitor<FnImagTrace, PrintTag> : public ParenPrinter<FnImagTrace>
{ 
  static void visit(FnImagTrace op, PrintTag t) 
    { t.os_m << "imagTrace"; }
};

// FnTraceColor
template <>
struct TagVisitor<FnTraceColor, PrintTag> : public ParenPrinter<FnTraceColor>
{ 
  static void visit(FnTraceColor op, PrintTag t) 
    { t.os_m << "traceColor"; }
};

// FnTraceSpin
template <>
struct TagVisitor<FnTraceSpin, PrintTag> : public ParenPrinter<FnTraceSpin>
{ 
  static void visit(FnTraceSpin op, PrintTag t) 
    { t.os_m << "traceSpin"; }
};

// FnTransposeSpin
template <>
struct TagVisitor<FnTransposeSpin, PrintTag> : public ParenPrinter<FnTransposeSpin>
{
  static void visit(FnTransposeSpin op, PrintTag t)
  {
    t.os_m << "transposeSpin"; }
};

// FnReal
template <>
struct TagVisitor<FnReal, PrintTag> : public ParenPrinter<FnReal>
{ 
  static void visit(FnReal op, PrintTag t) 
    { t.os_m << "real"; }
};

// FnImag
template <>
struct TagVisitor<FnImag, PrintTag> : public ParenPrinter<FnImag>
{ 
  static void visit(FnImag op, PrintTag t) 
    { t.os_m << "imag"; }
};

// FnLocalNorm2
template <>
struct TagVisitor<FnLocalNorm2, PrintTag> : public ParenPrinter<FnLocalNorm2>
{ 
  static void visit(FnLocalNorm2 op, PrintTag t) 
    { t.os_m << "localNorm2"; }
};

// FnTimesI
template <>
struct TagVisitor<FnTimesI, PrintTag> : public ParenPrinter<FnTimesI>
{ 
  static void visit(FnTimesI op, PrintTag t) 
    { t.os_m << "timesI"; }
};

// FnTimesMinusI
template <>
struct TagVisitor<FnTimesMinusI, PrintTag> : public ParenPrinter<FnTimesMinusI>
{ 
  static void visit(FnTimesMinusI op, PrintTag t) 
    { t.os_m << "timesMinusI"; }
};

// FnSeedToFloat
template <>
struct TagVisitor<FnSeedToFloat, PrintTag> : public ParenPrinter<FnSeedToFloat>
{ 
  static void visit(FnSeedToFloat op, PrintTag t) 
    { t.os_m << "seedToFloat"; }
};

// FnSpinProjectDir0Plus
template <>
struct TagVisitor<FnSpinProjectDir0Plus, PrintTag> : public ParenPrinter<FnSpinProjectDir0Plus>
{ 
  static void visit(FnSpinProjectDir0Plus op, PrintTag t) 
    { t.os_m << "spinProjectDir0Plus"; }
};

// FnSpinProjectDir1Plus
template <>
struct TagVisitor<FnSpinProjectDir1Plus, PrintTag> : public ParenPrinter<FnSpinProjectDir1Plus>
{ 
  static void visit(FnSpinProjectDir1Plus op, PrintTag t) 
    { t.os_m << "spinProjectDir1Plus"; }
};

// FnSpinProjectDir2Plus
template <>
struct TagVisitor<FnSpinProjectDir2Plus, PrintTag> : public ParenPrinter<FnSpinProjectDir2Plus>
{ 
  static void visit(FnSpinProjectDir2Plus op, PrintTag t) 
    { t.os_m << "spinProjectDir2Plus"; }
};

// FnSpinProjectDir3Plus
template <>
struct TagVisitor<FnSpinProjectDir3Plus, PrintTag> : public ParenPrinter<FnSpinProjectDir3Plus>
{ 
  static void visit(FnSpinProjectDir3Plus op, PrintTag t) 
    { t.os_m << "spinProjectDir3Plus"; }
};

// FnSpinProjectDir0Minus
template <>
struct TagVisitor<FnSpinProjectDir0Minus, PrintTag> : public ParenPrinter<FnSpinProjectDir0Minus>
{ 
  static void visit(FnSpinProjectDir0Minus op, PrintTag t) 
    { t.os_m << "spinProjectDir0Minus"; }
};

// FnSpinProjectDir1Minus
template <>
struct TagVisitor<FnSpinProjectDir1Minus, PrintTag> : public ParenPrinter<FnSpinProjectDir1Minus>
{ 
  static void visit(FnSpinProjectDir1Minus op, PrintTag t) 
    { t.os_m << "spinProjectDir1Minus"; }
};

// FnSpinProjectDir2Minus
template <>
struct TagVisitor<FnSpinProjectDir2Minus, PrintTag> : public ParenPrinter<FnSpinProjectDir2Minus>
{ 
  static void visit(FnSpinProjectDir2Minus op, PrintTag t) 
    { t.os_m << "spinProjectDir2Minus"; }
};

// FnSpinProjectDir3Minus
template <>
struct TagVisitor<FnSpinProjectDir3Minus, PrintTag> : public ParenPrinter<FnSpinProjectDir3Minus>
{ 
  static void visit(FnSpinProjectDir3Minus op, PrintTag t) 
    { t.os_m << "spinProjectDir3Minus"; }
};

// FnSpinReconstructDir0Plus
template <>
struct TagVisitor<FnSpinReconstructDir0Plus, PrintTag> : public ParenPrinter<FnSpinReconstructDir0Plus>
{ 
  static void visit(FnSpinReconstructDir0Plus op, PrintTag t) 
    { t.os_m << "spinReconstructDir0Plus"; }
};

// FnSpinReconstructDir1Plus
template <>
struct TagVisitor<FnSpinReconstructDir1Plus, PrintTag> : public ParenPrinter<FnSpinReconstructDir1Plus>
{ 
  static void visit(FnSpinReconstructDir1Plus op, PrintTag t) 
    { t.os_m << "spinReconstructDir1Plus"; }
};

// FnSpinReconstructDir2Plus
template <>
struct TagVisitor<FnSpinReconstructDir2Plus, PrintTag> : public ParenPrinter<FnSpinReconstructDir2Plus>
{ 
  static void visit(FnSpinReconstructDir2Plus op, PrintTag t) 
    { t.os_m << "spinReconstructDir2Plus"; }
};

// FnSpinReconstructDir3Plus
template <>
struct TagVisitor<FnSpinReconstructDir3Plus, PrintTag> : public ParenPrinter<FnSpinReconstructDir3Plus>
{ 
  static void visit(FnSpinReconstructDir3Plus op, PrintTag t) 
    { t.os_m << "spinReconstructDir3Plus"; }
};

// FnSpinReconstructDir0Minus
template <>
struct TagVisitor<FnSpinReconstructDir0Minus, PrintTag> : public ParenPrinter<FnSpinReconstructDir0Minus>
{ 
  static void visit(FnSpinReconstructDir0Minus op, PrintTag t) 
    { t.os_m << "spinReconstructDir0Minus"; }
};

// FnSpinReconstructDir1Minus
template <>
struct TagVisitor<FnSpinReconstructDir1Minus, PrintTag> : public ParenPrinter<FnSpinReconstructDir1Minus>
{ 
  static void visit(FnSpinReconstructDir1Minus op, PrintTag t) 
    { t.os_m << "spinReconstructDir1Minus"; }
};

// FnSpinReconstructDir2Minus
template <>
struct TagVisitor<FnSpinReconstructDir2Minus, PrintTag> : public ParenPrinter<FnSpinReconstructDir2Minus>
{ 
  static void visit(FnSpinReconstructDir2Minus op, PrintTag t) 
    { t.os_m << "spinReconstructDir2Minus"; }
};

// FnSpinReconstructDir3Minus
template <>
struct TagVisitor<FnSpinReconstructDir3Minus, PrintTag> : public ParenPrinter<FnSpinReconstructDir3Minus>
{ 
  static void visit(FnSpinReconstructDir3Minus op, PrintTag t) 
    { t.os_m << "spinReconstructDir3Minus"; }
};

// FnChiralProjectPlus
template <>
struct TagVisitor<FnChiralProjectPlus, PrintTag> : public ParenPrinter<FnChiralProjectPlus>
{ 
  static void visit(FnChiralProjectPlus op, PrintTag t) 
    { t.os_m << "chiralProjectPlus"; }
};

// FnChiralProjectMinus
template <>
struct TagVisitor<FnChiralProjectMinus, PrintTag> : public ParenPrinter<FnChiralProjectMinus>
{ 
  static void visit(FnChiralProjectMinus op, PrintTag t) 
    { t.os_m << "chiralProjectMinus"; }
};

// FnCmplx
template <>
struct TagVisitor<FnCmplx, PrintTag> : public ParenPrinter<FnCmplx>
{ 
  static void visit(FnCmplx op, PrintTag t) 
    { t.os_m << "cmplx"; }
};

// FnOuterProduct
template <>
struct TagVisitor<FnOuterProduct, PrintTag> : public ParenPrinter<FnOuterProduct>
{ 
  static void visit(FnOuterProduct op, PrintTag t) 
    { t.os_m << "outerProduct"; }
};

// FnLocalInnerProduct
template <>
struct TagVisitor<FnLocalInnerProduct, PrintTag> : public ParenPrinter<FnLocalInnerProduct>
{ 
  static void visit(FnLocalInnerProduct op, PrintTag t) 
    { t.os_m << "localInnerProduct"; }
};

// FnLocalInnerProductReal
template <>
struct TagVisitor<FnLocalInnerProductReal, PrintTag> : public ParenPrinter<FnLocalInnerProductReal>
{ 
  static void visit(FnLocalInnerProductReal op, PrintTag t) 
    { t.os_m << "localInnerProductReal"; }
};

// FnQuarkContract13
template <>
struct TagVisitor<FnQuarkContract13, PrintTag> : public ParenPrinter<FnQuarkContract13>
{ 
  static void visit(FnQuarkContract13 op, PrintTag t) 
    { t.os_m << "quarkContract13"; }
};

// FnQuarkContract14
template <>
struct TagVisitor<FnQuarkContract14, PrintTag> : public ParenPrinter<FnQuarkContract14>
{ 
  static void visit(FnQuarkContract14 op, PrintTag t) 
    { t.os_m << "quarkContract14"; }
};

// FnQuarkContract23
template <>
struct TagVisitor<FnQuarkContract23, PrintTag> : public ParenPrinter<FnQuarkContract23>
{ 
  static void visit(FnQuarkContract23 op, PrintTag t) 
    { t.os_m << "quarkContract23"; }
};

// FnQuarkContract24
template <>
struct TagVisitor<FnQuarkContract24, PrintTag> : public ParenPrinter<FnQuarkContract24>
{ 
  static void visit(FnQuarkContract24 op, PrintTag t) 
    { t.os_m << "quarkContract24"; }
};

// FnQuarkContract12
template <>
struct TagVisitor<FnQuarkContract12, PrintTag> : public ParenPrinter<FnQuarkContract12>
{ 
  static void visit(FnQuarkContract12 op, PrintTag t) 
    { t.os_m << "quarkContract12"; }
};

// FnQuarkContract34
template <>
struct TagVisitor<FnQuarkContract34, PrintTag> : public ParenPrinter<FnQuarkContract34>
{ 
  static void visit(FnQuarkContract34 op, PrintTag t) 
    { t.os_m << "quarkContract34"; }
};

// FnColorContract
template <>
struct TagVisitor<FnColorContract, PrintTag> : public ParenPrinter<FnColorContract>
{ 
  static void visit(FnColorContract op, PrintTag t) 
    { t.os_m << "colorContract"; }
};


//-----------------------------------------------------------------------------
// Additional operator tags 
//-----------------------------------------------------------------------------

template <>
struct TagVisitor<FnSum, PrintTag> : public ParenPrinter<FnSum>
{ 
  static void visit(FnSum op, PrintTag t) 
    { t.os_m << "sum"; }
};

template <>
struct TagVisitor<FnNorm2, PrintTag> : public ParenPrinter<FnNorm2>
{ 
  static void visit(FnNorm2 op, PrintTag t) 
    { t.os_m << "norm2"; }
};

template <>
struct TagVisitor<OpGammaConstMultiply, PrintTag> : public ParenPrinter<OpGammaConstMultiply>
{ 
  static void visit(OpGammaConstMultiply op, PrintTag t) 
    { t.os_m << "*"; }
};

template <>
struct TagVisitor<OpMultiplyGammaConst, PrintTag> : public ParenPrinter<OpMultiplyGammaConst>
{ 
  static void visit(OpMultiplyGammaConst op, PrintTag t) 
    { t.os_m << "*"; }
};

template <>
struct TagVisitor<OpGammaTypeMultiply, PrintTag> : public ParenPrinter<OpGammaTypeMultiply>
{ 
  static void visit(OpGammaTypeMultiply op, PrintTag t) 
    { t.os_m << "*"; }
};

template <>
struct TagVisitor<OpMultiplyGammaType, PrintTag> : public ParenPrinter<OpMultiplyGammaType>
{ 
  static void visit(OpMultiplyGammaType op, PrintTag t) 
    { t.os_m << "*"; }
};


//------------------------------------------------------------------------
// Special optimizations
//------------------------------------------------------------------------

template <>
struct TagVisitor<OpAdjMultiply, PrintTag> : public ParenPrinter<OpAdjMultiply>
{ 
  static void visit(OpAdjMultiply op, PrintTag t) 
    { t.os_m << "adjMultiply"; }
};

template <>
struct TagVisitor<OpMultiplyAdj, PrintTag> : public ParenPrinter<OpMultiplyAdj>
{ 
  static void visit(OpMultiplyAdj op, PrintTag t) 
    { t.os_m << "multiplyAdj"; }
};

template <>
struct TagVisitor<OpAdjMultiplyAdj, PrintTag> : public ParenPrinter<OpAdjMultiplyAdj>
{ 
  static void visit(OpAdjMultiplyAdj op, PrintTag t) 
    { t.os_m << "adjMultiplyAdj"; }
};

template <>
struct TagVisitor<FnTraceMultiply, PrintTag> : public ParenPrinter<FnTraceMultiply>
{ 
  static void visit(FnTraceMultiply op, PrintTag t) 
    { t.os_m << "traceMultiply"; }
};

template <>
struct TagVisitor<FnTraceColorMultiply, PrintTag> : public ParenPrinter<FnTraceColorMultiply>
{ 
  static void visit(FnTraceColorMultiply op, PrintTag t) 
    { t.os_m << "traceColorMultiply"; }
};

template <>
struct TagVisitor<FnTraceSpinMultiply, PrintTag> : public ParenPrinter<FnTraceSpinMultiply>
{ 
  static void visit(FnTraceSpinMultiply op, PrintTag t) 
    { t.os_m << "traceSpinMultiply"; }
};

template <>
struct TagVisitor<FnTraceSpinQuarkContract13, PrintTag> : public ParenPrinter<FnTraceSpinQuarkContract13>
{ 
  static void visit(FnTraceSpinQuarkContract13 op, PrintTag t) 
    { t.os_m << "traceSpinQuarkContract13"; }
};

#endif  // ! defined(QDP_USE_PROFILING)

} // namespace QDP

#endif  // QDP_PROFILE_INCLUDE
