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
// Class:
// UnaryFunction
// UnaryCastFunction
// BinaryFunction
// TrinaryFunction
// AssignFunction
//-----------------------------------------------------------------------------

#ifndef PETE_TOOLS_PRINTFUNCTIONS_H
#define PETE_TOOLS_PRINTFUNCTIONS_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
//
// The classes defined here are print functors that will print definitions
// of C++ functions for given operators and classes.  They have static member
// functions called print() that take an output stream and a number of
// operator and class definitions.  Values from the operator and class
// definition objects are used to fill in the definition that is printed.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Typedefs:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Includes:
//-----------------------------------------------------------------------------

#include <iostream>

using std::endl;

#include <string>

using std::string;

#include "Tools/ClassDescriptor.h"
#include "Tools/OperatorDescriptor.h"
#include "Tools/Join.h"

//-----------------------------------------------------------------------------
// Forward Declarations:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Full Description:
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// UnaryFunction
//
// Prints the definition of a unary function that creates a PETE object for a
// given operator from a given class.  First it checks if the class has any
// template arguments.  If the class has template args then the function must
// be templated on those arguments and we must use the keyword "typename"
// to define some types.
//-----------------------------------------------------------------------------

class UnaryFunction
{
public:
  template<class OSTR>
  void print(OSTR& ostr,const OperatorDescriptor& opdef,
	     const ClassDescriptor& class1) const 
  {
    ostr << endl;

    // Print a comment above the function
    if (class1.userClass())
    {
      cerr << "func= " << opdef.function() << "  comment1 = " << opdef.comment1() << endl;
      ostr << opdef.comment1() << endl;
      ostr << opdef.comment2() << endl;
    }

    string args = joinWithComma(opdef.argDef(), class1.argDef(1));
    bool temp = (args.size() > 0);
    string typenameString = temp ? "typename " : "";

    if (temp)
    {
      ostr << "template<" << args <<">" << endl;
    }
    ostr << "inline " << typenameString << "MakeReturn<UnaryNode<"
	 << opdef.tag() << "," << endl
	 << "  " << typenameString << "CreateLeaf<"
	 << class1.inputClass(1) << " >::Leaf_t>," << endl
         << "  typename UnaryReturn<C1," << opdef.tag()
         << " >::Type_t >::Expression_t" << endl
	 << "" << opdef.function() << "(const " << class1.inputClass(1) 
	 << " & l)"
	 << endl
	 << "{" << endl
	 << "  typedef UnaryNode<" << opdef.tag() << "," << endl
	 << "    " << typenameString << "CreateLeaf<"
	 << class1.inputClass(1) << " >::Leaf_t> Tree_t;" << endl
         << "  typedef typename UnaryReturn<C1," << opdef.tag()
         << " >::Type_t Container_t;" << endl
	 << "  return MakeReturn<Tree_t,Container_t>::make(Tree_t(" << endl
	 << "    CreateLeaf<"
	 << class1.inputClass(1) << " >::make(l)));" << endl
	 << "}" << endl;
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class UnaryCastFunction
{
public:
  template<class OSTR>
  void print(OSTR& ostr,const OperatorDescriptor& opdef,
	     const ClassDescriptor& class1) const
  {
    ostr << endl;

    // Print a comment above the function
    if (class1.userClass())
    {
      cerr << "func= " << opdef.function() << "  comment1 = " << opdef.comment1() << endl;
      ostr << opdef.comment1() << endl;
      ostr << opdef.comment2() << endl;
    }

    ostr << "template<" << joinWithComma("class T1",class1.argDef(2))
	 <<">" << endl
	 << "inline typename MakeReturn<UnaryNode<" << opdef.tag() << "<T1>,"
	 << endl
	 << "  typename CreateLeaf<" << class1.inputClass(2)
	 << " >::Leaf_t>," << endl
         << "  typename UnaryReturn<C2," << opdef.tag()
         << " >::Type_t >::Expression_t" << endl
	 << "" << opdef.function() << "(const T1&, const "
	 << class1.inputClass(2) << " & l)" << endl
	 << "{" << endl
	 << "  typedef UnaryNode<" << opdef.tag() << "<T1>," << endl
	 << "    typename CreateLeaf<" << class1.inputClass(2)
	 << " >::Leaf_t> Tree_t;" << endl
         << "  typedef typename UnaryReturn<C2," << opdef.tag()
         << " >::Type_t Container_t;" << endl
	 << "  return MakeReturn<Tree_t,Container_t>::make(Tree_t(" << endl
	 << "    CreateLeaf<"
	 << class1.inputClass(2) << " >::make(l)));" << endl
	 << "}" << endl;
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class BinaryFunction
{
public:
  template<class OSTR>
  void print(OSTR& ostr,const OperatorDescriptor& opdef,
	     const ClassDescriptor& class1,
	     const ClassDescriptor& class2) const
  {
    ostr << endl;

    // Print a comment above the function
    if (class1.userClass() && class2.userClass())
    {
      cerr << "func= " << opdef.function() << "  comment = " << opdef.comment1() << endl;
      ostr << opdef.comment1() << endl;
      ostr << opdef.comment2() << endl;
    }

    string args = joinWithComma(class1.argDef(1),class2.argDef(2));
    bool temp = (args.size() > 0);
    string typenameString = temp ? "typename " : "";

    string inp1 = (class1.scalarClass()) ? 
      "typename SimpleScalar<"+class1.inputClass(1)+">::Type_t" : 
      class1.inputClass(1);
    string inp2 = (class2.scalarClass()) ? 
      "typename SimpleScalar<"+class2.inputClass(2)+">::Type_t" : 
      class2.inputClass(2);

    string C1 = (class1.scalarClass()) ? inp1 : "C1";
    string C2 = (class2.scalarClass()) ? inp2 : "C2";
    string ret = "typename BinaryReturn<"+C1+","+C2+","+opdef.tag()+">::Type_t";

    string arg1 = (class1.scalarClass()) ? inp1+"(l)" : "l";
    string arg2 = (class2.scalarClass()) ? inp2+"(r)" : "r";

    if (temp)
    {
      ostr << "template<" << args <<">" << endl;
    }
    ostr << "inline " << typenameString << "MakeReturn<BinaryNode<"
	 << opdef.tag() << "," << endl
	 << "  " << typenameString << "CreateLeaf<"
         << inp1 << " >::Leaf_t," << endl
	 << "  " << typenameString << "CreateLeaf<"
         << inp2 << " >::Leaf_t>," << endl
         << "  " << ret << " >::Expression_t" << endl
	 << "" << opdef.function() << "(const " << class1.inputClass(1)
	 << " & l,const " << class2.inputClass(2) << " & r)" << endl
	 << "{" << endl;

#if 0
    if (class1.scalarClass())
      ostr << "  typedef " << inp1 << "  Scalar_t;" << endl;
    
    if (class2.scalarClass())
      ostr << "  typedef " << inp2 << "  Scalar_t;" << endl;
#endif

    ostr << "  typedef BinaryNode<" << opdef.tag() << "," << endl
	 << "    " << typenameString << "CreateLeaf<"
	 << inp1 << " >::Leaf_t," << endl
	 << "    " << typenameString << "CreateLeaf<"
	 << inp2 << " >::Leaf_t> Tree_t;" << endl
         << "  typedef " << ret << " Container_t;" << endl
	 << "  return MakeReturn<Tree_t,Container_t>::make(Tree_t(" << endl
	 << "    CreateLeaf<" << inp1 << " >::make(" << arg1 << "),"
	 << endl
	 << "    CreateLeaf<" << inp2 << " >::make(" << arg2 << ")));"
	 << endl
	 << "}" << endl;
  }
};

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

class TrinaryFunction
{
public:
  template<class OSTR>
  void print(OSTR& ostr,
	     const OperatorDescriptor& opdef,
	     const ClassDescriptor& class1,
	     const ClassDescriptor& class2,
	     const ClassDescriptor& class3) const
  {
    ostr << endl;

    // Print a comment above the function
    if (class1.userClass() && class2.userClass() && class3.userClass())
    {
      cerr << "func= " << opdef.function() << "  comment = " << opdef.comment1() << endl;
      ostr << opdef.comment1() << endl;
      ostr << opdef.comment2() << endl;
    }

    string args = joinWithComma(class1.argDef(1),
				class2.argDef(2),
				class3.argDef(3));
    bool temp = (args.size() > 0);
    string typenameString = temp ? "typename " : "";

    string inp1 = (class1.scalarClass()) ? 
      "typename SimpleScalar<"+class1.inputClass(1)+">::Type_t" : 
      class1.inputClass(1);
    string inp2 = (class2.scalarClass()) ? 
      "typename SimpleScalar<"+class2.inputClass(2)+">::Type_t" : 
      class2.inputClass(2);
    string inp3 = (class3.scalarClass()) ? 
      "typename SimpleScalar<"+class3.inputClass(3)+">::Type_t" : 
      class3.inputClass(3);

    string C1 = (class1.scalarClass()) ? inp1 : "C1";
    string C2 = (class2.scalarClass()) ? inp2 : "C2";
    string C3 = (class3.scalarClass()) ? inp3 : "C3";
    string ret = "typename TrinaryReturn<"+C1+","+C2+","+C3+","+opdef.tag()+">::Type_t";

    string arg1 = (class1.scalarClass()) ? inp1+"(a)" : "a";
    string arg2 = (class2.scalarClass()) ? inp2+"(b)" : "b";
    string arg3 = (class3.scalarClass()) ? inp3+"(c)" : "c";

    if (temp)
    {
      ostr << "template<" << args <<">" << endl;
    }
    ostr << "inline " << typenameString << "MakeReturn<TrinaryNode<"
	 << opdef.tag() << "," << endl
	 << "  " << typenameString << "CreateLeaf<"
         << inp1 << " >::Leaf_t," << endl
	 << "  " << typenameString << "CreateLeaf<"
         << inp2 << " >::Leaf_t," << endl
	 << "  " << typenameString << "CreateLeaf<"
         << inp3 << " >::Leaf_t>," << endl
         << "  " << ret << " >::Expression_t" << endl
	 << "" << opdef.function() << "(const " << class1.inputClass(1)
	 << " & a,const " << class2.inputClass(2) << " & b,const "
	 << class3.inputClass(3) << " & c)"
	 << endl
	 << "{" << endl;

    ostr << "  typedef TrinaryNode<" << opdef.tag() << "," << endl
	 << "    " << typenameString << "CreateLeaf<"
	 << inp1 << " >::Leaf_t," << endl
	 << "    " << typenameString << "CreateLeaf<"
	 << inp2 << " >::Leaf_t," << endl
	 << "    " << typenameString << "CreateLeaf<"
	 << inp3 << " >::Leaf_t> Tree_t;" << endl
         << "  typedef " << ret << " Container_t;" << endl
	 << "  return MakeReturn<Tree_t,Container_t>::make(Tree_t(" << endl
	 << "    CreateLeaf<" << inp1 << " >::make(" << arg1 << "),"
	 << endl
	 << "    CreateLeaf<" << inp2 << " >::make(" << arg2 << "),"
	 << endl
	 << "    CreateLeaf<" << inp3 << " >::make(" << arg3 << ")));"
	 << endl
	 << "}" << endl;
  }
};

//-----------------------------------------------------------------------------
// AssignFunctionForClass
//
// Print an operator function that takes LHS and RHS and calls the user-defined
// function evaluate(LHS,OperatorTag,RHS).
// This function allows us to define all the assignment operations (except for
// operator= which is replaced by assign) for classes that don't define them
// as member functions.
//-----------------------------------------------------------------------------

class AssignFunctionForClass
{
public:
  template<class OSTR>
  void print(OSTR& ostr,const OperatorDescriptor& opdef,
	     const ClassDescriptor& class1) const
  {
    ostr << endl;

    // Print a comment above the function
    if (class1.userClass())
    {
      cerr << "func= " << opdef.function() << "  comment1 = " << opdef.comment1() << endl;
      ostr << opdef.comment1() << endl;
      ostr << opdef.comment2() << endl;
    }

    ostr
      << "template<" << joinWithComma(class1.argDef(1),"class RHS")
      <<  ">" << endl
      << "inline" << endl
      << class1.inputClass(1) << "& " << opdef.function()
      << "(" << class1.inputClass(1) << "& lhs,const RHS& rhs)" << endl
      << "{" << endl
      << "  typedef typename CreateLeaf<RHS>::Leaf_t Leaf_t;" << endl
      << "  evaluate(lhs," << opdef.tag()
      << "(),MakeReturn<Leaf_t,C1>::make(CreateLeaf<RHS>::make(rhs)));" 
      << endl
      << "  return lhs;" << endl
      << "}" << endl;
  }
};

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

class AssignFunction
{
public:
  template<class OSTR>
  void print(OSTR& ostr,const OperatorDescriptor& opdef,
	     const ClassDescriptor& class1) const
  {
    ostr << endl;

    // Print a comment above the function
    if (class1.userClass())
    {
      cerr << "func= " << opdef.function() << "  comment1 = " << opdef.comment1() << endl;
      ostr << opdef.comment1() << endl;
      ostr << opdef.comment2() << endl;
    }

    ostr
      << "template<class LHS,class RHS>" << endl
      << "inline LHS&" << endl
      << opdef.function()
      << "(LHS& lhs,const RHS& rhs)" << endl
      << "{" << endl
      << "  typedef typename CreateLeaf<RHS>::Leaf_t Leaf_t;" << endl
      << "  evaluate(lhs," << opdef.tag()
      << "(),MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));" 
      << endl
      << "  return lhs;" << endl
      << "}" << endl;
  }
};

//////////////////////////////////////////////////////////////////////

#endif     // PETE_TOOLS_PRINTFUNCTIONS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PrintFunctions.h,v $   $Author: bjoo $
// $Revision: 1.4 $   $Date: 2005-07-20 11:06:53 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
