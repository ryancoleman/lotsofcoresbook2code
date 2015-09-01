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
// UnaryOp
// UnarySpecialOp
// UnaryBoolOp
// UnaryCastOp
// BinaryOp
// BinarySpecialOp
// BinaryBoolOp
// BinaryLeftOp
// BinaryAssignOp
// BinaryAssignBoolOp
// TrinaryOp
// InsertOp
//-----------------------------------------------------------------------------

#ifndef PETE_TOOLS_PRINTOPERATORS_H
#define PETE_TOOLS_PRINTOPERATORS_H

//////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
// Overview: 
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
// UnaryOp
//
// This class prints the basic unary operator tag for operations that pass
// the type of the argument through.
//-----------------------------------------------------------------------------

class UnaryOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    ostr << endl;
    if (opdef.templateArgs())
    {
      ostr << "template<" << opdef.argDef() << ">" << endl;
    }
    ostr << "struct " << opdef.tag(false) << endl
	 << "{" << endl;
    if (opdef.templateArgs())
    {
      ostr << "  PETE_EMPTY_CONSTRUCTORS_TEMPLATE("
	   << opdef.tag(false) << "," << opdef.argDef() << ")"
	   << endl;
    }
    else
    {
      ostr << "  PETE_EMPTY_CONSTRUCTORS(" << opdef.tag(false) << ")"
	   << endl;
    }
    ostr << "  template<class T>" << endl
	 << "  inline typename UnaryReturn<T, " << opdef.tag()
	 << " >::Type_t" << endl
	 << "  operator()(const T &a) const" << endl
	 << "  {" << endl
	 << "    " << opdef.expression() << endl
	 << "  }" << endl
	 << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// UnarySpecialOp
//
// This class prints the operator tag for unary operations that compute the
// return type by specializing the UnaryReturn struct. Tag is the same
// as a unary-op since a UnaryReturn<> struct must be provided in either
// case.
//-----------------------------------------------------------------------------

class UnarySpecialOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    UnaryOp().print(ostr, opdef);
  }
};

//-----------------------------------------------------------------------------
// UnaryBoolOp
//
// This class prints the operator tag for unary operations that return
// bools.
//-----------------------------------------------------------------------------

class UnaryBoolOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print Tag:

    UnaryOp().print(ostr, opdef);

    // Print UnaryReturn specialization:
    
    string args = joinWithComma(std::string("class T"), opdef.argDef());
    
    ostr << endl << "template<" << args << " >" << endl
	 << "struct UnaryReturn<T, " << opdef.tag() << " > {"
	 << endl
	 << "  typedef bool Type_t;" << endl
	 << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// UnaryCastOp
//
// This class prints the operator tag for unary operations that perform
// cast expressions.
//-----------------------------------------------------------------------------

class UnaryCastOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print Tag:

    ostr << endl
	 << "template <class T1>" << endl
	 << "struct " << opdef.tag() << "" << endl
	 << "{" << endl
	 << "  PETE_EMPTY_CONSTRUCTORS_TEMPLATE("
	 << opdef.tag(false) << ", T1)" << endl
	 << "  template<class T2>" << endl
	 << "  inline UnaryReturn<T2, " << opdef.tag() << "<T1> >" << endl
	 << "  operator()(const T2 &a) const" << endl
	 << "  {" << endl
	 << "    " << opdef.expression() << endl
	 << "  }" << endl
	 << "};" << endl;

    // Print UnaryReturn specialization:
    
    ostr << endl << "template<class T1, class T2>" << endl
	 << "struct UnaryReturn<T2, " << opdef.tag() << "<T1> > {"
	 << endl
	 << "  typedef T1 Type_t;" << endl
	 << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// BinaryOp
//
// This class prints the operator tag for binary operations that compute
// their return type in the default manner (by promotion).
//-----------------------------------------------------------------------------

class BinaryOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const 
  {
    ostr << endl;
    if (opdef.templateArgs())
      {
        ostr << "template<" << opdef.argDef() << ">" << endl;
      }
    ostr << "struct " << opdef.tag(false) << "" << endl
         << "{" << endl;
    if (opdef.templateArgs())
      {
	ostr << "  PETE_EMPTY_CONSTRUCTORS_TEMPLATE("
	     << opdef.tag(false) << "," << opdef.argDef() << ")"
	     << endl;
      }
    else
      {
	ostr << "  PETE_EMPTY_CONSTRUCTORS(" << opdef.tag(false) << ")"
	     << endl;
      }
    ostr << "  template<class T1, class T2>" << endl
         << "  inline typename BinaryReturn<T1, T2, "
         << opdef.tag() << " >::Type_t" << endl
         << "  operator()(const T1 &a, const T2 &b) const"
         << endl
         << "  {" << endl
         << "    " << opdef.expression() << endl
         << "  }" << endl
         << "};" << endl;
  }
private:
};

//-----------------------------------------------------------------------------
// BinarySpecialOp
//
// This class prints the operator tag for binary operations that compute
// their return type using a specialization of BinaryReturn.
//-----------------------------------------------------------------------------

class BinarySpecialOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const 
  {
    BinaryOp().print(ostr, opdef);
   }
private:
};

//-----------------------------------------------------------------------------
// BinaryBoolOp
//
// This class prints the operator tag for binary operations that return
// bools. Also generates the appropriate BinaryReturn specialization.
//-----------------------------------------------------------------------------

class BinaryBoolOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print tag:

    BinaryOp().print(ostr, opdef);
    
    // Print BinaryReturn specialization:

    string args = joinWithComma("class T1, class T2", 
				     opdef.argDef());
  
    ostr << endl << "template<" << args << " >" << endl
	 << "struct BinaryReturn<T1, T2, " 
	 << opdef.tag() << " > {"
	 << endl
	 << "  typedef bool Type_t;" << endl
	 << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// BinaryLeftOp
//
// This class prints the operator tag for binary operations that return
// the left operand. Also generates the appropriate BinaryReturn 
// specialization.
//-----------------------------------------------------------------------------

class BinaryLeftOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print tag:

    BinaryOp().print(ostr, opdef);

    // Print BinaryReturn specialization:

    string args = joinWithComma("class T1, class T2",
                                     opdef.argDef());

    ostr << endl << "template<" << args << " >" << endl
         << "struct BinaryReturn<T1, T2, "
         << opdef.tag() << " > {"
         << endl
         << "  typedef T1 Type_t;" << endl
         << "};" << endl;
   }
};

//-----------------------------------------------------------------------------
// BinaryAssignOp
//
// This class prints the operator tag for assignment operations.
//-----------------------------------------------------------------------------

class BinaryAssignOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print tag

    ostr << endl;
    if (opdef.templateArgs())
      {
        ostr << "template<" << opdef.argDef() << ">" << endl;
      }
    ostr << "struct " << opdef.tag(false) << "" << endl
         << "{" << endl;
    if (opdef.templateArgs())
      {
	ostr << "  PETE_EMPTY_CONSTRUCTORS_TEMPLATE("
	     << opdef.tag(false) << "," << opdef.argDef() << ")"
	     << endl;
      }
    else
      {
	ostr << "  PETE_EMPTY_CONSTRUCTORS(" << opdef.tag(false) << ")"
	     << endl;
      }
    ostr << "  template<class T1, class T2>" << endl
         << "  inline typename BinaryReturn<T1, T2, "
         << opdef.tag() << " >::Type_t" << endl
         << "  operator()(const T1 &a, const T2 &b) const"
         << endl
         << "  {" << endl
         << "    " << opdef.expression() << endl
         << "  }" << endl
         << "};" << endl;

    // Print BinaryReturn specialization:

    string args = joinWithComma("class T1, class T2",
                                     opdef.argDef());

    ostr << endl << "template<" << args << " >" << endl
         << "struct BinaryReturn<T1, T2, "
         << opdef.tag() << " > {"
         << endl
         << "  typedef T1 &Type_t;" << endl
         << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// BinaryAssignBoolOp
//
// This class prints the operator tag for assignment operations, which (for
// reasons that escape me -- SWH) might return a bool.
//-----------------------------------------------------------------------------

class BinaryAssignBoolOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print tag:

    ostr << endl;
    if (opdef.templateArgs())
      {
        ostr << "template<" << opdef.argDef() << ">" << endl;
      }
    ostr << "struct " << opdef.tag(false) << "" << endl
         << "{" << endl;
    if (opdef.templateArgs())
      {
        ostr << "  PETE_EMPTY_CONSTRUCTORS_TEMPLATE("
             << opdef.tag(false) << "," << opdef.argDef() << ")"
             << endl;
      }
    else
      {
        ostr << "  PETE_EMPTY_CONSTRUCTORS(" << opdef.tag(false) << ")"
             << endl;
      }
    ostr << "  template<class T1, class T2>" << endl
         << "  inline typename BinaryReturn<T1, T2, "
         << opdef.tag() << " >::Type_t" << endl
         << "  operator()(const T1 &a, const T2 &b) const"
         << endl
         << "  {" << endl
         << "    " << opdef.expression() << endl
         << "  }" << endl
         << "};" << endl;

    // Print BinaryReturn specialization:

    string args = joinWithComma("class T1, class T2",
                                     opdef.argDef());

    ostr << endl << "template<" << args << " >" << endl
         << "struct BinaryReturn<T1, T2, "
         << opdef.tag() << " > {"
         << endl
         << "  typedef bool Type_t;" << endl
         << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// TrinaryOp
//
// This class prints the operator tag for trinary operations.
//-----------------------------------------------------------------------------

class TrinaryOp
{
public:
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    // Print tag:

    ostr << endl;
    if (opdef.templateArgs())
      {
        ostr << "template<" << opdef.argDef() << ">" << endl;
      }
    ostr << "struct " << opdef.tag(false) << "" << endl
         << "{" << endl;
    if (opdef.templateArgs())
    {
      ostr << "  PETE_EMPTY_CONSTRUCTORS_TEMPLATE("
	   << opdef.tag(false) << "," << opdef.argDef() << ")"
	   << endl;
    }
    else
    {
      ostr << "  PETE_EMPTY_CONSTRUCTORS(" << opdef.tag(false) << ")"
	   << endl;
    }
    ostr << "  template<class T1, class T2, " 
	 << "class T3>" << endl
         << "  inline typename "
	 << "TrinaryReturn<T1, T2, T3, "
	 << opdef.tag() << " >" << endl
	 << "  ::Type_t" << endl
         << "  operator()(const T1 &a, const T2 &b, const T3 &c) const" << endl
         << "  {" << endl
         << "    " << opdef.expression() << endl
         << "  }" << endl
         << "};" << endl;
  }
};

//-----------------------------------------------------------------------------
// InsertOp
//
// This functor converts an operator descriptor into C++ code that can insert
// the operator into a map of vectors of operator descriptors.
// InsertOp should be constructed with the operator type string that says where
// in the map this type of operator is found.
//-----------------------------------------------------------------------------

class InsertOp
{
public:
  InsertOp(const string& optype) : optype_m(optype) { }
  template<class OSTR,class OPDEF>
  void print(OSTR& ostr,const OPDEF& opdef) const
  {
    ostr
      << "  m[\"" << optype_m << "\"].push_back(" << endl
      << "              OperatorDescriptor(\""
      << opdef.tag() << "\"," << endl
      << "                                 \""
      << opdef.function() << "\"," << endl
      << "                                 \""
      << opdef.expression() << "\"," << endl
      << "                                 \""
      << opdef.comment1() << "\"," << endl
      << "                                 \""
      << opdef.comment2() << "\"," << endl
      << "                                 \""
      << opdef.argDef() << "\"));"
      << endl;
  }
private:
  string optype_m;
};

#endif     // PETE_TOOLS_PRINTOPERATORS_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: PrintOperators.h,v $   $Author: edwards $
// $Revision: 1.3 $   $Date: 2003-02-28 03:41:21 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
