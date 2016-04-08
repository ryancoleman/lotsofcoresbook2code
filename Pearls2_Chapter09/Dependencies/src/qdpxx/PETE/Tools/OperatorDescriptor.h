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

#ifndef PETE_SRC_TOOLS_OPERATORDESCRIPTOR_H
#define PETE_SRC_TOOLS_OPERATORDESCRIPTOR_H

#include <iostream>

using std::ostream;
using std::endl;

#include <string>

using std::string;

#include "Tools/DescriptorBase.h"

class OperatorDescriptor: public DescriptorBase<6> {
public:

  //---------------------------------------------------------------------------
  // Constructors.

  OperatorDescriptor()
  { }

  OperatorDescriptor(const string &tag, const string &func,
		     const string &expr, const string &commnt1, 
		     const string &commnt2, const string &arg = "")
  {
    addData(0, tag);
    addData(1, func);
    addData(2, expr);
    addData(3, commnt1);
    addData(4, commnt2);
    addData(5, arg);
  }

  OperatorDescriptor(const OperatorDescriptor &model)
  : DescriptorBase<6>(model)
  { }

  //---------------------------------------------------------------------------
  // Trivial destructor. 
  
  ~OperatorDescriptor() { }

  //---------------------------------------------------------------------------
  // Copy-assignment operator: just copy members. 

  OperatorDescriptor &operator=(const OperatorDescriptor &rhs)
  {
    DescriptorBase<6>::operator=(rhs);
    
    return *this;
  }
  
  //---------------------------------------------------------------------------
  // Return strings/info. 

  const string tag(bool full = true) const
  {
    if (full)
      return str(0);
    else
      return str(0).substr(0, str(0).find('<'));
  }

  const string &function() const
  {
    return str(1);
  }

  const string &expression() const
  {
    return str(2);
  }

  const string &comment1() const
  {
    return str(3);
  }

  const string &comment2() const
  {
    return str(4);
  }

  const string &argDef() const
  {
    return str(5);
  }

  bool templateArgs() const
  {
    return argDef().size() != 0;
  }

};

inline ostream &operator<<(ostream &os, const OperatorDescriptor &o)
{
  os << "TAG  = " << o.tag() << endl;
  os << "FUNC = " << o.function() << endl;
  os << "EXPR = " << o.expression() << endl;
  os << "COMMENTA = " << o.comment1() << endl;
  os << "COMMENTB = " << o.comment2() << endl;
  os << "ARG  = " << o.argDef() << endl;
  
  return os;
}
  
#endif // PETE_SRC_TOOLS_OPERATORDESCRIPTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: OperatorDescriptor.h,v $   $Author: edwards $
// $Revision: 1.2 $   $Date: 2002-10-14 02:06:56 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
