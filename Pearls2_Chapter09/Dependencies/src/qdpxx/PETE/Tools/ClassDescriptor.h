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

#ifndef PETE_SRC_TOOLS_CLASSDESCRIPTOR_H
#define PETE_SRC_TOOLS_CLASSDESCRIPTOR_H

#include <iostream>
#include <stdio.h>
#include <string>

using std::string;
using std::ostream;
using std::endl;

#include "Tools/DescriptorBase.h"

class ClassDescriptor: public DescriptorBase<2> {
public:

  //---------------------------------------------------------------------------
  // Constructors.

  ClassDescriptor() {scalar_class = false; user_class = true;}

  ClassDescriptor(const string &ad, const string &ic, 
		  const bool _scalar_class, const bool _user_class)
  {
    addData(0, ad);
    addData(1, ic);
    scalar_class = _scalar_class;
    user_class = _user_class;
  }
  
  ClassDescriptor(const ClassDescriptor &model)
    : DescriptorBase<2>(model), scalar_class(model.scalar_class),
    user_class(model.user_class)
  { }

  //---------------------------------------------------------------------------
  // Trivial destructor. 
  
  ~ClassDescriptor() { }
  
  //---------------------------------------------------------------------------
  // Copy-assignment operator: just copy members. 

  ClassDescriptor &operator=(const ClassDescriptor &rhs)
  {
    DescriptorBase<2>::operator=(rhs);
    scalar_class = rhs.scalar_class;
    user_class = rhs.user_class;
    
    return *this;
  }
  
  //---------------------------------------------------------------------------
  // Return strings with numbers/args substituted. 

  string argDef(int i) const
  {
    return substituteNum(i, str(0));
  }

  string inputClass(int i) const
  {
    return substituteNum(i, str(1));
  }

  bool scalarClass() const
  {
    return scalar_class;
  }

  bool userClass() const
  {
    return user_class;
  }

private:
  bool scalar_class;   // indicates this is a number like class
  bool user_class;     // indicates this is read in from a file - a user class


  string substituteNum(int i, const string &s) const
  {
    char n[2];
    sprintf(n, "%d", i);
    string str(s), rep("[n]"), num(n);
    int pos;
    
    while ((pos = str.find(rep, 0)) < str.size())
      str.replace(pos, 3, num);
      
    return str;
  }

  // Currently substituteArg is unused.  Pooma used to convert
  // arguments without CreateLeaf, so this function was useful
  // for conversions like l -> Scalar<T1>(l) 

  string substituteArg(const string &arg, 
    const string &s) const
  {
    string str(s), rep("[arg]");
    int pos;
    
    while ((pos = str.find(rep, 0)) < str.size())
      str.replace(pos, 5, arg);
      
    return str;
  }
};

inline ostream &operator<<(ostream &os, const ClassDescriptor &o)
{
  os << "ARG   = " << o.argDef(1) << endl;
  os << "CLASS = " << o.inputClass(1) << endl;
  
  return os;
}
  
#endif // PETE_SRC_TOOLS_CLASSDESCRIPTOR_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: ClassDescriptor.h,v $   $Author: edwards $
// $Revision: 1.2 $   $Date: 2002-10-14 02:06:56 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
