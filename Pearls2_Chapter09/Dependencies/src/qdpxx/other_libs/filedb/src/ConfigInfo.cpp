/*----------------------------------------------------------------------------
 * Copyright (c) 2007      Jefferson Science Associates, LLC
 *                         Under U.S. DOE Contract No. DE-AC05-06OR23177
 *
 *                         Thomas Jefferson National Accelerator Facility
 *
 *                         Jefferson Lab
 *                         Scientific Computing Group,
 *                         12000 Jefferson Ave.,
 *                         Newport News, VA 23606
 *---------------------------------------------------------------------------- 
 *
 * Description:
 *     Configuration Information
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group 
 *     Jefferson Lab
 *
 * Revision History:
 *   $Log: ConfigInfo.cpp,v $
 *   Revision 1.3  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.2  2009/02/24 04:22:55  edwards
 *   Added include of time.h
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 */
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include "ConfigInfo.h"
#include <time.h>

using namespace std;

namespace FILEDB
{
  ConfigInfo::ConfigInfo (void)
    :config_ (0), index_ (0), inserted_ (0),
     type_ (FILEDB_FIX_DATA_SIZE), mtime_ (0), url_ ()
  {
    // empty
  }

  ConfigInfo::ConfigInfo (const ConfigInfo& data)
    :config_ (data.config_), index_ (data.index_),
     inserted_ (data.inserted_), type_ (data.type_), mtime_ (data.mtime_), 
     url_ (data.url_)
  {
    // empty
  }

  ConfigInfo&
  ConfigInfo::operator = (const ConfigInfo& data)
  {
    if (this != &data) {
      config_ = data.config_;
      index_ = data.index_;
      inserted_ = data.inserted_;
      type_ = data.type_;
      
      mtime_ = data.mtime_;
      url_ = data.url_;
    }
    return *this;
  }

  ConfigInfo::~ConfigInfo (void)
  {
    // empty
  }

  int
  ConfigInfo::configNumber (void) const
  {
    return config_;
  }

  void
  ConfigInfo::configNumber (int conf)
  {
    config_ = conf;
  }

  int
  ConfigInfo::index (void) const
  {
    return index_;
  }

  void
  ConfigInfo::index (int idx)
  {
    index_ = idx;
  }

  int
  ConfigInfo::inserted (void) const
  {
    return inserted_;
  }

  void
  ConfigInfo::insert (int flag)
  {
    inserted_ = flag;
  }

  int
  ConfigInfo::type (void) const
  {
    return type_;
  }

  void
  ConfigInfo::type (int t) 
  {
    type_ = t;
  }

  int
  ConfigInfo::modifiedTime (void) const
  {
    return mtime_;
  }

  void
  ConfigInfo::modifiedTime (int mtime)
  {
    mtime_ = mtime;
  }

  string
  ConfigInfo::urlname (void) const
  {
    return url_;
  }

  void
  ConfigInfo::urlname (const string& name)
  {
    url_ = name;
  }

  ostream&
  operator << (ostream& out, const ConfigInfo& info)
  {
    char mstr[128], *p;

    if (info.inserted_) {
      time_t tt = (time_t)info.mtime_;
      ::ctime_r (&tt, mstr);
      // convert the last \n into \0
      p = mstr;
      while (*p != '\n')
	p++;
      *p = '\0';

      out << info.config_ << " " << info.index_ << " "
	  << " " << info.inserted_  
	  << " " << mstr << " " << info.url_ << "\n";
    }
    else
      out << info.config_ << " " << info.index_ << " "
	  << " " << info.inserted_ 
	  << " " << "N/A" << " " << "N/A" << "\n";
    

    return out;
  }
}
