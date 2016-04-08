// -*- C++ -*-
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
 *     Data Analysis Configuration Information Data
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *   $Log: ConfigInfo.h,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 */
#ifndef _FILEDB_CONFIG_INFO_DATA_H
#define _FILEDB_CONFIG_INFO_DATA_H

#include <string>
#include <iostream>
#include "FileDB.h"

namespace FILEDB
{
  class ConfigInfo
  {
  private:    
    // config number
    int config_;
    // index number of this configuration
    int index_;
    // is this configuration inserted into the big database
    int inserted_;
    // type of this configuration: fixed or variable data size
    int type_;
    // last modified data
    int mtime_;
    // this configuration filename/or directory name
    std::string url_;
  public:
    /**
     * Constructor
     */
    ConfigInfo (void);
    ConfigInfo (const ConfigInfo& data);
    ConfigInfo& operator = (const ConfigInfo& data);

    /**
     * Destructor
     */
    ~ConfigInfo (void);

    /**
     * All the set and get functions
     * Configuration number
     */
    int configNumber (void) const;
    void configNumber (int number);

    /**
     * Index number of this configuration among all configurations
     */
    int index (void) const;
    void index (int idx);

    /**
     * Has this configuration been inserted into the final database
     */
    int inserted (void) const;
    void insert (int flag);

    /**
     * Check type of this configuration
     */
    int type (void) const;
    void type (int t);

    /**
     * Modified time of this configuration. The time is standard
     * unix time since 1/1/1970
     */
    int modifiedTime (void) const;
    void modifiedTime (int mtime);

    /**
     * Database or directory name associated with this configuration
     */
    std::string urlname (void) const;
    void urlname (const std::string& name);

    
    /**
     * overloaded output stream operator
     */
    friend std::ostream& operator << (std::ostream& out,
				      const ConfigInfo& info);
  };
}
#endif
