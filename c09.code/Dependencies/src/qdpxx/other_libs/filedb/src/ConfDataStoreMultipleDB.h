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
 *     A Base Class for storing one configuration of key objects 
 *     and associated vector in a hashed file
 * 
 *     This class is used by individual configuration generator to store
 *     data into a single database
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: ConfDataStoreMultipleDB.h,v $
 *   Revision 1.3  2009-10-11 01:50:13  edwards
 *   Changed the getBinary to actually call the internal getBinary. So far,
 *   this has been unused in this wrapper class.
 *
 *   Revision 1.2  2009/08/28 15:42:10  edwards
 *   Added a fileExists function. Removed unsupported code.
 *
 *   Revision 1.1  2009/08/22 12:32:30  edwards
 *   Added a wrapper of multiple DBs.
 *
 *
 *
 */
#ifndef _FILEDB_CONF_DATA_STORE_MULTIPLE_DB_H
#define _FILEDB_CONF_DATA_STORE_MULTIPLE_DB_H

#include "ConfDataStoreDB.h"
#include <vector>

namespace FILEDB
{
  /**
   * Class that holds multiple DBs. Can only be used in a read-only mode.
   */
  template <typename K, typename D>
  class ConfDataStoreMultipleDB 
  {
  protected:
    // The array of DBs
    std::vector< ConfDataStoreDB<K,D> > dbs_;

  public:

    /**
     * Empty constructor for a data store for multiple DBs
     */
    ConfDataStoreMultipleDB (void) {}

    /**
     * Destructor
     */
    virtual ~ConfDataStoreMultipleDB (void) {}


    /**
     * How much data and keys should be kept in memory in bytes
     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSize (const unsigned int size)
    {
      for(int i=0; i < dbs_.size(); ++i)
	dbs_[i].setCacheSize(size);
    }

    
    /**
     * Get maximum number of configurations
     */
    virtual unsigned int getMaxNumberConfigs (void) const
    {
      int ret = dbs_[0].getMaxNumberConfigs();

      return ret;
    }


    /**
     * Check if a DB file exists before opening.
     */
    virtual bool fileExists (const std::vector<std::string>& files) const
    {
      int ret = false;
      for(int i=0; i < dbs_.size(); ++i)
      {
	ret = dbs_[i].fileExists(files[i]);
	if (! ret)
	  return ret;
      }

      return ret;
    }


    /**
     * Open
     * @param files filenames holding all data and keys
     *
     * @return 0 on success, -1 on failure with proper errno set
     * 
     */
    virtual int open (const std::vector<std::string>& files)
    {
      dbs_.resize(files.size());

      int ret = 0;
      for(int i=0; i < dbs_.size(); ++i)
      {
	ret = dbs_[i].open(files[i], O_RDONLY, 0400);
	if (ret != 0)
	  return ret;
      }

      return ret;
    }


    virtual int close (void)
    {
      int ret = 0;
      for(int i=0; i < dbs_.size(); ++i)
      {
	ret = dbs_[i].close();
	if (ret != 0)
	  return ret;
      }

      return ret;
    }

    /**
     * Get data for a given key
     * @param key user supplied key
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int get (const K& key, D& data)
    {
      int ret = -1;

      try {
	for(int i=0; i < dbs_.size(); ++i) 
	{
	  ret = dbs_[i].get(key, data);
	  if (ret == 0)
	    break;
	}
      }
      catch (SerializeException& e) {
	std::cerr << "ConfDataStoreMultipleDB get error: " << e.what () << std::endl;
	ret = -1;
      }
      return ret;
    }


    /**
     * Get data for a given key in binary form
     *
     * @param key user supplied key in string format
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int getBinary (const std::string& key, std::string& data)
    {
      int ret = -1;
      for(int i=0; i < dbs_.size(); ++i)
      {
	ret = dbs_[i].getBinary(key, data);
	if (ret == 0)
	  break;
      }
      return ret;
    }

    /**
     * Does this key exist in the store
     * @param key a key object
     * @return true if the answer is yes
     */
    virtual bool exist(const K& key)
    {
      bool ret = false;
      
      try {
	for(int i=0; i < dbs_.size(); ++i)
	{
	  ret = dbs_[i].exist(key);
	  if (ret)
	    break;
	}
      }
      catch (SerializeException& e) {
	std::cerr << "Key check exist error: " << e.what () << std::endl;
	ret = false;
      }
      
      return ret;
    }


    /**
     * Get user user data from the metadata database
     *
     * @param user_data user supplied buffer to store user data
     * @return returns 0 if success. Otherwise failure.
     */
    virtual int getUserdata (std::string& user_data)
    {
      int ret = -1;

      if (dbs_.size() == 0)
	return ret;

      // Grab the first db
      ret = dbs_[0].getUserdata(user_data);

      return ret;
    }
  };
}
#endif

