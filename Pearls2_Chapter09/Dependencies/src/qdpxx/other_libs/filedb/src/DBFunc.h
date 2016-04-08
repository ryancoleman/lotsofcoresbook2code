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
 *     Template Functions for Database Storage and Retrieve
 *
 *     Functions are dealing with classes derived from DBKey and DBData
 * 
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: DBFunc.h,v $
 *   Revision 1.6  2009-08-28 15:42:22  edwards
 *   Added a fileExists function.
 *
 *   Revision 1.5  2009/03/05 00:40:05  edwards
 *   Changed include path of filehash files back to no relative path.
 *
 *   Revision 1.4  2009/03/04 19:13:05  edwards
 *   Changed some include guards and paths to filehash to be relative.
 *
 *   Revision 1.3  2009/03/04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.2  2009/03/02 23:58:21  chen
 *   Change implementation on keys iterator which get keys only
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#ifndef _FILEDB_DB_FUNC_H
#define _FILEDB_DB_FUNC_H

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <sstream>
#include "FileDB.h"
#include "ffdb_db.h"


namespace FILEDB
{
  /**
   * Return all keys to vectors in binary form of strings
   *
   */
  extern void binaryAllKeys (FFDB_DB* dbh, 
			     std::vector<std::string>& keys)
    throw (FileHashDBException);

  /**
   * Return all keys and data to vectors in binary form of strings
   *
   */
  extern void binaryAllPairs (FFDB_DB* dbh, 
			      std::vector<std::string>& keys, 
			      std::vector<std::string>& data) 
    throw (FileHashDBException);

  /**
   * get key and data pair from a database pointed by pointer dbh
   *
   * @param dbh database pointer
   * @key key associated with this data. This key must be string form
   * @data data to be stored into the database. It is in string form
   *
   * @return 0 on success. Otherwise failure
   */
  extern int getBinaryData (FFDB_DB* dbh, const std::string& key, std::string& data);


  /**
   * Insert key and data pair in string format into the database
   *
   * @param dbh database pointer
   * @key key associated with this data. This key must be string form
   * @data data to be stored into the database. It is in string form
   *
   * @return 0 on success. Otherwise failure
   */
  extern int insertBinaryData (FFDB_DB* dbh, const std::string& key, const std::string& data);
  
  /**
   * Open a database with name dbase. If this database does not exist,
   * a new database will be created.
   *
   * @param dbase     database file name
   * @param open_flags database flags for opening: it is can be regular unix
   * open flags: such as O_RDONLY, O_RDWR, O_TRUNC and so on
   * @param mode file creation mode
   * @return Database pointer. If something wrong, a NULL pointer is returned
   */
  template <typename K>
  FFDB_DB* openDatabase (const std::string& dbase, 
			 int open_flags, int mode, FFDB_HASHINFO* info)
  {
    K tkey;
    FFDB_DB* dbh;

    // open database
    if (tkey.hasHashFunc())
      info->hash = (K::hash);

    if (tkey.hasCompareFunc())
      info->cmp = (K::compare);

    // other information elements have been set before calling this one
    // open database
    std::cerr << "Open store database " << dbase << std::endl;

    // turn off umask
    ::umask (0);

    dbh = ffdb_dbopen (dbase.c_str(), open_flags, mode, (const void *)info);

    return dbh;
  }


  /**
   * Insert key and data pair into a database pointed by pointer dbh
   *
   * @param dbh database pointer
   * @key key associated with this data. This key must be subclass of DBKey
   * @data data to be stored into the database. This data must be subclass of
   * DBData
   * @param flag database put flag: FFDB_NOOVERWRITE or 0
   *
   * @return 0 on success. Otherwise failure
   */
  template <typename K, typename D>
  int insertData (FFDB_DB* dbh, const K& key, const D& data,
		  unsigned int flag = 0)
    throw (SerializeException)
  {
    int ret;
       
    // convert key into its binary form
    std::string keyObj;
    try {
      key.writeObject (keyObj);
    }
    catch (SerializeException& e) {
      throw;
    }
    
    // create key
    FFDB_DBT dbkey;
    dbkey.data = &keyObj[0];
    dbkey.size = keyObj.size();
          
    // Convert data into binary form
    std::string dataObj;
    try {
      data.writeObject (dataObj);
    }catch (SerializeException& e) {
      throw;
    }

    // create DBt object
    FFDB_DBT dbdata;
    dbdata.data = &dataObj[0];
    dbdata.size = dataObj.size();

    // now it is time to insert
    ret = dbh->put (dbh, &dbkey, &dbdata, flag);

    return ret;
  }


  /**
   * Insert key and data pair into a database pointed by pointer dbh
   *
   * @param dbh database pointer
   * @key key associated with this data. This key must be subclass of DBKey
   * @data data to be stored into the database. This data must be subclass of
   * DBData
   * @param flag database put flag: FFDB_NOOVERWRITE or 0
   *
   * @return 0 on success. Otherwise failure
   */
  template <typename K>
  int insertData (FFDB_DB* dbh, const K& key, const std::string& data,
		  unsigned int flag = 0)
    throw (SerializeException)
  {
    int ret;
       
    // convert key into its binary form
    std::string keyObj;
    try {
      key.writeObject (keyObj);
    }
    catch (SerializeException& e) {
      throw;
    }
    
    // create key
    FFDB_DBT dbkey;
    dbkey.data = const_cast<char*>(keyObj.c_str());
    dbkey.size = keyObj.size();
          
    // create DBt object
    FFDB_DBT dbdata;
    dbdata.data = const_cast<char*>(data.c_str());
    dbdata.size = data.size();

    // now it is time to insert
    ret = dbh->put (dbh, &dbkey, &dbdata, flag);

    return ret;
  }


  /**
   * get key and data pair from a database pointed by pointer dbh
   *
   * @param dbh database pointer
   * @key key associated with this data. This key must be subclass of DBKey
   * @data data to be stored into the database. This data must be subclass of
   * DBData
   *
   * @return 0 on success. Otherwise failure
   */
  template <typename K, typename D>
  int getData (FFDB_DB* dbh, const K& key, D& data) throw (SerializeException)
  {
    int ret = 0;

    // first convert key into binary buffer
    std::string keyObj;
    try {
      key.writeObject (keyObj);
    }
    catch (SerializeException& e) {
      ret = -1;
      throw;
    }

    // create key
    FFDB_DBT dbkey;
    dbkey.data = &keyObj[0];
    dbkey.size = keyObj.size();

    // create and empty dbt data object
    FFDB_DBT dbdata;
    dbdata.data = 0;
    dbdata.size = 0;

    // now retrieve data from database
    ret = dbh->get (dbh, &dbkey, &dbdata, 0);
    if (ret == 0) {
      try {
	// convert object into a string
	std::string dataObj;
	dataObj.assign((char*)dbdata.data, dbdata.size);
	data.readObject (dataObj);
	// I have to use free since I use malloc in c code
	free(dbdata.data);
      }
      catch (SerializeException& e) {
	ret = -1;
	throw;
      }
    }
    return ret;
  }

  /**
   * Get key and data pair from a database pointed by pointer dbh
   * The data item is in binary string form
   *
   * @param dbh database pointer
   * @key key associated with this data. This key must be subclass of DBKey
   * @data data to be stored into the database. This data must be subclass of
   * DBData
   *
   * @return 0 on success. Otherwise failure
   */
  template <typename K, typename D>
  int getData (FFDB_DB* dbh, const K& key, std::string& data) throw (SerializeException)
  {
    int ret = 0;

    // first convert key into binary buffer
    std::string keyObj;
    try {
      key.writeObject (keyObj);
    }
    catch (SerializeException& e) {
      ret = -1;
      throw;
    }

    // create key
    FFDB_DBT dbkey;
    dbkey.data = &keyObj[0];
    dbkey.size = keyObj.size();

    // create and empty dbt data object
    FFDB_DBT dbdata;
    dbdata.data = 0;
    dbdata.size = 0;

    // now retrieve data from database
    ret = dbh->get (dbh, &dbkey, &dbdata, 0);
    if (ret == 0) {
      // convert object into a string
      data.assign((char*)dbdata.data, dbdata.size);
      // I have to use free since I use malloc in c code
      free(dbdata.data);
    }
    return ret;
  }


  /**
   * Return all keys to a vector provided by an application
   *
   */
  template <typename K, typename D>
  void allKeys (FFDB_DB* dbh, std::vector<K>& keys) 
    throw (SerializeException, FileHashDBException)
  {
    FFDB_DBT  dbkey;
    ffdb_cursor_t *crp;
    K    arg;
    int  ret;
    
    try {
      // create cursor
      ret = dbh->cursor (dbh, &crp, FFDB_KEY_CURSOR);
      if (ret != 0) 
	throw FileHashDBException ("DBFunc allKeys", "Cannot create cursor");

      // get everything from meta dat
      dbkey.data = 0;
      dbkey.size = 0;
      while ((ret = crp->get (crp, &dbkey, 0, FFDB_NEXT)) == 0) {
	// convert into key object
	std::string keyObj;
	keyObj.assign((char*)dbkey.data, dbkey.size);
	arg.readObject (keyObj);

	// put this new key into the vector
	keys.push_back (arg);

	// free memory
	free (dbkey.data);

	dbkey.data = 0;
	dbkey.size = 0;
      }
      if (ret != FFDB_NOT_FOUND) 
	throw FileHashDBException ("DBFunc AllKeys", "Cursor next error");
    }
    catch (SerializeException& e) {
      throw;
    } 
    
    // close cursor
    if (crp)
      crp->close(crp);
  }



  /**
   * Return all keys and data to vectors provided by an application
   *
   */
  template <typename K, typename D>
  void allPairs (FFDB_DB* dbh, std::vector<K>& keys, std::vector<D>& data) 
    throw (SerializeException, FileHashDBException)
  {
    FFDB_DBT  dbkey, dbdata;
    ffdb_cursor_t* crp;
    K    arg;
    D    d;
    int  ret;
    
    try {
      // create cursor
      ret = dbh->cursor (dbh, &crp, FFDB_KEY_CURSOR);
      if (ret != 0) 
	throw FileHashDBException ("DBFunc AllPairs", "Create Cursor Error");

      // get everything from meta data
      dbkey.data = dbdata.data = 0;
      dbkey.size = dbdata.size = 0;
      while ((ret = crp->get (crp, &dbkey, &dbdata, FFDB_NEXT)) == 0) {
	// convert into key object
	std::string keyObj;
	keyObj.assign((char*)dbkey.data, dbkey.size);
	arg.readObject (keyObj);

	// convert into data object
	std::string dataObj;
	dataObj.assign((char*)dbdata.data, dbdata.size);
	d.readObject (dataObj);

	// put this new key into the vector
	keys.push_back (arg);
	
	// put this data into the other vector
	data.push_back (d);

	// free memory
	free (dbkey.data); free (dbdata.data);
	dbkey.data = dbdata.data = 0;
	dbkey.size = dbdata.size = 0;
      }
      if (ret != FFDB_NOT_FOUND) 
	throw FileHashDBException ("DBFunc AllPairs", "Cursor Next Error");
    }
    catch (SerializeException& e) {
      throw;
    }
    
    // close cursor
    if (crp != NULL)
      crp->close(crp);
  }


  /**
   * Check whether this database is empty or not
   *
   */
  template <typename K>
  int isDatabaseEmpty (FFDB_DB* dbh)
    throw (SerializeException, FileHashDBException)
  {
    FFDB_DBT dbkey, dbdata;
    ffdb_cursor_t* crp;
    K    arg;
    int  ret;
    
    try {
      // create cursor
      ret = dbh->cursor (dbh, &crp, FFDB_KEY_CURSOR);
      if (ret != 0) 
	throw FileHashDBException ("DBFunc isDatabaseEmpty", "cursor creation error");

      // get first thing from the database 
      dbkey.data = dbdata.data = 0;
      dbkey.size = dbdata.size = 0;
      ret = crp->get (crp, &dbkey, &dbdata, FFDB_NEXT);
      if (ret == 0) {
	free (dbkey.data);
	free (dbdata.data);
	crp->close(crp);
	return 0;
      }

      if (ret != FFDB_NOT_FOUND) 
	throw FileHashDBException ("DBFunc isDatabaseEmpty", "cursor first error");
    }
    catch (SerializeException& e) {
      throw;
    }
    
    // close cursor
    if (crp != NULL)
      crp->close(crp);

    return 1;
  }

  /**
   * Check whether this key exists in the database
   *
   * @param dbh database pointer
   * @param arg key we are interested
   *
   * @return 1 if this key exists. return 0 not there
   */
  template <typename K>
  int keyExist (FFDB_DB* dbh, const K& arg) 
    throw (SerializeException)
  {
    int ret;
    
    std::string keyObj;
    try {
      arg.writeObject (keyObj);
    }
    catch (SerializeException& e) {
      throw;
    }

    // create Db key object
    FFDB_DBT dbkey;
    dbkey.data = &keyObj[0];
    dbkey.size = keyObj.size();

    // dbdata is not used
    FFDB_DBT dbdata;
    dbdata.data = 0;
    dbdata.size = 0;
    
    try {
      ret = dbh->get (dbh, &dbkey, &dbdata, 0);
      if (ret == 0) {
	free (dbdata.data);
	return 1;
      }
    }
    catch (SerializeException& e) {
      throw;
    }
    return 0;
  }
  

  /**
   * flush data base to disk
   */
  void flushDatabase (FFDB_DB* dbh);


  /**
   * From filename convert to base directory name
   */
  std::string unixDirectoryName (const std::string& filename);


  /**
   * Check if a file exists
   */
  bool fileExists (const std::string& filename);
}
#endif
