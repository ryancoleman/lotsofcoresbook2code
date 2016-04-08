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
 *     Implementation for Simple String Class as DBKey and DBData
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: DBString.h,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#ifndef _FILEDB_DB_STRING_H
#define _FILEDB_DB_STRING_H

#include "DBKey.h"
#include "DBData.h"
#include "FileDB.h"

namespace FILEDB
{
  /**
   * Binary form of StringKey Class
   */
  class StringKey : public DBKey
  {
  public:
    /**
     * Default Constructor
     */
    StringKey (void);

    /**
     * Constructor
     * Using a std::string to construct a StringKey
     */
    StringKey (const std::string& str);
    StringKey (const char* str);

    /**
     * Copy Constructor
     */
    StringKey (const StringKey& key);

    /**
     * Assignment Operator
     */
    StringKey& operator = (const StringKey& key);


    /**
     * Assignment Operator
     */
    StringKey& operator = (const std::string& str);

    /**
     * Equal and not equal operators
     */
    int operator == (const StringKey& key2);
    int operator != (const StringKey& key2); 

    /**
     * Destructor
     */
    ~StringKey (void);

    /**
     * Conversion Operator
     */
    operator std::string (void);

    /**
     * Serial ID or class id 
     */
    const unsigned short serialID (void) const {return 102;}

    /**
     * I am going to use default hash string function
     */
    int hasHashFunc (void) const {return 0;}

    /**
     * I am going to use default compare function
     */
    int hasCompareFunc (void) const {return 0;}

    /**
     * Convert this object into binary form stored in the internal buffer
     *
     */
    void writeObject (std::string&) const throw (SerializeException);


    /**
     * Convert a buffer into an object
     *
     */
    void readObject (const std::string&)  throw (SerializeException);


    /**
     * Empty hash and compare functions. We are using default functions
     * provided by Berkeley DB
     */
    static unsigned int hash (const void* bytes, unsigned int len) {return 0;}
    static int compare (const FFDB_DBT* k1, const FFDB_DBT* k2) {return 0;}
   
  private:
    static const unsigned short MAGIC= (unsigned short)0x11ff;
    static const unsigned short SWAPPED_MAGIC=(unsigned short)0xff11;
    
    // internal string
    std::string str_;
  };


  /**
   * Binary form of userData Class
   */
  class UserData : public DBData
  {
  public:
    /**
     * Default Constructor
     */
    UserData (void);

    /**
     * Constructor
     * Using a std::string to construct a UserData
     */
    UserData (const std::string& str);
    UserData (const char* str);
    UserData (const char* str, unsigned int size);

    /**
     * Copy Constructor
     */
    UserData (const UserData& data);

    /**
     * Assignment Operator
     */
    UserData& operator = (const UserData& data);


    /**
     * Assignment Operator
     */
    UserData& operator = (const std::string& str);


    /**
     * Equal and not equal operators
     */
    int operator == (const UserData& data2);
    int operator != (const UserData& data2); 



    /**
     * Destructor
     */
    ~UserData (void);

    /**
     * Conversion Operator
     */
    operator std::string (void);

    /**
     * Serial ID or class id 
     */
    const unsigned short serialID (void) const {return 0x1144;}

    /**
     * Convert this object into binary form stored in the internal buffer
     *
     */
    void writeObject (std::string&) const throw (SerializeException);


    /**
     * Convert a buffer into an object
     *
     */
    void readObject (const std::string&)  throw (SerializeException);
   
  private:
    static const unsigned short MAGIC= (unsigned short)0x11ff;
    static const unsigned short SWAPPED_MAGIC=(unsigned short)0xff11;
    // internal string
    std::string str_;
  };

}
#endif
