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
 *     Interface for database key class
 *     Subclasses should at least provide "=" and "==" operators
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: DBKey.h,v $
 *   Revision 1.4  2009-03-05 00:40:05  edwards
 *   Changed include path of filehash files back to no relative path.
 *
 *   Revision 1.3  2009/03/04 19:13:05  edwards
 *   Changed some include guards and paths to filehash to be relative.
 *
 *   Revision 1.2  2009/03/04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#ifndef _FILEDB_DB_KEY_H
#define _FILEDB_DB_KEY_H

#include "Serializable.h"
#include "ffdb_db.h"

namespace FILEDB
{
  /**
   * Database Hash function definition
   */
  typedef unsigned int (*ffdb_key_hash) (const void* bytes, 
					 unsigned int len); 

  /**
   * Database Key compare function definition
   */
  typedef int (*ffdb_key_compare) (const FFDB_DBT* k1, const FFDB_DBT* k2);

  class DBKey : public Serializable
  {
  public:
    /**
     * Destructor
     */
    ~DBKey (void) {;}

    /**
     * Does this key provide its own hash function
     * If this class is going to provide the hash function, use the
     * above hash function definition to implement a static function
     * with name hash
     *
     * @return 1 this class provide hash function, 0 otherwise
     */
    virtual int hasHashFunc (void) const = 0;

    /**
     * Does this key provide its own btree key compare function
     * If this class is going to provide the compare function, use the
     * above compare function definition to implement a static function
     * with name compare
     *
     * @return 1 this class provide compare function, 0 otherwise
     */
    virtual int hasCompareFunc (void) const = 0;

  protected:
    /**
     * Constructor
     */
    DBKey (void) {;}

  };
}
#endif
