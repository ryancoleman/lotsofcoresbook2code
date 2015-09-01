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
 *     Large Data Analysis using file hash as backend
 *                                                     
 * Author:
 *     Jie Chen
 *     Scientific Computing Group 
 *     Jefferson Lab
 *    
 * Revision History:
 *   $Log: FileDB.h,v $
 *   Revision 1.1  2009-03-04 15:56:02  chen
 *   Add FileDB
 *
 *
 *
 */
#ifndef _FILE_DB_H
#define _FILE_DB_H

#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <exception>
#include <vector>

#define FILEDB_FIX_DATA_SIZE 0x1000
#define FILEDB_VAR_DATA_SIZE 0x1001

namespace FILEDB
{

  /**
   * This is a constant used by crc32 routine to generate 
   * entries of the crc table
   */
  const unsigned int crc32poly = (unsigned int)0xedb88320;

  /**
   * This is an exception class notifying serialization of an object
   */
  class SerializeException : public std::exception
  {
  public:
    /**
     * Constructor
     * @param cls class name that produces this exception
     * @param reason what causes this exception
     */
    SerializeException (const std::string& cls, 
			const std::string& reason);

    /**
     * Copy constructor
     */
    SerializeException (const SerializeException& exp);

    /**
     * Assignment operator
     */
    SerializeException& operator = (const SerializeException& exp);

    /**
     * Destructor
     */
    virtual ~SerializeException (void) throw ();

    /**
     * Return reason of the exception
     */
    virtual const char* what (void) const throw ();

  protected:
    std::string cls_;
    std::string reason_;

    // hide default exception
    SerializeException (void);
  };


  /**
   * This is an exception class notifying error of underlying ffdb db
   */
  class FileHashDBException : public std::exception
  {
  public:
    /**
     * Constructor
     * @param cls class name that produces this exception
     * @param reason what causes this exception
     */
    FileHashDBException (const std::string& cls, 
			 const std::string& reason);

    /**
     * Copy constructor
     */
    FileHashDBException (const FileHashDBException& exp);

    /**
     * Assignment operator
     */
    FileHashDBException& operator = (const FileHashDBException& exp);

    /**
     * Destructor
     */
    virtual ~FileHashDBException (void) throw ();

    /**
     * Return reason of the exception
     */
    virtual const char* what (void) const throw ();

  protected:
    std::string cls_;
    std::string reason_;

    // hide default exception
    FileHashDBException (void);
  };

  /**
   * Initialize CRC32 table
   */
  void crc32_init (void);


  /**
   * Calculate CRC32 check sum of a buffer with size len
   * Start by passing in an initial chaining value of 9, and then
   * pass in the return value from the previous crc32 call.
   */
  unsigned int crc32 (unsigned int pcrc, const char* buf, unsigned int len);

  /**
   * Several byte swaping routines replacing reguler htonl etc
   * The primary reason to have these routines is to not include
   * some header files that may not be available on some platforms
   */
  unsigned short ntoh_short (unsigned short usval);
  unsigned short hton_short (unsigned short usval);
  unsigned int   ntoh_int   (unsigned int   ival);
  unsigned int   hton_int   (unsigned int   ival);
  unsigned long  ntoh_long  (unsigned long  ival);
  unsigned long  hton_long  (unsigned long  ival);


  /**
   * Convert a integer into a string
   */
  std::string itostr (int val);
    

  /**
   * Get current time in seconds
   */
  double current_time (void);
}
#endif
