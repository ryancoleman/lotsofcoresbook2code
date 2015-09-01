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
 *     Interface adopted from java to serialize a data and key for data
 *     analysis
 *
 *     This class uses an import and export buffer semantic. Namely,
 *     an internal object is serialized into a binary string, and 
 *     imported from a binary string. Management of the input and output
 *     buffers is by the caller.
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: Serializable.h,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 */
#ifndef _FILEDB_SERIALIZABLE_H
#define _FILEDB_SERIALIZABLE_H

#include "FileDB.h"
#include <string>

namespace FILEDB
{
  class Serializable
  {
  public:
    /**
     * Destructor
     */
    virtual ~Serializable (void) {;}

    /**
     * Get the serial id of this class
     */
    virtual const unsigned short serialID (void) const = 0;

    /**
     * Return this object into a binary form
     */
    virtual void writeObject (std::string& output) const throw (SerializeException) = 0;


    /**
     * Convert input object retrieved from database or network into an object
     */
    virtual void readObject (const std::string& input) throw (SerializeException) = 0;


  protected:
    /**
     * Constructor
     */
    Serializable (void) {;}
  };
}
#endif
