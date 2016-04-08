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
 *     Large Data Analysis using database as backend
 *                                                     
 * Author:
 *     Jie Chen
 *     Scientific Computing Group 
 *     Jefferson Lab
 *    
 * Revision History:
 *   $Log: FileDB.cpp,v $
 *   Revision 1.1  2009-03-04 15:56:02  chen
 *   Add FileDB
 *
 *
 */
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include "FileDB.h"

using namespace std;

namespace FILEDB {
  /**
   * Implementation of Serialization Exception Class
   */
  SerializeException::SerializeException (const std::string& cls,
					  const std::string& reason)
    :exception(), cls_ (cls), reason_(reason)
  {
    // empty;
  }

  /**
   * Copy constructor
   */
  SerializeException::SerializeException (const SerializeException& exp)
    :cls_ (exp.cls_), reason_ (exp.reason_)
  {
    // empty
  }


  /**
   * Assignment operator
   */
  SerializeException &
  SerializeException::operator = (const SerializeException& exp)
  {
    if (this != &exp) {
      cls_ = exp.cls_;
      reason_ = exp.reason_;
    }
    return *this;
  }

  /**
   * Desstructor
   */
  SerializeException::~SerializeException (void) throw ()
  {
    // empty
  }

  /**
   * Return the reason of this exception
   */
  const char* 
  SerializeException::what (void) const throw ()
  {
    return (cls_ + ": " + reason_).c_str();
  }


  /**
   * Implementation of FileHashDB Exception Class
   */
  FileHashDBException::FileHashDBException (const std::string& cls,
					    const std::string& reason)
    :exception(), cls_ (cls), reason_(reason)
  {
    // empty;
  }

  /**
   * Copy constructor
   */
  FileHashDBException::FileHashDBException (const FileHashDBException& exp)
    :cls_ (exp.cls_), reason_ (exp.reason_)
  {
    // empty
  }


  /**
   * Assignment operator
   */
  FileHashDBException &
  FileHashDBException::operator = (const FileHashDBException& exp)
  {
    if (this != &exp) {
      cls_ = exp.cls_;
      reason_ = exp.reason_;
    }
    return *this;
  }

  /**
   * Desstructor
   */
  FileHashDBException::~FileHashDBException (void) throw ()
  {
    // empty
  }

  /**
   * Return the reason of this exception
   */
  const char* 
  FileHashDBException::what (void) const throw ()
  {
    return (cls_ + ": " + reason_).c_str();
  }



  double current_time (void)
  {
    struct timeval ct;

    gettimeofday (&ct, 0);

    return (double)(ct.tv_sec + ct.tv_usec/1000000.0);
  } 


  /**
   * this part of code is from GNUnet
   */
  static unsigned int crc_table[256];
  static int crc32_inited = 0;

  /**
   * This routine initialize the CRC table above
   */
  void crc32_init (void)
  {
    if (!crc32_inited) {
      unsigned int i, j;
      unsigned int  h = 1;
      crc_table[0] = 0;
      for (i = 128; i; i >>= 1) {
	h = (h >> 1) ^ ((h & 1) ? crc32poly : 0);
	/* h is now crc_table[i] */
	for (j = 0; j < 256; j += 2 * i)
	  crc_table[i + j] = crc_table[j] ^ h;
      }
      crc32_inited = 1;
    }
  }


  /**
   * Caculate crc32 checksum for a buffer with length len
   */
  unsigned int crc32 (unsigned int crc, const char* buf, unsigned int len)
  {
    crc ^= 0xffffffff;
    while (len--)
      crc = (crc >> 8) ^ crc_table[(crc ^ *buf++) & 0xff];
    return crc ^ 0xffffffff;
  }


  /**
   * Convert integer to string
   */
  std::string itostr (int val)
  {
    std::ostringstream o;
    if (!(o << val)) {
      cerr << "Integer " << val << " conversion to string error. " << endl;
      ::exit (1);
    }
    return o.str ();
  }

  /**
   * Byte Swapping Routines replacing ntohl etc
   */
#if defined(_FILEDB_BIG_ENDIAN)
  unsigned short ntoh_short (unsigned short sval)
  {
    return sval;
  }

  unsigned short hton_short (unsigned short sval)
  {
    return sval;
  }

  unsigned int   ntoh_int   (unsigned int ival)
  {
    return ival;
  }

  unsigned int   hton_int   (unsigned int ival)
  {
    return ival;
  }


  unsigned long  ntoh_long  (unsigned long lval)
  {
    return lval;
  }

  unsigned long  hton_long  (unsigned long lval)
  {
    return lval;
  }
#else

  /**
   * Several used macros
   */
#define __FILEDB_BSWAP_16(x) \
  (( ( (x) >> 8) & 0xffu) | ( ( (x) & 0xffu) << 8))

#define __FILEDB_BSWAP_32(x) \
  (( ( (x) & 0xff000000u) >> 24) | ( ( (x) & 0x00ff0000u) >> 8) | \
   ((  (x) & 0x0000ff00u) <<  8) | ( ( (x) & 0x000000ffu) << 24))

#define __FILEDB_BSWAP_64(x) \
  (( ( (x) & 0xff00000000000000ull) >> 56) \
   | (((x) & 0x00ff000000000000ull) >> 40) \
   | (((x) & 0x0000ff0000000000ull) >> 24) \
   | (((x) & 0x000000ff00000000ull) >> 8 ) \
   | (((x) & 0x00000000ff000000ull) << 8 ) \
   | (((x) & 0x0000000000ff0000ull) << 24) \
   | (((x) & 0x000000000000ff00ull) << 40) \
   | (((x) & 0x00000000000000ffull) << 56))

  unsigned short ntoh_short (unsigned short sval)
  {
    return __FILEDB_BSWAP_16(sval);
  }

  unsigned short hton_short (unsigned short sval)
  {
    return __FILEDB_BSWAP_16(sval);    
  }

  unsigned int   ntoh_int   (unsigned int ival)
  {
    return __FILEDB_BSWAP_32(ival);
  }

  unsigned int   hton_int   (unsigned int ival)
  {
    return __FILEDB_BSWAP_32(ival);
  }


  unsigned long  ntoh_long  (unsigned long lval)
  {
    if (sizeof(long) == sizeof(int))
      return __FILEDB_BSWAP_32(lval);
    else
      return __FILEDB_BSWAP_64(lval);
  }

  unsigned long  hton_long  (unsigned long lval)
  {
    if (sizeof(long) == sizeof(int))
      return __FILEDB_BSWAP_32(lval);
    else
      return __FILEDB_BSWAP_64(lval);
  }
#endif
}
