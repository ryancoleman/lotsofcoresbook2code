// -*- C++ -*-
/*
 **************************************************************************
 *                                                                        *
 *          General Purpose Hash Function Algorithms Library              *
 *                                                                        *
 * Author: Arash Partow - 2002                                            *
 * URL: http://www.partow.net                                             *
 * URL: http://www.partow.net/programming/hashfunctions/index.html        *
 *                                                                        *
 * Copyright notice:                                                      *
 * Free use of the General Purpose Hash Function Algorithms Library is    *
 * permitted under the guidelines and in accordance with the most current *
 * version of the Common Public License.                                  *
 * http://www.opensource.org/licenses/cpl.php                             *
 *                                                                        *
 **************************************************************************
*/



#ifndef _FILEDB_INCLUDE_GENERALHASHFUNCTION_C_H
#define _FILEDB_INCLUDE_GENERALHASHFUNCTION_C_H


#include <cstdio>

namespace FILEDB {

  typedef unsigned int (*hash_func) (char* str, unsigned int len);

  unsigned int RSHash  (char* str, unsigned int len);
  unsigned int JSHash  (char* str, unsigned int len);
  unsigned int PJWHash (char* str, unsigned int len);
  unsigned int ELFHash (char* str, unsigned int len);
  unsigned int BKDRHash(char* str, unsigned int len);
  unsigned int SDBMHash(char* str, unsigned int len);
  unsigned int DJBHash (char* str, unsigned int len);
  unsigned int DEKHash (char* str, unsigned int len);
  unsigned int BPHash  (char* str, unsigned int len);
  unsigned int FNVHash (char* str, unsigned int len);
  unsigned int APHash  (char* str, unsigned int len);
  unsigned int db_hash_func (char* str, unsigned int len);

  const int NUM_HASHS=12;
  const int HASH_BUCKETS=53027;
};
			     
#endif
