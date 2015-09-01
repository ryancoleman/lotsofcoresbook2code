/**
 * Copyright (C) <2008> Jefferson Science Associates, LLC
 *                      Under U.S. DOE Contract No. DE-AC05-06OR23177
 *
 *                      Thomas Jefferson National Accelerator Facility
 *
 *                      Jefferson Lab
 *                      Scientific Computing Group,
 *                      12000 Jefferson Ave.,      
 *                      Newport News, VA 23606 
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ----------------------------------------------------------------------------
 * Description:
 *     Hash functions used in the hash based database
 *     This is taken from Berkeley DB version 4
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_hash_func.c,v $
 *     Revision 1.1  2009-02-20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
/*-
 * See the file LICENSE for redistribution information.
 *
 * Copyright (c) 1996,2007 Oracle.  All rights reserved.
 */
/*
 * Copyright (c) 1990, 1993
 *	Margo Seltzer.  All rights reserved.
 */
/*
 * Copyright (c) 1990, 1993
 *	The Regents of the University of California.  All rights reserved.
 *
 * This code is derived from software contributed to Berkeley by
 * Margo Seltzer.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * $Id: ffdb_hash_func.c,v 1.1 2009-02-20 20:44:47 chen Exp $
 */

/*
 * __ham_func2 --
 *	Phong Vo's linear congruential hash.
 *
 * PUBLIC: u_int32_t __ham_func2 __P((DB *, const void *, u_int32_t));
 */
#include <stdio.h>
#include "ffdb_db.h"
#include "ffdb_hash_func.h"

#define	DCHARHASH(h, c)	((h) = 0x63c63cd9*(h) + 0x9c39c33d + (c))

unsigned int
__ham_func2(key, len)
     const void *key;
     unsigned int len;
{
  const unsigned char *e, *k;
  unsigned int h;
  unsigned char c;

  k = key;
  e = k + len;
  for (h = 0; k != e;) {
    c = *k++;
    if (!c && k > e)
      break;
    DCHARHASH(h, c);
  }
  return (h);
}

/*
 * __ham_func3 --
 *	Ozan Yigit's original sdbm hash.
 *
 * Ugly, but fast.  Break the string up into 8 byte units.  On the first time
 * through the loop get the "leftover bytes" (strlen % 8).  On every other
 * iteration, perform 8 HASHC's so we handle all 8 bytes.  Essentially, this
 * saves us 7 cmp & branch instructions.
 *
 * PUBLIC: unsigned int __ham_func3 __P((DB *, const void *, unsigned int));
 */
unsigned int
__ham_func3(key, len)
     const void *key;
     unsigned int len;
{
  const unsigned char *k;
  unsigned int n, loop;

  if (len == 0)
    return (0);

#define	HASHC	n = *k++ + 65599 * n
  n = 0;
  k = key;

  loop = (len + 8 - 1) >> 3;
  switch (len & (8 - 1)) {
  case 0:
    do {
      HASHC;
    case 7:
      HASHC;
    case 6:
      HASHC;
    case 5:
      HASHC;
    case 4:
      HASHC;
    case 3:
      HASHC;
    case 2:
      HASHC;
    case 1:
      HASHC;
    } while (--loop);
  }
  return (n);
}

/*
 * __ham_func4 --
 *	Chris Torek's hash function.  Although this function performs only
 *	slightly worse than __ham_func5 on strings, it performs horribly on
 *	numbers.
 *
 * PUBLIC: unsigned int __ham_func4 __P((DB *, const void *, unsigned int));
 */
unsigned int
__ham_func4(key, len)
     const void *key;
     unsigned int len;
{
  const unsigned char *k;
  unsigned int h, loop;

  if (len == 0)
    return (0);

#define	HASH4a	h = (h << 5) - h + *k++;
#define	HASH4b	h = (h << 5) + h + *k++;
#define	HASH4	HASH4b
  h = 0;
  k = key;
  
  loop = (len + 8 - 1) >> 3;
  switch (len & (8 - 1)) {
  case 0:
    do {
      HASH4;
    case 7:
      HASH4;
    case 6:
      HASH4;
    case 5:
      HASH4;
    case 4:
      HASH4;
    case 3:
      HASH4;
    case 2:
      HASH4;
    case 1:
      HASH4;
    } while (--loop);
  }
  return (h);
}

/*
 * Fowler/Noll/Vo hash
 *
 * The basis of the hash algorithm was taken from an idea sent by email to the
 * IEEE Posix P1003.2 mailing list from Phong Vo (kpv@research.att.com) and
 * Glenn Fowler (gsf@research.att.com).  Landon Curt Noll (chongo@toad.com)
 * later improved on their algorithm.
 *
 * The magic is in the interesting relationship between the special prime
 * 16777619 (2^24 + 403) and 2^32 and 2^8.
 *
 * This hash produces the fewest collisions of any function that we've seen so
 * far, and works well on both numbers and strings.
 *
 * PUBLIC: unsigned int __ham_func5 __P((DB *, const void *, unsigned int));
 */
unsigned int
__ham_func5(key, len)
     const void *key;
     unsigned int len;
{
  const unsigned char *k, *e;
  unsigned int h;

  
  k = key;
  e = k + len;
  for (h = 0; k < e; ++k) {
    h *= 16777619;
    h ^= *k;
  }
  return (h);
}

/*
 * __ham_test --
 *
 * PUBLIC: unsigned int __ham_test __P((DB *, const void *, unsigned int));
 */
unsigned int
__ham_test(key, len)
     const void *key;
     unsigned int len;
{
  return ((unsigned int)*(char *)key);
}


/**
 * A simple implementation of log2 on an integer
 */
unsigned int
__ffdb_log2 (unsigned int num)
{
  unsigned int i, limit;

  limit = 1;
  for (i = 0; limit < num; limit = limit << 1, i++)
    ;
  return i;
}

/**
 * Default key compare function
 */
int
__ham_defcmp(const FFDB_DBT* a, const FFDB_DBT* b)
{
  size_t len;
  unsigned char *p1, *p2;

  /*
   * Returns:
   *	< 0 if a is < b
   *	= 0 if a is = b
   *	> 0 if a is > b
   *
   * XXX
   * If a size_t doesn't fit into a long, or if the difference between
   * any two characters doesn't fit into an int, this routine can lose.
   * What we need is a signed integral type that's guaranteed to be at
   * least as large as a size_t, and there is no such thing.
   */
  len = a->size > b->size ? b->size : a->size;
  for (p1 = a->data, p2 = b->data; len--; ++p1, ++p2)
    if (*p1 != *p2)
      return ((long)*p1 - (long)*p2);
  return ((long)a->size - (long)b->size);
}

/**
 * Setting default hash function
 */
unsigned int (*__ffdb_default_hash)(const void* key, unsigned int len) = __ham_func5;

/**
 * Setting default compare function
 */
int (*__ffdb_default_cmp)(const FFDB_DBT *a, const FFDB_DBT *b)= __ham_defcmp;




/**
 * CRC32 Checksum routine
 */

/**
 * this part of code is from GNUnet
 */
static unsigned int _ffdb_crc_table[256];
static int _ffdb_crc32_inited = 0;

#define _CRC32POLY (unsigned int)0xedb88320

/**
 * This routine initialize the CRC table above
 */
void __ffdb_crc32_init (void)
{
  if (!_ffdb_crc32_inited) {
    unsigned int i, j;
    unsigned int  h = 1;
    _ffdb_crc_table[0] = 0;
    for (i = 128; i; i >>= 1) {
      h = (h >> 1) ^ ((h & 1) ? _CRC32POLY : 0);
      /* h is now crc_table[i] */
      for (j = 0; j < 256; j += 2 * i)
	_ffdb_crc_table[i + j] = _ffdb_crc_table[j] ^ h;
    }
    _ffdb_crc32_inited = 1;
  }
}


  /**
   * Caculate crc32 checksum for a buffer with length len
   */
unsigned int 
__ffdb_crc32_checksum (unsigned int crc, const unsigned char* buf, 
		       unsigned int len)
{
  crc ^= 0xffffffff;
  while (len--)
    crc = (crc >> 8) ^ _ffdb_crc_table[(crc ^ *buf++) & 0xff];
  return crc ^ 0xffffffff;
}








