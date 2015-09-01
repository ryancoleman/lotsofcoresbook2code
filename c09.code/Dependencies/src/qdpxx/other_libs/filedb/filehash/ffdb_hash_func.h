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
 *     
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_hash_func.h,v $
 *     Revision 1.1  2009-02-20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#ifndef _FFDB_HASH_FUNC_H
#define _FFDB_HASH_FUNC_H

extern unsigned int __ham_func2(const void* key, unsigned int len);
extern unsigned int __ham_func3(const void* key, unsigned int len);
extern unsigned int __ham_func4(const void* key, unsigned int len);
extern unsigned int __ham_func5(const void* key, unsigned int len);
extern unsigned int __ffdb_log2(unsigned int num);
extern int          __ham_defcmp(const FFDB_DBT* a, const FFDB_DBT* b);

extern unsigned int (*__ffdb_default_hash)(const void* key, unsigned int len);
extern int (*__ffdb_default_cmp)(const FFDB_DBT *a, const FFDB_DBT *b);

extern void __ffdb_crc32_init (void);
extern unsigned int  __ffdb_crc32_checksum (unsigned int crc, 
					    const unsigned char* buffer,
					    unsigned int len);
#endif
