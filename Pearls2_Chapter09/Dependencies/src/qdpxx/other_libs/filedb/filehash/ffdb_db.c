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
 *     Pure File Based Hash Database
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_db.c,v $
 *     Revision 1.1  2009-02-20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>

#include "ffdb_db.h"
#include "ffdb_pagepool.h"

/**
 * Decleration of hash open function
 */
extern FFDB_DB*
__ffdb_hash_open (const char* fname, int flags, int node,
		  const FFDB_HASHINFO* info, int dflags);

FFDB_DB *
ffdb_dbopen(const char* fname, int flags, int mode, const void* openinfo)
{
#define	USE_OPEN_FLAGS							\
  (O_CREAT | O_EXCL | O_EXLOCK | O_NONBLOCK | O_RDONLY |		\
   O_RDWR | O_SHLOCK | O_TRUNC)

  /* check whether a wrong flag is set */
  if (FFDB_FLAG_ISSET(flags, (~USE_OPEN_FLAGS))) {
    errno = EINVAL;
    return 0;
  }
  return __ffdb_hash_open(fname, flags & USE_OPEN_FLAGS, mode,
			  openinfo, flags);
}

static int
__ffdb_dberr (void)
{
  return (FFDB_ERROR);
}

/*
 * __DBPANIC -- Stop.
 *
 * Parameters:
 *	dbp:	pointer to the DB structure.
 */
void
ffdb_dbpanic(FFDB_DB* dbp)
{
  /* The only thing that can succeed is a close. */
  dbp->del = (int (*)())__ffdb_dberr;
  dbp->fd = (int (*)())__ffdb_dberr;
  dbp->get = (int (*)())__ffdb_dberr;
  dbp->put = (int (*)())__ffdb_dberr;
  dbp->cursor = (int (*)())__ffdb_dberr;
  dbp->sync = (int (*)())__ffdb_dberr;
}

