/*
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
 *     $Log: ffdb_db.h,v $
 *     Revision 1.3  2009-03-04 19:12:28  edwards
 *     Renamed DB_HASH and __db to avoid name collisions with Berkeley DB.
 *
 *     Revision 1.2  2009/03/02 23:58:21  chen
 *     Change implementation on keys iterator which get keys only
 *
 *     Revision 1.1  2009/02/20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#ifndef _FFDB_FILEDB_H
#define _FFDB_FILEDB_H

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

typedef unsigned int	pgno_t;
typedef unsigned short	indx_t;


/*
 * Typical return status
 */
#define FFDB_ERROR -1
#define FFDB_SUCCESS 0
#define FFDB_NOT_FOUND 1
#define FFDB_SPECIAL 2

/*
 *Little endian and big endien
 */
#ifndef BYTE_ORDER
#define	LITTLE_ENDIAN	1234		/* LSB first: i386, vax */
#define	BIG_ENDIAN	4321		/* MSB first: 68000, ibm, net */
#endif

/*
 * Some file open flags may not be available for some machines
 */
#ifndef O_EXLOCK
#define O_EXLOCK 0
#endif

#ifndef O_SHLOCK		       
#define	O_SHLOCK 0
#endif

/*
 * Errno definition for invalid file type
 */
#ifndef EFTYPE
#define EFTYPE EINVAL
#endif


/*
 * We only have DB_HASH in this package
 */
typedef enum {FFDB_HASH = 1} FFDB_DBTYPE;

/*
 * Still use the key data pair structure
 */
typedef struct _ffdb_dbt_{
  void *data;                          /* data                 */
  unsigned int size;                   /* data length in bytes */
}FFDB_DBT;

/*
 * DB access method and cursor operation values.  Each value is an operation
 * code to which additional bit flags are added.
 *
 * Most are not used
 */
#define	FFDB_AFTER		 1	/* Dbc.put */
#define	FFDB_APPEND		 2	/* Db.put */
#define	FFDB_BEFORE		 3	/* Dbc.put */
#define	FFDB_CONSUME		 4	/* Db.get */
#define	FFDB_CONSUME_WAIT	 5	/* Db.get */
#define	FFDB_CURRENT		 6	/* Dbc.get, Dbc.put, DbLogc.get */
#define	FFDB_FIRST		 7	/* Dbc.get, DbLogc->get */
#define	FFDB_GET_BOTH		 8	/* Db.get, Dbc.get */
#define	FFDB_GET_BOTHC		 9	/* Dbc.get (internal) */
#define	FFDB_GET_BOTH_RANGE	10	/* Db.get, Dbc.get */
#define	FFDB_GET_RECNO		11	/* Dbc.get */
#define	FFDB_JOIN_ITEM		12	/* Dbc.get; don't do primary lookup */
#define	FFDB_KEYFIRST		13	/* Dbc.put */
#define	FFDB_KEYLAST		14	/* Dbc.put */
#define	FFDB_LAST		15	/* Dbc.get, DbLogc->get */
#define	FFDB_NEXT		16	/* Dbc.get, DbLogc->get */
#define	FFDB_NEXT_DUP		17	/* Dbc.get */
#define	FFDB_NEXT_NODUP		18	/* Dbc.get */
#define	FFDB_NODUPDATA		19	/* Db.put, Dbc.put */
#define	FFDB_NOOVERWRITE	20	/* Db.put */
#define	FFDB_NOSYNC		21	/* Db.close */
#define	FFDB_POSITION		22	/* Dbc.dup */
#define	FFDB_PREV		23	/* Dbc.get, DbLogc->get */
#define	FFDB_PREV_DUP		24	/* Dbc.get */
#define	FFDB_PREV_NODUP		25	/* Dbc.get */
#define	FFDB_SET		26	/* Dbc.get, DbLogc->get */
#define	FFDB_SET_RANGE		27	/* Dbc.get */
#define	FFDB_SET_RECNO		28	/* Db.get, Dbc.get */
#define	FFDB_UPDATE_SECONDARY	29	/* Dbc.get, Dbc.del (internal) */
#define	FFDB_WRITECURSOR	30	/* Db.cursor */
#define	FFDB_WRITELOCK		31	/* Db.cursor (internal) */

/**
 * Two different cursor type: one traverse key the other traverse data
 */
#define FFDB_KEY_CURSOR  0x1001
#define FFDB_DATA_CURSOR 0x1004 

/**
 * Forward decleration of cursor
 */
typedef struct _ffdb_cursor_
{
  /* Get routine */
  /* If data is null (0), caller is not interested in data */
  int (*get) (struct _ffdb_cursor_ *c, 
	      FFDB_DBT* key,  FFDB_DBT* data, unsigned int flags);
  /* Close this cursor */
  int (*close)(struct _ffdb_cursor_ *c);

  /* type of this cursor */
  int type;

  /* internal pointer */
  void* internal;

}ffdb_cursor_t;




/* Access method description structure. */
typedef struct __ffdb {
  FFDB_DBTYPE   type;			/* Underlying db type. */
  int (*close)	(struct __ffdb *);
  int (*del)	(const struct __ffdb *, const FFDB_DBT *, unsigned int);
  int (*get)	(const struct __ffdb *, const FFDB_DBT *, FFDB_DBT *, unsigned int);
  int (*put)	(const struct __ffdb *, FFDB_DBT *, const FFDB_DBT *, unsigned int);
  int (*sync)	(const struct __ffdb *, unsigned int);
  int (*cursor) (const struct __ffdb *, ffdb_cursor_t **, unsigned int type);
  void *internal;			/* Access method private. */
  int (*fd)	(const struct __ffdb *);
} FFDB_DB;

/*
 * Hash database magic number and version
 */
#define FFDB_HASHMAGIC 0xcece3434
#define FFDB_HASHVERSION 5

/*
 * How do we store key and data on a page
 * 1) key and data try to be on the primary page
 * 2) key points to pageno and offset where data are
 */
#define FFDB_STORE_EMBED    0x00ffddee
#define FFDB_STORE_INDIRECT 0x00ff1100
 
/*
 * Structure used to pass parameters to the hashing routines. 
 */
typedef struct {
  unsigned int	bsize;		 /* bucket size */
  unsigned int	nbuckets;	 /* number of buckets */
  unsigned long	cachesize;	 /* bytes to cache */
  int           rearrangepages;  /* to rearrange page on open/close to save
				  * space
				  */
  unsigned int   userinfolen;    /* how many bytes for user information */
  unsigned int   numconfigs;     /* number of configurations */
  unsigned int  (*hash) (const void *, unsigned int); /* hash function */
                                /* key compare func */
  int           (*cmp) (const FFDB_DBT *, const FFDB_DBT *); 
} FFDB_HASHINFO;


/*
 * Internal byte swapping code if we are using little endian
 */
/*
 * Little endian <==> big endian 32-bit swap macros.
 *	M_32_SWAP	swap a memory location
 *	P_32_SWAP	swap a referenced memory location
 *	P_32_COPY	swap from one location to another
 */
#define	M_32_SWAP(a) {							\
	unsigned int _tmp = a;						\
	((char *)&a)[0] = ((char *)&_tmp)[3];				\
	((char *)&a)[1] = ((char *)&_tmp)[2];				\
	((char *)&a)[2] = ((char *)&_tmp)[1];				\
	((char *)&a)[3] = ((char *)&_tmp)[0];				\
}
#define	P_32_SWAP(a) {							\
	unsigned int_tmp = *(unsigned int *)a;				\
	((char *)a)[0] = ((char *)&_tmp)[3];				\
	((char *)a)[1] = ((char *)&_tmp)[2];				\
	((char *)a)[2] = ((char *)&_tmp)[1];				\
	((char *)a)[3] = ((char *)&_tmp)[0];				\
}
#define	P_32_COPY(a, b) {						\
	((char *)&(b))[0] = ((char *)&(a))[3];				\
	((char *)&(b))[1] = ((char *)&(a))[2];				\
	((char *)&(b))[2] = ((char *)&(a))[1];				\
	((char *)&(b))[3] = ((char *)&(a))[0];				\
}

/*
 * Little endian <==> big endian 16-bit swap macros.
 *	M_16_SWAP	swap a memory location
 *	P_16_SWAP	swap a referenced memory location
 *	P_16_COPY	swap from one location to another
 */
#define	M_16_SWAP(a) {							\
	unsigned short _tmp = a;					\
	((char *)&a)[0] = ((char *)&_tmp)[1];				\
	((char *)&a)[1] = ((char *)&_tmp)[0];				\
}
#define	P_16_SWAP(a) {							\
         unsigned short _tmp = *(unsigned short *)a;			\
	((char *)a)[0] = ((char *)&_tmp)[1];				\
	((char *)a)[1] = ((char *)&_tmp)[0];				\
}
#define	P_16_COPY(a, b) {						\
	((char *)&(b))[0] = ((char *)&(a))[1];				\
	((char *)&(b))[1] = ((char *)&(a))[0];				\
}


#define FFDB_DEFAULT_UINFO_LEN 4000
/**
 * The file contains user provided information right after
 * the header page
 */
typedef struct _user_info_ 
{
  unsigned int len;
  unsigned char* uinfo;
}ffdb_user_info_t;

/**
 * The file contains configuration information right after the above
 * user information 
 */
typedef struct _config_info_
{
  int config;                   /* configuration number          */
  int index;                    /* index into all configurations */
  int inserted;                 /* configuration inserted        */
  int type;                     /* type of configuration (fixed) */
  int mtime;                    /* modified time of this config  */
#define _FFDB_MAX_FNAME 128
  char fname[_FFDB_MAX_FNAME];
}ffdb_config_info_t;


/**
 * All configuration information 
 */
typedef struct _all_config_info_
{
  int numconfigs;
  ffdb_config_info_t *allconfigs;
}ffdb_all_config_info_t;


#ifdef __cplusplus
extern "C"
{
#endif
/*
 * Open a database handle
 * @param fname database filename
 * @param flags database open flags
 * @param mode typical file onwership mode
 * @openinfo user supplied information for opening a database
 * @return a pointer to FFDB_DB structure. return 0 if something wrong
 */
extern FFDB_DB*
ffdb_dbopen (const char* fname, int flags, int mode, const void* openinfo);


/**
 * Set a paticular configuration information
 * 
 * @param db pointer to underlying database
 * @param config a configuration structure to be set
 *
 * @return 0 on success. -1 on failure with a proper errno set
 */
extern int
ffdb_set_config (FFDB_DB* db, ffdb_config_info_t* config);

/**
 * Get a paticular configuration information
 *
 * @param db pointer to underlying database
 * @param confignum the configuration number
 * @param config retrieved configuration information will be here
 *
 * @return 0 on success, -1 on failure with a proper errno set
 */
extern int
ffdb_get_config (const FFDB_DB* db, unsigned int confignum,
		 ffdb_config_info_t* config);   



/**
 * Set all configurations
 *
 * @param db pointer to underlying database
 * @param configs all configuration information
 *
 * @return 0 on success -1 on failure with a proper errno set
 */
extern int
ffdb_set_all_configs (FFDB_DB* db, ffdb_all_config_info_t* configs);



/**
 * Get all configurations
 * caller should free memory of configs->allconfigs
 *
 * @param db pointer to underlying database
 * @param configs all configuration information
 *
 * @return 0 on success -1 on failure with a proper errno set
 */
extern int
ffdb_get_all_configs (const FFDB_DB* db, ffdb_all_config_info_t* configs);


/**
 * Get number of configurations information
 *
 * @param db pointer to underlying database
 *
 * @return number of configurations allocated
 */
extern unsigned int
ffdb_num_configs (const FFDB_DB* db);


/**
 * Set user information for the database
 *
 * @param db pointer to underlying database
 * @param data user data
 * @param len user data len
 *
 * @return 0 on success. -1 on failure with a proper errno
 */
extern int
ffdb_set_user_info (FFDB_DB* db, unsigned char* data, unsigned int len);

/**
 * Get user information for the database
 *
 * @param db pointer to underlying database
 * @param data user data
 * @param len user data len. Caller allocate space for data and pass 
 * initial data length. On return, the actual data length will be stored
 * in len.
 *
 * @return 0 on success. -1 on failure with a proper errno
 * 
 */
extern int
ffdb_get_user_info (const FFDB_DB* db, unsigned char data[], 
		    unsigned int* len);


/**
 * Get maximum user information length in bytes allocated
 *
 * @param db pointer to underlying database
 * @return number of bytes allocated for user information
 */
extern unsigned int
ffdb_max_user_info_len (const FFDB_DB* db);


/*
 * A routine which reset the database handle under panic mode
 */
extern void ffdb_dbpanic(FFDB_DB* dbp);

#ifdef __cplusplus
};
#endif

#endif
