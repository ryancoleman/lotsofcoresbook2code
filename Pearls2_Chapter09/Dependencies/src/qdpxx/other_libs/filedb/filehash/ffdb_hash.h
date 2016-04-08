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
 *     Hash  based  database implementation
 *     This implementation is based on Berkeley DB and 
 *     P. Larson "Dynamic Hash Tables" Communication of ACM, 1988 Vol 31, Num 4
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_hash.h,v $
 *     Revision 1.2  2009-04-21 18:51:20  chen
 *     Fix bugs related to number of pages upon moving pages in addition to clean pages on disk when the pages has been moved
 *
 *     Revision 1.1  2009/02/20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#ifndef _FFDB_HASH_H
#define _FFDB_HASH_H

#include "ffdb_db.h"
#include "ffdb_pagepool.h"

/* Operations */
typedef enum {
  HASH_GET, HASH_PUT, HASH_PUTNEW, HASH_DELETE, HASH_FIRST, HASH_NEXT
} ACTION;

/**
 * forward decleration
 */
struct _ffdb_crs_;
typedef struct _ffdb_crs_ ffdb_crs_t;

/**
 * Hash Table Information 
 * This is stored in the first page of a hash based file
 *
 */
typedef struct hashhdr {	/* Disk resident portion */
  int	magic;		        /* Magic NO for hash tables */
  int	version;	        /* Version ID */
  int   lorder;                 /* byte order stored in the database */
  int	bsize;		        /* Bucket/Page Size */
  int	bshift;		        /* Bucket shift */
  unsigned int ovfl_point;	/* Where overflow pages are being allocated */
  unsigned int max_bucket;      /* ID of Maximum bucket in use */
  unsigned int high_mask;	/* Mask to modulo into entire table */
  unsigned int low_mask;	/* Mask to modulo into lower half of table */
  unsigned int ffactor;	        /* Fill factor */
  unsigned int nkeys;		/* Number of keys in hash table */
  unsigned int hdrpages;	/* Size of table header */
  unsigned int uinfolen;        /* user information length in bytes */
  unsigned short   uinfo_page;  /* User information start at this page */
  unsigned short   uinfo_npages;/* User information number of pages    */
  unsigned short   cfig_page;   /* config information start at this page */
  unsigned short   cfig_npages; /* config information number of pages    */
  unsigned int     num_cfigs;   /* number of configurations              */
  unsigned int	h_charkey;      /* value of hash(CHARKEY) */
  unsigned int  num_moved_pages;/* number of moved pages at the last level */
#define NCACHED	32		/* number of spare points */
  pgno_t spares[NCACHED];       /* indicating starting page number at this 
				 * spliting stage
				 */
  pgno_t free_pages[NCACHED];   /* pages containing free pages information */
  unsigned int chksum;          /* header crc checksum value */
} ffdb_hashhdr_t;


/**
 * Hash table definition
 */
typedef struct htab {		/* Memory resident data structure */
  FFDB_TAILQ_HEAD(_ffdb_cursor_queue, _ffdb_crs_) curs_queue;
  ffdb_hashhdr_t hdr;		/* Header */
  unsigned int (*hash) (const void *, unsigned int); /* Hash Function */
  int      (*h_compare) (const FFDB_DBT *, const FFDB_DBT *);
  int   mborder;		/* machine byte order */
  int	flags;	                /* Flag values */
  int	fp;		        /* File pointer */
  char *fname;        	        /* File path */
  char *bigdata_buf;	        /* Temporary Buffer for BIG data */
  int bigdata_len; 	        /* Length of bigdata_buf */
  char *bigkey_buf;	        /* Temporary Buffer for BIG keys */
  int bigkey_len;	        /* Length of bigkey_buf */
  unsigned short  *split_buf;	/* Temporary buffer for splits */
  ffdb_crs_t      *seq_cursor;	/* Cursor used for hash_seq */
  int	db_errno;               /* Error Number -- for DBM compatability */
  int	new_file;	        /* Indicates if fd is backing store or no */
  int	save_file;	        /* Indicates whether we need to flush file at
				 * exit */
  pgno_t curr_dpage;            /* current data page number */
  int   rearrange_pages;        /* rearrange pages to save disk space */
  ffdb_pagepool_t *mp;		/* mpool for buffer management */
  pthread_mutex_t lock;		/* lock */
} ffdb_htab_t;


/* hash item for key and data pointer */
typedef struct hash_ent_ {
  pgno_t 		pgno;                  /* page number      */
  void*                 pagep;                 /* memory of page   */
  pgno_t		bucket;                /* bucket number    */
  indx_t		pgndx;                 /* index within the page */
  int    	        status;                /* data found or not */
  int    	        seek_size;             /* combined size */
  pgno_t		key_off;               /* key offset    */
  pgno_t                key_len;               /* key length    */
  pgno_t		data_off;              /* data offset   */
  unsigned int          data_chksum;           /* data checksum */
  unsigned int   	caused_expand;         /* cause expand  */
} ffdb_hent_t;


#define	ITEM_ERROR	-1
#define ITEM_CLEAN      0
#define	ITEM_OK		1
#define	ITEM_NO_MORE	2

#define	UNKNOWN		0xffffffff		/* for num_items */
#define	NO_EXPAND	0xfffffffe 


/**
 * Cursor of database supporting only get for now
 */
struct _ffdb_crs_
{
  /* Entry in the cursor queue managed by hash table */
  FFDB_TAILQ_ENTRY(_ffdb_crs_) queue;
  /* hash table pointer */
  ffdb_htab_t* hashp;
  /* pointer to public cursor */
  ffdb_cursor_t *pcursor;
  /* Current key or data depending type of the cursor */
  ffdb_hent_t item;
  /* internal lock for the cursor */
  pthread_mutex_t lock;	
};


/**
 * Constants
 */
#define	MAX_BSIZE		262144		/* 2^18 */
#define MIN_BUFFERS		6
#define MINHDRSIZE		512
#define DEF_CACHESIZE	        134217728       /* 2^27 default cache */
#define DEF_BUCKET_SIZE		4096
#define DEF_BUCKET_SHIFT	12		/* log2(BUCKET) */
#define DEF_SEGSIZE		256
#define DEF_SEGSIZE_SHIFT	8		/* log2(SEGSIZE)	 */
#define DEF_DIRSIZE		256
#define DEF_FFACTOR		65536
#define MIN_FFACTOR		4
#define SPLTMAX			8
#define CHARKEY			"%$sniglet^&"
#define NUMKEY			1038583
#define BYTE_SHIFT		3
#define INT32_T_TO_BYTE		2    /* Convert int size to byte << 2  = x4 */
#define INT32_T_BYTE_SHIFT	5    /* number of bits >> 5 (2^5) = 32 bits   
                                      * There are total of 32 level of
				      * splits */
#define ALL_SET			((u_int32_t)0xFFFFFFFF)
#define ALL_CLEAR		0

#define PTROF(X)	((BUFHEAD *)((ptr_t)(X)&~0x3))
#define ISMOD(X)	((ptr_t)(X)&0x1)
#define DOMOD(X)	((X) = (int8_t *)((ptr_t)(X)|0x1))
#define ISDISK(X)	((ptr_t)(X)&0x2)
#define DODISK(X)	((X) = (int8_t *)((ptr_t)(X)|0x2))

#define BITS_PER_MAP	32

/* Given the address of the beginning of a big map, clear/set the nth bit */
#define CLRBIT(A, N)	((A)[(N)/BITS_PER_MAP] &= ~(1<<((N)%BITS_PER_MAP)))
#define SETBIT(A, N)	((A)[(N)/BITS_PER_MAP] |= (1<<((N)%BITS_PER_MAP)))
#define ISSET(A, N)	((A)[(N)/BITS_PER_MAP] & (1<<((N)%BITS_PER_MAP)))

/* Overflow management */
/**
 * Overflow page numbers are allocated per split point.  At each doubling of
 * the table, we can allocate extra pages.  So, an overflow page number has
 * the top 5 bits indicate which split point and the lower 11 bits indicate
 * which page at that split point is indicated (pages within split points are
 * numberered starting with 1).
 */

#define SPLITSHIFT	11
#define SPLITMASK	0x7FF /* 11111111111  11 bit of all 1 */
#define SPLITNUM(N)	(((unsigned int)(N)) >> SPLITSHIFT)
#define OPAGENUM(N)	((N) & SPLITMASK)
#define	OADDR_OF(S,O)	((unsigned int)((unsigned int)(S) << SPLITSHIFT) + (O))

#define POW2(N)  (1 << (N))

/** 
 * Bucket address to page number: bucket address is obtained
 * from __call_hash code.
 * The bucket number + number of header page + number of over flow page
 */
#if 0
#define BUCKET_TO_PAGE(B) \
	((B) + hashp->hdr.spares[__ffdb_log2(B + 1)])
#endif

#if 1
#define BUCKET_TO_PAGE(B,page)					 \
  do {                                                           \
    unsigned int slevel = __ffdb_log2((B) + 1);			 \
    (page) = ((B) > 0) ? ((B) - POW2(slevel - 1) + hashp->hdr.spares[slevel]) : (hashp->hdr.spares[slevel]); \
  }while(0)
#endif

#define MAX_PAGES(H) (0xFFFFFFFF)

/* Shorthands for accessing structure */
#define METADATA_PGNO 0
#define SPLIT_PGNO 0xFFFF

/**
 * One of the most important function: hash function
 * return bucket number
 */
extern unsigned int _ffdb_call_hash (ffdb_htab_t* hashp, const void* k, 
				     unsigned int len);

/**
 * Get a new page
 * This routine is always called before get_page
 *
 * @param hashp the hash table pointer
 * @addr  address for this page: can be a bucket address or others
 * @addrtype can be any of the above page type
 *
 * @return 0 on success. return -1 otherwise
 */
extern int ffdb_new_page (ffdb_htab_t* hashp, pgno_t addr, 
			  unsigned int addrtype);

/**
 * Get a page and return page address
 *
 * @param hashp the hash table pointer
 * @addr  address for this page: can be a bucket address or others
 * @addrtype can be any of the above page type
 * @page returned allocated page number
 *
 * @return memory address of this page on success. Otherwise NULL (0) returned
 */
extern void* ffdb_get_page (ffdb_htab_t* hashp, pgno_t addr, 
			    unsigned int addrtype, unsigned int flags,
			    pgno_t* page);


/**
 * Put a page back to pagepool
 *
 * @param hashp the hash table pointer
 * @mem   memory address of this page
 * @addrtype can be any of the above page type
 * @dirty is this page dirty
 *
 * @return 0 on success. Otherwise failure
 */
extern int ffdb_put_page (ffdb_htab_t* hashp, void* mem, 
			  unsigned int addrtype, int dirty);


/**
 * Release item and the related page associated with this item
 *
 * This routine is used when the item is no longer needed, which
 * usually happened when hash_get did not get anything
 *
 * @param hashp the hash table pointer
 * @param item hash entry item
 *
 * @return 0 on success. -1 on failure
 */
extern int ffdb_release_item (ffdb_htab_t* hashp,
			      ffdb_hent_t* item);


/**
 * Find a page where a key and data pair should reside (or already reside)
 *
 * @param hashp the hash table pointer
 * @param key   the key
 * @param val   the data
 * @param item  the information for the hash entry
 * 
 * @return 0 success. check item->status to see whether the item is here. 
 * return -1 if there are some failure.
 */
extern int ffdb_find_item (ffdb_htab_t* hashp,
			   FFDB_DBT* key, FFDB_DBT* val,
			   ffdb_hent_t* item);


/**
 * Get item from database. The item contains page and index 
 * information obtained from ffdb_find_item call
 *
 * @param hashp the hash table pointer
 * @param key   the key
 * @param val   the data pointer
 * @param item  the information for the hash entry
 * @param freepage whether to put page back pointed by item
 *
 * @return return 0 (always since we know the key is there). return -1 otherwise
 */
extern int ffdb_get_item (ffdb_htab_t* hashp,
			  const FFDB_DBT* key, FFDB_DBT* val,
			  ffdb_hent_t* item, int freepage);

/**
 * Add a pair of key and data into the hash database
 *
 * @param hashp the hash table pointer
 * @param key the pointer to the key
 * @param val pointer to the data
 * @param item the information for this pair of hash entry
 * @param replace to replace existing one or just add new one
 *
 * @return 0 on success. return -1 on failure
 */
extern int ffdb_add_pair (ffdb_htab_t* hashp,
			  FFDB_DBT* key, const FFDB_DBT* val,
			  ffdb_hent_t* item,
			  int replace);

/**
 * Add an over flow page to a page identified by item->pgno. The reason of
 * doing this because a new item cannot fit into the page designated by hash
 * bucket directly.
 *
 * @param hashp the hash table pointer
 * @param key the pointer to the key
 * @param val pointer to the data
 * @param item the information for this pair of hash entry
 *
 * @return 0 on success. return -1 on failure
 */
extern int ffdb_add_ovflpage (ffdb_htab_t* hashp, 
			      FFDB_DBT* key, const FFDB_DBT* val, 
			      ffdb_hent_t* item);


/**
 * Split a bucket: this happens when a bucket is full. This bucket may not be 
 * splitted right away (overflow pages needed), but it will eventually 
 * will be splitted
 * 
 * @param hashp the pointer to hash table
 * @param oldbucket the bucket number for the old bucket to be split
 * @param newbucket the bucket number for the new bucket
 * @param isdoubling whether this is doubling level
 * 
 * @return returns 0 on success, -1 on failure
 *
 */
extern int ffdb_split_bucket (ffdb_htab_t* hashp, unsigned int oldbucket,
			      unsigned int newbucket, int isdoubling);



/**
 * Rearrange pages when we finished writing the database to fill empty
 * pages left by not used hash buckets at the last level
 *
 * @param hashp the pointer to hash table
 *
 * @return always return 0
 */
extern int ffdb_rearrage_pages_on_close (ffdb_htab_t* hashp);



/**
 * Rearrange pages when we open the database to move data pages
 * that are using the pages for regular bucket pages
 *
 * @param hashp the pointer to hash table
 *
 * @return always return 0
 */
extern int ffdb_rearrage_pages_on_open (ffdb_htab_t* hashp);


/**
 * Reduce file size if possible. For example, the first time
 * a file was not doing rearrange pages and second time this
 * file is doing page rearrange
 *
 * @param hashp the pointer to hash table
 */
extern void ffdb_reduce_filesize (ffdb_htab_t* hashp);


/**
 * Find the last data page from a given position
 *
 * @param hashp the usual hash table pointer
 * @param start we are looking for the data page from start backward
 *
 * @return last data page number
 */
extern pgno_t ffdb_last_data_page (ffdb_htab_t* hashp, pgno_t start);


/**
 * Set configuration information
 * 
 * @param hashp the pointer to the hash table
 * @paran configs all configuration information
 *
 * @return 0 on success -1 on failure with a proper errno set
 */
extern int ffdb_set_configs (ffdb_htab_t* hashp, ffdb_all_config_info_t* configs);

/**
 * Get configuration information
 * 
 * @param hashp the pointer to the hash table
 * @paran configs all configuration information
 *
 * @return 0 on success -1 on failure with a proper errno set
 */
extern int ffdb_get_configs (ffdb_htab_t* hashp, ffdb_all_config_info_t* configs);


/**
 * Set a particular configuration information
 * @param hashp the pointer to the hash table
 * @paran config a configuration information
 *
 * @return 0 on success -1 on failure with a proper errno set 
 */
extern int ffdb_set_config_info (ffdb_htab_t* hashp, ffdb_config_info_t* config);

/**
 * Get a particular configuration information
 * @param hashp the pointer to the hash table
 * @paran confignum a particuler config number
 * @paran config a configuration information
 *
 * @return 0 on success -1 on failure with a proper errno set 
 */
extern int ffdb_get_config_info (ffdb_htab_t* hashp, 
				 unsigned int confignum, 
				 ffdb_config_info_t* config);

/**
 * Set user information
 * @param hashp the pointer to the hash table
 * @param data caller provided data
 * @param len length of caller provided data
 *
 * @return return 0 on success and -1 otherwise
 */
extern int ffdb_set_uinfo (ffdb_htab_t* hashp, unsigned char* data,
			   unsigned int len);


/**
 * Get user information
 * @param hashp the pointer to the hash table
 * @param data caller provided data
 * @param len length of caller provided data space
 *
 * @return return 0 on success and -1 otherwise
 */
extern int ffdb_get_uinfo (ffdb_htab_t* hashp, unsigned char data[],
			   unsigned int *len);


/**
 * Cursor Get routine 
 */
extern int ffdb_cursor_find_by_key (ffdb_htab_t* hashp, ffdb_crs_t* cursor,
				    FFDB_DBT* key, FFDB_DBT* data,
				    unsigned int flags);



/**
 * Get all page information and display them. This is for debug purpose
 */
extern void ffdb_disp_all_page_info (ffdb_htab_t* hashp);


#endif
