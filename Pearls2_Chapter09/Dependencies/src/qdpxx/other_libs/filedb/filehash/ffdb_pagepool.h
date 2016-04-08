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
 *     Memory Cache Page Pool Header File
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_pagepool.h,v $
 *     Revision 1.4  2009-07-18 21:02:36  chen
 *     Fix a bug when cache size is not big enough by changing pgp->npages to be number of pages in the file and maxpgno to be maximum page number currently in use
 *
 *     Revision 1.3  2009/04/17 00:42:14  chen
 *     add dump stack routine
 *
 *     Revision 1.2  2009/03/04 01:44:26  chen
 *     Combine multiple writes
 *
 *     Revision 1.1  2009/02/20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#ifndef _FFDB_PAGE_POOL_H
#define _FFDB_PAGE_POOL_H

#include <stdio.h>
#include <pthread.h>

#include "ffdb_cq.h"

/**
 * Some commonly used macros
 */
#define FFDB_FLAG_CLR(flags, f)   ((flags) &= ~(f))
#define FFDB_FLAG_ISSET(flags, f) ((flags) & (f))
#define FFDB_FLAG_SET(flags, f)   ((flags) |= (f))
#define FFDB_POWER_OF_TWO(x)      (( (x) & (x - 1)) == 0)

#define FFDB_LOCK_INIT(lock)      (pthread_mutex_init ( &(lock), 0))
#define FFDB_LOCK_FINI(lock)      (pthread_mutex_destroy (&(lock)))
#define FFDB_LOCK(lock)           (pthread_mutex_lock(&(lock)))
#define FFDB_UNLOCK(lock)         (pthread_mutex_unlock(&(lock)))
#define FFDB_COND_INIT(cond)      (pthread_cond_init(&(cond), 0))
#define FFDB_COND_FINI(cond)      (pthread_cond_destroy(&(cond)))
#define FFDB_COND_WAIT(cond,lock) (pthread_cond_wait(&(cond), &(lock)))
#define FFDB_COND_SIGNAL(cond)    (pthread_cond_signal(&(cond)))
#define FFDB_COND_BROADCAST(cond) (pthread_cond_broadcast(&(cond)))
#define FFDB_THREAD_NULL(id)      (memset(&id, 0, sizeof(pthread_t)))
#define FFDB_THREAD_ID            (pthread_self())
#define FFDB_THREAD_SAME(id1,id2) (pthread_equal(id1,id2))


/**
 * Clear first a few bytes on a new page
 */
#define FFDB_CLEANHDR_SIZE        32


/**
 * What fractions of pages inside cache to flush out when
 * The maximum of cache is reached
 *
 * This is the inverse of the above value
 */
#define FFDB_WRITE_FRAC           5


/*
 * Common flags --
 *	Interfaces which use any of these common flags should never have
 *	interface specific flags in this range.
 */
#define	FFDB_CREATE	      0x00000001    /* Create file as necessary. */
#define	FFDB_NOMMAP	      0x00000010    /* Don't mmap underlying file. */
#define	FFDB_RDONLY	      0x00000020    /* Read-only (O_RDONLY). */
#define	FFDB_DIRECT	      0x00000040    /* No Buffering (OS)  */
#define	FFDB_THREAD	      0x00000080    /* Applications are threaded. */
#define	FFDB_TRUNCATE	      0x00000100    /* Discard existing DB (O_TRUNC). */

/**
 * cache page operation flag
 */
#define	FFDB_PAGE_CREATE 0x00000001	/* page needs to be written */
#define	FFDB_PAGE_DIRTY	 0x00000002	/* page needs to be written */
#define	FFDB_PAGE_EDIT	 0x00000004	/* Modify without copying. */
#define	FFDB_PAGE_FREE	 0x00000008	/* Free page if present. */
#define	FFDB_PAGE_LAST   0x00000010	/* Return the last page. */
#define	FFDB_PAGE_NEW	 0x00000020	/* Create a new page. */

#define	FFDB_PAGE_PINNED 0x00000100	/* page is pinned into memory */
#define	FFDB_PAGE_VALID	 0x00000200	/* page address is valid */
#define	FFDB_PAGE_LOCKED 0x00000400	/* page should stay in memory */

#define	FFDB_PAGE_IGNOREPIN 0x00001000	/* Ignore if the page is pinned.*/
#define FFDB_PAGE_REQUEST   0x00002000  /* Allocate a new page with a
					   specific page number. */
#define	FFDB_PAGE_NEXT	    0x00004000  /* Allocate a new page with
					   next page number. */




#define	FFDB_MAX_PAGE_NUMBER	0xffffffff	/* >= # of pages in a file */
#define	FFDB_MAX_REC_NUMBER	0xffffffff	/* >= # of records in a tree */

/**
 * Maximum and minimum page size
*/
#define	FFDB_MIN_PGSIZE	0x000200   /* Minimum page size (512). */
#define	FFDB_MAX_PGSIZE	0x040000   /* Maximum page size (262144). */


/*
 * The memory pool scheme is a simple one.  Each in-memory page is referenced
 * by a bucket which is threaded in up to two of three ways.  All active pages
 * are threaded on a hash chain (hashed by page number) and an lru chain.
 * Inactive pages are threaded on a free chain.  Each reference to a memory
 * pool is handed an opaque MPOOL cookie which stores all of this information.
 */
#define	FFDB_HASHSIZE	16384
#define	FFDB_HASHKEY(pgno)	((pgno - 1 + FFDB_HASHSIZE) % FFDB_HASHSIZE)

/**
 * Forward decleration of structure
 */
struct _ffdb_bkt;

/**
 * The waiters of a bucket defined in the following
 */
typedef struct _ffdb_bkt_waiter {
  FFDB_CIRCLEQ_ENTRY(_ffdb_bkt_waiter) q; /* pointer inside waiter queue */
  int wakeup;                             /* flag for conditional variable */
  struct _ffdb_bkt* bp;                   /* which bucket we are waiting  */
  pthread_cond_t cv;                      /* conditional variable         */
}ffdb_bkt_waiter_t;

/**
 * The BKT structures are the elements of the queues.
 */
typedef struct _ffdb_bkt {
  FFDB_CIRCLEQ_ENTRY(_ffdb_bkt) hq;                    /* hash queue */
  FFDB_CIRCLEQ_ENTRY(_ffdb_bkt) lq;                    /* LRU queue */
  FFDB_CIRCLEQ_HEAD(_ffdb_wqh, _ffdb_bkt_waiter) wqh;  /* waiter queue head */  
  void    *page;		                       /* page */
  pgno_t   pgno;		                       /* page number */
  unsigned int ref;                                    /* how many using it */
  unsigned int waiters; 		               /* number of waiters */
  unsigned int flags;		                       /* flags (state)*/
  pthread_t    owner;			               /* owner of this page */
} ffdb_bkt_t;

/**
 * The bucket structure for sorting and flushing to disk purpose
 */
typedef struct _ffdb_sbkt {
  FFDB_SLIST_ENTRY(_ffdb_sbkt) sl;
  ffdb_bkt_t* bp;			       
}ffdb_sbkt_t;

/**
 * Define a head pointer pointing to the head of the single list
 */
typedef FFDB_SLIST_HEAD(_ffdb_slh, _ffdb_sbkt) ffdb_slh_t;

/**
 * Define user supplied pgio function
 */
typedef void (*ffdb_pgiofunc_t) (void* arg, pgno_t pgno, void* mem);



/*
 * The memory page pool structure keeping track of number pages and so on
 */
typedef struct _ffdb_pagepool_
{
  FFDB_CIRCLEQ_HEAD(_ffdb_lqh, _ffdb_bkt) lqh; /* lru queue head */
  FFDB_CIRCLEQ_HEAD(_ffdb_hqh, _ffdb_bkt) hqh[FFDB_HASHSIZE]; /* hash queue array */
  pgno_t	curcache;		/* current number of cached pages */
  pgno_t	maxcache;		/* max number of cached pages */
  pgno_t	npages;			/* number of pages in the file */
  pgno_t	maxpgno;		/* maximum pages number in use */
  unsigned int	pagesize;		/* file page size */
  unsigned int  fileflags;              /* file creation flag */
  int	        fd;		        /* file descriptor */
  int           close_fd;   		/* do i close fd on exit */
  /* page in conversion routine */
  ffdb_pgiofunc_t pgin;
  /* page out conversion routine */
  ffdb_pgiofunc_t pgout;
  void	*pgcookie;		       /* cookie for page in/out routines */
#ifdef _FFDB_STATISTICS
  unsigned int	cachehit;
  unsigned int	cachemiss;
  unsigned int	pagealloc;
  unsigned int	pagereuse;
  unsigned int	pageflush;
  unsigned int	pageswap;
  unsigned int	pageget;
  unsigned int	pagenew;
  unsigned int	pageput;
  unsigned int	pagechange;
  unsigned int	pageread;
  unsigned int	pagewrite;
  unsigned int  pagewait;
#endif  
  pthread_mutex_t lock;
}ffdb_pagepool_t;

#ifdef _cplusplus
extern "C" {
#endif

/**
 * Create pagepool memory handle
 *
 * @param pagepool this is a pointer for ffdb_pagepool_t
 * @param flags not used: has to be zero
 * @return 0 on success, otherwise errno will be returned
 */
extern int
ffdb_pagepool_create (ffdb_pagepool_t** pagepool,
		      unsigned int flags);

/**
 * Open a page cache poll
 *
 * @param filename a backend filename
 * @param flags must be zero or by bitwise OR'ing together one or
 * more of the following value: FFDB_CREATE, FFDB_DIRECT, FFDB_RDONLY
 * @param mode the same as chmod call 
 * @param pagesize a fixed page size to use. This must be power of two
 * @param maxcache how many pages to keep in the page cache poll
 * @param a pointer to ffdb_pagepool_t
 *
 * @return 0 if everything is ok. Otherwise errno is returned
 */
extern int
ffdb_pagepool_open (ffdb_pagepool_t* pagepool,
		    const char* filename,
		    unsigned int flags,
		    int mode,
		    const unsigned int pagesize,
		    const unsigned int maxcache);


/**
 * Open a page cache poll
 *
 * @param fd an opened file descriptor for the backend file
 * @param pagesize a fixed page size to use. This must be power of two
 * @param maxcache how many pages to keep in the page cache poll
 * @param a pointer to ffdb_pagepool_t
 *
 * @return 0 if everything is ok. Otherwise errno is returned
 */
extern int
ffdb_pagepool_open_fd (ffdb_pagepool_t* pagepool,
		       int fd,
		       const unsigned int pagesize,
		       const unsigned int maxcache);


/**
 * Create a new page not backed by any file
 * 
 * @param pgp cache page pool pointer
 * @param pageno returned page number if this routine is success
 * @param flags if flags = FFDB_PAGE_REQUEST, page will be created using
 * page number stored in pageno. If flags = 0, page will be created and new 
 * pagenumber is returned.
 * @param mem returned memory address of this page.
 *
 */
extern int
ffdb_pagepool_new_page (ffdb_pagepool_t* pgp, pgno_t* pageno, 
			unsigned int flags, void** mem);



/**
 * Get a page from the cache poll backed by the file. 
 * If the file is not yet created, a new page
 * will be allocated and eventually written back to the file. Once this
 * page is claimed by a thread, other threads cannot access this thread.
 *
 * @param pgp cache page pool pointer
 * @param pageno requested page number
 * @param flags this flag can be either 0, or by bitwise inclusively OR'ing
 * the following values: 
 * FFDB_PAGE_CREATE if the specified page does not exist, create it.
 * FFDB_PAGE_DIRTY  this page will be modified before leaving the cache
 * FFDB_PAGE_LAST return the last page of the source file and copy it page
 * number to the memory location of the pagno
 * FFDB_NEW create a new page in the file, and copy its page number into the
 * the memory localtion of the pageno.
 * @param mem returned memory address of this page.
 * @return 0 on success. Otherwise return errno
 *
 */
extern int
ffdb_pagepool_get_page (ffdb_pagepool_t* pgp, pgno_t* pageno,
			unsigned int flags, void** mem);


/**
 * Put a page back into the cache so that other threads can access this page
 *
 * 
 * @param pgp cache page pool pointer
 * @param mem user cached page memory will be returned
 * @param flags this flags can be FFDB_PAGE_DIRTY if this page has been modified
 * 
 * @return 0 if this page is successfully be put back into the cache. Otherwise
 * return -1 with an appropriate errno set
 */
extern int
ffdb_pagepool_put_page (ffdb_pagepool_t* pgp, void* mem, 
			unsigned int flags);




/**
 * Put a page back into the cache so that other threads can access this page
 * with new page number for this page
 *
 * 
 * @param pgp cache page pool pointer
 * @param mem user cached page memory will be returned
 * @param newpgnum a new pagenumber associated with this page
 * 
 * @return 0 if this page is successfully be put back into the cache. Otherwise
 * return -1 with an appropriate errno set
 */
extern int
ffdb_pagepool_change_page (ffdb_pagepool_t* pgp, void* mem, 
			   pgno_t newpagenum);



/**
 * Flush all dirty pages back to the back end file. However, if any modified 
 * pages are in use. They will be ignored
 *
 * @param  pgp cache page pool pointer
 */
extern int
ffdb_pagepool_sync (ffdb_pagepool_t* pgp);

/**
 * Close the page poll pointer and any resource associated with this file
 * This implies all dirty pages are flushed out, 
 *
 * @param pgp cache page poll pointer
 * @return 0 the file is closed successfully and all pages are written back
 * to th back file. If there are still pinned pages by other threads, return -1
 * and errno is set to EAGAIN.
 */
extern int
ffdb_pagepool_close (ffdb_pagepool_t* pgp);


/**
 * Really delete page from cache. This is used either by internal
 * other routines or user knows this page is not used anymore
 * @param pgp cache page poll pointer
 * @param mem memory pointer of this page
 * @return 0 this page is deleted and memory is freed. Otherwise return -1 or
 * appropriate errno: EAGAIN: there are waiters on this page.
 */
extern int
ffdb_pagepool_delete (ffdb_pagepool_t* pgp, void* mem);


  /**
   * Add user callback to page in/out from a backend file
   * @param pgp a pagepool pointer
   * @param pgin a function for a pgae read in from backend file
   * @param pgoit a function for a page flush out to disk
   * @param cookie a void* user supplied argument
   */
extern void
ffdb_pagepool_filter (ffdb_pagepool_t* pgp, ffdb_pgiofunc_t pgin,
		      ffdb_pgiofunc_t pgout, void* cookie);


/**
 * Simple utility to dump stack trace
 * Now it only has implementation on linux
 */
extern void 
ffdb_dump_stack (void);


#ifdef _FFDB_STATISTICS
/** 
 *Print out statistics information of this page pool
 * @param pgp the pointer to pagepool
 */
extern void
ffdb_pagepool_stat (ffdb_pagepool_t* pgp);
#endif


#ifdef _cplusplus
};
#endif

#endif
