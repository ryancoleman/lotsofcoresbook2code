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
 *     Memory Cache Page Pool Implementation based on both old and new 
 *     Berkeley DB
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_pagepool.c,v $
 *     Revision 1.7  2009-07-18 21:02:36  chen
 *     Fix a bug when cache size is not big enough by changing pgp->npages to be number of pages in the file and maxpgno to be maximum page number currently in use
 *
 *     Revision 1.6  2009/05/08 17:37:31  chen
 *     Fix a major bug (not clean out moved page inside page pool cache)
 *
 *     Revision 1.5  2009/04/21 18:51:20  chen
 *     Fix bugs related to number of pages upon moving pages in addition to clean pages on disk when the pages has been moved
 *
 *     Revision 1.4  2009/04/17 00:42:14  chen
 *     add dump stack routine
 *
 *     Revision 1.3  2009/03/04 01:44:26  chen
 *     Combine multiple writes
 *
 *     Revision 1.2  2009/02/24 04:25:05  edwards
 *     Check if O_LARGEFILE is defined before attempting to add it to the lflags
 *     used in the "open" call.
 *
 *     Revision 1.1  2009/02/20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef __USE_LARGEFILE64
#define __USE_LARGEFILE64
#endif

#ifndef __LARGEFILE64_SOURCE
#define __LARGEFILE64_SOURCE
#endif

// BJ: Added these to fix the compiler which complained when 
// I tried to compile this file with -std=c99
// See also: http://www.openwall.com/lists/owl-dev/2012/02/14/1
#ifndef __USE_XOPEN2K8
#define __USE_XOPEN2K8
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif

#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <ffdb_db.h>
#include "ffdb_pagepool.h"

#ifdef __linux

#include <execinfo.h>
#define SBUF_SIZE 400

void ffdb_dump_stack(void)
{
  int i, nptrs;
  void *buffer[SBUF_SIZE];
  char **strings;

  nptrs = backtrace(buffer, SBUF_SIZE);
  printf("backtrace() returned %d addresses\n", nptrs);

  strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL) {
    perror("backtrace_symbols");
    exit(EXIT_FAILURE);
  }

  fprintf (stderr, "[bt] Execution path: \n");
  for (i = 0; i < nptrs; i++)
    fprintf(stderr, "[bt] %s\n", strings[i]);

  free (strings);
}
#else
void ffdb_dump_stack (void)
{
  /* empty */
}
#endif


/**
 * Flush out some pages in order of page numbers
 */
static int 
_ffdb_pagepool_sync_i (ffdb_pagepool_t* pgp, unsigned int numpages);

/* Test for valid page sizes. */
#define	IS_VALID_PAGESIZE(x)						\
	(FFDB_POWER_OF_TWO(x) && (x) >= FFDB_MIN_PGSIZE && ((x) <= FFDB_MAX_PGSIZE))

/**
 * Do a shallow copy of bucket from real bucket to a simple bucket
 */
static void
_ffdb_shallow_copy_bk (ffdb_sbkt_t* sbp, ffdb_bkt_t* bp)
{
  sbp->bp = bp;
}

/**
 * Comparision function for two simple buckets
 */
static int
_ffdb_sbkt_cmp (ffdb_sbkt_t* b1, ffdb_sbkt_t* b2)
{
  return b1->bp->pgno - b2->bp->pgno;
}


/**
 * Non recursive merge sort of a single linked list according to
 * page number
 */
static void
_ffdb_slist_merge_sort (ffdb_slh_t* head)
{
  ffdb_sbkt_t *p, *q, *e, *tail;
  int insize, nmerges, psize, qsize, i;

  if (FFDB_SLIST_EMPTY(head))
    return;

  /* First time merge two elements */
  insize = 1;

  while (1) {
    p = FFDB_SLIST_FIRST(head);
    FFDB_SLIST_FIRST(head) = 0;
    tail = 0;

    /* count number of merges we do in this pass */
    nmerges = 0; 
    while (p) {
      nmerges++;
      /* there merges to be done */
      q = p;
      psize = 0;
      for (i = 0; i < insize; i++) {
	psize++;
	q = FFDB_SLIST_NEXT(q,sl);
	if (!q) break;
      }

      /* if q has not fallen off end, we have two lists to merge */
      qsize = insize;
      
      /* now we have two lists: merge them */
      while (psize > 0 || (qsize > 0 && q)) {
	  
	  /* decide whether next element of merge comes from q or q */
	  if (psize == 0) {
	    /* p is empty, e must come from q */
	    e = q;
	    q = FFDB_SLIST_NEXT(q,sl);
	    qsize--;
	  }
	  else if (qsize == 0 || !q) {
	    /* q is empty, e must come from p */
	    e = p;
	    p = FFDB_SLIST_NEXT(p,sl);
	    psize--;
	  }
	  else if (_ffdb_sbkt_cmp (p, q) <= 0) {
	    /* First element of p is lower (or same)
	     * e must be from p */
	    e = p;
	    p = FFDB_SLIST_NEXT(p,sl);
	    psize--;
	  }
	  else {
	    /* First element of q is lower (or same)
	     * e must be from q */
	    e = q;
	    q = FFDB_SLIST_NEXT(q,sl);
	    qsize--;
	  }
	  
	  /* add the next element to the merged list */
	  if (tail) {
	    FFDB_SLIST_INSERT_AFTER(tail,e,sl);
	  }
	  else {
	    FFDB_SLIST_FIRST(head) = e;
	  }
	  
	  /* assign tail to the new one */
	  tail = e;
	}

	/* Now p has stepped 'insize' places along, and q has too */
	p = q;
    }
    
    /* Now we finished one pass of the list using merge size of insize
     * We have to set the tail next pointer to zero */
    FFDB_SLIST_NULL_NEXT(tail,sl);

    /* If we have done only one merge, we are done */
    if (nmerges <= 1) /* allowing for nmerges == 0, the empty list case */
      break;

    /* Otherwise repeat, merging list twice the size */
    insize *= 2;
  }
    
  /* When we get here, we are finished */
}


/*
 * _ffdb_pagepool_write
 *	Write a dirty page to disk.
 * This routine is called with page being locked
 */
static int
_ffdb_pagepool_write(ffdb_pagepool_t* pgp, 
		     ffdb_bkt_t* bp)
{
  off_t offset;
  int nbytes;
  int ret = 0;

#ifdef _FFDB_STATISTICS
  ++pgp->pagewrite;
#endif

  /* Run through the user's filter. */
  if (pgp->pgout)
    (pgp->pgout)(pgp->pgcookie, bp->pgno, bp->page);

  offset =  (off_t)pgp->pagesize * bp->pgno;

  if (lseek(pgp->fd, offset, SEEK_SET) != offset) 
    ret = -1;
  else {
    if ((nbytes = write(pgp->fd, bp->page, pgp->pagesize)) != pgp->pagesize) 
      ret = -1;
  }

  if (ret == 0)
    FFDB_FLAG_CLR(bp->flags, FFDB_PAGE_DIRTY);

  /* Update how many pages this file holds now */
  if (bp->pgno >= pgp->npages) {
    pgp->npages = bp->pgno + 1;
  }

  return ret;
}


/*
 * _ffdb_clean_page_ondisk
 *	Clean a page on disk
 * 
 */
static int
_ffdb_clean_page_ondisk (ffdb_pagepool_t* pgp, pgno_t num)
{
  off_t offset;
  int nbytes;
  int ret = 0;
  char *cleanbuf;

#ifdef _FFDB_STATISTICS
  ++pgp->pagewrite;
#endif

  /* allocate clean memory */
  cleanbuf = (char *)calloc (pgp->pagesize, sizeof(char));
  if (!cleanbuf) {
    fprintf (stderr, "Cannot allocate a clean buffer for page %d\n", num);
    exit (1);
  }

  offset =  (off_t)pgp->pagesize * num;

  if (lseek(pgp->fd, offset, SEEK_SET) != offset) 
    ret = -1;
  else {
    if ((nbytes = write(pgp->fd, cleanbuf, pgp->pagesize)) != pgp->pagesize) 
      ret = -1;
  }

  /* free memory */
  free (cleanbuf);

  return ret;
}

/**
 * Get a page from cache when the page is not used by a thread
 * @param pgp pagepool pointer
 * @param bkt a new pointer to a bucket
 * @return 0 on success, -1 return no bucket can be reused. otherwise
 * error on flushing the unused page
 * Upon returning of this routine, the reused bucket's pinned flag is set
 *
 * This routine is called when pgp->lock is held
 */
static int
_ffdb_pagepool_reuse_bkt (ffdb_pagepool_t* pgp, ffdb_bkt_t** retbp)
{
  struct _ffdb_hqh *head;
  int ret = 0;
  int needwrite = 0;
  int reusepage = 0;
  ffdb_bkt_t* bp = 0;

  *retbp = 0;
  /**
   * Sanity check: LRU queues should not be empty
   */
  if (FFDB_CIRCLEQ_EMPTY(&pgp->lqh)) {
    fprintf (stderr, "_ffdb_pagepool_bkt: LRU queue is empty. Quit!\n");
    abort ();
  }

  /**
   * Walk the LRU queue now
   */
  FFDB_CIRCLEQ_FOREACH(bp, &pgp->lqh, lq) {

    reusepage = 0;
      
    if (!(FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_LOCKED)) && 
	!(FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_PINNED)) &&
	bp->waiters == 0) {
      /* This page is not locked and pinned, so it can be reused */
      /* And no one is waiting for it */
      reusepage = 1;
      /**
       * Flush if is dirty. 
       * This may not be efficient
       */
      if (FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_DIRTY)) {
	if (_ffdb_pagepool_write(pgp, bp) == -1) {
	  fprintf (stderr, "_ffdb_pagepool_bkt: page %d flush error\n",
		   bp->pgno);
	  ret = errno;
	  reusepage = 0;
	}
	else
	  needwrite = 1;
      }
      
      /**
       * If reusepage flag is set, this page can be reused
       */
      if (reusepage) {
#ifdef _FFDB_STATISTICS
	++pgp->pageswap;
#endif
	/* Remove from the hash and lru queues. */
	head = &pgp->hqh[FFDB_HASHKEY(bp->pgno)];
	FFDB_CIRCLEQ_REMOVE(head, bp, hq);
	FFDB_CIRCLEQ_REMOVE(&pgp->lqh, bp, lq);
#if 0
	fprintf (stderr, "Reuse remove page number %d\n", bp->pgno);
#endif

	bp->ref = 0;
	bp->waiters = 0;
	bp->flags = 0;
	bp->owner = FFDB_THREAD_ID;

	/* Now I need to set flags before unlock */
	bp->flags = FFDB_PAGE_PINNED | FFDB_PAGE_VALID;

#ifdef _FFDB_STATISTICS
	++pgp->pagereuse;
#endif
	*retbp = bp;
      }
    }
    if (reusepage)
      break;
  }
  /* not found any reuseable page */
  if ((*retbp) == 0)
    ret = -1;

  /* Now flush out some fraction of pages to speed up performance */
  if (needwrite) {
#ifdef _FFDB_DEBUG
    fprintf (stderr, "Flush %d pages out\n", pgp->maxcache/FFDB_WRITE_FRAC);
#endif
    _ffdb_pagepool_sync_i (pgp, pgp->maxcache/FFDB_WRITE_FRAC);
  }

  return ret;
}


/**
 * Get a page from the cache or create a new one
 * There is no lock on this routine. Higher level routines have to do it
 *
 * When this routine is one, the newly found (created) bucket should have
 * its pinned flag set
 *
 * This routine is called when the pgp->lock is held
 */
static ffdb_bkt_t *
_ffdb_pagepool_new_bkt (ffdb_pagepool_t* pgp)
{
  ffdb_bkt_t *bp = 0;

  /* If under cache limit, or there are no pages can be flushed
   * we always create new page
   *
   * valgrind keeps complaining about uninitialized memory. 
   */
  if ((bp = (ffdb_bkt_t *)malloc(sizeof(ffdb_bkt_t) + pgp->pagesize)) == 0)
    return 0;

#ifdef _FFDB_STATISTICS
  ++pgp->pagealloc;
#endif
  ++pgp->curcache;

  bp->page = (char *)bp + sizeof(ffdb_bkt_t);
  bp->ref = 0;
  bp->waiters = 0;
  bp->flags = 0;
  bp->owner = FFDB_THREAD_ID;

  /* change the bp flags */
  bp->flags = FFDB_PAGE_PINNED | FFDB_PAGE_VALID;

  /* initialize waiter queue */
  FFDB_CIRCLEQ_INIT(&(bp->wqh));

  return (bp);
}

/**
 * Create ffdb_pagepool handle used by all threads of a process
 */
int
ffdb_pagepool_create (ffdb_pagepool_t** pgp, unsigned int flags)
{
  int ret, i;

  ffdb_pagepool_t *p = (ffdb_pagepool_t *)malloc(sizeof(ffdb_pagepool_t));
  if (!p) {
    fprintf (stderr, "cannot allocate space for ffdb_pagepool_t structure.\n");
    return errno;
  }
  /**
   * initialize each element of this structure
   */
  memset (p, 0, sizeof(ffdb_pagepool_t));
  p->fd = -1;

  /**
   * Initialize LRU and hash table
   */
  FFDB_CIRCLEQ_INIT (&(p->lqh));
  for (i = 0; i < FFDB_HASHSIZE; i++) 
    FFDB_CIRCLEQ_INIT (&(p->hqh[i]));  

  /**
   * Create a pthread mutex lock
   */
  if ((ret = FFDB_LOCK_INIT (p->lock)) != 0) {
    free (p);
    return ret;
  }

  *pgp = p;
  return 0;
}

/**
 * Open a page cache poll object using a given file
 */
int
ffdb_pagepool_open (ffdb_pagepool_t* pgp,
		    const char* filename,
		    unsigned int flags,
		    int mode,
		    const unsigned int pagesize,
		    const unsigned int maxcache)
{
  int lflags, fd;
  struct stat sb;

  /* check whether some other threads already opened this file */
  FFDB_LOCK(pgp->lock);
  if (pgp->fd != -1) {
    FFDB_UNLOCK(pgp->lock);
    return 0;
  }

  /**
   * first check open flags
   */
  unsigned int ok_flags = FFDB_CREATE | FFDB_RDONLY | FFDB_DIRECT;
  if (FFDB_FLAG_ISSET(flags, ~ok_flags)) {
    fprintf (stderr, "ffdb_pagepool_open wrong flags specification\n");
    goto openerr;
  }
  
  /**
   * Check page size
   */
  if (!IS_VALID_PAGESIZE(pagesize)) {
    fprintf (stderr, "ffdb_pagepool_open wrong pagesize\n");
    errno = EINVAL;
    goto openerr;
  }
   
  pgp->fileflags = flags;

  /**
   * Let us open the file
   */
  lflags = 0;
  if (FFDB_FLAG_ISSET(flags, FFDB_CREATE))
    lflags |= O_CREAT;

  if (FFDB_FLAG_ISSET(flags, FFDB_RDONLY))
    lflags |= O_RDONLY;
  else
    lflags |= O_RDWR;

#ifdef O_DIRECT
  if (FFDB_FLAG_ISSET(flags, FFDB_DIRECT))
    lflags |= O_DIRECT;
#endif
  

  /**
   * Always support large file if I am on a 32bit linux machine
   */
#if defined(O_LARGEFILE)
  if (sizeof(long) == sizeof(int))
    lflags |= O_LARGEFILE;
#endif

  fd = open (filename, lflags, mode);
  if (fd <= 0) {
    fprintf (stderr, "ffdb_pagepool_open cannot open file %s\n",
	     filename);
    goto openerr;
  }


  /* Now check a few more things about this file */
  if (fstat (fd,  &sb) != 0) 
    goto openerr;

  if (!S_ISREG(sb.st_mode)) {
    errno = EINVAL;
    goto openerr;
  }
  /* File size should be multiple of pagesize */
  if (sb.st_size % pagesize != 0) {
    fprintf (stderr, "File %s filesize %ld is not multiple of pagesize %d\n",
	     filename, sb.st_size, pagesize);
    errno = EINVAL;
    goto openerr;
  }

  /* Set up some attribute of ffdb_pagepool structure */
  pgp->fd = fd;
  pgp->close_fd = 1;
  pgp->maxcache = maxcache;
  pgp->pagesize = pagesize;
  
  /* number of pages I am holding */
  pgp->npages = sb.st_size / pagesize;

  /* maximum page number this pagepool has now */
  pgp->maxpgno = pgp->npages;

  /* unlock the code */
  FFDB_UNLOCK(pgp->lock);
  return 0;

 openerr:
  FFDB_UNLOCK(pgp->lock);
  return errno;
}


/**
 * Open a page cache poll object using a given file
 */
int
ffdb_pagepool_open_fd (ffdb_pagepool_t* pgp,
		       int fd,
		       const unsigned int pagesize,
		       const unsigned int maxcache)
{
  unsigned int oflags = 0;
  unsigned int flags = 0;
  struct stat sb;

  /* check whether some other threads already opened this file */
  FFDB_LOCK(pgp->lock);
  if (pgp->fd != -1) {
    FFDB_UNLOCK(pgp->lock);
    return 0;
  }

  /* get flags of this opened file descriptor 
   * this flag has the same values defined in open(2) system call
   */
  oflags = fcntl (fd, F_GETFL);
  if ((oflags & O_ACCMODE) == O_RDONLY)
    FFDB_FLAG_SET (flags, FFDB_RDONLY);

  /**
   * first check open flags
   */
  unsigned int ok_flags = FFDB_CREATE | FFDB_RDONLY | FFDB_DIRECT;
  if (FFDB_FLAG_ISSET(flags, ~ok_flags)) {
    fprintf (stderr, "ffdb_pagepool_open wrong flags specification\n");
    goto openerr;
  }
  pgp->fileflags = flags;

  /**
   * Check page size
   */
  if (!IS_VALID_PAGESIZE(pagesize)) {
    fprintf (stderr, "ffdb_pagepool_open wrong pagesize\n");
    errno = EINVAL;
    goto openerr;
  }
   
  /* Now check a few more things about this file */
  if (fstat (fd,  &sb) != 0) 
    goto openerr;

  if (!S_ISREG(sb.st_mode)) {
    errno = EINVAL;
    goto openerr;
  }

  /* File size should be multiple of pagesize */
  if (sb.st_size % pagesize != 0) {
    fprintf (stderr, "filesize %ld is not multiple of pagesize %d\n",
	     sb.st_size, pagesize);
    errno = EINVAL;
    goto openerr;
  }

  /* Set up some attribute of ffdb_pagepool structure */
  pgp->fd = fd;
  pgp->close_fd = 0;
  pgp->maxcache = maxcache;
  pgp->pagesize = pagesize;

  /* number of pages I am holding */
  pgp->npages = sb.st_size / pagesize;

  /* maximum page number the pool has now */
  pgp->maxpgno = pgp->npages;

  /* unlock the code */
  FFDB_UNLOCK(pgp->lock);
  return 0;

 openerr:
  FFDB_UNLOCK(pgp->lock);
  return errno;
}

/**
 * Get a new page from the back source file
 */
static int
_ffdb_pagepool_load_new_page (ffdb_pagepool_t* pgp, pgno_t pageno,
			      unsigned int flags, void** mem)
{
  int status, nbytes;
  off_t off;
  struct _ffdb_hqh *head;
  ffdb_bkt_t* bp = 0;

  /**
   * Get a BKT from the cache.  Assign a new page number, attach
   * it to the head of the hash chain, the tail of the lru chain,
   * and return.
   *
   */
  if (pgp->curcache > pgp->maxcache) {
    /**
     * If the cache is max'd out, walk the lru list for a buffer we
     * can flush.  If we find one, write it (if necessary) and take it
     * off any lists.  If we don't find anything we grow the cache anyway.
     * The cache never shrinks.
     */
    status = _ffdb_pagepool_reuse_bkt (pgp, &bp);
    if (status == -1) {
      /* cannot find page to reuse */
      bp = _ffdb_pagepool_new_bkt (pgp);
    }
  }
  else 
    bp = _ffdb_pagepool_new_bkt (pgp);


  if (!bp) {
    /* This has to be successful. This is a new page */
    fprintf (stderr, "ffdb_pagepool_load_new_page: cannot get new page of page number %d\n", pageno);
    abort ();
  }
  
  /**
   * The obtained bucket has pinned flag set, we own this page.
   * It is time to populate this page using back file
   */
  off = (off_t)pgp->pagesize * (pageno);
  if (lseek (pgp->fd, off, SEEK_SET) != off) {
    fprintf (stderr, "ffdb_pagepool_load_new_page: cannot seek to a right position.\n");
    return errno;
  }
  nbytes = read (pgp->fd, bp->page, pgp->pagesize);
  if (nbytes != pgp->pagesize && nbytes > 0) {
    fprintf (stderr, "ffdb_pagepool_load_new_page: cannot read back end file\n");
    return errno;
  }
  else if (nbytes == 0) 
    memset (bp->page, 0, pgp->pagesize);

#ifdef _FFDB_STATISTICS
  ++pgp->pageread;
#endif

  /* Set page number */
  bp->pgno = pageno;
  bp->owner = FFDB_THREAD_ID;
  bp->ref = 1;
  bp->waiters = 0;

  /* Change flags of this page since I own this page now */
  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_DIRTY) ||
      FFDB_FLAG_ISSET(flags, FFDB_PAGE_EDIT))
    FFDB_FLAG_SET(bp->flags, FFDB_PAGE_DIRTY);
  
  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_LOCKED))
    FFDB_FLAG_SET(bp->flags, FFDB_PAGE_LOCKED);
  
  *mem = bp->page; 

  /* insert this page into LRU and hash bucket */
  head = &pgp->hqh[FFDB_HASHKEY(bp->pgno)];
#if 0
  {
    ffdb_bkt_t* bk;
    FFDB_CIRCLEQ_FOREACH(bk, head, hq) {
      if (bk->pgno == bp->pgno) {
	fprintf (stderr, "Hash has this page %d alreay\n", bp->pgno);
	pause ();
      }
    }

    FFDB_CIRCLEQ_FOREACH(bk, &pgp->lqh, lq) {
      if (bk->pgno == bp->pgno) {
	fprintf (stderr, "LRU has this page %d alreay\n", bp->pgno);
	pause ();
      }
    } 
  }
#endif

  FFDB_CIRCLEQ_INSERT_HEAD(head, bp, hq);
  FFDB_CIRCLEQ_INSERT_TAIL(&pgp->lqh, bp, lq);
#if 0
  fprintf (stderr, "Load page insert pageno %d\n", bp->pgno);
#endif

  /**
   * Check whether page in callback routine 
   */
  if (pgp->pgin) {
    (pgp->pgin)(pgp->pgcookie, bp->pgno, bp->page);
  }

  return 0;
}
				
				

/**
 * Create a new page not from the back source file
 * @param pgp pagepool pointer
 * @param pageno address of either requested page number or retured page number
 * @param flags request flags.  if flags = FFDB_PAGE_REQUEST, page 
 * will be created using  page number stored in pageno. If flags = 0, 
 * page will be created and new pagenumber is returned.
 * @param mem returned memory address of this page
 * @return 0 on success, otherwise return either errno or -1.
 *
 * This code should be called with lock held
 */
static int
_ffdb_pagepool_new_page_i (ffdb_pagepool_t* pgp, pgno_t* pageno,
			   unsigned int flags,  void** mem)
{
  int status;
  struct _ffdb_hqh *head;
  ffdb_bkt_t *bp = 0;

  if (pgp->maxpgno == FFDB_MAX_PAGE_NUMBER) {
    (void)fprintf(stderr, "ffdb_pagepool_new_page: page allocation overflow.\n");
    abort();
  }
#ifdef _FFDB_STATISTICS
  ++pgp->pagenew;
#endif

  /*
   * Get a BKT from the cache.  Assign a new page number, attach
   * it to the head of the hash chain, the tail of the lru chain,
   * and return.
   *
   */
  if (pgp->curcache > pgp->maxcache) {
    /**
     * If the cache is max'd out, walk the lru list for a buffer we
     * can flush.  If we find one, write it (if necessary) and take it
     * off any lists.  If we don't find anything we grow the cache anyway.
     * The cache never shrinks.
     */
    status = _ffdb_pagepool_reuse_bkt (pgp, &bp);
    if (status == -1) {
      /* cannot find page to reuse */
      bp = _ffdb_pagepool_new_bkt (pgp);
    }
  }
  else
    bp = _ffdb_pagepool_new_bkt (pgp);

  /**
   * If we do not have a page, we return -1
   */
  if (!bp) {
    fprintf (stderr, "_ffdb_pagepool_new_page_i: cannot find a new page\n");
    return -1;
  }
  if (FFDB_FLAG_ISSET (flags, FFDB_PAGE_REQUEST)) {
    bp->pgno = *pageno;
    /* new pages can be less than last page because there may be holes */
    if (bp->pgno > pgp->maxpgno) 
      pgp->maxpgno = bp->pgno;
  } else {
    pgp->maxpgno++;
    bp->pgno = *pageno = pgp->maxpgno;
  }

  /* Now we have one thread holding this page */
  bp->ref = 1;
  /* now assign my thread to owner */
  bp->owner = FFDB_THREAD_ID;

  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_DIRTY) ||
      FFDB_FLAG_ISSET(flags, FFDB_PAGE_EDIT))
    FFDB_FLAG_SET(bp->flags, FFDB_PAGE_DIRTY);
  
  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_LOCKED))
    FFDB_FLAG_SET(bp->flags, FFDB_PAGE_LOCKED);


  head = &pgp->hqh[FFDB_HASHKEY(bp->pgno)];
#if 0
  {
    ffdb_bkt_t* bk;
    FFDB_CIRCLEQ_FOREACH(bk, head, hq) {
      if (bk->pgno == bp->pgno) {
	fprintf (stderr, "Hash has this page %d alreay\n", bp->pgno);
	pause ();
      }
    }

    FFDB_CIRCLEQ_FOREACH(bk, &pgp->lqh, lq) {
      if (bk->pgno == bp->pgno) {
	fprintf (stderr, "LRU has this page %d alreay\n", bp->pgno);
	pause ();
      }
    } 
  }
#endif
  FFDB_CIRCLEQ_INSERT_HEAD(head, bp, hq);
  FFDB_CIRCLEQ_INSERT_TAIL(&pgp->lqh, bp, lq);

  *mem = bp->page;

  /**
   * since this is a new page. clear first few bytes 
   * to tell applications that this is a clean page 
   */
  memset (bp->page, 0, FFDB_CLEANHDR_SIZE);

  return 0;
}


/**
 * Create a new page not from the back source file
 */
int
ffdb_pagepool_new_page (ffdb_pagepool_t* pgp, pgno_t* pageno,
			unsigned int flags, void** mem)
{
  int status;

  *mem = 0;
  FFDB_LOCK(pgp->lock);
  status = _ffdb_pagepool_new_page_i (pgp, pageno, flags, mem);
  FFDB_UNLOCK(pgp->lock);
  
  return status;
}

/**
 * Get a cached page from the page poll
 * There will be no difference of treatment on the thread getting
 * this page, i.e. readonly or write
 * flags can be 0 or OR'ING the following values
 * FFDB_PAGE_CREATE if the specified page does not exist, create it.
 * FFDB_PAGE_DIRTY  this page will be modified before leaving the cache
 * FFDB_PAGE_LOCK   this page will stay in cache
 * number to the memory location of the pagno
 * FFDB_NEW create a new page in the file, and copy its page number into the
 * the memory localtion of the pageno.
 *
 * Since each page is locked by checking whether pinned flag is set,
 * so as long as pinned flag is changed, one is ok to modify other
 * attribute of the page
 */
int
ffdb_pagepool_get_page (ffdb_pagepool_t* pgp, pgno_t* pageno,
			unsigned int flags, void** mem)
{
  int ret, found;
  ffdb_bkt_t* bp;
  struct _ffdb_hqh *head;

  /* Set memory pointer to NULL */
  *mem = 0;

  FFDB_LOCK (pgp->lock);
#ifdef _FFDB_STATISTICS
  pgp->pageget++;
#endif

  /**
   * Check flag for consistence
   */
  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_DIRTY)) {
    if (FFDB_FLAG_ISSET(pgp->fileflags, FFDB_RDONLY)) {
	fprintf (stderr, "ffdb_pagepool_get: DIRTY_PAGE flag cannot be used on readonly file.\n");
	errno = EINVAL;
	FFDB_UNLOCK (pgp->lock);
	return errno;
      }
  }
  
  /**
   * Handle the case of new page
   * No need to lock page here because each page is a new alloced page
   */
  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_NEW)) {
    if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_CREATE)) {
	fprintf (stderr, "ffdb_pagepool_get: PAGE_NEW flag cannot be used with PAGE_CREATE flag\n");
	errno = EINVAL;
	FFDB_UNLOCK (pgp->lock);
	return errno;
    }
    fprintf (stderr, "calling new page i for page number %d\n", pageno);
    ret = _ffdb_pagepool_new_page_i (pgp, pageno, flags, mem);
    FFDB_UNLOCK (pgp->lock);    
    return ret;
  }

#if 0
  {
    ffdb_bkt_t* bk;

    head = &pgp->hqh[FFDB_HASHKEY(*pageno)];
    FFDB_CIRCLEQ_FOREACH(bk, head, hq) {
      if (bk->pgno == *pageno ) {
	fprintf (stderr, "Hash has this page %d alreay 0x%x\n", *pageno,
		 bk->page);
      }
    }

    FFDB_CIRCLEQ_FOREACH(bk, &pgp->lqh, lq) {
      if (bk->pgno == *pageno) {
	fprintf (stderr, "LRU has this page %d alreay 0x%x\n", *pageno,
		 bk->page);
      }
    } 
  }
#endif

  /**
   * Try to find a page from existing cache. This page has to be
   * not pinned by other threads
   */
  found = 0;
  head = &pgp->hqh[FFDB_HASHKEY(*pageno)];
  FFDB_CIRCLEQ_FOREACH(bp, head, hq) {
    if (bp->pgno == *pageno) {
      found = 1;
      break;
    }
  }

#ifdef _FFDB_STATISTICS  
  if (found)
    pgp->cachehit++;
  else
    pgp->cachemiss++;
#endif

  if (found) { /* Now I am still holding the lock */
#ifdef _FFDB_DEBUG
    fprintf (stderr, "Found page %d in pagepool at memory 0x%x\n", *pageno,
	     bp->page);
#endif
    /**
     * If I am the owner, the bp will be returned
     */
    if (FFDB_THREAD_SAME(bp->owner, FFDB_THREAD_ID)) {
      FFDB_FLAG_SET(bp->flags, FFDB_PAGE_PINNED);
      bp->ref++;
    }
    else {
      /**
       * A different thread try to access this page
       */
      if (bp->waiters > 0 || FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_PINNED)) {
	/* Create a waiter for this page */
	ffdb_bkt_waiter_t *waiter;
	ffdb_bkt_waiter_t* fw = 0;

	waiter = (ffdb_bkt_waiter_t *)malloc(sizeof(ffdb_bkt_waiter_t));
	if (!waiter) {
	  fprintf (stderr, "ffdb_pagepool_get: cannot allocate space for waiter object\n");
	  abort ();
	}
	waiter->wakeup = 0;
	waiter->bp = bp;
	FFDB_COND_INIT(waiter->cv);
	bp->waiters++;
	FFDB_CIRCLEQ_INSERT_HEAD(&bp->wqh, waiter, q);
	
	/* Now this waiter is waiting for the page */
	while (waiter->wakeup == 0) 
	  FFDB_COND_WAIT(waiter->cv, pgp->lock);
	
	/* Now waiter is done, we should have the page now */
	fw = FFDB_CIRCLEQ_LAST(&bp->wqh);
	if (fw != waiter) {
	  fprintf (stderr, "ffdb_pagepool_get_page: waiter wakeup with wrong waiter pointer\n");
	  abort ();
	}
	/* remove this from queue */
	FFDB_CIRCLEQ_REMOVE(&bp->wqh, waiter, q);

	bp->waiters--;

      	/* Destroy the conditional variable */
	FFDB_COND_FINI(waiter->cv);

	free (waiter);
      }
      /* Now I have grabed the page */
      bp->ref++;

      /* Now set owner of this page */
      bp->owner = FFDB_THREAD_ID;

      /* Now I get hold of this page */      
      FFDB_FLAG_SET(bp->flags, FFDB_PAGE_PINNED);
    }

    /* Change flags of this page since I own this page now */
    if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_DIRTY) ||
	FFDB_FLAG_ISSET(flags, FFDB_PAGE_EDIT))
      FFDB_FLAG_SET(bp->flags, FFDB_PAGE_DIRTY);
  
    if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_LOCKED))
      FFDB_FLAG_SET(bp->flags, FFDB_PAGE_LOCKED);
  
    *mem = bp->page; 
    
    /**
     * The following removal and reinsert must be done at the same 
     * time.
     * If the item is removed first, the lock is given up when 
     * a thread is waiting on the conditional variable.
     * Another thread come in to request the same page, it will not
     * find the page and will create a new page with the same page
     * number. One will have two entries in the hash and LRU with
     * the same page number
     */

    /* remove this page from hash and LRU */
    FFDB_CIRCLEQ_REMOVE(head, bp, hq);
    FFDB_CIRCLEQ_REMOVE(&pgp->lqh, bp, lq);
    /* We found this page in the cache so we have to 
     * move this page to the head of the hash chain and the tail of the
     * lru chain
     */
    FFDB_CIRCLEQ_INSERT_HEAD(head, bp, hq);
    FFDB_CIRCLEQ_INSERT_TAIL(&pgp->lqh, bp, lq);
#if 0
    fprintf (stderr, "Insert pageno %d\n", bp->pgno);
#endif

    FFDB_UNLOCK(pgp->lock);
    return 0;
  }
  else {
    /* We have not found a page with this page number, therefore
     * we have to create new pages using the pagenumber
     */
    /**
     * We cannot find this page with provided page number
     * if flag FFDB_PAGE_CREATE is set, we have to create this page
     */
    if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_CREATE) &&
	*pageno >= pgp->npages ) {
      ret = _ffdb_pagepool_new_page_i (pgp, pageno, 
				       flags | FFDB_PAGE_REQUEST, mem);
    }
    else {
      ret = _ffdb_pagepool_load_new_page (pgp, *pageno, flags, mem);
    }


    FFDB_UNLOCK (pgp->lock);

    return ret;
  }
  /* Should never get here */
  fprintf (stderr, "ffdb_pagepool_get_page: should not get here\n");
  return 0;
}

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
int
ffdb_pagepool_put_page (ffdb_pagepool_t* pgp, void* mem, 
			unsigned int flags)
{
  ffdb_bkt_t* bp;
  ffdb_bkt_waiter_t* sleeper = 0;

  FFDB_LOCK(pgp->lock);
#ifdef _FFDB_STATISTICS
  pgp->pageput++;
#endif
  
  bp = (ffdb_bkt_t *)((char *)mem - sizeof (ffdb_bkt_t));

  if (!FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_PINNED)) {
    fprintf (stderr, "ffdb_pagepool_put_page: page %d is not pinned.\n",
	     bp->pgno);
#ifdef _FFDB_STATISTICS
    ffdb_pagepool_stat (pgp);
#endif
    ffdb_dump_stack ();
    abort();
  }

  /**
   * Here is a scary thing, if the flags has DIRTY bit set,
   * but bp->flags does not have the DIRTY bit set,
   * how do we control multiple threads simultaneous writes
   * suggestion: always set flag when you create pages
   */
  if (FFDB_FLAG_ISSET(flags, FFDB_PAGE_DIRTY))
    FFDB_FLAG_SET(bp->flags, FFDB_PAGE_DIRTY);

  /**
   * Derefence the page
   */
  bp->ref--;
  /*
   * I am giving up the ownership
   */
  FFDB_THREAD_NULL(bp->owner);
  
  /** 
   * Now Unpin the page and wake up other threads waiting 
   */
  FFDB_FLAG_CLR(bp->flags, FFDB_PAGE_PINNED);

  if (bp->waiters > 0) {
    sleeper = FFDB_CIRCLEQ_LAST(&bp->wqh);
    sleeper->wakeup = 0xdeafbeaf;
  }

  FFDB_UNLOCK(pgp->lock);

  /* bp->waiters will not be changed until other threads are waken up */
#if 0
  fprintf (stderr, "Signaling page %d free \n", bp->pgno);
#endif
  if (sleeper)
    FFDB_COND_SIGNAL(sleeper->cv);
  return 0;
}


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
int
ffdb_pagepool_change_page (ffdb_pagepool_t* pgp, void* mem, 
			   pgno_t newpagenum)
{
  struct _ffdb_hqh* head;
  ffdb_bkt_t* bp;
  unsigned int oldpagenum;

  FFDB_LOCK(pgp->lock);
#ifdef _FFDB_STATISTICS
  pgp->pagechange++;
#endif
  
  bp = (ffdb_bkt_t *)((char *)mem - sizeof (ffdb_bkt_t));

  if (!FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_PINNED)) {
    fprintf (stderr, "ffdb_pagepool_put_page: page %d is not pinned.\n",
	     bp->pgno);
#ifdef _FFDB_STATISTICS
    ffdb_pagepool_stat (pgp);
#endif
    abort();
  }

  /**
   * if there are waiters, we cannot change this page 
   * just do a regular put and come back later to change
   */
  if (bp->waiters > 0) {
    fprintf (stderr, "ffdb_pagepool_change: page %d has waiters\n",bp->pgno);
    FFDB_UNLOCK(pgp->lock);
    
    ffdb_pagepool_put_page (pgp, mem, 0);
    return EAGAIN;
  }

  /**
   * Here is a scary thing, if the flags has DIRTY bit set,
   * but bp->flags does not have the DIRTY bit set,
   * how do we control multiple threads simultaneous writes
   * suggestion: always set flag when you create pages
   */
  FFDB_FLAG_SET(bp->flags, FFDB_PAGE_DIRTY);

  /**
   * Derefence the page
   */
  bp->ref--;

  /* only thread holding this page can delete this page, and no other
   * threads waiting on this page 
   */
  head = &pgp->hqh[FFDB_HASHKEY(bp->pgno)];

  /* Remove from the hash and lru queues. */
  FFDB_CIRCLEQ_REMOVE(head, bp, hq);
  FFDB_CIRCLEQ_REMOVE(&pgp->lqh, bp, lq);

  /**
   * Remember old pagenumber
   */
  oldpagenum = bp->pgno;

  /**
   * Change page number
   */
  bp->pgno = newpagenum;

  /*
   * I am giving up the ownership
   */
  FFDB_THREAD_NULL(bp->owner);
  
  /** 
   * Now Unpin the page and wake up other threads waiting 
   */
  FFDB_FLAG_CLR(bp->flags, FFDB_PAGE_PINNED);

  /**
   * Now add this entry back to the lists
   */
  head = &pgp->hqh[FFDB_HASHKEY(bp->pgno)];
  FFDB_CIRCLEQ_INSERT_HEAD(head, bp, hq);
  FFDB_CIRCLEQ_INSERT_TAIL(&pgp->lqh, bp, lq);

  /**
   * Change number of pages if pages are moved back
   */
  if (newpagenum >= pgp->npages) {
    pgp->npages = newpagenum + 1;
#if 0
    fprintf (stderr, "reset number of pages = %d\n", pgp->npages);
#endif
  }

  /**
   * Clean out old disk content
   */
  _ffdb_clean_page_ondisk (pgp, oldpagenum);

  FFDB_UNLOCK(pgp->lock);
  return 0;
}

/**
 * Flush number of pages to disk
 * If numpages == 0, flush all dirty pages to disk
 */
static int
_ffdb_pagepool_sync_i (ffdb_pagepool_t* pgp, unsigned int numpages)
{
  unsigned int num;
  ffdb_bkt_t* bp;
  ffdb_sbkt_t* sbp;
  ffdb_sbkt_t* next;
  ffdb_slh_t slh;
  FFDB_SLIST_INIT (&slh);

  /* Walk through every bucket and check whether it is pinned
   * If it is not pinned and it is dirty, I will sort these pages
   * according to page number
   */
  num = 0;
  FFDB_CIRCLEQ_FOREACH(bp, &pgp->lqh, lq){
    if (!FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_PINNED) &&
	bp->waiters == 0 && 
	FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_DIRTY)) {
      /* insert this bucket into a single linked list */
      sbp = (ffdb_sbkt_t *)malloc(sizeof(ffdb_sbkt_t));
      if (!sbp) {
	fprintf (stderr, "ffdb_pagepool_sync: cannot allocate space for single list element.\n");
	abort ();
      }
      _ffdb_shallow_copy_bk (sbp, bp);
      /* add this to the list */
      FFDB_SLIST_INSERT_HEAD (&slh, sbp, sl);
      num++;
      if (numpages > 0 && num >= numpages)
	break;
    }
  }

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Flushed %d pages out\n", num);
#endif

  /* Do a merge sort on the list slh according to pageno */
  _ffdb_slist_merge_sort (&slh);

  /* Now walk through the sorted list, and dump pages to the back end file */
  sbp = FFDB_SLIST_FIRST(&slh);
  next = 0;
  while (sbp) {
    next = FFDB_SLIST_NEXT(sbp, sl);

    if (FFDB_FLAG_ISSET(sbp->bp->flags, FFDB_PAGE_DIRTY)) {

      if (_ffdb_pagepool_write (pgp, sbp->bp) != 0) {
	fprintf (stderr, "ffdb_pagepool_sync: writing page %d error.\n",
		 sbp->bp->pgno);
	return -1;
      }
    }
#ifdef _FFDB_STATISTICS
    ++pgp->pageflush;
#endif
    
    /* Free memory of each simple bucket */
    free (sbp);
    sbp = next;
  }

  return 0;
}


/**
 * Flush all dirty pages back to the back end file. However, if any modified 
 * pages are in use. They will be ignored
 *
 * @param  pgp cache page pool pointer
 */
int
ffdb_pagepool_sync (ffdb_pagepool_t* pgp)
{
  int ret;

  FFDB_LOCK (pgp->lock);

  ret = _ffdb_pagepool_sync_i (pgp, 0);

  FFDB_UNLOCK (pgp->lock);

  return ret;
}



/**
 * Close the page poll pointer and any resource associated with this file
 * This implies all dirty pages are flushed out, 
 *
 * @param pgp cache page poll pointer
 * @return 0 the file is closed successfully and all pages are written back
 * to th back file. If there are still pinned pages by other threads, return -1
 * and errno is set to EAGAIN.
 */
int
ffdb_pagepool_close (ffdb_pagepool_t* pgp)
{
  ffdb_bkt_t* bp;

  /* First Sync Everything to disk */
  ffdb_pagepool_sync (pgp);

  FFDB_LOCK(pgp->lock);
  
  /* Free Every BUCKET */
  bp = FFDB_CIRCLEQ_FIRST(&pgp->lqh);  
  while (!FFDB_CIRCLEQ_EMPTY(&pgp->lqh)) {
    FFDB_CIRCLEQ_REMOVE(&pgp->lqh, bp, lq);

    free (bp);
    bp = FFDB_CIRCLEQ_FIRST(&pgp->lqh);  
  }

  /* close file descriptor */
  if (pgp->close_fd)
    close (pgp->fd);
  
  FFDB_UNLOCK(pgp->lock);  

  /* destroy lock */
  FFDB_LOCK_FINI(pgp->lock);

  free (pgp);
  return 0;
}

/**
 * Really delete page from cache. This is used either by internal
 * other routines or user knows this page is not used anymore
 * @param pgp cache page poll pointer
 * @param mem memory pointer of this page
 * @return 0 this page is deleted and memory is freed. Otherwise return -1
 */
int
ffdb_pagepool_delete (ffdb_pagepool_t* pgp, void* mem)
{
  struct _ffdb_hqh* head;
  ffdb_bkt_t* bp;
  int ret = 0;

  /* I have to lock this routine to prevent race condition
   * to ffdb_pagepool_find
   */
  FFDB_LOCK(pgp->lock);

  /* first get the page bucket pointer of this memory */
  bp = (ffdb_bkt_t *)((char *)mem - sizeof(ffdb_bkt_t));

  /* only thread holding this page can delete this page, and no other
   * threads waiting on this page 
   */
  head = &pgp->hqh[FFDB_HASHKEY(bp->pgno)];

  if (bp->ref == 1) {
    /* sanity check: page pin flag must be set */
    if (!FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_PINNED)) {
      fprintf (stderr, "ffdb_pagepool_delete: this page is not pinned\n");
      abort ();
    }

    /* if there are waiters, we cannot delete this page we need put
     */
    if (bp->waiters > 0) {
      fprintf (stderr, "ffdb_pagepool_delete: page %d has waiters\n",bp->pgno);
      FFDB_UNLOCK(pgp->lock);
      return EAGAIN;
    }
    
    /* Check whether this page is dirty. If it is , flush to disk */
    if (FFDB_FLAG_ISSET(bp->flags, FFDB_PAGE_DIRTY)) {
      if (_ffdb_pagepool_write (pgp, bp) != 0) {
	fprintf (stderr, "ffdb_pagepool_delete: page %d is dirty and flush it to disk encountered error.\n", bp->pgno);
	FFDB_UNLOCK(pgp->lock);
	return errno;
      }
    }
    
    /* Remove from the hash and lru queues. */

    FFDB_CIRCLEQ_REMOVE(head, bp, hq);
    FFDB_CIRCLEQ_REMOVE(&pgp->lqh, bp, lq);

    /* Decrease number of pages in the cache */
    --pgp->curcache;

    FFDB_UNLOCK(pgp->lock);
    /* free memory */
    free(bp); 

  }
  else {
    fprintf (stderr, "ffdb_pagepool_delete: other threads are still using this page %d ref = %d\n", bp->pgno, bp->ref);

    FFDB_UNLOCK(pgp->lock);

    ret = -1;
  }
  return ret;
}

/**
 * User supplied filter code
 */
void
ffdb_pagepool_filter (ffdb_pagepool_t* pgp, 
		      ffdb_pgiofunc_t pgin, ffdb_pgiofunc_t pgout,
		      void* arg)
{
  FFDB_LOCK(pgp->lock);

  pgp->pgin = pgin;
  pgp->pgout = pgout;
  pgp->pgcookie = arg;

  FFDB_UNLOCK(pgp->lock);
}


#ifdef _FFDB_STATISTICS
void
ffdb_pagepool_stat (ffdb_pagepool_t* pgp)
{
  ffdb_bkt_t *bp;
  int cnt;
  char *sep;
  ffdb_sbkt_t* sbp;
  ffdb_sbkt_t* next;
  ffdb_slh_t slh;
  FFDB_SLIST_INIT (&slh);

  fprintf(stderr, "%u pages in the file\n", pgp->npages);
  fprintf(stderr,
		"page size %u, cacheing %u pages of %u page max cache\n",
		pgp->pagesize, pgp->curcache, pgp->maxcache);
  fprintf(stderr, "%u page puts, %u page gets, %u page new\n",
		pgp->pageput, pgp->pageget, pgp->pagenew);
  fprintf(stderr, "%u page allocs, %u page reuse, %u page swap, %u page flushes\n",
	  pgp->pagealloc, pgp->pagereuse, pgp->pageswap, pgp->pageflush);
  if (pgp->cachehit + pgp->cachemiss)
    fprintf(stderr,
		  "%.0f%% cache hit rate (%u hits, %u misses)\n", 
		  ((double)pgp->cachehit / (pgp->cachehit + pgp->cachemiss))
		  * 100, pgp->cachehit, pgp->cachemiss);
  fprintf(stderr, "%u page reads, %u page writes\n",
	  pgp->pageread, pgp->pagewrite);
  
  sep = "";
  cnt = 0;
  FFDB_CIRCLEQ_FOREACH(bp, &pgp->lqh, lq) {
    /* insert this bucket into a single linked list */
    sbp = (ffdb_sbkt_t *)malloc(sizeof(ffdb_sbkt_t));
    if (!sbp) {
      fprintf (stderr, "ffdb_pagepool_sync: cannot allocate space for single list element.\n");
      abort ();
    }
    _ffdb_shallow_copy_bk (sbp, bp);
    /* add this to the list */
    FFDB_SLIST_INSERT_HEAD (&slh, sbp, sl);
  }
  /* Do a merge sort on the list slh according to pageno */
  _ffdb_slist_merge_sort (&slh);

  /* Now walk through the sorted list, and dump pages to the back end file */
  sbp = FFDB_SLIST_FIRST(&slh);
  next = 0;
  while (sbp) {
    next = FFDB_SLIST_NEXT(sbp, sl);

    fprintf(stderr, "%s%d", sep, sbp->bp->pgno);
    if (sbp->bp->flags & FFDB_PAGE_DIRTY)
      fprintf(stderr, "d");
    if (sbp->bp->flags & FFDB_PAGE_PINNED)
      fprintf(stderr, "P");
    if (sbp->bp->flags & FFDB_PAGE_LOCKED)
      fprintf(stderr, "L");
    if (++cnt == 10) {
      sep = "\n";
      cnt = 0;
    } else
      sep = ", ";
    
    free (sbp);
    sbp = next;
  }
  fprintf(stderr, "\n");
}
#endif

