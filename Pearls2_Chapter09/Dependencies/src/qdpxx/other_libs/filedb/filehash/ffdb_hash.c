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
 *     $Log: ffdb_hash.c,v $
 *     Revision 1.9  2009-06-03 02:26:39  edwards
 *     Turn off some debugging.
 *
 *     Revision 1.8  2009/05/08 17:37:31  chen
 *     Fix a major bug (not clean out moved page inside page pool cache)
 *
 *     Revision 1.7  2009/05/05 17:55:07  chen
 *     fix hdr->nkeys which over count when replace an existing keys
 *
 *     Revision 1.6  2009/04/21 18:51:19  chen
 *     Fix bugs related to number of pages upon moving pages in addition to clean pages on disk when the pages has been moved
 *
 *     Revision 1.5  2009/03/14 15:41:36  edwards
 *     Turned off some debugging output.
 *
 *     Revision 1.4  2009/03/04 19:12:28  edwards
 *     Renamed DB_HASH and __db to avoid name collisions with Berkeley DB.
 *
 *     Revision 1.3  2009/03/04 18:03:13  chen
 *     Add flush disk when doubling
 *
 *     Revision 1.2  2009/03/02 23:27:26  chen
 *     Test DBMerge Code
 *
 *     Revision 1.1  2009/02/20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _FFDB_DEBUG
#include <assert.h>
#endif

#include "ffdb_db.h"
#include "ffdb_pagepool.h"
#include "ffdb_page.h"
#include "ffdb_hash_func.h"
#include "ffdb_hash.h"


#ifdef _FFDB_STATISTICS
unsigned int hash_accesses, hash_collisions, hash_expansions, hash_overflows, hash_bigpages;
#endif


/**
 * Forward decleration of various routines needed for DB
 */
static int _ffdb_hash_close  (FFDB_DB *);
static int _ffdb_hash_delete (const FFDB_DB *, const FFDB_DBT *, 
			      unsigned int);
static int _ffdb_hash_fd     (const FFDB_DB *);
static int _ffdb_hash_get    (const FFDB_DB *, const FFDB_DBT *, FFDB_DBT *, 
			      unsigned int);
static int _ffdb_hash_put    (const FFDB_DB *, FFDB_DBT *, const FFDB_DBT *, 
			      unsigned int);
static int _ffdb_hash_sync   (const FFDB_DB *, unsigned int);
static int _ffdb_hash_cursor (const FFDB_DB *, ffdb_cursor_t**, unsigned int);
static int _ffdb_cursor_get  (ffdb_cursor_t* cursor, FFDB_DBT* key,
			      FFDB_DBT* data, unsigned int flags);
static int _ffdb_cursor_close (ffdb_cursor_t* cursor);

/**
 * Test big endian or little endian
 */
static int
_ffdb_test_endian (void)
{
  union {
    int l;
    char c[sizeof(int)];
  }u;

  u.l = 1;
  if (u.c[sizeof(int) - 1] == 1)
    return BIG_ENDIAN;

  return LITTLE_ENDIAN;
}

/**
 * Swap headers only if the byte order is little endian
 */
static void
_ffdb_swap_header (ffdb_htab_t* hashp)
{
  ffdb_hashhdr_t *hdrp;
  int i;

  hdrp = &hashp->hdr;

  M_32_SWAP(hdrp->magic);
  M_32_SWAP(hdrp->version);
  M_32_SWAP(hdrp->lorder);
  M_32_SWAP(hdrp->bsize);
  M_32_SWAP(hdrp->bshift);
  M_32_SWAP(hdrp->ovfl_point);
  M_32_SWAP(hdrp->max_bucket);
  M_32_SWAP(hdrp->high_mask);
  M_32_SWAP(hdrp->low_mask);
  M_32_SWAP(hdrp->ffactor);
  M_32_SWAP(hdrp->nkeys);
  M_32_SWAP(hdrp->hdrpages);
  M_32_SWAP(hdrp->uinfolen);
  M_16_SWAP(hdrp->uinfo_page);
  M_16_SWAP(hdrp->uinfo_npages);
  M_16_SWAP(hdrp->cfig_page);
  M_16_SWAP(hdrp->cfig_npages);
  M_32_SWAP(hdrp->num_cfigs);
  M_32_SWAP(hdrp->h_charkey);
  M_32_SWAP(hdrp->num_moved_pages);
  for (i = 0; i < NCACHED; i++) 
    M_32_SWAP(hdrp->spares[i]);
  for (i = 0; i < NCACHED; i++) 
    M_32_SWAP(hdrp->free_pages[i]);
  M_32_SWAP(hdrp->chksum);
}


/*
 * Hashp->hdr needs to be byteswapped.
 */
static void
_ffdb_swap_header_copy(ffdb_hashhdr_t* srcp, ffdb_hashhdr_t* destp)
{
  int i;

  P_32_COPY(srcp->magic, destp->magic);
  P_32_COPY(srcp->version, destp->version);
  P_32_COPY(srcp->lorder, destp->lorder);
  P_32_COPY(srcp->bsize, destp->bsize);
  P_32_COPY(srcp->bshift, destp->bshift);
  P_32_COPY(srcp->ovfl_point, destp->ovfl_point);
  P_32_COPY(srcp->max_bucket, destp->max_bucket);
  P_32_COPY(srcp->high_mask, destp->high_mask);
  P_32_COPY(srcp->low_mask, destp->low_mask);
  P_32_COPY(srcp->ffactor, destp->ffactor);
  P_32_COPY(srcp->nkeys, destp->nkeys);
  P_32_COPY(srcp->hdrpages, destp->hdrpages);
  P_32_COPY(srcp->uinfolen, destp->uinfolen);
  P_16_COPY(srcp->uinfo_page, destp->uinfo_page);
  P_16_COPY(srcp->uinfo_npages, destp->uinfo_npages);
  P_16_COPY(srcp->cfig_page, destp->cfig_page);
  P_16_COPY(srcp->cfig_npages, destp->cfig_npages);
  P_32_COPY(srcp->num_cfigs, destp->num_cfigs);
  P_32_COPY(srcp->h_charkey, destp->h_charkey);
  P_32_COPY(srcp->num_moved_pages, destp->num_moved_pages);
  for (i = 0; i < NCACHED; i++) 
    P_32_COPY(srcp->spares[i], destp->spares[i]);
  for (i = 0; i < NCACHED; i++) 
    P_32_COPY(srcp->free_pages[i], destp->free_pages[i]);
  P_32_COPY(srcp->chksum, destp->chksum);
}

/**
 * Flush out header onto disk
 */
static int
_ffdb_hput_header (ffdb_htab_t* hashp)
{
  ffdb_hashhdr_t *whdrp;
  ffdb_hashhdr_t whdr;
  unsigned int num_copied = 0;
  unsigned int chksum = 0;

  whdrp = &hashp->hdr;

  if (hashp->mborder == LITTLE_ENDIAN) {
    whdrp = &whdr;
    _ffdb_swap_header_copy(&hashp->hdr, whdrp);
  }

  /* calculate checksum value */
  chksum = __ffdb_crc32_checksum (chksum, (const unsigned char *)whdrp,
				  sizeof(ffdb_hashhdr_t) - sizeof(unsigned int));
  /* hashp->hdr.chksum = chksum; */
  if (hashp->mborder == LITTLE_ENDIAN) 
    M_32_SWAP(chksum);
  whdrp->chksum = chksum;


  /* write the header */
  lseek(hashp->fp, 0, SEEK_SET);
  num_copied = write(hashp->fp, whdrp, sizeof(ffdb_hashhdr_t));
  if (num_copied != sizeof(ffdb_hashhdr_t)) {
    fprintf(stderr, "hash: could not write hash header");
    return -1;
  }

  return 0;
}


/**
 * Flush meta header information to backend file
 */
static int
_ffdb_flush_meta (ffdb_htab_t* hashp)
{
  /* if this file does not need to be saved, do not do anything */
  if (!hashp->save_file)
    return 0;

  /* just rewrite the first three items */
  hashp->hdr.magic = FFDB_HASHMAGIC;
  hashp->hdr.version = FFDB_HASHVERSION;
  hashp->hdr.h_charkey = hashp->hash(CHARKEY, sizeof(CHARKEY));

  return _ffdb_hput_header (hashp);
}


/**
 * This is the routine called by hash close
 * lock is held by other routines
 */
static int
_ffdb_hdestroy (ffdb_htab_t* hashp)
{
  int save_errno = 0;

  if (hashp->rearrange_pages && hashp->save_file)
    ffdb_rearrage_pages_on_close (hashp);

#ifdef _FFDB_STATISTICS
  { 
    int i;
    fprintf(stderr, "hdestroy: accesses %u collisions %u\n",
	    hash_accesses, hash_collisions);
    fprintf(stderr,
	    "hdestroy: expansions %u\n", hash_expansions);
    fprintf(stderr,
	    "hdestroy: overflows %u\n", hash_overflows);
    fprintf(stderr,
	    "hdestroy: big key/data pages %u\n", hash_bigpages);
    fprintf(stderr,
	    "keys %u max bucket %d\n", hashp->hdr.nkeys, hashp->hdr.max_bucket);

    fprintf (stderr,
	     " Number of pages moved %d\n", hashp->hdr.num_moved_pages);
    
    for (i = 0; i < NCACHED; i++)
      fprintf(stderr,
		    "spares[%d] = %d\n", i, hashp->hdr.spares[i]);

    for (i = 0; i < NCACHED; i++)
      fprintf(stderr,
		    "freepages[%d] = %d\n", i, hashp->hdr.free_pages[i]);

    fprintf (stderr, "Average ffactor = %d\n", hashp->hdr.nkeys/(hashp->hdr.max_bucket + 1));
  }
#endif

  /* flush meta information header to disk */
  if (_ffdb_flush_meta (hashp) != 0)
    if (save_errno == 0)
      save_errno = errno;

  /* close the pagepool */
#ifdef _FFDB_STATISTICS
  ffdb_pagepool_stat (hashp->mp);
#endif
  ffdb_pagepool_sync (hashp->mp);
  ffdb_pagepool_close (hashp->mp);

  /* Reduce file size if possible */
  if (hashp->rearrange_pages)
    ffdb_reduce_filesize (hashp);

  if (hashp->fp != -1)
    close (hashp->fp);

  if (hashp->fname)
    free (hashp->fname);

  /* free allocated memory */
  if (hashp->split_buf)
    free(hashp->split_buf);
  if (hashp->bigdata_buf)
    free(hashp->bigdata_buf);
  if (hashp->bigkey_buf)
    free(hashp->bigkey_buf);

  if (save_errno) {
    errno = save_errno;
    return -1;
  }
  return 0;
}

/**
 * Initialize hashtable values
 * return 0 on sucess, errno or -1 on failure
 */
static int
_ffdb_init_hash(ffdb_htab_t* hashp, const char* file, 
		FFDB_HASHINFO* info)
{
  struct stat statbuf;

  /* This is the machine we are writing the database on
   * and we are using the native format
   */
  hashp->hdr.lorder = _ffdb_test_endian ();
#if 0
  if (hashp->hdr.lorder == LITTLE_ENDIAN) 
    fprintf (stderr, "Info: Data will be stored in little endian format.\n");
  else
    fprintf (stderr, "Info: Data will be stored in big endian format.\n");
#endif

  hashp->hdr.nkeys = 0;
  hashp->hdr.bsize = DEF_BUCKET_SIZE;
  hashp->hdr.bshift = DEF_BUCKET_SHIFT;
  hashp->hdr.ffactor = DEF_FFACTOR;
  hashp->hash = __ffdb_default_hash;
  hashp->h_compare = __ffdb_default_cmp;
  memset(hashp->hdr.spares, 0, sizeof(hashp->hdr.spares));

  /* Fix bucket size to be optimal for file system */
  /* This is not very reliable for NFS based file system */
  if (stat(file, &statbuf) != 0)
    return errno;
  hashp->hdr.bsize = statbuf.st_blksize;
  hashp->hdr.bshift = __ffdb_log2(hashp->hdr.bsize);
	
  if (info) {
    if (info->bsize) {
      /* Round pagesize up to power of 2 */
      hashp->hdr.bshift = __ffdb_log2(info->bsize);
      hashp->hdr.bsize = 1 << hashp->hdr.bshift;
      if (hashp->hdr.bsize > MAX_BSIZE) {
	errno = EINVAL;
	return errno;
      }
      /* the header occupies one page    */
      if (hashp->hdr.bsize < sizeof(ffdb_hashhdr_t)) {
	fprintf (stderr, "Hash page is too small to handle meta info.\n");
	errno = EINVAL;
	return errno;
      }
    }

    /* Set rearrange page flag */
    hashp->rearrange_pages = info->rearrangepages;
    
    if (info->hash)
      hashp->hash = info->hash;
    if (info->cmp)
      hashp->h_compare = info->cmp;
  }

  return 0;
}


/**
 * Functions to get/put hash header.  We access the file directly.
 */
unsigned int
_ffdb_hget_header(ffdb_htab_t *hashp, unsigned int page_size)
{
  unsigned num_copied, i, newchksum;
  unsigned char *hdr_dest;

  num_copied = 0;
  i = 0;

  hdr_dest = (unsigned char *)&hashp->hdr;
  
  /* 
   * XXX
   * This should not be printing to stderr on a "normal" error case.
   */
  lseek(hashp->fp, 0, SEEK_SET);
  num_copied = read(hashp->fp, hdr_dest, sizeof(ffdb_hashhdr_t));
  if (num_copied != sizeof(ffdb_hashhdr_t)) {
    fprintf(stderr, "Fatal error : hash could not retrieve header");
    return 0;
  }

  /* Now I need compare checksum value */
  newchksum = 0;
  newchksum = __ffdb_crc32_checksum (newchksum, hdr_dest,
				     sizeof(ffdb_hashhdr_t) - sizeof(unsigned int));

  if (hashp->mborder == LITTLE_ENDIAN)
    _ffdb_swap_header(hashp);

  if (newchksum != hashp->hdr.chksum) {
    fprintf (stderr, "Hash Meta Data Checksum error!!!!!! 0x%x != 0x%x (retrieved)\n", newchksum, hashp->hdr.chksum);
    exit (1);
  }

#if 0
  if (hashp->hdr.lorder == LITTLE_ENDIAN) 
    fprintf (stderr, "Info: Data are stored in little endian format.\n");
  else
    fprintf (stderr, "Info: Data are stored in big endian format.\n");
#endif

  return num_copied;
}


/**
 * Read in pages containing user information and configuration information
 */
static int
_ffdb_read_user_info (ffdb_htab_t *hashp)
{
  unsigned short i;
  void *mem;
  pgno_t tpage = INVALID_PGNO;
  int   dirty = 0;
  unsigned int flags = 0;
  pgno_t start = hashp->hdr.uinfo_page;
  unsigned short npages = hashp->hdr.uinfo_npages;

  if (hashp->new_file) {
    for (i = 0; i < npages; i++)
      ffdb_new_page (hashp, start + i, HASH_UINFO_PAGE);
  }

  if (hashp->save_file) {
    flags = FFDB_PAGE_DIRTY | FFDB_PAGE_LOCKED;
    dirty = 1;
  }

  for (i = 0; i < npages; i++) {
    mem = ffdb_get_page (hashp, start + i, HASH_UINFO_PAGE,
			 flags, &tpage);
    ffdb_put_page (hashp, mem, HASH_UINFO_PAGE, dirty);
  }
  return 0;
}


static int
_ffdb_read_config_info (ffdb_htab_t *hashp)
{
  unsigned short i;
  void *mem;
  pgno_t tpage = INVALID_PGNO;
  pgno_t start = hashp->hdr.cfig_page;
  unsigned short npages = hashp->hdr.cfig_npages;
  int   dirty = 0;
  unsigned int flags = 0;


  if (hashp->new_file) {
    for (i = 0; i < npages; i++)
      ffdb_new_page (hashp, start + i, HASH_CONFIG_PAGE);
  }

  if (hashp->save_file) {
    flags = FFDB_PAGE_DIRTY | FFDB_PAGE_LOCKED;
    dirty = 1;
  }


  for (i = 0; i < npages; i++) {
    mem = ffdb_get_page (hashp, start + i, HASH_CONFIG_PAGE,
			 flags, &tpage);
    ffdb_put_page (hashp, mem, HASH_CONFIG_PAGE, dirty);
  }
  return 0;
}

/**
 * Calculate how many pages configuration information
 * occupies
 */
static unsigned int
_ffdb_num_config_info_pages (ffdb_htab_t* hashp, unsigned int num)
{
  unsigned int npages, numelems;

  npages = numelems = 0;
  while (num > 0) {
    /* number of configurations on this page */
    numelems = (hashp->hdr.bsize - PAGE_OVERHEAD)/(sizeof(ffdb_config_info_t));
    if (numelems > num)
      numelems = num;
    num -= numelems;
    npages++;
  } 
  return npages;
}


/**
 * Initialize the hash table if this is a new file
 * Initialize hash table with 'nb' number of buckets
 */
static int
_ffdb_init_htab (ffdb_htab_t* hashp, int nb, FFDB_HASHINFO* info)
{
  int l2, nbuckets, nhdr_pages;
  unsigned int uinfosize;
  pgno_t i, pg;

  /**
   * log2 of elements
   */
  l2 = __ffdb_log2 (nb);

  /**
   * Number of buckets = 2^l2
   */
  nbuckets = 1 << l2;


  /* Calculate number of header pages */
  /* real header occupies one */
  pg = nhdr_pages = 1;

  /* Check how many user informantion length */
  if (info && info->userinfolen > 0) 
    hashp->hdr.uinfolen = info->userinfolen;
  else
    hashp->hdr.uinfolen = FFDB_DEFAULT_UINFO_LEN;
  uinfosize = hashp->hdr.uinfolen + sizeof(unsigned int);

  /* This is user information */
  hashp->hdr.uinfo_page = pg;

  if (uinfosize % (hashp->hdr.bsize - PAGE_OVERHEAD) == 0)
    nhdr_pages += (uinfosize/(hashp->hdr.bsize - PAGE_OVERHEAD));
  else
    nhdr_pages += (uinfosize/(hashp->hdr.bsize - PAGE_OVERHEAD) + 1);
  hashp->hdr.uinfo_npages = nhdr_pages - pg;

  /* Calculate configuration number of pages */
  pg = nhdr_pages;
  hashp->hdr.cfig_page = pg;
  if (info && info->numconfigs > 0) 
    hashp->hdr.num_cfigs = info->numconfigs;
  else
    hashp->hdr.num_cfigs = 1;

  /* total number of configuration size */
  nhdr_pages += (_ffdb_num_config_info_pages (hashp, hashp->hdr.num_cfigs));
  hashp->hdr.cfig_npages = nhdr_pages - pg;

  /* Assign the number of header pages */
  hashp->hdr.hdrpages = nhdr_pages;

  /* populate spares */
  hashp->hdr.spares[0] = nhdr_pages;

  if (l2 >= 1)
    hashp->hdr.spares[1] = nhdr_pages + 1;
  for (i = 2; i <= l2; i++)
    hashp->hdr.spares[i] = hashp->hdr.spares[i - 1] + POW2(i - 2);

  for (; i < NCACHED; i++)
    hashp->hdr.spares[i] = 0;

  /* initialize free pages */
  for (i = 0; i < NCACHED; i++)
    hashp->hdr.free_pages[i] = INVALID_PGNO;

  /* current spliting level */
  hashp->hdr.ovfl_point = l2;

  /* high mask and low mask value */
  hashp->hdr.max_bucket = hashp->hdr.high_mask = nbuckets - 1;
  hashp->hdr.low_mask = (nbuckets >> 1) - 1;

  return 0;
}



/**
 * Interface to outside DB calls
 */
FFDB_DB *
__ffdb_hash_open (const char* fname, int flags, int mode, const void* arg)
{
  struct stat statbuf;
  FFDB_DB *dbp;
  ffdb_htab_t *hashp;
  unsigned int csize;
  unsigned long tcsize;
  int ret, new_table;
  FFDB_HASHINFO *info = (FFDB_HASHINFO *)arg;

  /**
   * We are not going to have write only system :-)
   */
  if ((flags & O_ACCMODE) == O_WRONLY) {
    errno = EINVAL;
    return 0;
  }

  /**
   * We are not going to allow null filename
   */
  if (!fname) {
    fprintf (stderr, "A NULL filename is not allowed.\n");
    errno = EINVAL;
    return 0;
  }

  /* initialize crc32 check sum */
  __ffdb_crc32_init ();

  /**
   * Create the memory resident hash table
   */
  hashp = (ffdb_htab_t *)malloc(sizeof(ffdb_htab_t));
  if (!hashp) {
    fprintf (stderr, "Cannot allocate space to hash table memory object.\n");
    return 0;
  }
  memset (hashp, 0, sizeof(ffdb_htab_t));
  hashp->fp = -1;
  hashp->curr_dpage = INVALID_PGNO;
  hashp->rearrange_pages = 1;

  /* check this machine byte order */
  hashp->mborder = _ffdb_test_endian();

  /*
   * Even if user wants write only, we need to be able to read
   * the actual file, so we need to open it read/write. But, the
   * field in the hashp structure needs to be accurate so that
   * we can check accesses.
   */
  hashp->flags = flags;
  hashp->save_file = (hashp->flags & O_RDWR);

  /**
   * Check whether this is a new database
   */
  new_table = 0;
  if ((flags & O_TRUNC) ||(stat(fname, &statbuf) != 0 && (errno == ENOENT))) {
    if (errno == ENOENT)
      errno = 0;	/* In case someone looks at errno. */
    new_table = 1;
  }

  /**
   * Now open file according to flags and mode
   */
  if ((hashp->fp = open(fname, flags, mode)) == -1) {
    free (hashp);
    return 0;
  }
  
  /* File will be closed on execve call */
  (void)fcntl(hashp->fp, F_SETFD, 1);
  
  /* Process arguments to set up hash table header. */
  if (new_table) {
    if (_ffdb_init_hash(hashp, fname, info) != 0) {
      close (hashp->fp);
      free (hashp);
      return 0;
    }
  } 
  else {
    /* Table already exists */
    if (info && info->hash)
      hashp->hash = info->hash;
    else
      hashp->hash = __ffdb_default_hash;

    if (info && info->cmp)
      hashp->h_compare = info->cmp;
    else
      hashp->h_compare = __ffdb_default_cmp;

    if (info)
      hashp->rearrange_pages = info->rearrangepages;

    /* copy metadata from page into header */
    if (_ffdb_hget_header(hashp,(info && info->bsize ? info->bsize : DEF_BUCKET_SIZE)) != sizeof(ffdb_hashhdr_t)) {
      close (hashp->fp);
      free (hashp);
      errno = EFTYPE;
      return 0;
    }


    /* Verify file type, versions and hash function */
    if (hashp->hdr.magic != FFDB_HASHMAGIC) {
      close (hashp->fp);
      free (hashp);
      errno = EFTYPE;
      return 0;
    }
    
    /* check hash version */
    if (hashp->hdr.version != FFDB_HASHVERSION) {
      close (hashp->fp);
      free (hashp);
      errno = EFTYPE;
      return 0;
    }

    /* compare the calculated hash value and stored hash value */
    if (hashp->hash(CHARKEY, sizeof(CHARKEY)) != hashp->hdr.h_charkey) {
      close (hashp->fp);
      free (hashp);
      errno = EFTYPE;
      return 0;
    }
  }

  /**
   * Get ready to open pagepool
   * cachesize from info is in bytes. We converted into number of pages
   */
  if (info && info->cachesize) {
    // need to check the number of cache pages not exceeding the 
    // maximum number of integers
    tcsize = info->cachesize / hashp->hdr.bsize;
    if (tcsize > 0x7FFFFFFF) {
      fprintf (stderr, "Number of desired cache pages exceeds the maximum of integer due to rediculous large cache memory request of %lu MBytes.\n",
	       (info->cachesize >> 20));
      exit (1);
    }
    csize = tcsize;
  }
  else 
    csize = DEF_CACHESIZE / hashp->hdr.bsize;


  /**
   * Create page pool first
   */
  if ((ret = ffdb_pagepool_create (&hashp->mp, 0)) != 0) {
    fprintf (stderr, "Cannot open pagepool for %s\n", fname);
    close (hashp->fp);
    free (hashp);
    errno = ret;
    return 0;
  }
  
  /**
   * Open memory page pool
   */
  if ((ret = ffdb_pagepool_open_fd(hashp->mp, hashp->fp, 
				   hashp->hdr.bsize, csize)) != 0) {
    fprintf (stderr, "Cannot open pagepool for %s\n", fname);
    close (hashp->fp);
    free (hashp);
    errno = ret;
    return 0;
  }

  /**
   * Now add page in and out filter
   */
  ffdb_pagepool_filter(hashp->mp, ffdb_pgin_routine, ffdb_pgout_routine, hashp);

  /*
   * For a new table, set up the appropriate hashtable information
   */
  if (new_table) {
    if (_ffdb_init_htab(hashp, info && info->nbuckets ? info->nbuckets : 1, info) != 0) {
      fprintf (stderr, "Cannot initialize hash header information.\n");
      ffdb_pagepool_close (hashp->mp);
      close (hashp->fp);
      free (hashp);
      errno = ENOMEM;
      return 0;
    }
  }

  /**
   * Now save the file name 
   */
  hashp->fname = (char *)malloc (strlen (fname) + 1);
  if (!hashp->fname) {
    close (hashp->fp);
    free (hashp);
    return 0;
  }
  strcpy (hashp->fname, fname);

  /* initialize the cursor queue */
  FFDB_TAILQ_INIT(&hashp->curs_queue);
  hashp->seq_cursor = NULL;


  /* get a chunk of memory for our split buffer */
  hashp->split_buf = (unsigned short *)malloc(hashp->hdr.bsize);
  if (!hashp->split_buf) {
    ffdb_pagepool_close (hashp->mp);
    close (hashp->fp);
    free (hashp->fname);
    free (hashp);
    return 0;
  }

  hashp->new_file = new_table;


  /* always read in user information page and configuration page */
  _ffdb_read_user_info (hashp);
  _ffdb_read_config_info (hashp);

  if (!(dbp = (FFDB_DB *)malloc(sizeof(FFDB_DB)))) {
    ffdb_pagepool_close (hashp->mp);
    close (hashp->fp);
    free (hashp->fname);
    free (hashp->split_buf);
    free (hashp);
    return 0;
  }

  /* assign all appropriate routines */
  dbp->internal = hashp;
  dbp->close = _ffdb_hash_close;
  dbp->del = _ffdb_hash_delete;
  dbp->fd = _ffdb_hash_fd;
  dbp->get = _ffdb_hash_get;
  dbp->put = _ffdb_hash_put;
  dbp->sync = _ffdb_hash_sync;
  dbp->cursor = _ffdb_hash_cursor;
  dbp->type = FFDB_HASH;

#ifdef _FFDB_DEBUG
  (void)fprintf(stderr,
		"%s\n%s%p\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%x\n%s%x\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n",
		"init_htab:",
		"TABLE POINTER   ", (void *)hashp,
		"BUCKET SIZE     ", hashp->hdr.bsize,
		"BUCKET SHIFT    ", hashp->hdr.bshift,
		"FILL FACTOR     ", hashp->hdr.ffactor,
		"MAX BUCKET      ", hashp->hdr.max_bucket,
		"OVFL POINT      ", hashp->hdr.ovfl_point,
		"HIGH MASK       ", hashp->hdr.high_mask,
		"LOW  MASK       ", hashp->hdr.low_mask,
		"NKEYS           ", hashp->hdr.nkeys,
		"Num moved pages ", hashp->hdr.num_moved_pages,
		"Header Pages    ", hashp->hdr.hdrpages,
		"Uinfo start page ", hashp->hdr.uinfo_page,
		"Uinfo Num pages ",  hashp->hdr.uinfo_npages,
		"Config start at ",  hashp->hdr.cfig_page,
		"Config Num pages ", hashp->hdr.cfig_npages,
		"Number of configs ", hashp->hdr.num_cfigs,
		"Data Page Start ", hashp->curr_dpage);
  {
    unsigned int i;
    for (i = 0; i < NCACHED; i++)
      fprintf(stderr,
		    "spares[%d] = %d\n", i, hashp->hdr.spares[i]);

    for (i = 0; i < NCACHED; i++)
      fprintf(stderr,
		    "freepages[%d] = %d\n", i, hashp->hdr.free_pages[i]);
  }
#endif

#ifdef _FFDB_STATISTICS
  hash_overflows = hash_accesses = hash_collisions = hash_expansions = 0;
  hash_bigpages = 0;
#endif

  /* Check whether I need to move data pages back to original places if
   * new insert is expected
   */
  if (!new_table && hashp->save_file && hashp->rearrange_pages) 
    ffdb_rearrage_pages_on_open (hashp);

  /* Find out current data page which is the last data page we are using */
  /* Otherwise the datapage assigned to new insert will start at wrong place */
  if (hashp->hdr.spares[hashp->hdr.ovfl_point + 1] > 0) {
    hashp->curr_dpage = ffdb_last_data_page (hashp, hashp->hdr.spares[hashp->hdr.ovfl_point + 1] - 1);

#if 0
    fprintf (stderr, "Info: datapage will start at %d for level %d\n",
	     hashp->curr_dpage, hashp->hdr.ovfl_point);
#endif
  }


#if 0
  /**
   * Display all page information
   */
  ffdb_disp_all_page_info (hashp);
#endif


#if 0
  (void)fprintf(stderr,
		"%s\n%s%p\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%x\n%s%x\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n%s%d\n",
		"init_htab:",
		"TABLE POINTER   ", (void *)hashp,
		"BUCKET SIZE     ", hashp->hdr.bsize,
		"BUCKET SHIFT    ", hashp->hdr.bshift,
		"FILL FACTOR     ", hashp->hdr.ffactor,
		"MAX BUCKET      ", hashp->hdr.max_bucket,
		"OVFL POINT      ", hashp->hdr.ovfl_point,
		"HIGH MASK       ", hashp->hdr.high_mask,
		"LOW  MASK       ", hashp->hdr.low_mask,
		"NKEYS           ", hashp->hdr.nkeys,
		"Num moved pages ", hashp->hdr.num_moved_pages,
		"Header Pages    ", hashp->hdr.hdrpages,
		"Uinfo start page ", hashp->hdr.uinfo_page,
		"Uinfo Num pages ",  hashp->hdr.uinfo_npages,
		"Config start at ",  hashp->hdr.cfig_page,
		"Config Num pages ", hashp->hdr.cfig_npages,
		"Number of configs ", hashp->hdr.num_cfigs,
		"Data Page Start ", hashp->curr_dpage);
  {
    unsigned int i;
    for (i = 0; i < NCACHED; i++)
      fprintf(stderr,
		    "spares[%d] = %d\n", i, hashp->hdr.spares[i]);

    for (i = 0; i < NCACHED; i++)
      fprintf(stderr,
		    "freepages[%d] = %d\n", i, hashp->hdr.free_pages[i]);
  }
#endif

  /* finally intialize lock */
  FFDB_LOCK_INIT (hashp->lock);
  return dbp;
}


/*****************************************************************************
 * The following are routines handling i/o operations of database            *
 *****************************************************************************/

/**
 * Close the database
 */
static int
_ffdb_hash_close (FFDB_DB* dbp)
{
  ffdb_htab_t* hashp;
  int retval;
  
  if (!dbp) 
    return -1;
  hashp = (ffdb_htab_t *)dbp->internal;

  FFDB_LOCK(hashp->lock);

  retval = _ffdb_hdestroy (hashp);

  /* free lock */
  FFDB_LOCK_FINI(hashp->lock);

  /* finally free hashp itself */
  free (hashp);

  /* free memory */
  free (dbp);

  return retval;
}

/**
 * One of the most important function: hash function
 * return bucket number
 */
unsigned int
_ffdb_call_hash (ffdb_htab_t* hashp, const void* k, unsigned int len)
{
  unsigned int n, bucket;

  n = hashp->hash (k, len);
  bucket = (n & hashp->hdr.high_mask); /* n mod 2^(i + 1) */

  if (bucket > hashp->hdr.max_bucket) {
    bucket = (n & hashp->hdr.low_mask); /* n mod 2^i */
  }

  return bucket;
}

/**
 * Expand hash table
 * this happens when a bucket cannot hold any more keys
 *
 * We are not checking fill factor because we are trying to save
 * disk space instead of trying to speed up the access
 */
static int
_ffdb_expand_table (ffdb_htab_t* hashp)
{
  unsigned int old_bucket, new_bucket;
  int spare_indx, isdoubling, ret;
  pgno_t p;

#ifdef _FFDB_STATISTICS
  hash_expansions++;
#endif
  isdoubling = 0;

  /* The number of buckets is increased by one, obviously */
  new_bucket = ++hashp->hdr.max_bucket;
  
  /* which bucket to split */
  old_bucket = (hashp->hdr.max_bucket & hashp->hdr.low_mask);

  /* If new bucket is greater than high mask, we do doubling again */
  if (new_bucket > hashp->hdr.high_mask) {
    hashp->hdr.low_mask = hashp->hdr.high_mask;
    hashp->hdr.high_mask = new_bucket | hashp->hdr.low_mask;
  }

  /*
   * If the split point is increasing (hdr.max_bucket's log base 2
   * increases), we need to copy the current contents of the spare
   * split bucket to the next bucket.
   */
  spare_indx = __ffdb_log2(hashp->hdr.max_bucket + 1);

  if (spare_indx > hashp->hdr.ovfl_point) {
    /* update current data page value */
    hashp->curr_dpage = INVALID_PGNO;
    hashp->hdr.ovfl_point = spare_indx;
    isdoubling = 1;
  }


  BUCKET_TO_PAGE(new_bucket, p);
  if (p == MAX_PAGES(hashp)) {
    fprintf (stderr, "Reach maximum number of pages %d\n",
	     MAX_PAGES(hashp));
    exit (1);
  }
  
  /* Get a new bucket */
  if (ffdb_new_page (hashp, new_bucket, HASH_BUCKET_PAGE) != 0)
    return -1;

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Get a new expanded page at bucket %d split old bucket %d\n", new_bucket, old_bucket);
#endif

  if (isdoubling) {
    /* Write out meta header here */
    _ffdb_flush_meta (hashp);
    ffdb_pagepool_sync (hashp->mp); 
  }
  
  ret = ffdb_split_bucket (hashp, old_bucket, new_bucket, isdoubling);
  return ret;
}


/**
 * Get key from the database
 * returns 0: on success
 * returns 1: where the key not found
 * returns -1: there is an error
 * currently there is no flag is used
 */
static int
_ffdb_hash_get (const FFDB_DB* dbp, const FFDB_DBT* key,
		FFDB_DBT* data, unsigned int flag)
{
  ffdb_htab_t* hashp;
  ffdb_hent_t item;
  unsigned int bucket;
  int status;

#ifdef _FFDB_DEBUG
  fprintf (stderr, "ffdb_hash_get key %s\n", (char *)key->data);
#endif

  hashp = (ffdb_htab_t *)dbp->internal;

#ifdef _FFDB_STATISTICS
  FFDB_LOCK (hashp->lock);
  hash_accesses++;
  FFDB_UNLOCK (hashp->lock);
#endif

  /* initialize item */
  memset (&item, 0, sizeof (ffdb_hent_t));
  /* Calculate the hash item size */
  item.seek_size = PAIRSIZE(key, data);

  /* calculate hash value for this key */
  bucket = _ffdb_call_hash (hashp, key->data, key->size);
  item.bucket = bucket;
  
#ifdef _FFDB_DEBUG
  fprintf (stderr, "Hash get Key %s hash bucket = %d\n", (char *)key->data, bucket);
#endif

  /* Now I need find a page on which this key may reside */
  status = ffdb_find_item (hashp, (FFDB_DBT *)key, 0, &item);
  if (status != 0){ /* Something is really wrong */
    return -1;
  }

  if (item.status == ITEM_NO_MORE) {
    ffdb_release_item (hashp, &item);
    return FFDB_NOT_FOUND;
  }

  /* Now the item is found, item contains page information */
  /* page os released after the call */
  status = ffdb_get_item (hashp, key, data, &item, 1);

  return status;
}


/**
 * Put a key and data pair into the database
 */
static int
_ffdb_hash_put (const FFDB_DB* dbp, FFDB_DBT* key, const FFDB_DBT* data,
		unsigned int flag)
{
  ffdb_htab_t* hashp;
  ffdb_hent_t item;
  unsigned int bucket;
  int status, newkey;

#if 0
#ifdef _FFDB_DEBUG
  fprintf (stderr, "ffdb_hash_put key %s data %s\n", (char *)key->data,
	   (char *)data->data);
#endif
#endif

  hashp = (ffdb_htab_t *)dbp->internal;

#ifdef _FFDB_STATISTICS
  FFDB_LOCK (hashp->lock);
  hash_accesses++;
  FFDB_UNLOCK (hashp->lock);
#endif

  /* This is a new key */
  newkey = 1;

  /* check flags: if there is a flag it must be FFDB_NOOVERWRITE */
  if (flag && flag != FFDB_NOOVERWRITE) {
    FFDB_LOCK (hashp->lock);
    hashp->db_errno = errno = EINVAL;
    FFDB_UNLOCK (hashp->lock);
    return -1;
  }
  
  /* check file permission, if this is a read only file, cannot do it */
  if ((hashp->flags & O_ACCMODE) == O_RDONLY) {
    FFDB_LOCK (hashp->lock);
    hashp->db_errno = errno = EPERM;
    FFDB_UNLOCK (hashp->lock);
    return -1;
  }

  /* initialize item */
  memset (&item, 0, sizeof (ffdb_hent_t));
  /* Calculate the hash item size */
  item.seek_size = PAIRSIZE(key, data);

  /* calculate hash value for this key */
  bucket = _ffdb_call_hash (hashp, key->data, key->size);
  item.bucket = bucket;

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Key %s hash bucket = %d\n", (char *)key->data, bucket);
#endif

  /* Now I need find a page on which this key may reside */
  /* Find item is thread safe */
  status = ffdb_find_item (hashp, key, (FFDB_DBT *)data, &item);
  if (status != 0)  /* Something is really wrong */
    return status;


  FFDB_LOCK(hashp->lock);
  if (item.status == ITEM_NO_MORE) {
    /* There is no item found, we need to insert this item */
    /* Find out whether there is space on this page to fit this pair */
    if (PAIRFITS(item.pagep, key, data)) {
#ifdef _FFDB_DEBUG
      fprintf (stderr, "This data item bucket %d fit with page %d\n", bucket, item.pgno);
#endif
      if ((status = ffdb_add_pair (hashp, key, data, &item, 0)) != 0) {
	FFDB_UNLOCK (hashp->lock);  
	return status;
      }
    }
    else {
      /* Data will not fit on the page */
#ifdef _FFDB_DEBUG
      fprintf (stderr, "This data item will not fit in bucket %d with page %d: overflow page will be added\n", bucket, item.pgno);
#endif

      /* First chain an overflow page */
      if ((status = ffdb_add_ovflpage (hashp, key, data, &item)) != 0) {
	FFDB_UNLOCK (hashp->lock);  
	return status;
      }

      /* Now I need to expand the table even though the current bucket
       * may not be splited at this moment
       */
      if ((status = _ffdb_expand_table (hashp)) != 0) {
	hashp->hdr.nkeys++;
	FFDB_UNLOCK (hashp->lock);
	return status;
      }
    }
    /* Now I should have the data on the page */
  }
  else if (item.status == ITEM_OK) {
    newkey = 0;
    /* This item already exists. Check flag to make sure it is 
     * a replace flag
     */
    if (flag && flag == FFDB_NOOVERWRITE) {
      FFDB_UNLOCK (hashp->lock);  
      return -1;
    }

    if ((status = ffdb_add_pair (hashp, key, data, &item, 1)) != 0) {
      FFDB_UNLOCK (hashp->lock);  
      return -1;
    }
  }      	

  /* update number key information */
  if (newkey)
    hashp->hdr.nkeys++;

  FFDB_UNLOCK (hashp->lock);  
  return 0;
}


/**
 * Delete a key from the database
 */
static int
_ffdb_hash_delete (const FFDB_DB* dbp, const FFDB_DBT* key, unsigned int flag)
{
  return 0;
}


/**
 * Return file descriptor for the opened database
 */
static int
_ffdb_hash_fd (const FFDB_DB* dbp)
{
  ffdb_htab_t *hashp;

  if (!dbp) {
    return -1;
  }

  hashp = (ffdb_htab_t *)dbp->internal;
  if (hashp->fp == -1) {
    errno = ENOENT;
    return -1;
  }
  return hashp->fp;
}

/**
 * Sync the database to file (disk)
 */
static int
_ffdb_hash_sync (const FFDB_DB *dbp, unsigned int flags)
{
  ffdb_htab_t* hashp;
  
  if (!dbp) 
    return -1;
  hashp = (ffdb_htab_t *)dbp->internal;

  FFDB_LOCK(hashp->lock);

  /* flush meta information header to disk */
  if (_ffdb_flush_meta (hashp) != 0) {
    FFDB_UNLOCK(hashp->lock);
    return -1;
  }

  ffdb_pagepool_sync (hashp->mp);

  FFDB_UNLOCK(hashp->lock);
  return 0;
}

/**
 * Return a new cursor
 */
static int
_ffdb_hash_cursor (const FFDB_DB* dbp, ffdb_cursor_t** c, unsigned int type)
{
  ffdb_cursor_t *nc;
  ffdb_crs_t* inc;

  if (type != FFDB_KEY_CURSOR && type != FFDB_DATA_CURSOR) {
    fprintf (stderr, "provided cursor type is wrong.\n");
    errno = EINVAL;
    return -1;
  }
  nc = (ffdb_cursor_t *)malloc(sizeof(ffdb_cursor_t));
  if (!nc) {
    fprintf (stderr, "Cannot allocate space for a new cursor.\n");
    *c = 0;
    return -1;
  }
  /* Create internal cursor */
  inc = (ffdb_crs_t *)malloc(sizeof(ffdb_crs_t));
  if (!inc) {
    fprintf (stderr, "Cannot allocate space for cursor internal structure.\n");
    free (nc);
    *c = 0;
    return -1;
  }
  inc->hashp = (ffdb_htab_t *)dbp->internal;
  inc->pcursor = nc;
  memset (&inc->item, 0, sizeof(ffdb_hent_t));  
  inc->item.status = ITEM_CLEAN;
  FFDB_TAILQ_INSERT_TAIL(&(inc->hashp->curs_queue), inc, queue);

  /* Set internal pointer */
  nc->internal = inc;
  /* type of this cursor */
  nc->type = type;
  /* Set up close and get function for the cursor */
  nc->get = _ffdb_cursor_get;
  nc->close = _ffdb_cursor_close;

  *c = nc;
  FFDB_LOCK_INIT(inc->lock);

  return 0;
}


/**
 * Code related to configuration and user information
 * which are particular to QCD code
 */
int
ffdb_set_all_configs (FFDB_DB* db, ffdb_all_config_info_t* configs)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  /* Check whether configs have the same number of configs the database
   * thinks it should have
   */
  if (configs->numconfigs != hashp->hdr.num_cfigs) {
    fprintf (stderr, "Number of configurations mismatch: %d (dbase thinks) != %d (requested). \n", hashp->hdr.num_cfigs, configs->numconfigs);
    errno = EINVAL;
    return -1;
  }
  
  /* No need to lock the whole process because each page is locked */
  return ffdb_set_configs (hashp, configs);
}

int
ffdb_get_all_configs (const FFDB_DB* db, ffdb_all_config_info_t* configs)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  configs->numconfigs = hashp->hdr.num_cfigs;
  configs->allconfigs = (ffdb_config_info_t *)malloc(configs->numconfigs * sizeof(ffdb_config_info_t));
  if (!configs->allconfigs) {
    fprintf (stderr, "Cannot allocate memory space for all configurations.\n");
    errno = ENOMEM;
    return -1;
  }
  
  /* No need to lock the whole process because each page is locked */
  return ffdb_get_configs (hashp, configs);
}

unsigned int
ffdb_num_configs (const FFDB_DB* db)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  return hashp->hdr.num_cfigs;
}

int
ffdb_set_config (FFDB_DB* db, ffdb_config_info_t* config)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  if (config->index >= hashp->hdr.num_cfigs || config->index < 0) {
    fprintf (stderr, "ffdb_set_config provides a wrong configuration index number %d ( shoud be less than %d) \n",
	     config->index, hashp->hdr.cfig_npages);
    errno = EINVAL;
    return -1;
  }

  return ffdb_set_config_info (hashp, config);
}


int
ffdb_get_config (const FFDB_DB* db, unsigned int confignum, 
		 ffdb_config_info_t* config)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  if (config->index >= hashp->hdr.cfig_npages || config->index < 0) {
    fprintf (stderr, "ffdb_set_config provides a wrong configuration index number %d ( shoud be less than %d) \n",
	     config->index, hashp->hdr.cfig_npages);
    errno = EINVAL;
    return -1;
  }

  return ffdb_get_config_info (hashp, confignum, config);
}

/**
 * Set user information
 */
int
ffdb_set_user_info (FFDB_DB* db, unsigned char* data, unsigned int len)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  if (len > hashp->hdr.uinfolen) {
    fprintf (stderr, "%d exceeds maximum allowed user data is %d\n", 
	     len, hashp->hdr.uinfolen);
    errno = EINVAL;
    return -1;
  }
  
  return ffdb_set_uinfo (hashp, data, len);
}

/**
 * Set user information
 */
int
ffdb_get_user_info (const FFDB_DB* db, unsigned char data[], unsigned int* len)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  return ffdb_get_uinfo (hashp, data, len);
}


unsigned int 
ffdb_max_user_info_len (const FFDB_DB* db)
{
  ffdb_htab_t* hashp;

  hashp = (ffdb_htab_t *)db->internal;

  return hashp->hdr.uinfolen;
}


/************************************************************************
 * Cursor related routines                                              *
 ************************************************************************/
static int
_ffdb_cursor_get  (ffdb_cursor_t* cursor, FFDB_DBT* key,
		   FFDB_DBT* data, unsigned int flags)
{
  int status;
  ffdb_crs_t*  icrs = (ffdb_crs_t *)cursor->internal;
  ffdb_htab_t* hashp = icrs->hashp;
  unsigned int realflags;

  /* Check flag to see whether it is valid */
  if (flags != FFDB_FIRST && flags != FFDB_LAST &&
      flags != FFDB_NEXT && flags != FFDB_PREV) {
    fprintf (stderr, "Unsupported cursor get flag %d\n", flags);
    return -1;
  }
    
  realflags = flags;
  if (cursor->type == FFDB_KEY_CURSOR) {
    /* If the cursor has not been initialized */
    if (flags == FFDB_NEXT && icrs->item.status == ITEM_CLEAN)
      realflags = FFDB_FIRST;
    
    if (flags == FFDB_PREV && icrs->item.status == ITEM_CLEAN)
      realflags = FFDB_LAST;

    FFDB_LOCK(icrs->lock);
    status = ffdb_cursor_find_by_key (hashp, icrs, key, data, realflags);
    FFDB_UNLOCK(icrs->lock);
  }
  else {
    status = 0;
  }
  return status;
}


static int
_ffdb_cursor_close (ffdb_cursor_t* cursor)
{
  ffdb_crs_t* inc = (ffdb_crs_t *)cursor->internal;
  ffdb_htab_t* hashp = inc->hashp;

  FFDB_TAILQ_REMOVE(&(hashp->curs_queue), inc, queue);

  /* Check whether there is a page need to be released */
  if (inc->item.pagep) {
    ffdb_put_page (hashp, inc->item.pagep, TYPE(inc->item.pagep), 0);
    inc->item.pagep = 0;
  }

  FFDB_LOCK_FINI(inc->lock);
  /* Free all memory */
  free (inc);
  free (cursor);
  return 0;
}
