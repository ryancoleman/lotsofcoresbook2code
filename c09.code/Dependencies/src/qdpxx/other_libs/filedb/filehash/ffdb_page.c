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
 *     Hash database disk page layout definition
 *     All meta information are stored on disk in network byte order
 *
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *     $Log: ffdb_page.c,v $
 *     Revision 1.7  2009-10-08 18:19:54  edwards
 *     Fixed printf message to be "%u".
 *
 *     Revision 1.6  2009/07/18 21:02:52  chen
 *     minor test code change
 *
 *     Revision 1.5  2009/05/08 17:37:31  chen
 *     Fix a major bug (not clean out moved page inside page pool cache)
 *
 *     Revision 1.4  2009/04/21 18:51:20  chen
 *     Fix bugs related to number of pages upon moving pages in addition to clean pages on disk when the pages has been moved
 *
 *     Revision 1.3  2009/03/05 20:33:12  chen
 *     release bucket page too early fixed
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include "ffdb_db.h"
#include "ffdb_page.h"
#include "ffdb_hash.h"
#include "ffdb_hash_func.h"

/**
 * get next data page number either a new or reuse from a free page
 */
static pgno_t _ffdb_data_page (ffdb_htab_t* hashp, int new_page, int* reuse);
static pgno_t _ffdb_ovfl_page (ffdb_htab_t* hashp, int* reuse);

#ifdef _FFDB_STATISTICS
extern unsigned int hash_accesses, hash_collisions, hash_expansions, hash_overflows, hash_bigpages;
#endif

/**
 * Swap data pointer
 */
#define M_DATAP_SWAP(dp) do {					\
    M_32_SWAP((dp)->first);					\
    M_32_SWAP((dp)->offset);                                    \
    M_32_SWAP((dp)->len);                                       \
    M_32_SWAP((dp)->chksum);					\
  }while(0)

/**
 * Swap data header
 */
#define M_DATA_HEADER_SWAP(header)                              \
  do {                                                          \
    M_32_SWAP((header)->len);					\
    M_32_SWAP((header)->status);				\
    M_32_SWAP((header)->next);					\
    M_32_SWAP((header)->key_page);				\
    M_32_SWAP((header)->key_idx);				\
  }while(0)


/**
 * Swap configuration information structure
 */
#define M_CONFIG_INFO_SWAP(cinfo)                              \
  do {                                                         \
    M_32_SWAP((cinfo)->config);                                \
    M_32_SWAP((cinfo)->index);                                 \
    M_32_SWAP((cinfo)->inserted);                              \
    M_32_SWAP((cinfo)->type);                                  \
    M_32_SWAP((cinfo)->mtime);                                 \
  }while(0)


/**
 * Swap page header in depending on what type of page
 */
static void
_ffdb_swap_page_metainfo_in (void *p)
{
  unsigned int i, next, nelems;
  ffdb_data_header_t* header;

  M_32_SWAP(CURR_PGNO(p));
  M_32_SWAP(PREV_PGNO(p));
  M_32_SWAP(NEXT_PGNO(p));
  M_32_SWAP(PAGE_SIGN(p));
  M_16_SWAP(NUM_ENT(p));
  M_16_SWAP(TYPE(p));
  M_32_SWAP(CHKSUM(p));
  M_32_SWAP(OFFSET(p));

#ifdef _FFDB_DEBUG
  fprintf (stderr, "SWAP Header Read page %d \n", CURR_PGNO(p));
  fprintf (stderr, "Page Information: PGNO %d PREV PGNO %d NEXT PGNO %d PAGE SIGN 0x%x\n",
	   CURR_PGNO(p), PREV_PGNO(p),NEXT_PGNO(p), PAGE_SIGN(p));
  fprintf (stderr, "Page TYPE 0x%x OFFSET %d\n", TYPE(p), OFFSET(p));
#endif    
  
  switch (TYPE(p)) {
  case HASH_BUCKET_PAGE:
    for (i = 0; i < NUM_ENT(p); i++) {
      /* swap the offset */
      M_32_SWAP(KEY_OFF(p, i));
      M_32_SWAP(KEY_LEN(p, i));
      M_32_SWAP(DATAP_OFF(p, i));
      /* SWAP Data Pointer */
      M_DATAP_SWAP(DATAP(p, i));
    }
    break;
  case HASH_DATA_PAGE:
    nelems = NUM_ENT(p);
    M_32_SWAP(FIRST_DATA_POS(p));
    next = FIRST_DATA_POS(p);
    for (i = 0; i < nelems; i++) {
      /* There are multiple data on this page */
      header = BIG_DATA_HEADER(p, next);
      /* Swap header */
      M_DATA_HEADER_SWAP(header);
#ifdef _FFDB_DEBUG
      fprintf (stderr, "data header len = %d, status = 0x%x next = %d key_page = %d key_idx = %d\n",
	       header->len, header->status, header->next, header->key_page, header->key_idx);
#endif
      next = header->next;
    }
    break;
  case HASH_FREE_PAGE:
    for (i = 0; i < NUM_ENT(p); i++) 
      /* swap the offset */
      M_32_SWAP(FREE_PAGE(p, i));
    break;
  case HASH_UINFO_PAGE:
    M_32_SWAP(USER_DATA_LEN(p));
    break;
  case HASH_CONFIG_PAGE:
    for (i = 0; i < NUM_ENT(p); i++) 
      M_CONFIG_INFO_SWAP(CONFIG_INFO(p, i));
    break;
  default:
    break;
  }
}

/**
 * Calculate check sum for a page
 */
static unsigned int
_ffdb_page_checksum (ffdb_htab_t* hashp, void* pagep)
{
  unsigned int chksum = 0;
  
  /* check sum everything before checksum entry */
  chksum = __ffdb_crc32_checksum (chksum, pagep, CHKSUM_POS);
  /* check sum everything after the checksum entry */
  chksum = __ffdb_crc32_checksum (chksum, 
				  (unsigned char *)pagep + CHKSUM_OFFSET,
				  hashp->hdr.bsize - CHKSUM_OFFSET);

  return chksum;
}

/**
 * Initialize a new page
 */
static int
_ffdb_init_page (ffdb_htab_t* hashp, void* p, pgno_t page, 
		 unsigned int addrtype)
{
  /* Page number */
  CURR_PGNO(p) = page;
  PREV_PGNO(p) = NEXT_PGNO(p) = INVALID_PGNO;
  PAGE_SIGN(p) = FFDB_PAGE_MAGIC;
  /* Number of items */
  NUM_ENT(p) = 0;
  /* Page Type */
  TYPE(p) = addrtype;
  CHKSUM(p) = 0;
  switch(addrtype) {
  case HASH_BUCKET_PAGE:
  case HASH_OVFL_PAGE:
    OFFSET(p) = hashp->hdr.bsize - 1;
    /* I will treat overflow pages as bucket page */
    TYPE(p) = HASH_BUCKET_PAGE;
    break;
  case HASH_BIG_PAGE:
    OFFSET(p) = BIG_PAGE_OVERHEAD;
    FIRST_DATA_POS(p) = BIG_PAGE_OVERHEAD;
    break;
  case HASH_UINFO_PAGE:
    USER_DATA_LEN(p) = 0;
    break;
  case HASH_CONFIG_PAGE:
    OFFSET(p) = BIG_PAGE_OVERHEAD;
    break;
  case HASH_FREE_PAGE:
    OFFSET(p) = 0;
    /* Not using offset: actually just pad */
    break;
  default:
    break;
  }
  return 0;
}    

/**
 * This is the routine called right after a page is read
 */
void
ffdb_pgin_routine (void* arg, pgno_t pgno, void* page)
{
  ffdb_htab_t* hashp;
  unsigned int chksum = 0;

  hashp = (ffdb_htab_t *)arg;

  /* calculate checksum before byte swapped */
  chksum = _ffdb_page_checksum (hashp, page);

#if 0
  /* First swap the header to disk */
  if (hashp->mborder == LITTLE_ENDIAN) 
    _ffdb_swap_page_metainfo_in (page);
#endif

#if 1
  if (hashp->hdr.lorder != hashp->mborder)
    _ffdb_swap_page_metainfo_in (page);  
#endif

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Read page %d in at memory %p\n", pgno, page);
  fprintf (stderr, "Page Information: PGNO %d PREV PGNO %d NEXT PGNO %d PAGE SIGN 0x%x\n",
	   CURR_PGNO(page), PREV_PGNO(page),NEXT_PGNO(page), PAGE_SIGN(page));
  fprintf (stderr, "Page TYPE 0x%x OFFSET %d\n", TYPE(page), OFFSET(page));
#endif  
  
  if (PAGE_SIGN(page) != FFDB_PAGE_MAGIC) {
    /* This is a unintialized page */
    _ffdb_init_page (hashp, page, pgno, HASH_DELETED_PAGE);
  }
  else {
    if (chksum != CHKSUM(page)) {
      /* The check sum chould 0 because this page is empty page */
      fprintf (stderr, "Reading page %d checksum mismatch 0x%x (new) != 0x%x (on disk)\n", pgno, chksum, CHKSUM(page));
      exit (123);
    }
  }
}

/**
 * Swap page meta information out before this page is written to disk
 */
static void
_ffdb_swap_page_metainfo_out (void *p)
{
  unsigned int i, next;
  unsigned int type = TYPE(p);
  unsigned int num = NUM_ENT(p);
  ffdb_data_header_t* header;

#ifdef _FFDB_DEBUG
  fprintf (stderr, "SWAP Header Write page %d \n", CURR_PGNO(p));
  fprintf (stderr, "Page Information: PGNO %d PREV PGNO %d NEXT PGNO %d PAGE SIGN 0x%x\n",
	   CURR_PGNO(p), PREV_PGNO(p),NEXT_PGNO(p), PAGE_SIGN(p));
  fprintf (stderr, "Page TYPE 0x%x OFFSET %d\n", TYPE(p), OFFSET(p));
#endif     

  M_32_SWAP(CURR_PGNO(p));
  M_32_SWAP(PREV_PGNO(p));
  M_32_SWAP(NEXT_PGNO(p));
  M_32_SWAP(PAGE_SIGN(p));
  M_16_SWAP(NUM_ENT(p));
  M_16_SWAP(TYPE(p));
  M_32_SWAP(CHKSUM(p));
  M_32_SWAP(OFFSET(p));

  switch (type) {
  case HASH_BUCKET_PAGE:
    for (i = 0; i < num; i++) {
      M_DATAP_SWAP(DATAP(p,i));
      /* Swap other information */
      M_32_SWAP(KEY_OFF(p, i));
      M_32_SWAP(KEY_LEN(p, i));
      M_32_SWAP(DATAP_OFF(p, i));
    }
    break;
  case HASH_DATA_PAGE:
    /* get next data item header position */
    next = FIRST_DATA_POS(p);
    M_32_SWAP(FIRST_DATA_POS(p));
    for (i = 0; i < num; i++) {
      /* get data header and next item before swapping */
      header = BIG_DATA_HEADER(p, next);

#ifdef _FFDB_DEBUG
      fprintf (stderr, "data header len = %d, status = 0x%x next = %d key_page = %d key_idx = %d\n",
	       header->len, header->status, header->next, header->key_page, header->key_idx);
#endif

      next = header->next;
      M_DATA_HEADER_SWAP(header);
    }
    break;
  case HASH_FREE_PAGE:
    for (i = 0; i < num; i++) 
      /* swap the offset */
      M_32_SWAP(FREE_PAGE(p, i));
    break;
  case HASH_UINFO_PAGE:
    M_32_SWAP(USER_DATA_LEN(p));
    break;
  case HASH_CONFIG_PAGE:
    for (i = 0; i < num; i++) 
      M_CONFIG_INFO_SWAP(CONFIG_INFO(p, i));
    break;
  default:
    break;
  }
}

void
ffdb_pgout_routine (void* arg, pgno_t pgno, void* page)
{
  ffdb_htab_t* hashp;
  unsigned int chksumval;

  hashp = (ffdb_htab_t *)arg;

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Write page %d memory %p\n", pgno, page);
  fprintf (stderr, "Page Information: PGNO %d PREV PGNO %d NEXT PGNO %d NUM_ENT %d\n",
	   CURR_PGNO(page), PREV_PGNO(page),NEXT_PGNO(page), NUM_ENT(page));
  fprintf (stderr, "Page TYPE 0x%x OFFSET %d PAGE_SIGN 0x%x \n", TYPE(page), OFFSET(page), PAGE_SIGN(page));
#endif

  /* First swap the header to disk */
#if 0
  if (hashp->mborder == LITTLE_ENDIAN) 
    _ffdb_swap_page_metainfo_out (page);
#endif

#if 1
  if (hashp->hdr.lorder != hashp->mborder)
    _ffdb_swap_page_metainfo_out (page);
#endif

  /* Need to calculte check sum here for different type of page */
  /* calculate checksum after byte swapped */
  chksumval = _ffdb_page_checksum (hashp, page);
  if (hashp->hdr.lorder != hashp->mborder)
    M_32_SWAP(chksumval);
  CHKSUM(page) = chksumval;
  
}


/**
 * Get a real data item from data pages pointed by the data pointer
 * 
 * @param hashp the pointer to hash table
 * @param item  hash entry item pointer
 * @param val   the pointer to data storage: if val->data is null, 
 * a new malloc is used, caller has to free. If val->data is not null,
 * a quick memcpy using val->size as the length.
 * @param datap data pointer containg information about data
 *
 * @return 0 on success, -1 otherwise
 */
static int
_ffdb_get_data (ffdb_htab_t* hashp, ffdb_hent_t* item,
		FFDB_DBT* val, ffdb_datap_t* datap)
{
  pgno_t next, tp;
  void* pagep;
  ffdb_data_header_t* header;
  unsigned int start, rlen, idx, copylen, newchksum;
  int needfree = 0;

  /* Get first page where the data item resides */
  pagep = ffdb_get_page (hashp, datap->first, HASH_DATA_PAGE, 0, &tp);
  if (!pagep) {
    fprintf (stderr, "Cannot get data page at %d \n", datap->first);
    return -1;
  }
  
  /* now do a quick sanity check */
  assert (datap->first == CURR_PGNO(pagep));
  assert (NUM_ENT(pagep) >= 1);

  /* Jump to the offset pointed by datap */
  header = BIG_DATA_HEADER(pagep,datap->offset);
 
  assert (header->len == datap->len);
  assert (header->status == DATA_VALID);
  assert (header->key_page == item->pgno);
  assert (header->key_idx == item->pgndx);

  /* Now check whether I have allocated space to data */
  if (val->data && val->size > 0) {
    if (val->size < header->len) {
      fprintf (stderr, "Warning: application allocated space %d < data size %d\n",
	       val->size, header->len);
      return -1;
    }
    else
      val->size = header->len;
  }
  else {
    val->data = (char *)malloc(header->len * sizeof (char));
    val->size = header->len;
    needfree = 1;
  }
  
  /* Now I am ready to copy */
  rlen = val->size;
  start = datap->offset + sizeof(ffdb_data_header_t);
  next = NEXT_PGNO(pagep);
  while (rlen > 0) {
    /* where copy starts in val->data */
    idx = val->size - rlen;
    /* find out how many bytes of data to copy */
    if (start + rlen <= hashp->hdr.bsize) 
      copylen = rlen;
    else 
      copylen = hashp->hdr.bsize - start;

    memcpy ((unsigned char *)val->data + idx, pagep + start, copylen);
    /* release this page */
    ffdb_put_page (hashp, pagep, HASH_DATA_PAGE, 0);

    rlen -= copylen;
    if (rlen > 0) { /* multiple pages */
      /* get next page */
      pagep = ffdb_get_page (hashp, next, HASH_DATA_PAGE, 0, &tp);
      if (!pagep) {
	fprintf (stderr, "Cannot get data page at %d\n", next);
	val->size = 0;
	if (needfree) {
	  free (val->data);
	  val->data = 0;
	}
	return -1;
      }
      start = BIG_PAGE_OVERHEAD;
      next = NEXT_PGNO(pagep);
    }
  }
  /* Now item is copied, run check sum */
  newchksum = 0;
  newchksum = __ffdb_crc32_checksum (newchksum, val->data,
				     val->size);

  if (val->size == header->len) {
    if (newchksum != datap->chksum) {
      fprintf (stderr, "Val = %s size = %d\n", (char *)val->data, val->size);
      fprintf (stderr, "Get data checksum mismatch 0x%x (calculated) != 0x%x (stored)\n", newchksum, datap->chksum);
      if (needfree) {
	free (val->data);
	val->data = 0;
      }
      val->size = 0;
      return -1;
    }
  }
  return 0;
}

/**
 * Add a datum into a data page
 *
 * @param hashp the pointer to hash table
 * @param key_page what page key is stored on
 * @param key_index index of this key on the key page
 * @param val   the pointer to the datum
 * @param pagep memory pointer of the beginning of the data page
 * @param dpage data page number
 * @datap datap regular hash data pointer residing on hash page
 *
 * @return returns 0 on success, -1 system failure
 */
static int
_ffdb_add_data (ffdb_htab_t *hashp, pgno_t key_page, pgno_t key_index,
		const FFDB_DBT* val, void* mem, pgno_t pnum,
		ffdb_datap_t* datap)
{
  unsigned int start, fspace, npages, idx, copylen;
  ffdb_data_header_t header;
  pgno_t tp, cpage, currp, prevp, fp;
  void *cpagep, *currpagep;
  int reuse, rlen;

  /* Go to start of the free byte which is 4 byte aligned 
   * The Highest_Free result can be 0 because the previous write 
   * only leaves either no space or not enough space for a
   * data header
   *
   * -------------Important--------------
   * When data is put on a page, unless there are space at least for
   * a data header, this page is marked full
   * 
   */
  reuse = 0;
  start = HIGHEST_FREE(mem);
  if (start == 0) {
    /* I need to allocate another page */
    reuse = 0;
    cpage = _ffdb_data_page (hashp, 1, &reuse);
    cpagep = ffdb_get_page (hashp, cpage, HASH_DATA_PAGE, FFDB_CREATE, &tp);
    if (!cpagep) {
      fprintf (stderr, "Cannot allocate next page at %d\n", cpage);
      ffdb_put_page (hashp, mem, HASH_DATA_PAGE, 0);
      return -1;
    }

    /* this has to be success */
    assert (cpagep != 0);

    /* this page could be reuse page from overflow pages */
    if (reuse) 
      _ffdb_init_page (hashp, cpagep, cpage, HASH_DATA_PAGE);      

    NEXT_PGNO(mem) = cpage;
    PREV_PGNO(cpagep) = pnum;

    /* Now free the old data page */
    ffdb_put_page (hashp, mem, HASH_DATA_PAGE, 0);

    /* update start value */
    start = HIGHEST_FREE(cpagep);
  }
  else {
    cpage = pnum;
    cpagep = mem;
  }

  /* data header stays in the same page */
  /* find out how much space left on this page */
  fspace = hashp->hdr.bsize - start;

  /* The following has to be true */
  assert (fspace >= BIG_DATA_OVERHEAD);

  header.len = val->size;
  header.status = DATA_VALID;
  header.key_page = key_page;
  header.key_idx = key_index;
  header.next = 0;

  if (BIG_DATA_TOTAL_SIZE(val) <= fspace) {
    header.next = start + BIG_DATA_TOTAL_SIZE(val);
    /* align the address */
    ALIGN_ADDR(header.next);

    /* If there is no space for header, jump to next page */
    if (header.next + BIG_DATA_OVERHEAD >= hashp->hdr.bsize)
      header.next = 0;

    /* now copy header on to memory */
    memcpy (cpagep + start, &header, BIG_DATA_OVERHEAD);
    /* Now copy data on to memory */
    memmove (cpagep + start + BIG_DATA_OVERHEAD, val->data, val->size);

    
    /* Update the page header */
    HIGHEST_FREE(cpagep) = header.next;
    NUM_ENT(cpagep)++;

    /* Put this page back */
    ffdb_put_page (hashp, cpagep, HASH_DATA_PAGE, 1);
  }
  else {
    /* this data expands into multiple pages */
    header.next = 0;

    /* how many bytes are left for other pages */
    rlen = val->size;

    if (fspace - BIG_DATA_OVERHEAD > 0) {
    /* copy part of data onto this page */
      memmove (cpagep + start + BIG_DATA_OVERHEAD, val->data, 
	       fspace - BIG_DATA_OVERHEAD);
      rlen = val->size - (fspace - BIG_DATA_OVERHEAD); 
    }
    /* Now copy data to each page */
    npages = 0;
    
    /* get a free page or a new page */
    /* fp is the first page of the chain */
    reuse = 0;
    fp = currp = _ffdb_data_page (hashp, 1, &reuse);
    prevp = cpage;
    
    while (rlen > 0) {
      currpagep = ffdb_get_page (hashp, currp, HASH_DATA_PAGE,FFDB_CREATE, &tp);
      if (!currpagep) {
	fprintf (stderr, "Cannot allocate data page at %d\n", currp);
	return -1;
      }
      if (reuse)
	_ffdb_init_page (hashp, currpagep, currp, HASH_DATA_PAGE);
      
      PREV_PGNO(currpagep) = prevp;

      /* We do not need to increase number of items, since this is part
       * of data on this page, header is on another page 
       */
      /* check whether this page is enough for this datum */
      if (rlen <= hashp->hdr.bsize - BIG_PAGE_OVERHEAD - BIG_DATA_OVERHEAD) {
	/* Free Memory and First Data Position must be aligned */
	HIGHEST_FREE(currpagep) = BIG_PAGE_OVERHEAD + rlen;
	ALIGN_ADDR(HIGHEST_FREE(currpagep));
	FIRST_DATA_POS(currpagep) = (BIG_PAGE_OVERHEAD + rlen);
	ALIGN_ADDR(FIRST_DATA_POS(currpagep));
	
	/* we copy remaining part of data */
	copylen = rlen;
      }
      else {
	/* there is no space for another data */
	HIGHEST_FREE(currpagep) = 0;
	FIRST_DATA_POS(currpagep) = 0; 
	/* we copy most of the rest of page */
	copylen = hashp->hdr.bsize - BIG_PAGE_OVERHEAD;
      }

      /* where to start copy the data */
      idx = val->size - rlen;
      memmove (currpagep + BIG_PAGE_OVERHEAD, ((unsigned char *)val->data + idx),
	       copylen);
      
      /* reduce number of bytes */
      rlen -= copylen;

      if (rlen > 0) {
	/* remember the previous page number */
	prevp = currp;

	/* get next page number */
	reuse = 0;
	currp = _ffdb_data_page (hashp, 1, &reuse);

	/* update next page number */
	NEXT_PGNO(currpagep) = currp;
      }

      /* put this page out */
      ffdb_put_page (hashp, currpagep, HASH_DATA_PAGE, 1);

      /* update number of pages */
      npages++;
    }

    
    /* Now data copoied, I need to update header 
     * information on the first page
     */

    /* now copy header on to memory */
    memcpy (cpagep + start, &header, BIG_DATA_OVERHEAD);

    /* Update the page header */
    NEXT_PGNO(cpagep) = fp;
    HIGHEST_FREE(cpagep) = header.next;
    NUM_ENT(cpagep)++;

    /* Put this page back */
    ffdb_put_page (hashp, cpagep, HASH_DATA_PAGE, 1); 

  }
  
  /* Update data header */
  datap->first = cpage;
  datap->offset = start;
  
  return 0;
}

/**
 * Replace a data inside database data page
 *
 * @param hashp the usual pointer to the hash table
 * @param val the new value to replace old value
 * @param mem the data page memory pointer (this is the first page)
 * @param pnum the data page page number
 * @param datap data pointer from ket page with all updated information
 *
 * @return returns 0 on success, -1 on failure
 */
static int
_ffdb_replace_data (ffdb_htab_t* hashp, const FFDB_DBT* val,
		    void* mem, pgno_t pnum,
		    ffdb_datap_t* datap)
{
  ffdb_data_header_t* header;
  unsigned int fspace, copylen, idx;
  pgno_t tp, currp;
  void *currpagep;
  int rlen;

  /* Get old data header. The new data header is at the same place */
  header = BIG_DATA_HEADER(mem, datap->offset);

  /* data header length need to be changed and header next stays the same */
  header->len = datap->len;

  /* copy data on to data pages */
  fspace = hashp->hdr.bsize - datap->offset;
  /* The following has to be true */
  assert (fspace >= BIG_DATA_OVERHEAD);

  /* Data can fit in the page */
  if (BIG_DATA_TOTAL_SIZE(val) <= fspace) {
    /* Copy data onto the page, header is changed already in memory */
    memmove (mem + datap->offset + BIG_DATA_OVERHEAD, val->data,
	     val->size);

    /* Put data page away */
    ffdb_put_page (hashp, mem, TYPE(mem), 1);
  }
  else { /* This data expands multiple pages */
    /* how many bytes are left for othe pages */
    rlen = val->size;

    /* Copy part of data to this page */
    if (fspace - BIG_DATA_OVERHEAD > 0) {
    /* copy part of data onto this page */
      memmove (mem + datap->offset + BIG_DATA_OVERHEAD, val->data, 
	       fspace - BIG_DATA_OVERHEAD);
      rlen = val->size - (fspace - BIG_DATA_OVERHEAD); 
    }
    
    /* Now copy remaining data to multiple pages */
    currp = NEXT_PGNO(mem);

    while (rlen > 0) {
      /* load this page */
      currpagep = ffdb_get_page (hashp, currp, HASH_DATA_PAGE, 0, &tp);
      if (!currpagep) {
	fprintf (stderr, "Cannot allocate page %d to overwrite data\n",
		 currp);
	return -1;
      }

      /* check whether data will fit this page */
      if (rlen <= hashp->hdr.bsize - BIG_PAGE_OVERHEAD - BIG_DATA_OVERHEAD) 
	copylen = rlen;
      else 
	copylen = hashp->hdr.bsize - BIG_PAGE_OVERHEAD;

      /* where to start copy the data */
      idx = val->size - rlen;
      memmove (currpagep + BIG_PAGE_OVERHEAD, 
	       ((unsigned char *)val->data + idx), copylen);
      
      /* reduce number of remaining bytes */
      rlen -= copylen;

      if (rlen > 0) 
	currp = NEXT_PGNO (currpagep);

      /* release current page */
      ffdb_put_page (hashp, currpagep, TYPE(currpagep), 1);
    }

    /* Now put the first page out */
    ffdb_put_page (hashp, mem, TYPE(mem), 1);
  }
  return 0;
}


/* When a page split, the hash items on the page will be split into
 * the old page and a new page. The data items pointed by the hash items
 * have to be updated to reflect data header information about key page
 * and key index inside the key page
 */
static int
_ffdb_update_data_info (ffdb_htab_t* hashp, ffdb_datap_t* datap,
			pgno_t kpage, unsigned int kpgindx)
{
  void* dpagep;
  pgno_t dpage;
  ffdb_data_header_t* data_header;

  /* Get data page given from datap */
  dpagep = ffdb_get_page (hashp, datap->first, HASH_DATA_PAGE, 0, &dpage);
  if (!dpagep) {
    fprintf (stderr, "Cannot retrieve data page at %u\n", datap->first);
    return -1;
  }

  /* Access this data item and update header information */
  data_header = BIG_DATA_HEADER (dpagep, datap->offset);
  data_header->key_page = kpage;
  data_header->key_idx = kpgindx;

  /* Commit change */
  ffdb_put_page (hashp, dpagep, HASH_DATA_PAGE, 1);

  return 0;
}

/**
 * Get a free page from free map page if there is one
 * Return 0 when there is no free page
 */
static pgno_t
_ffdb_reuse_free_ovflpage (ffdb_htab_t* hashp)
{
  pgno_t fpage, tp, num, np, nextp;
  void *fpagep;
  unsigned int clevel = hashp->hdr.ovfl_point;

  /* check current free page number at this level */
  num = 0;
  if ((fpage = hashp->hdr.free_pages[clevel]) != INVALID_PGNO) {
    fpagep = ffdb_get_page (hashp, fpage, HASH_FREE_PAGE, 0, &tp);
    if (!fpagep) {
      fprintf (stderr, "Fatal error: free map page at %d could not be found\n",
	       fpage);
      abort ();
    }

    if ((np = NUM_FREE_PAGES(fpagep)) > 0) {
      num = FREE_PAGE(fpagep, np - 1);
      /* reduce number of free pages by one */
      NUM_FREE_PAGES(fpagep)--;
      ffdb_put_page (hashp, fpagep, HASH_FREE_PAGE, 1);
    }
    else {
      nextp = NEXT_PGNO(fpagep);

      /* There is no page free anymore, so I can free this 
       * page for a normal use
       *
       * No need to clean this page, we will reset this page upon reloading
       */
      ffdb_put_page (hashp, fpagep, HASH_OVFL_PAGE, 1);
#ifdef _FFDB_DEBUG
      fprintf (stderr, "Free page %d is not longer needed, put it back\n",
	       fpage);
#endif
      num = fpage;
      /* nextp could be INVALID_PGNO */
      hashp->hdr.free_pages[clevel] = nextp;
#ifdef _FFDB_DEBUG
      fprintf (stderr, "Free map page %d contains no free pages, reuse it\n",
	       fpage);
#endif
    }
  }
  return num;
}



/**
 * Find out what is next data page number given current page number
 * We need first to check freed overflow pages
 *
 * If there are somthing really wrong, the page released by this call
 * cannot be reclaimed. (We will live with the consequence)
 */
static pgno_t
_ffdb_data_page (ffdb_htab_t* hashp, int new_page, int* reuse)
{
  pgno_t num = 0;
  pgno_t maxp = 0;
  /* get next level of overflow point */
  unsigned int level = hashp->hdr.ovfl_point + 1;

  *reuse = 0;
  if (hashp->curr_dpage == INVALID_PGNO) {
    /* The data page and overflow page starts at the following page number */
    BUCKET_TO_PAGE(hashp->hdr.high_mask, hashp->curr_dpage);
    hashp->curr_dpage++;

    BUCKET_TO_PAGE(hashp->hdr.max_bucket, maxp);
    fprintf (stderr, "Doubling at level %d max_bucket at %d data page starts at %d\n", hashp->hdr.ovfl_point, maxp, hashp->curr_dpage);

#ifdef _FFDB_DEBUG
    fprintf (stderr, "data page starts at %d for level %d\n",
	     hashp->curr_dpage, hashp->hdr.ovfl_point);
#endif
  }

  if (hashp->hdr.spares[level] == 0) 
    hashp->hdr.spares[level] = hashp->curr_dpage + 1;

  if (!new_page) 
    num = hashp->curr_dpage;
  else {
    num = _ffdb_reuse_free_ovflpage (hashp);
    if (num > 0) {
#ifdef _FFDB_DEBUG
      fprintf (stderr, "Reuse previously freed overflow page %d\n", num);
#endif
      *reuse = 1;
    }
    else {
      num = hashp->hdr.spares[level];
      hashp->hdr.spares[level]++;
    }
    hashp->curr_dpage = num;
  }
  return num;
}

/**
 * Find out what is next overflow page number given current page number
 * We need first to check freed overflow pages
 *
 * If there are somthing really wrong, the page released by this call
 * cannot be reclaimed. (We will live with the consequence)
 */
static pgno_t
_ffdb_ovfl_page (ffdb_htab_t* hashp, int* reuse)
{
  pgno_t num = 0;
  pgno_t maxp = 0;
  /* get next level of overflow point */
  unsigned int level = hashp->hdr.ovfl_point + 1;

  *reuse = 0;
  if (hashp->curr_dpage == INVALID_PGNO) {
    /* The data page and overflow page starts at the following page number */
    BUCKET_TO_PAGE(hashp->hdr.high_mask, hashp->curr_dpage);
    hashp->curr_dpage++;

    BUCKET_TO_PAGE(hashp->hdr.max_bucket, maxp);
    fprintf (stderr, "Doubling at level %d max_bucket at %d over flow page starts at %d\n", hashp->hdr.ovfl_point, maxp, hashp->curr_dpage);
#ifdef _FFDB_DEBUG
    fprintf (stderr, "overflow page starts at %d for level %d\n",
	     hashp->curr_dpage, hashp->hdr.ovfl_point);
#endif
  }

  if (hashp->hdr.spares[level] == 0) 
    hashp->hdr.spares[level] = hashp->curr_dpage + 1;

  num = _ffdb_reuse_free_ovflpage (hashp);
  if (num > 0) {
#ifdef _FFDB_DEBUG
    fprintf (stderr, "Overflow reuse previously freed overflow page %d\n", num);
#endif
    *reuse = 1;
  }
  else {
    num = hashp->hdr.spares[level];
    hashp->hdr.spares[level]++;
  }
  return num;
}


/**
 * Search for the last data page on this level starting from the last
 * going backward to search
 *
 * This routine has to success since there is at least one data page
 * at any level
 */
pgno_t
ffdb_last_data_page (ffdb_htab_t* hashp, pgno_t start)
{
  pgno_t page, ret, tp;
  void* pagep;
  int  done = 0;

  /* We start from the first page above reguler hash pages */
  ret = INVALID_PGNO;
  page = start;
  while (!done) {
    pagep = ffdb_get_page (hashp, page, HASH_RAW_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get page %d while looking for the last data page.\n",
	       page);
      abort ();
    }
    if (TYPE(pagep) == HASH_DATA_PAGE && NEXT_PGNO(pagep) == INVALID_PGNO) {
      done = 1;
      ret = page;
    }
    ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
    page--;
  }

  return ret;
}

/**
 * Put a page back into pagepool
 */
int
ffdb_put_page (ffdb_htab_t* hashp, void* mem,
	       unsigned int addrtype, int dirty)
{
  unsigned int flags = 0;

  if (dirty) {
    flags |= FFDB_PAGE_DIRTY;
  }
  return ffdb_pagepool_put_page (hashp->mp, mem, flags);
}

/**
 * get a brand new page
 */
int
ffdb_new_page (ffdb_htab_t* hashp, pgno_t addr, 
	       unsigned int addrtype)
{
  pgno_t page;
  void* mem;
  int   status;

  switch (addrtype) {
  case HASH_BUCKET_PAGE:
    BUCKET_TO_PAGE(addr,page);
#ifdef _FFDB_DEBUG
    fprintf (stderr, "New: Bucket %d mapped to page %d\n", addr, page);
#endif
    break;
  default:
    page = addr;
    break;
  }
  status = ffdb_pagepool_new_page (hashp->mp, &page,
				   FFDB_PAGE_REQUEST, &mem);
  if (status != 0)
    return status;
  
  /* We have to initialize the page */
  _ffdb_init_page (hashp, mem, page, addrtype);

  /* Flush back the page so a get page will have the correct content */
  ffdb_put_page (hashp, mem, addrtype, 1);

  return 0;
}

/**
 * Add a free overflow page the the free information page
 *
 * free page is opage
 * fpagep is page pointer for the first page of the free information page
 */
static void
_ffdb_update_free_page (ffdb_htab_t* hashp, void* fpagep, pgno_t opage)
{
  unsigned int nf = NUM_FREE_PAGES(fpagep);

  FREE_PAGE(fpagep, nf) = opage;
  
  NUM_FREE_PAGES(fpagep)++;
}

/**
 * Free an overflow page: this happens when a bucket is split. 
 * We need to update free page in the previous level
 *
 * Return whether this page should be deleted
 *
 * The freed overflow page is represented by 'memp'
 */
static void
_ffdb_free_ovflpage (ffdb_htab_t* hashp, void* memp, int isdoubling,
		     int* deleteit)
{
  unsigned int level, nf;
  pgno_t fpage, tp, opage, nextpage, xpage;
  void *fpagep, *fnext, *xpagep;
  int xtra_page, reuse;

  /* Set delete flag to true */
  *deleteit = 1;

  /* get overflow page number */
  opage = CURR_PGNO(memp);

  /* check what level this page belongs to */
  for (level = 0; level <= hashp->hdr.ovfl_point; level++) {
    if (opage < hashp->hdr.spares[level]) {
      break;
    }
  }
  /* level is one level down */
  level--;

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Free overflow page %d for level %d isdoubling = %d\n",
	   opage, level, isdoubling);
#endif

  reuse = 0;
  if (hashp->hdr.free_pages[level] == INVALID_PGNO) {
    /* I just use this page to keep track free pages */
    fpage = opage;
    hashp->hdr.free_pages[level] = fpage;
    *deleteit = 0;

    fnext = fpagep = memp;

    /* Clean this page out into a free page */
    _ffdb_init_page (hashp, fpagep, fpage, HASH_FREE_PAGE);
#ifdef _FFDB_DEBUG
    fprintf (stderr, "Reuse %d to be free page \n", fpage);
#endif
  }
  else {
    fpage = hashp->hdr.free_pages[level];

    /* get this page */
    fnext = fpagep = ffdb_get_page (hashp, fpage, HASH_FREE_PAGE, 
				    FFDB_CREATE, &tp);
    if (!fpagep) {
      fprintf (stderr, "cannot allocate %d for free overflow pages\n", fpage);
      abort ();
    }
    if (reuse)
      _ffdb_init_page (hashp, fpagep, fpage, HASH_FREE_PAGE);

#ifdef _FFDB_DEBUG
    fprintf (stderr, "Get a page %d keeping free page at level %d\n",
	     tp, level);
#endif

    xtra_page = 1;
    /* I have to check chain of pages */
    while (fnext) {
      /* get current number of free pages contained inside this page */
      nf = NUM_FREE_PAGES(fpagep);
      if (nf < MAX_NUM_FREE_PAGES - 1){ /* There is space on this page */
	/* Now write this page number to the free page */
	_ffdb_update_free_page (hashp, fpagep, opage);
	xtra_page = 0;
	/* put this page out */
	ffdb_put_page (hashp, fpagep, TYPE(fpagep), 1);
	/* We are done */
	fnext = 0;
	xtra_page = 0;
      }
      else {
	nextpage = NEXT_PGNO(fpagep);
	if (nextpage == INVALID_PGNO)
	  fnext = 0;
	else {
	  fnext = ffdb_get_page (hashp, nextpage, HASH_FREE_PAGE, 0, &tp);
	  if (!fnext) {
	    fprintf (stderr, "Cannot get free bitmap page %d\n",
		     nextpage);
	    /* release the previous page */
	    ffdb_put_page (hashp, fpagep, TYPE(fpagep), 0);
	    return;
	  }
	  /* release the previous page */
	  ffdb_put_page (hashp, fpagep, TYPE(fpagep), 0);

	  /* assign new pointers */
	  fpagep = fnext;
	  fpage = nextpage;
	}
      }
    }
  
    /* If we need extra page */
    reuse = 0;
    if (xtra_page) {
      /* We need an extra page to store free page. We just use this page */
      xpage = opage;
#ifdef _FFDB_DEBUG
      fprintf (stderr, "reuse overflow page %d as chained free page %d\n",
	       xpage, fpage);
#endif
      *deleteit = 0;

      /* initialize this page */
      xpagep = memp;
      _ffdb_init_page (hashp, xpagep, xpage, HASH_FREE_PAGE);

#ifdef _FFDB_STATISTICS
      hash_overflows++;
#endif

      /* set up correct link */
      NEXT_PGNO (fpagep) = xpage;
      PREV_PGNO (xpagep) = fpage;

      /* release the previous page */
      ffdb_put_page (hashp, fpagep, TYPE(fpagep), 0); 
    }
  }
}

/**
 * Delete a page: usually an overflow page
 */
static int
ffdb_delete_page (ffdb_htab_t* hashp, void* pagep, 
		  unsigned int addrtype, int isdoubling)
{
  int deleteit = 0;
  if (addrtype == HASH_OVFL_PAGE) {
    _ffdb_free_ovflpage (hashp, pagep, isdoubling, &deleteit);
    if (deleteit) {
      /* if this page need to be deleted */
      _ffdb_init_page (hashp, pagep, CURR_PGNO(pagep), HASH_DELETED_PAGE);
      return ffdb_pagepool_delete (hashp->mp, pagep); 
    }
    else 
      return ffdb_put_page (hashp, pagep, TYPE(pagep), 1);
  }
  else 
    return ffdb_pagepool_delete (hashp->mp, pagep);
}

/**
 * Get a page
 */
void *
ffdb_get_page (ffdb_htab_t* hashp, pgno_t addr, 
	       unsigned int addrtype, unsigned flags, pgno_t* page)
{
  int   status;
  void* mem = 0;

  switch (addrtype) {
  case HASH_BUCKET_PAGE:
    BUCKET_TO_PAGE(addr, *page);
#ifdef _FFDB_DEBUG
    fprintf (stderr, "Get: Bucket %d mapped to page %d\n", addr, *page);
#endif
    break;
  default:
    *page = addr;
    break;
  }
  status = ffdb_pagepool_get_page (hashp->mp, page, flags, &mem);

  if (status != 0)
    return 0;

  /* If this page is a new page, we need to initialize it.
   * A new page is a page with page number 0 since zero page number
   * will never occur
   */
  if (FFDB_FLAG_ISSET(flags, FFDB_CREATE) && CURR_PGNO(mem) == 0) 
    _ffdb_init_page (hashp, mem, *page, addrtype);

  return mem;
}

/**
 * Release item and the related page associated with this item
 *
 * This routine is used when the item is no longer needed, which
 * usually happened when hash_get did not get anything
 */
int ffdb_release_item (ffdb_htab_t* hashp,
		       ffdb_hent_t* item)
{
  int status;
  /* Now I do not need this page */
  status = ffdb_put_page (hashp, item->pagep, HASH_RAW_PAGE, 0);

  item->pagep = 0;
  item->pgno = INVALID_PGNO;
  item->bucket = INVALID_PGNO;
  item->status = ITEM_ERROR;

  return status;
}


/**
 * Routines to find a page for key and data pair
 */
int ffdb_find_item (ffdb_htab_t* hashp,
		    FFDB_DBT* key, FFDB_DBT* val,
		    ffdb_hent_t* item)
{
  unsigned int i, found, done;
  pgno_t nextp;
  unsigned int chksum = 0;
  FFDB_DBT ekey;
  unsigned char *ekdata = 0;
  ffdb_datap_t* datap = 0;


  /* first get page for this bucket */
  item->pagep = ffdb_get_page (hashp, item->bucket, HASH_BUCKET_PAGE,
			       FFDB_PAGE_CREATE, &item->pgno);
  if (item->pagep == 0) {
    fprintf (stderr, "Cannot get page for bucket %d\n", item->bucket);
    item->status = ITEM_ERROR;
    return -1;
  }
#ifdef _FFDB_DEBUG
  fprintf (stderr, "Allocate page %d for bucket %d num entry %d\n", item->pgno, 
	   item->bucket, NUM_ENT(item->pagep));
#endif

  /* do a quick checksum on data */
  if (val) {
    chksum = 0;
    chksum = __ffdb_crc32_checksum (chksum, val->data, val->size);
    item->data_chksum = chksum;
  }
  else
    item->data_chksum = 0;

  if (NUM_ENT(item->pagep) == 0) { /* new page or nothing on it */
    item->pgndx = 0;
    item->status = ITEM_NO_MORE;
    _ffdb_init_page (hashp, item->pagep, item->pgno, HASH_BUCKET_PAGE);
    return 0;
  }

#ifdef _FFDB_STATISTICS
  FFDB_LOCK (hashp->lock);
  hash_collisions++;
  FFDB_UNLOCK (hashp->lock);
#endif

  /* set item page number information */
  item->pgno = CURR_PGNO (item->pagep);
  item->status = ITEM_NO_MORE;
  found = 0;
  done = 0;
  while (!done) {
    for (i = 0; i < NUM_ENT(item->pagep); i++) {
      /* get key for this index */
      ekdata = KEY(item->pagep, i);
      ekey.data = ekdata;
      ekey.size = KEY_LEN(item->pagep, i);
      /* get data for this index */
      datap = DATAP(item->pagep, i);

      /* We do not allow duplicated keys */
      if (hashp->h_compare(key, &ekey) == 0) {
	found = 1;
	break;
      }
    }
    if (found) /* We have found an item */
      done = 1;
    else {
      /* This is the page */
      if (NEXT_PGNO(item->pagep) == INVALID_PGNO)
	done = 1;
      else {
	nextp = NEXT_PGNO(item->pagep);

	/* Put this page back */
	ffdb_put_page (hashp, item->pagep, HASH_BUCKET_PAGE, 0);

	/* get the next page */
	item->pagep = ffdb_get_page (hashp, nextp, HASH_OVFL_PAGE,
				     0, &item->pgno);
	if (item->pagep == 0) {
	  fprintf (stderr, "Cannot get next page for bucket %d at page %d\n", 
		   item->bucket, item->pgno);
	  item->status = ITEM_ERROR;
	  return -1;
	}
#ifdef _FFDB_DEBUG
	fprintf (stderr, "Could not find item on primary page, check next page at %d\n", nextp);
#endif
      }
    }
  }
  /* no matter whether we found the item or not, the index = i */
  item->pgndx = i;

  if (found) {
    /* We found the key */
    item->status = ITEM_OK;
    item->key_off = KEY_OFF(item->pagep, i);
    item->key_len = KEY_LEN(item->pagep, i);
    item->data_off = DATAP_OFF(item->pagep, i);
  }
  return 0;
}


int ffdb_get_item (ffdb_htab_t* hashp,
		   const FFDB_DBT* key, FFDB_DBT* val,
		   ffdb_hent_t* item, int freepage)
{
  int status;
  ffdb_datap_t* datap;

  /* first get data pointer */
  datap = DATAP(item->pagep, item->pgndx);
#ifdef _FFDB_DEBUG
  fprintf (stderr, "Found Key %s information: \n", (char *)key->data);
  fprintf (stderr, "At first page %d\n", datap->first);
  fprintf (stderr, "At offset %d\n", datap->offset);
  fprintf (stderr, "Length of data %d\n", datap->len);
  fprintf (stderr, "Checksum of data is 0x%x\n", datap->chksum);
  fprintf (stderr, "freepage is = %d\n", freepage);
#endif

  /* Now I have to hop to data page to get this data item */
  status = _ffdb_get_data (hashp, item, val, datap);
  if (status != 0) {
    fprintf (stderr, "Cannot get data on page %d at offset %d\n",
	     datap->first, datap->offset);
    return status;
  }

 /* Now I do not need this page */
  if (freepage)
    ffdb_put_page (hashp, item->pagep, HASH_RAW_PAGE, 0);

  return 0;
}

/**
 * Add a pair of key and data onto a page (hash page) represented by
 * page address and page number
 */
static int
_ffdb_add_item_on_page (ffdb_htab_t* hashp, void* pagep, pgno_t page,
			FFDB_DBT* key, const FFDB_DBT* val,
			unsigned int data_chksum)
{
  unsigned int n, off, soff;
  ffdb_datap_t datap;
  pgno_t dpage, fpage;
  void* memp;
  int status, reuse;
  
  n = NUM_ENT(pagep);
  /* Find place to put the key */
  off = OFFSET(pagep) - key->size + 1;
  memmove (pagep + off, key->data, key->size);

  /* Set Key Offset Value */
  KEY_OFF(pagep, n) = off;
  KEY_LEN(pagep, n) = key->size;

  /*  Find place to put data pointer value */
  off -= sizeof(ffdb_datap_t);
  soff = off;
  /* I have to align this pointer to 4 byte boundary */
  ALIGN_DATAP_OFFSET_VAL(off);

  /* I will fill 0 to the gap of data and key */
  /* ------off---soff-----000key */
  if (off != soff) 
    memset (pagep + off + sizeof(ffdb_datap_t), 0, soff - off);


  /* Here I have to figure out where to put the data */
  reuse = 0;
  fpage = _ffdb_data_page(hashp, 0, &reuse);

#ifdef _FFDB_DEBUG
  fprintf (stderr, "Get data page %d to store data from hash page %d at level %d\n",
	   fpage, page, hashp->hdr.ovfl_point);
#endif

  /* Now I need check this page: get this page first */
  memp = ffdb_get_page (hashp, fpage, HASH_DATA_PAGE, FFDB_CREATE,
			&dpage);
  if (!memp) {
    fprintf (stderr, "cannot get data page for at page number %d\n",
	     dpage);
    return -1;
  }
  

  /* Now Put data on this page: this data can expand multiple pages */
  datap.first = fpage;
  datap.offset = 0;
  datap.len = val->size;
  datap.chksum = data_chksum;

  /* add data to data page provided key page and key index in the page
   * datap offset and first page is updated in the add_data call 
   */
  status = _ffdb_add_data (hashp, page, n, val, memp, dpage, &datap);
  if (status != 0) {
    fprintf (stderr, "cannot put data into data page at page number %d\n",
	     dpage);
    return -1;
  }
#ifdef _FFDB_DEBUG
  fprintf (stderr, "Data len %d is stored at page %d offset %d checksum 0x%x\n",
	   datap.len, datap.first, datap.offset, datap.chksum);
#endif
  memmove (pagep + off, &datap, sizeof(ffdb_datap_t));
  DATAP_OFF(pagep, n) = off;

  /* Update this page information */
  NUM_ENT(pagep) = n + 1;
  OFFSET(pagep) = off - 1;

  return 0;
}

/**
 * Replace a hash item on the page identified by item structure
 * Data can be only replaced when the new data has size <= the existing
 * data size
 */
static int
_ffdb_replace_item_on_page (ffdb_htab_t* hashp,
			    FFDB_DBT* key, const FFDB_DBT* val,
			    ffdb_hent_t* item)
{
  int status;
  ffdb_datap_t* datap;
  pgno_t dpage;
  void* memp;

  /* Get current data pointer information of this key */
  datap = DATAP (item->pagep, item->pgndx);
  
#ifdef _FFDB_DEBUG
  fprintf (stderr, "Replace Key %s information: \n", (char *)key->data);
  fprintf (stderr, "At first page %d\n", datap->first);
  fprintf (stderr, "At offset %d\n", datap->offset);
  fprintf (stderr, "Length of data %d\n", datap->len);
  fprintf (stderr, "Checksum of data is 0x%x\n", datap->chksum);
  fprintf (stderr, "With new data size of %d\n", val->size);
  fprintf (stderr, "new check sum = 0x%x\n", item->data_chksum);
#endif
  if (val->size > datap->len) {
    fprintf (stderr, "Replacement data size %d is larger than the existing data size %d\n", val->size, datap->len);
    return -1;
  }
  
  /* Now I can put data on the page pointed by data pointer */
  memp = ffdb_get_page (hashp, datap->first, HASH_DATA_PAGE, 0,
			&dpage);
  if (!memp) {
    fprintf (stderr, "Cannot get data page %d to replace data at offset %d\n",	     datap->first, datap->offset);
    return -1;
  }

  /* Change data pointer value */
  datap->len = val->size;
  datap->chksum = item->data_chksum;

  status = _ffdb_replace_data (hashp, val, memp, dpage, datap);

  if (status != 0) {
    fprintf (stderr, "Cannot replace data of len %d on page %d at offset %d\n",
	     datap->len, dpage, datap->offset);
    return -1;
  }

  return 0;
}


/**
 * Add a pair of key and data into the hash database
 *
 * The page associated with this pair should be cached
 *
 */
int ffdb_add_pair (ffdb_htab_t* hashp,
		   FFDB_DBT* key, const FFDB_DBT* val,
		   ffdb_hent_t* item, int replace)
{
  int status;
  if (!replace)
    status = _ffdb_add_item_on_page (hashp, item->pagep, item->pgno,
				     key, val, item->data_chksum);
  else
    status = _ffdb_replace_item_on_page (hashp, key, val, item);

  if (status != 0) 
    /* Now I do not need this page */
    ffdb_put_page (hashp, item->pagep, HASH_BUCKET_PAGE, 0);
  else
    ffdb_put_page (hashp, item->pagep, HASH_BUCKET_PAGE, 1);

  return status;
}


/**
 * Add a pair to an overflow page
 */
int ffdb_add_ovflpage (ffdb_htab_t* hashp,
		       FFDB_DBT* key, const FFDB_DBT* val,
		       ffdb_hent_t* item)
{
  pgno_t ovflpage, tp;
  void*  opagep;
  int    status, reuse;

  /* Find out next overflow page number */
  reuse = 0;
  ovflpage = _ffdb_ovfl_page (hashp, &reuse);
#ifdef _FFDB_DEBUG
  fprintf (stderr, "Allocate an overflow page %d for primary page %d at level %d\n",
	   ovflpage, item->pgno, hashp->hdr.ovfl_point);
#endif
  
  /* Get this page and this page is initialized */
  opagep = ffdb_get_page (hashp, ovflpage, HASH_OVFL_PAGE, FFDB_CREATE, &tp);
  if (!opagep) {
    fprintf (stderr, "Cannot get an over flow page %d for primary page %d\n",
	     ovflpage, item->pgno);
    /* Now I do not need this page */
    ffdb_put_page (hashp, item->pagep, HASH_BUCKET_PAGE, 0);
    return -1; 
  }
  if (reuse)
    _ffdb_init_page (hashp, opagep, ovflpage, HASH_OVFL_PAGE);

#ifdef _FFDB_STATISTICS
  hash_overflows++;
#endif
  
  /* Set up page link */
  NEXT_PGNO(item->pagep) = ovflpage;
  PREV_PGNO(opagep) = item->pgno;

  /* Add this pair to the new page */
  status = _ffdb_add_item_on_page (hashp, opagep, ovflpage,
				   key, val, item->data_chksum);

  if (status != 0) {
    ffdb_put_page (hashp, opagep, HASH_OVFL_PAGE, 0);
    ffdb_put_page (hashp, item->pagep, HASH_BUCKET_PAGE, 0);
  }
  else {
    ffdb_put_page (hashp, opagep, HASH_OVFL_PAGE, 1);
    ffdb_put_page (hashp, item->pagep, HASH_BUCKET_PAGE, 1);
  }
  return 0;
}

/**
 * Just copy content of key and its associated data pointer to a page
 * at the right location
 * We already checked this page can fit this key and data pointer
 */
static int
_ffdb_write_key_datap_to_page (ffdb_htab_t* hashp, FFDB_DBT* key, 
			       ffdb_datap_t *datap,
			       void* pagep, pgno_t page)
{
  unsigned int n, off, soff;
  int status;
  
  n = NUM_ENT(pagep);
  /* Find place to put the key */
  off = OFFSET(pagep) - key->size + 1;
  memmove (pagep + off, key->data, key->size);

  /* Set Key Offset Value */
  KEY_OFF(pagep, n) = off;
  KEY_LEN(pagep, n) = key->size;

  /*  Find place to put data pointer value */
  off -= sizeof(ffdb_datap_t);
  soff = off;
  /* I have to align this pointer to 4 byte boundary */
  ALIGN_DATAP_OFFSET_VAL(off);

  /* I will fill 0 to the gap of data and key */
  /* ------off---soff-----000key */
  if (off != soff) 
    memset (pagep + off + sizeof(ffdb_datap_t), 0, soff - off);

  /* write data pointer value */
  memmove (pagep + off, datap, sizeof(ffdb_datap_t));
  DATAP_OFF(pagep, n) = off;

  /* Need update data item information about key page and key index */
  /* The data items are on data pages, this could be slow           */
  status = _ffdb_update_data_info (hashp, datap, page, n);

  /* Update this page information */
  NUM_ENT(pagep) = n + 1;
  OFFSET(pagep) = off - 1;

  return status;
}

/**
 * Add a pair of key and data pointer value to a bucket page
 * which could be the first of a few overflow pages
 *
 * The data pointer points to data pages that could be allocated
 * in previous levels
 */
static int
_ffdb_add_key_datap_to_bucket (ffdb_htab_t* hashp, FFDB_DBT* key, 
			       ffdb_datap_t *datap, 
			       unsigned int bucket)
{
  pgno_t page, ovflpage, nextpage, tp;
  void *pagep, *opagep, *npagep;
  int needovfl, reuse;

  /* First grab a page related to this bucket */
  npagep = pagep = ffdb_get_page (hashp, bucket, HASH_BUCKET_PAGE, 
				  FFDB_CREATE, &page);
  if (!pagep) {
    fprintf (stderr, "Cannot get page for bucket %d\n", bucket);
    return -1;
  }

  needovfl = 1;
  nextpage = page;
  while (npagep) {
    /* check whether this pair should fit on this page */
    /* the following macro does not care the data part */
    if (PAIRFITS (npagep, key, dumb)) {
      _ffdb_write_key_datap_to_page (hashp, key, datap, pagep, page);
      needovfl = 0;
      /* release this page */
      ffdb_put_page (hashp, pagep, HASH_BUCKET_PAGE, 1);
      /* I am done here */
      npagep = 0;
      needovfl = 0;
    }
    else {
      nextpage = NEXT_PGNO (pagep);
      if (nextpage == INVALID_PGNO)
	npagep = 0;
      else {
	npagep = ffdb_get_page (hashp, nextpage, HASH_OVFL_PAGE, 0, &tp);
	if (!npagep) {
	  fprintf (stderr, "Cannot get overflow (expanded) page %d \n",
		   nextpage);
	  /* release the previous page */
	  ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
	  return -1;
	}
	/* release the previous page */
	ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
	
	/* assign pointers */
	pagep = npagep;
	page = nextpage;
      }
    }
  }
      
  if (needovfl) {
    reuse = 0;
    ovflpage = _ffdb_ovfl_page (hashp, &reuse);
#ifdef _FFDB_DEBUG
    fprintf (stderr, "Get an overflow page (expanded) %d for page %d\n",
	     ovflpage, page);
#endif
    /* Get this page and this page is initialized */
    opagep = ffdb_get_page (hashp, ovflpage, HASH_OVFL_PAGE,
			    FFDB_CREATE, &tp);
    if (!opagep) {
      fprintf (stderr, "Cannot allocate expanded overflow page %d for page %d\n",
	       ovflpage, page);
      /* release the previous page */
      ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
      return -1;
    }
    if (reuse) 
      _ffdb_init_page (hashp, opagep, ovflpage, HASH_OVFL_PAGE);

#ifdef _FFDB_STATISTICS
    hash_overflows++;
#endif

    /* set up correct link */
    NEXT_PGNO (pagep) = ovflpage;
    PREV_PGNO (opagep) = page;

    /* release the previous page , but the head is changed so 
     * we have to write it to disk
     */
    ffdb_put_page (hashp, pagep, TYPE(pagep), 1); 

    /* add key and data pointer to this page */
    _ffdb_write_key_datap_to_page (hashp, key, datap, opagep, ovflpage);

    /* release this page */
    ffdb_put_page (hashp, opagep, HASH_OVFL_PAGE, 1);
  }
  return 0;
}

/**
 * Split a bucket: this happens when a bucket is full. This bucket may not be 
 * splitted right away (overflow pages needed), but it will eventually 
 * will be splitted
 */
int 
ffdb_split_bucket (ffdb_htab_t* hashp, unsigned int oldbucket,
		   unsigned int newbucket, int isdoubling)
{
  unsigned int i;
  FFDB_DBT key;
  void *oldpagep, *temp_pagep;
  pgno_t oldpage, nextpage, tp;
  unsigned char *kdata = 0;
  ffdb_datap_t* datap = 0;
  int base_page = 1;

  /* This page may not be in memory or in use so FFDB_CREATE flag has to be 
   * used
   */
  oldpagep = ffdb_get_page (hashp, oldbucket, HASH_BUCKET_PAGE, FFDB_CREATE, &oldpage);
  if (!oldpagep) {
    fprintf (stderr, "Cannot find page for bucket %d\n", oldbucket);
    return -1;
  }

  /* get split buffer pointer */
  temp_pagep = hashp->split_buf;
  /* Copy the old page content to the temp page */
  memcpy (temp_pagep, oldpagep, hashp->hdr.bsize);

  /* reset the old page */
  _ffdb_init_page (hashp, oldpagep, oldpage, HASH_BUCKET_PAGE);
  ffdb_put_page (hashp, oldpagep, HASH_BUCKET_PAGE, 1);

  while (temp_pagep != 0) {
    /* get key and data pointer value */
    for (i = 0; i < NUM_ENT(temp_pagep); i++) {
      kdata = KEY(temp_pagep, i);
      key.data = kdata;
      key.size = KEY_LEN(temp_pagep, i);
      datap = DATAP(temp_pagep, i);

      /* Now we need to put this key and data pointer pair */
      if (_ffdb_call_hash (hashp, key.data, key.size) == oldbucket) 
	/* this stays with old page without changing data pointer value */
	_ffdb_add_key_datap_to_bucket (hashp, &key, datap, oldbucket);
      else 
	_ffdb_add_key_datap_to_bucket (hashp, &key, datap, newbucket);
    }
    
    /* get next page number */
    nextpage = NEXT_PGNO(temp_pagep);
    
    /* if this is the base page (regular hash page) */
    if (base_page) 
      base_page = 0;
    else  /* free overflow page */ 
      ffdb_delete_page (hashp, temp_pagep, HASH_OVFL_PAGE, isdoubling);
    
    if (nextpage != INVALID_PGNO) 
      temp_pagep = ffdb_get_page (hashp, nextpage, HASH_OVFL_PAGE, 0, &tp);
    else
      break;
  }
  return 0;
}


/**
 * Update Data Page Content when a key page are rearranged to fill gaps.
 * Inside the key page (overflow), there are datap value pointing to
 * data pages which contain key_page information for each data. These
 * key_page value need to be updated to the new location.
 * 
 */
static int
_ffdb_update_data_pages_content (ffdb_htab_t* hashp, void *pagep, 
				 pgno_t oldpage, pgno_t newpage,
				 pgno_t oldfirst, pgno_t oldlast,
				 pgno_t newfirst, pgno_t newlast)
{
  unsigned int i;
  void* dpagep;
  pgno_t dp;
  ffdb_data_header_t *header = 0;
  ffdb_datap_t* datap = 0;

  for (i = 0; i < NUM_ENT(pagep); i++) {
    datap = DATAP(pagep, i);

    /* First check whether this data page pointed by datap->first
     * is a page to be moved. If yes, we are doing nothing, since
     * its key_page content is already updated when the data page
     * is moved
     */
    if (datap->first >= newfirst && datap->first <= newlast)
      continue;

    /* get data page pointed back by this pointer */
    dpagep = ffdb_get_page (hashp, datap->first, HASH_DATA_PAGE, 0, &dp);
    if (!dpagep) {
      fprintf (stderr, "Cannot get data page %d pointed by overflow page %d item %d\n",
	       datap->first, oldpage, i);
      return -1;
    }
    
    /* change information on the key page */
    header = BIG_DATA_HEADER(dpagep, datap->offset);

    assert (header->key_page == oldpage);

    /* key_page data item entry points to this new page */
    header->key_page = newpage;

    /* Put back the page */
    ffdb_put_page (hashp, dpagep, TYPE(dpagep), 1);
  }
  return 0;
}


/**
 * Update all Key Pages (either primary or overflow) content when a data page
 * is moved. The data page contains one or more key page information through
 * key_page field. A new hash key item is accessed through key_page and 
 * key_indx pair. The data pointer value (datap->first) of this hash item
 * has to be updated to the new page.
 *
 * If the key_page is also among pages that are moved, the data header's 
 * key_page will be updated to its new location.
 */
static int
_ffdb_update_key_pages_content (ffdb_htab_t* hashp, void *pagep, 
				pgno_t oldpage, pgno_t newpage,
				pgno_t oldfirst, pgno_t oldlast,
				pgno_t newfirst, pgno_t newlast)
{
  unsigned int i, next;
  void* kpagep;
  pgno_t kp, newkp;
  ffdb_data_header_t *header = 0;
  ffdb_datap_t* datap = 0;

  next = FIRST_DATA_POS(pagep);
  for (i = 0; i < NUM_ENT(pagep); i++) {
    header = BIG_DATA_HEADER(pagep, next);

    /* get key page pointed back by this header */
    kpagep = ffdb_get_page (hashp, header->key_page, HASH_RAW_PAGE, 0, &kp);
    if (!kpagep) {
      fprintf (stderr, "Cannot get key page %d from data page %d item %d\n",
	       header->key_page, oldpage, i);
      return -1;
    }
    
    /* change information on the key page */
    datap = DATAP(kpagep, header->key_idx);

#ifdef _FFDB_DEBUG
    fprintf (stderr, "page %d Datap checksum = 0x%x header key = %d\n", oldpage,
	     datap->chksum, header->key_page);
#endif

    assert (datap->first == oldpage);

    /* Key page data item entry points to this new page */
    datap->first = newpage;

    /* Check whether the key page actually is going to be moved */
    /* this works based on assumption data pages are moved first */
    /* the data page entry points back to the new hash page that
     * is going to be moved to newkp
     */
    if (header->key_page >= oldfirst && header->key_page <= oldlast) {
      newkp = newfirst + (header->key_page - oldfirst);
#ifdef _FFDB_DEBUG
      fprintf (stderr, "Data Page %d contains hash key entry %d that moves to %d\n", oldpage, header->key_page, newkp);
#endif
      /* Now my header has to be changed to point to new page as well */
      header->key_page = newkp;
    }
    
    next = header->next;

    /* Put back the page */
    ffdb_put_page (hashp, kpagep, TYPE(kpagep), 1);
  }
  return 0;
}

/**
 * Clean out range of pages from pagepool cache
 * @param hashp hashtable pointer
 * @param first the first page to be removed
 * @param last the last page to be removed
 */
static void _ffdb_remove_pages (ffdb_htab_t* hashp, 
				pgno_t first, pgno_t last)
{
  unsigned int numpages, i;
  pgno_t page, tp;
  void *pagep;

  /* total number of pages to be removed */
  numpages = last - first + 1;

  page = first;
  for (i = 0; i < numpages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_RAW_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Warning: cannot find moved old page %d \n",
	       page);
      page++;
      continue;
    }
    if (ffdb_pagepool_delete (hashp->mp, pagep) != 0) 
      fprintf (stderr, "Cannot delete a moved page %d\n", page);
    page++;
  }
}
  


/**
 * Move pages when the file closes
 *
 * @param hashp hash table pointer
 * @param type this type of pages to be moved
 * @param oldfirst first page of old pages to be moved
 * @param oldlast the last page of old page to be moved
 * @param newfirst first page of new pages to move into
 * @param newlast last page of new pages to move into
 *
 * @return last page number on the linked list after successful move. 
 * If there is no data pages being moved, INVALID_PGNO is returned
 */
static pgno_t
_ffdb_move_pages (ffdb_htab_t* hashp, unsigned short type,
		  pgno_t oldfirst, pgno_t oldlast,
		  pgno_t newfirst, pgno_t newlast)
{
  pgno_t dlast, page, rpage, tp;
  pgno_t prevp, nextp;
  unsigned int numpages, i;
  void *pagep, *prevpagep, *nextpagep;

  dlast = INVALID_PGNO;
  /* Create page header array to keep track old page information */
  numpages = oldlast - oldfirst + 1;

  /* Now we can move pages */
  page = oldfirst;
  rpage = newfirst;
  for (i = 0; i < numpages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_RAW_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot grab old page %d to move to new page %d.\n",
	       page, rpage);
      return dlast;
    }

    if (TYPE(pagep) != type) {
      page++;
      rpage++;
      ffdb_put_page (hashp, pagep, HASH_RAW_PAGE, 0);
      continue;
    }

    /* Get its original previous and next page number */
    nextp = NEXT_PGNO(pagep);
    prevp = PREV_PGNO(pagep);

    if (nextp != INVALID_PGNO) {
      if (nextp < oldfirst || nextp > oldlast) {
	/* The next page is not going to be moved, its previous link 
	 * Should be updated
	 */
	nextpagep = ffdb_get_page (hashp, nextp, HASH_RAW_PAGE, 0, &tp);
	if (!nextpagep) {
	  fprintf (stderr, "Cannot get next page %d not in the moving list for page %d\n", nextp, page);
	  ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
	  return dlast;
	}

	PREV_PGNO(nextpagep) = rpage;
	ffdb_put_page (hashp, nextpagep, TYPE(nextpagep), 1);
      }
      else {
	/* I have to update this page next link to a new page number */
	NEXT_PGNO(pagep) = newfirst + (nextp - oldfirst);
      }
    }
    else {
      /* This is last data page in these moved page.
       * update last new location of this page
       */
      if (type == HASH_DATA_PAGE)
	dlast = rpage;
    }


    if (prevp != INVALID_PGNO) {
      if (prevp < oldfirst || prevp > oldlast) {
	/* The previous page is not going to be moved, its previous link 
	 * Should be updated
	 */
	prevpagep = ffdb_get_page (hashp, prevp, HASH_RAW_PAGE, 0, &tp);
	if (!prevpagep) {
	  fprintf (stderr, "Cannot get previous page %d not in the moving list for page %d\n", prevp, page);
	  ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
	  return dlast;
	}

	NEXT_PGNO(prevpagep) = rpage;
	ffdb_put_page (hashp, prevpagep, TYPE(prevpagep), 1);
      }
      else {
	/* The previous link of this page has to be updated */
	PREV_PGNO(pagep) = newfirst + (prevp - oldfirst);
      }
    }

    /* Now update current page number */
    CURR_PGNO(pagep) = rpage;

#ifdef _FFDB_DEBUG
    fprintf (stderr, "Move old page type 0x%x %d to new page %d\n", type, page, rpage);
#endif

    /* Now update the hash data content on the primary pages
     * pointed by key field of data on this page
     */
    if (type == HASH_DATA_PAGE) {
      if (_ffdb_update_key_pages_content (hashp, pagep, page, rpage, 
					  oldfirst, oldlast,
					  newfirst, newlast) != 0) {
	/* Release this page and update page number */
	ffdb_put_page (hashp, pagep, HASH_DATA_PAGE, 0);
	return dlast;
      }
#ifdef _FFDB_DEBUG
      fprintf (stderr, "Update primary page for data page %d %d page = %d page sign = 0x%x\n", page, rpage, CURR_PGNO(pagep), PAGE_SIGN(pagep));
#endif
    }
    else if (type == HASH_BUCKET_PAGE) {
      /* Update the key page entry inside each data page pointed 
       * by the overflow page. The data pages that have been moved
       *will not be updated since they are updated by the above routine
       */
      if (_ffdb_update_data_pages_content (hashp, pagep, page, rpage,
					   oldfirst, oldlast,
					   newfirst, newlast) != 0) {
	/* Release this page and update page number */
	ffdb_put_page (hashp, pagep, HASH_DATA_PAGE, 0);
	return dlast;
      }
#ifdef _FFDB_DEBUG
      fprintf (stderr, "Update data page for primary page %d %d page = %d page sign = 0x%x\n", page, rpage, CURR_PGNO(pagep), PAGE_SIGN(pagep));
#endif
    }
      

    /* Release this page and update page number */
    ffdb_pagepool_change_page (hashp->mp, pagep, CURR_PGNO(pagep));

    /* update local variables */
    rpage++;
    page++;    
  }
    
  return dlast;
}

/**
 * Optional to reorganize the pages to fill holes left by 
 * not fully using all buckets pages at a particular level
 */
int
ffdb_rearrage_pages_on_close (ffdb_htab_t* hashp)
{
  pgno_t funused, lunused;
  pgno_t dlast, dfirst;
  unsigned int num_unused, num_pages;
  unsigned int num_page_moved;

  /* Find out how many pages are not used by reguler hash pages */
  num_unused = hashp->hdr.high_mask - hashp->hdr.max_bucket;
  if (num_unused == 0)
    return 0;

  /* First and last unused page */
  BUCKET_TO_PAGE(hashp->hdr.max_bucket + 1, funused);
  BUCKET_TO_PAGE(hashp->hdr.high_mask, lunused);
#if 1
  fprintf (stderr, "There are %d pages (buckets) not used starting at page %d and end at %d on last level %d\n",
	   num_unused, funused, lunused, hashp->hdr.ovfl_point);
#endif

  /* This is the last page used so far */
  dlast = hashp->hdr.spares[hashp->hdr.ovfl_point + 1] - 1;
  assert (dlast > 0);
  assert (dlast > lunused);

#if 0
  fprintf (stderr, "Last data page = %d\n", dlast);
#endif

  /* Now we need to figure out what pages to move
   * These pages may contain overflow pages and data pages
   * which require different treatment. We move data pages
   * first
   */
  num_pages = dlast - lunused;
  num_page_moved = (num_pages < num_unused) ? num_pages : num_unused;

  /* First page to be moved */
  dfirst = dlast - num_page_moved + 1;
  /* Update last free page to be used */
  lunused = funused + num_page_moved  - 1;

#if 1
  fprintf (stderr, "We are going to move %d pages from %d to %d to new pages from %d to %d\n", num_page_moved, dfirst, dlast, funused, lunused);
#endif

  /* move data pages first */
  _ffdb_move_pages (hashp, HASH_DATA_PAGE, dfirst, dlast,
		    funused, lunused);

  /* move overflow page next */
  _ffdb_move_pages (hashp, HASH_BUCKET_PAGE, dfirst, dlast,
		    funused, lunused);

  _ffdb_move_pages (hashp, HASH_FREE_PAGE, dfirst, dlast,
		    funused, lunused);

  hashp->hdr.num_moved_pages = num_page_moved;

  return 0;
}


/**
 * Optional to reorganize the pages to move pages occuping
 * reguler hash pages when the database expects new inserts
 */
int
ffdb_rearrage_pages_on_open (ffdb_htab_t* hashp)
{
  pgno_t oldfirst, oldlast, newfirst, newlast;
  pgno_t lastdpage;
  unsigned int num_page_moved;

  /* Check whether any pages has been moved when the database is closed */
  if (hashp->hdr.num_moved_pages == 0)
    return 0;
  
  if (hashp->hdr.max_bucket == hashp->hdr.high_mask) 
    return 0;

  /* We have not used every single hash bucket at the last level */
  BUCKET_TO_PAGE(hashp->hdr.max_bucket + 1, oldfirst);

  /* We will get number of moved page from meta header which
   * records how many pages are moved when this file was closed
   */
  num_page_moved = hashp->hdr.num_moved_pages;
  
  oldlast = oldfirst + num_page_moved - 1;


  /* We figure out where these page are going to */
  newlast = hashp->hdr.spares[hashp->hdr.ovfl_point + 1] - 1;
  newfirst = newlast - num_page_moved + 1;

#if 1
  fprintf (stderr, "Page can be moved starting at page %d and end at %d on last level %d to page %d to %d\n",
	   oldfirst, oldlast, hashp->hdr.ovfl_point, newfirst, newlast);
#endif


  /* Move data pages first */
  lastdpage = _ffdb_move_pages (hashp, HASH_DATA_PAGE, oldfirst, oldlast,
				newfirst, newlast);

  /* move overflow page next */
  _ffdb_move_pages (hashp, HASH_BUCKET_PAGE, oldfirst, oldlast,
		    newfirst, newlast);
  
  _ffdb_move_pages (hashp, HASH_FREE_PAGE, oldfirst, oldlast,
		    newfirst, newlast);

  /* If there are no data pages moved, unlikely but possible */
  if (lastdpage == INVALID_PGNO) 
    hashp->curr_dpage = ffdb_last_data_page (hashp, newfirst - 1);

  /* Reset the number of moved pages */
  hashp->hdr.num_moved_pages = 0;

  /* Clean out those pages from page cache to make memory
   * the same as if there are no pages have been here
   */
  _ffdb_remove_pages (hashp, oldfirst, oldlast);

  return 0;
}


void 
ffdb_reduce_filesize (ffdb_htab_t* hashp)
{
  pgno_t last;
  off_t length;

  /* This is the last page before all pages are moved */
  last = hashp->hdr.spares[hashp->hdr.ovfl_point + 1] - 1;

  /* This is the last page after the pages are moved */
  last = last - hashp->hdr.num_moved_pages;

  /* Now we need to shrink the file size */
  if (hashp->hdr.num_moved_pages > 0) {
    length = (off_t)hashp->hdr.bsize * (last + 1);

    if (ftruncate (hashp->fp, length) != 0) {
      fprintf (stderr, "Cannot truncate database %s to length %ld\n",
	       hashp->fname, length);
    }
  }
}

/**
 * Code related to configurations
 */
int 
ffdb_set_configs (ffdb_htab_t* hashp, ffdb_all_config_info_t* configs)
{
  unsigned int i, k, cnum, numelems, tnums;
  pgno_t page, tp;
  void* pagep;
  ffdb_config_info_t* cfig;

  /* index into all configs */
  cnum = 0;
  /* total number of configurations */
  tnums = hashp->hdr.num_cfigs;

  page = hashp->hdr.cfig_page;
  for (i = 0; i < hashp->hdr.cfig_npages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_CONFIG_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get configuration page %d.\n", page);
      return -1;
    }
    /* number of configurations on this page */
    numelems = (hashp->hdr.bsize - PAGE_OVERHEAD)/(sizeof(ffdb_config_info_t));
    if (numelems > tnums)
      numelems = tnums;
    NUM_ENT(pagep) = numelems;

    for (k = 0; k < NUM_ENT(pagep); k++) {
      cfig = CONFIG_INFO(pagep, k);
      memcpy (cfig, &(configs->allconfigs[cnum++]),
	      sizeof(ffdb_config_info_t));
    }

    /* release page */
    ffdb_put_page (hashp, pagep, TYPE(pagep), 1);
    page++;
    tnums -= numelems;
  }
  return 0;
}


int 
ffdb_get_configs (ffdb_htab_t* hashp, ffdb_all_config_info_t* configs)
{
  unsigned int i, k, cnum;
  pgno_t page, tp;
  void* pagep;
  ffdb_config_info_t* cfig;

  cnum = 0;
  page = hashp->hdr.cfig_page;
  for (i = 0; i < hashp->hdr.cfig_npages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_CONFIG_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get configuration page %d.\n", page);
      return -1;
    }
    
    for (k = 0; k < NUM_ENT(pagep); k++) {
      cfig = CONFIG_INFO(pagep, k);
      memcpy (&(configs->allconfigs[cnum++]), cfig,
	      sizeof(ffdb_config_info_t));
    }

    /* release page */
    ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
    page++;
  }
  return 0;
}


int 
ffdb_set_config_info (ffdb_htab_t* hashp, ffdb_config_info_t* config)
{
  unsigned int i, k;
  pgno_t page, tp;
  void* pagep;
  ffdb_config_info_t* cfig;

  page = hashp->hdr.cfig_page;
  for (i = 0; i < hashp->hdr.cfig_npages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_CONFIG_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get configuration page %d.\n", page);
      return -1;
    }
    
    for (k = 0; k < NUM_ENT(pagep); k++) {
      cfig = CONFIG_INFO(pagep, k);
      if (cfig->index == config->index && cfig->config == config->config) {
	memcpy (cfig, config, sizeof(ffdb_config_info_t));
	ffdb_put_page (hashp, pagep, TYPE(pagep), 1);
	return 0;
      }
    }
    ffdb_put_page (hashp, pagep, TYPE(pagep), 0);    /* release page */
    page++;
  }
  return -1;
}


int 
ffdb_get_config_info (ffdb_htab_t* hashp, unsigned int confignum,
		      ffdb_config_info_t* config)
{
  unsigned int i, k;
  pgno_t page, tp;
  void* pagep;
  ffdb_config_info_t* cfig;

  page = hashp->hdr.cfig_page;
  for (i = 0; i < hashp->hdr.cfig_npages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_CONFIG_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get configuration page %d.\n", page);
      return -1;
    }
    
    for (k = 0; k < NUM_ENT(pagep); k++) {
      cfig = CONFIG_INFO(pagep, k);
      if (cfig->config == confignum) {
	memcpy (config, cfig, sizeof(ffdb_config_info_t));
	ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
	return 0;
      }
    }
    ffdb_put_page (hashp, pagep, TYPE(pagep), 0);    /* release page */
    page++;
  }
  return -1;
}


/**
 * Set user meta information in the 2nd page
 */
int 
ffdb_set_uinfo (ffdb_htab_t* hashp, unsigned char* data,
		unsigned int len)
{
  unsigned int i, copylen, offset, idx, rlen;
  pgno_t page, tp;
  void* pagep;
  
  idx = 0;
  rlen = len;
  page = hashp->hdr.uinfo_page;
  for (i = 0; i < hashp->hdr.uinfo_npages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_UINFO_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get page %d to store user information.\n",
	       page);
      return -1;
    }
    if (i == 0) {
      /* Set number of item to 1 on the first page */
      NUM_ENT(pagep) = 1;
      USER_DATA_LEN(pagep) = rlen;
      if (rlen > hashp->hdr.bsize - PAGE_OVERHEAD - sizeof(unsigned int)) 
	copylen = hashp->hdr.bsize - PAGE_OVERHEAD - sizeof(unsigned int);
      else 
	copylen = rlen;
      offset = PAGE_OVERHEAD + sizeof(unsigned int);
    }
    else {
      /* Set number of item to 0 on the other pages */
      NUM_ENT(pagep) = 0;
      if (rlen > hashp->hdr.bsize - PAGE_OVERHEAD) 
	copylen = hashp->hdr.bsize - PAGE_OVERHEAD;
      else
	copylen = rlen;
      offset = PAGE_OVERHEAD;
    }
    /* where to start copy from data */
    idx = len - rlen;

    memcpy ((unsigned char *)pagep + offset, &data[idx], copylen);

    /* put page out */
    ffdb_put_page (hashp, pagep, TYPE(pagep), 1);

    rlen -= copylen;

    page++;
  }
  return 0;
}


/**
 * Get user meta information in the 2nd page
 */
int 
ffdb_get_uinfo (ffdb_htab_t* hashp, unsigned char data[],
		unsigned int *len)
{
  unsigned int i, copylen, offset, idx, rlen;
  pgno_t page, tp;
  void* pagep;
  
  idx = 0;
  rlen = *len;
  page = hashp->hdr.uinfo_page;
  for (i = 0; i < hashp->hdr.uinfo_npages; i++) {
    pagep = ffdb_get_page (hashp, page, HASH_UINFO_PAGE, 0, &tp);
    if (!pagep) {
      fprintf (stderr, "Cannot get page %d to retrieve user information.\n",
	       page);
      return -1;
    }
    if (i == 0) {
      if (*len < USER_DATA_LEN(pagep)) {
	fprintf (stderr, "Error: provided user buffer size %d < stored data length %d.\n", *len, USER_DATA_LEN(pagep));
	return -1;
      }
      *len = USER_DATA_LEN(pagep);
      rlen = *len;

      if (rlen > hashp->hdr.bsize - PAGE_OVERHEAD - sizeof(unsigned int)) 
	copylen = hashp->hdr.bsize - PAGE_OVERHEAD - sizeof(unsigned int);
      else 
	copylen = rlen;
      offset = PAGE_OVERHEAD + sizeof(unsigned int);
    }
    else {
      if (rlen > hashp->hdr.bsize - PAGE_OVERHEAD) 
	copylen = hashp->hdr.bsize - PAGE_OVERHEAD;
      else
	copylen = rlen;
      offset = PAGE_OVERHEAD;
    }
    /* where to start copy from data */
    idx = *len - rlen;

    memcpy (&data[idx], (unsigned char *)pagep + offset, copylen);

    /* put page out */
    ffdb_put_page (hashp, pagep, TYPE(pagep), 0);

    rlen -= copylen;

    page++;
  }
  return 0;
}

/***************************************************************************
 *         Cursor related routines                                         *
 ***************************************************************************/
int 
ffdb_cursor_find_by_key (ffdb_htab_t* hashp, ffdb_crs_t* cursor,
			 FFDB_DBT* key, FFDB_DBT* data,
			 unsigned int flags)
{
  pgno_t bucket, tp, nextp;
  unsigned char* ekdata;
  unsigned int   eksize;

  if (flags == FFDB_FIRST) {
    if (cursor->item.pagep) {
      /* need to release the previous page */
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      cursor->item.pagep = 0;
    }
      
    bucket = 0;
    cursor->item.pagep = ffdb_get_page (hashp, bucket,
					HASH_BUCKET_PAGE, 0, &tp);
    if (!(cursor->item.pagep)) {
      fprintf (stderr, "Cannot get page for the first bucket\n");
      cursor->item.status = ITEM_ERROR;
      return -1;
    }
    /* Skip empty buckets */
    while (NUM_ENT(cursor->item.pagep) == 0 && 
	   bucket <= hashp->hdr.max_bucket) {
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      /* Get next bucket */
      cursor->item.pagep = 0;
      bucket++;

      cursor->item.pagep = ffdb_get_page (hashp, bucket, HASH_BUCKET_PAGE, 
					  0, &tp);
      if (!(cursor->item.pagep)) {
	fprintf (stderr, "Cannot get page for bucket %d for cursor.\n",
		 bucket);
	cursor->item.status = ITEM_ERROR;
	return -1;
      }
    }
    if (bucket > hashp->hdr.max_bucket) {
      /* empty database */
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      cursor->item.status = ITEM_NO_MORE;
      cursor->item.pagep = 0;
      return FFDB_NOT_FOUND;
    }
    cursor->item.bucket = bucket;
    cursor->item.pgno = tp;
    cursor->item.pgndx = 0;
    cursor->item.status = ITEM_OK;
  }
  else if (flags == FFDB_LAST) {
    if (cursor->item.pagep) {
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      cursor->item.pagep = 0;
    }
    bucket = hashp->hdr.max_bucket;
    
    cursor->item.pagep = ffdb_get_page (hashp, bucket,
					HASH_BUCKET_PAGE, 0, &tp);
    if (!(cursor->item.pagep)) {
      fprintf (stderr, "Cannot get page for the last bucket\n");
      cursor->item.status = ITEM_ERROR;
      return -1;
    }
    /* Skip empty buckets */
    while (NUM_ENT(cursor->item.pagep) == 0) {
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      /* Get next bucket */
      cursor->item.pagep = 0;
      bucket--;

      cursor->item.pagep = ffdb_get_page (hashp, bucket, HASH_BUCKET_PAGE, 
					  0, &tp);
      if (!(cursor->item.pagep)) {
	fprintf (stderr, "Cannot get page for bucket %d for cursor.\n",
		 bucket);
	cursor->item.status = ITEM_ERROR;
	return -1;
      }
      if (bucket == 0 && NUM_ENT(cursor->item.pagep) == 0) {
	/* this is the last one and it is empty */
	ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
	cursor->item.status = ITEM_NO_MORE;
	cursor->item.pagep = 0;
	return FFDB_NOT_FOUND;
      }
    }

    cursor->item.pgno = tp;
    cursor->item.bucket = bucket;    
    cursor->item.pgndx = NUM_ENT(cursor->item.pagep) - 1;
    cursor->item.status = ITEM_OK;
  }
  else if (flags == FFDB_NEXT) {
    /* We have reached the last data on the key page */
    if (cursor->item.pgndx == NUM_ENT(cursor->item.pagep) - 1) {
      /* Get next page (overplow page) */
      nextp = NEXT_PGNO(cursor->item.pagep);
      /* put back this page */
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      if (cursor->item.bucket >= hashp->hdr.max_bucket &&
	  nextp == INVALID_PGNO) {
	/* We are done */
	cursor->item.status = ITEM_NO_MORE;
	cursor->item.pagep = 0;
	return FFDB_NOT_FOUND;
      }
      if (nextp == INVALID_PGNO) {
	/* Increase bucket by one: next bucket */
	cursor->item.bucket++;
	/* Get new page */
	cursor->item.pagep = ffdb_get_page (hashp, cursor->item.bucket,
					    HASH_BUCKET_PAGE, 0, &tp);

	if (!(cursor->item.pagep)) {
	  fprintf (stderr, "Cannot get page for bucket %d for cursor.\n",
		   cursor->item.bucket);
	  cursor->item.status = ITEM_ERROR;
	  return -1;
	}

	/* Skip empty buckets */
	while (NUM_ENT(cursor->item.pagep) == 0) {
	  ffdb_put_page (hashp, cursor->item.pagep, 
			 TYPE(cursor->item.pagep), 0);
	  if (cursor->item.bucket == hashp->hdr.max_bucket) {
	    cursor->item.pagep = 0;
	    /* End of buckets and we are done */
	    cursor->item.status = ITEM_NO_MORE;
	    return FFDB_NOT_FOUND;
	  }
	  /* Get next bucket */
	  cursor->item.pagep = 0;
	  cursor->item.bucket++;

	  cursor->item.pagep = ffdb_get_page (hashp, cursor->item.bucket, 
					      HASH_BUCKET_PAGE, 0, &tp);
	  if (!(cursor->item.pagep)) {
	    fprintf (stderr, "Cannot get page for bucket %d for cursor.\n",
		     cursor->item.bucket);
	    cursor->item.status = ITEM_ERROR;
	    return -1;
	  }
	}
      }
      else {
	/* Get new page */
	cursor->item.pagep = ffdb_get_page (hashp, nextp,
					    HASH_RAW_PAGE, 0, &tp);
	if (!cursor->item.pagep) {
	  fprintf (stderr, "Cannot get page for next cursor bucket %d\n",
		   cursor->item.bucket);
	  cursor->item.status = ITEM_ERROR;
	  return -1;
	}
      }
      cursor->item.pgno = tp;
      cursor->item.pgndx = 0;
      cursor->item.status = ITEM_OK;
    }
    else {
      /* Increase page index by one */
      cursor->item.pgndx++;
    }
  }
  else if (flags == FFDB_PREV) {
    /* We have reached the beginning of the key page */
    if (cursor->item.pgndx == 0) {
      /* Remember next page since we always go to the bucket page */
      nextp = NEXT_PGNO(cursor->item.pagep);
      /* put this page out */
      ffdb_put_page (hashp, cursor->item.pagep, TYPE(cursor->item.pagep), 0);
      if (cursor->item.bucket == 0 && nextp == INVALID_PGNO) {
	/* We are done */
	cursor->item.status = ITEM_NO_MORE;
	cursor->item.pagep = 0;
	return FFDB_NOT_FOUND;
      }
      if (nextp == INVALID_PGNO) {
	/* decrease bucket by one */
	cursor->item.bucket--;
	/* Get new page */
	cursor->item.pagep = ffdb_get_page (hashp, cursor->item.bucket,
					    HASH_BUCKET_PAGE, 0, &tp);


	/* Skip empty buckets */
	while (NUM_ENT(cursor->item.pagep) == 0) {
	  ffdb_put_page (hashp, cursor->item.pagep, 
			 TYPE(cursor->item.pagep), 0);
	  if (cursor->item.bucket == 0) {
	    /* this is the last one and it is empty */
	    cursor->item.status = ITEM_NO_MORE;
	    cursor->item.pagep = 0;
	    return FFDB_NOT_FOUND;
	  }

	  /* Get next bucket */
	  cursor->item.pagep = 0;
	  cursor->item.bucket--;

	  cursor->item.pagep = ffdb_get_page (hashp, cursor->item.bucket, 
					      HASH_BUCKET_PAGE, 0, &tp);
	  if (!(cursor->item.pagep)) {
	    fprintf (stderr, "Cannot get page for bucket %d for cursor.\n",
		     cursor->item.bucket);
	    cursor->item.status = ITEM_ERROR;
	    return -1;
	  }

	}
      }
      else {
	/* Get next page */
	cursor->item.pagep = ffdb_get_page (hashp, nextp,
					    HASH_RAW_PAGE, 0, &tp);
	if (!cursor->item.pagep) {
	  fprintf (stderr, "Cannot get page for prev cursor bucket %d\n",
		   cursor->item.bucket);
	  cursor->item.status = ITEM_ERROR;
	  return -1;
	}
      }
      cursor->item.pgno = tp;
      cursor->item.status = ITEM_OK;
      cursor->item.pgndx = NUM_ENT(cursor->item.pagep) - 1;
    }
    else {
      /* decrease page index by one */
      cursor->item.pgndx--;
    }
  }
  else {
    fprintf (stderr, "Unsupported cursor flag %d\n", flags);
    return -1;
  }

  /* Get Key data and size */
  ekdata = KEY(cursor->item.pagep, cursor->item.pgndx);
  eksize = KEY_LEN(cursor->item.pagep, cursor->item.pgndx);

  if (key->data && key->size > 0) {
    /* User supplied space */
    if (key->size < eksize) {
      fprintf (stderr, "Warning: application provided key space %d < key stored on disk with size %d\n",
	       key->size, eksize);
      return -1;
    }
    else
      key->size = eksize;
  }
  else {
    key->data = (unsigned char *)malloc(eksize * sizeof(unsigned char));
    key->size = eksize;
  }
  memcpy (key->data, ekdata, eksize);

  if (data) {
    /* Get data for this key */
    cursor->item.status = ITEM_OK;
    cursor->item.key_off = KEY_OFF(cursor->item.pagep, cursor->item.pgndx);
    cursor->item.key_len = eksize;
    cursor->item.data_off = DATAP_OFF(cursor->item.pagep, cursor->item.pgndx);

    return ffdb_get_item (hashp, key, data, &cursor->item, 0);
  }
  return 0;
}

/**
 * Dump out all page information for debug purpose
 */
void
ffdb_disp_all_page_info (ffdb_htab_t* hashp)
{
  unsigned int k;
  void *pagep;
  pgno_t page, tp;
  unsigned int numpages = hashp->mp->npages;

  fprintf (stderr, "There are total of %d pages \n", numpages);

  /**
   * Dump header page first
   */
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

  for (k = 1; k < numpages; k++) {
    page = k;
    pagep = ffdb_get_page (hashp, page,
			   HASH_RAW_PAGE, 0, &tp);
    if (pagep) {
      fprintf (stderr, "Page Information: PGNO %d PREV PGNO %d NEXT PGNO %d PAGE SIGN 0x%x ",
	       CURR_PGNO(pagep), PREV_PGNO(pagep),NEXT_PGNO(pagep), PAGE_SIGN(pagep));
      fprintf (stderr, "Page TYPE 0x%x OFFSET %d\n", TYPE(pagep), OFFSET(pagep));   
      ffdb_put_page (hashp, pagep, TYPE(pagep), 0);
    }
    else
      fprintf (stderr, "cannot get page %d\n", page);
  }
}

