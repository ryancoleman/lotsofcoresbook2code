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
 *     $Log: ffdb_page.h,v $
 *     Revision 1.1  2009-02-20 20:44:47  chen
 *     initial import
 *
 *
 *
 */
#ifndef _FFDB_PAGE_H
#define _FFDB_PAGE_H

/**
 * Invalid page number
 */
#define INVALID_PGNO 0xffffffff

/**
 * A single data is pointed by the following structure
 * On the page where the data is stored, the data len is the first
 * 4 bytes stored.
 */
typedef struct _ffdb_datap_
{
  pgno_t first;          /* First page where the data resides */
  pgno_t offset;         /* Offset within the page            */
  pgno_t len;            /* total length of data              */
  unsigned int chksum;   /* crc checksum                      */
}ffdb_datap_t;


/**
 * Alignment for the data pointer value
 */
#define DATAP_ALIGNMENT (unsigned int)4
#define ALIGN_DATAP_OFFSET_VAL(val)                                 \
  do {                                                              \
    unsigned int mask = DATAP_ALIGNMENT - 1;                        \
    if (val & mask)                                                 \
      val = (val) & ~mask;					    \
  }while(0)

/**
 * Every page has this signature so that a uninitialized page
 * can be distinguished
 */
#define FFDB_PAGE_MAGIC (unsigned int)0xe0f1a2b9

/**
 * Regular page layout format (overflow page has the same format)
 *
 * Hash pages store meta-data beginning at the top of the page (offset 0) 
 * and key/data values beginning at the bottom of the page (offset pagesize).
 * Fields are always accessed via macros so that we can change the page
 * format without too much pain. The reason to keep the length of the key
 * is the data pointer has to be aligned on 4 byte boundary. Therefore there
 * may be a gap between the data and the key. In the accessor
 * macros below, P means a pointer to the page, I means an index of the
 * particular entry being accessed. 
 *
 * Key and Data pair are from the end of a page arranged like
 *
 *       ............   Data2   Key2    Data1   Key1
 *
 * Hash base page format
 * BYTE ITEM                    NBYTES  TYPE            ACCESSOR MACRO
 * ---- ------------------      ------  --------        --------------
 * 0    current page number     4       pgno_t          CURR_PGNO(p)
 * 4    previous page number    4       pgno_t          PREV_PGNO(P)
 * 8    next page number        4       pgno_t          NEXT_PGNO(P)
 * 12   page signature          4       pgno_t          PAGE_SIGN(P)
 * 16   # pairs on page         2       indx_t          NUM_ENT(P)
 * 18   page type               2       indx_t          TYPE(P)
 * 20   check sum (crc)         4       pgno_t          CHKSUM(P)
 * 24   highest free byte       4       pgno_t          OFFSET(P)
 * 28   key offset 0            4       pgno_t          KEY_OFF(P, I)
 * 32   key len 0               4       pgno_t          KEY_LEN(P, I)
 * 36   data offset 0           4       pgno_t          DATA_OFF(P, I)
 * 40   key  offset 1           4       pgno_t          KEY_OFF(P, I)
 * 44   key  len 1              4       pgno_t          KEY_LEN(P, I)
 * 48   data offset 1           4       pgno_t          DATA_OFF(P, I)
 * ...etc...
 */
/* Indices (in bytes) of the beginning of each of these entries */
#define I_CURR_PGNO      0
#define I_PREV_PGNO	 4
#define I_NEXT_PGNO	 8
#define I_PAGE_SIGN      12
#define I_ENTRIES	 16
#define I_TYPE		 18
#define I_CHKSUM         20
#define I_HF_OFFSET	 24

/* Checksum location for this page */
#define CHKSUM_POS I_CHKSUM
#define CHKSUM_OFFSET I_HF_OFFSET

/* Overhead is everything prior to the first key/data pair. */
#define PAGE_OVERHEAD	(I_HF_OFFSET + sizeof(pgno_t))

/* overhead for one key and one data: two offsets + one key len */
/* Here is the reason limiting page size: the offset is a unsigned short */
#define PAIR_OVERHEAD (3*sizeof(pgno_t))

/* macro to retrive a value of Type Y from page P at offset F */
#define FIND_VALUE(P, T, F) (((T *)((unsigned char *)(P) + F))[0])

/**
 * The following macros to access different field on a page
 */
#define NUM_ENT(P)      (FIND_VALUE((P), indx_t, I_ENTRIES))
#define CURR_PGNO(P)	(FIND_VALUE((P), pgno_t, I_CURR_PGNO))
#define PREV_PGNO(P)	(FIND_VALUE((P), pgno_t, I_PREV_PGNO))
#define NEXT_PGNO(P)	(FIND_VALUE((P), pgno_t, I_NEXT_PGNO))
#define PAGE_SIGN(P)	(FIND_VALUE((P), pgno_t, I_PAGE_SIGN))
#define TYPE(P)		(FIND_VALUE((P), indx_t, I_TYPE))
#define OFFSET(P)	(FIND_VALUE((P), pgno_t, I_HF_OFFSET))
#define CHKSUM(P)       (FIND_VALUE((P), pgno_t, I_CHKSUM))
#define HIGHEST_FREE(p) (OFFSET((p)))

/* we keep current page no on page */
#define ADDR(P)         (CURR_PGNO(P))

/* retrive key/data pair offsets and data for a given index */
/* value of data pointer offset has to be aligned to 4 byte boundary */
#define DATAP_OFF(P, N) \
  FIND_VALUE(P,pgno_t, PAGE_OVERHEAD + N * PAIR_OVERHEAD + 2*sizeof(pgno_t))
#define KEY_OFF(P, N) \
  FIND_VALUE(P,pgno_t, PAGE_OVERHEAD + N * PAIR_OVERHEAD)
#define KEY_LEN(P, N) \
  FIND_VALUE(P,pgno_t, PAGE_OVERHEAD + N * PAIR_OVERHEAD + sizeof(pgno_t))


/* Key value with index N on page P */
#define KEY(P, N)   (((unsigned char *)(P) + KEY_OFF((P), (N))))
/* Data pointer value with index N on page P */
#define DATAP(P, N) ((ffdb_datap_t *)(((unsigned char *)(P) + DATAP_OFF((P), (N)))))

/* calculate various sizes on a page */
/* This is over estimate: since we have to align data pointer on 4 byte boundary
 * therefore the worst case is sizeof(ffdb_datap_t) + sizeof(int) - 1
 */
#define PAIRSIZE(K,D) (PAIR_OVERHEAD + (K)->size + sizeof(ffdb_datap_t) + sizeof(int) - 1)

/*
 * Since these are all unsigned, we need to guarantee that we never go
 * negative.  Offset values are 0-based and overheads are one based (i.e.
 * one byte of overhead is 1, not 0), so we need to convert OFFSETs to
 * 1-based counting before subtraction.
 */
#define FREESPACE(P) \
  ((OFFSET((P)) + 1 - PAGE_OVERHEAD - (NUM_ENT((P)) * PAIR_OVERHEAD)))

/**
 * Whether a regular key and pair would fit
 */
#define PAIRFITS(P,K,D)	((PAIRSIZE((K),(D))) <= FREESPACE((P)))

/**
 * Page allocation definition
 */
#define A_BUCKET	1000
#define A_OVFL		1100
#define A_DATA          1100
#define A_RAW		4000

/**
 * Data address alignment on 4 byte boundary
 */
#define ADDR_ALIGNMENT (unsigned int)4
#define ALIGN_ADDR(val)                                             \
  do {                                                              \
    unsigned int mask = ADDR_ALIGNMENT - 1;                         \
    if (val & mask)                                                 \
      val = (val + mask) & ~mask;                                   \
  }while(0)


/**
 * big page (data page) format
 *
 * Currently the checksum is not used because the data is checksumed
 * Highest free byte is also 4 byte aligned
 * Highest free byte is the offset location of next free byte from 
 * the beginning of the page
 *
 * The number of data equals the number of data with headers on the page.
 * If a data item has its header in another page, this data item is not
 * counted in this number. But if a data item has its header in this page
 * and rest data in another page, it will be counted here.
 *
 * The first data offset is the offset location where the first data starts
 *
 * Each data item starts at 4 byte boundary
 * 0    current page number     4       pgno_t          CURR_PGNO(p)
 * 4    previous page number    4       pgno_t          PREV_PGNO(P)
 * 8    next page number        4       pgno_t          NEXT_PGNO(P)
 * 12   page signature          4       pgno_t          PAGE_SIGN(P)
 * 16   # data  on page         2       indx_t          NUM_ENT(P)
 * 18   page type               2       indx_t          TYPE(P)
`* 20   crc checksum            4       pgno_t          CKSUM(P)
 * 24   highest free byte       4       pgno_t          OFFSET(P)
 * 28   first data offset       4       pgno_t          FIRSTDATA(P)
 * 
 * ............................................................
 * free byte  data length       4       pgno_t
 *            data_status       4       deleted or not
 *            next              4       pgno_t
 *            key_page          4       pgno_t
 *            key_idx           4       pgno_t
 *      data
 *
 * When a page has been used, highest free byte = 0
 * When a single data span multiple pages, firstdata points to next
 * data on the page if the data can be fit in the page. Otherwise
 * firstdata points to 0.
 */

/**
 * This is the definition of the header of each data item to allow 
 * easy reversal lookup of the key
 */
#define DATA_VALID    0xffeeffee
#define DATA_INVALID  0x11221122

typedef struct _ffdb_data_header_
{
  pgno_t  len;                /* Length of this data item          */
  pgno_t  status;             /* data is deleted or not            */
  pgno_t  next;               /* next data item on this page       */
  pgno_t  key_page;           /* page number where the key resides */
  pgno_t  key_idx;            /* index within the key page to find key */
}ffdb_data_header_t;

#define I_FIRST_DATA_POS   28

#define FIRST_DATA_POS(P)	(FIND_VALUE((P), pgno_t, I_FIRST_DATA_POS))

/**
 * Big page overhead
 */
#define BIG_PAGE_OVERHEAD	(I_FIRST_DATA_POS + sizeof(pgno_t))

/**
 * Big data header
 */
#define BIG_DATA_OVERHEAD       (sizeof(ffdb_data_header_t))

/**
 * Total big data size including header information
 */
#define BIG_DATA_TOTAL_SIZE(val)(BIG_DATA_OVERHEAD + (val)->size)

/**
 * Get data header pointed by offset on a page
 */
#define BIG_DATA_HEADER(P,offset) ((ffdb_data_header_t *)(((unsigned char *)(P) + (offset))))

/**
 * A single key must be fit into a page
 * we are not dealing with the case that a single key cannot be fit
 * into a page
 */
#define FFDB_MAX_KEYSIZE(h) ((h->hdr.bsize) - PAGE_OVERHEAD - PAIROVERHEAD - sizeof(ffdb_datap_t))

#define PAIRFITS(P,K,D)	((PAIRSIZE((K),(D))) <= FREESPACE((P)))

/**
 * Free Map Page: This is the pages descrbing which pages are free.
 * The free pages are the result of pages being deleted (most likely
 * data pages.
 * The free map pages are in for each spliting level. If there are
 * no free pages, there will be no free map pages
 *
 * 0    current page number     4       pgno_t          CURR_PGNO(p)
 * 4    previous page number    4       pgno_t          PREV_PGNO(P)
 * 8    next page number        4       pgno_t          NEXT_PGNO(P)
 * 12   page signature          4       pgno_t          PAGE_SIGN(P)
 * 16   # data  on page         2       indx_t          NUM_ENT(P)
 * 18   page type               2       indx_t          TYPE(P)
`* 20   crc checksum            4       pgno_t          CKSUM(P)
 * 24   pad                     4       pgno_t          not used
 * 28   free page number1       4       pgno_t          
 * 32   free page number2       4       pgno_t          
 */

/**
 * Get free overflow/datapage number
 */
#define FREE_PAGE_START    28
#define NUM_FREE_PAGES(P) (NUM_ENT((P)))
#define FREE_PAGE(P,IDX) (FIND_VALUE((P), pgno_t, FREE_PAGE_START + (IDX)*sizeof(pgno_t)))
#define MAX_NUM_FREE_PAGES ((hashp->hdr.bsize - FREE_PAGE_START)/sizeof(pgno_t))

/**
 * Page Type Definition
 */
#define HASH_BUCKET_PAGE     0x1001
#define HASH_BIG_PAGE        0x2001
#define HASH_DATA_PAGE       HASH_BIG_PAGE
#define HASH_OVFL_PAGE       0x3001
#define HASH_UINFO_PAGE      0x4001
#define HASH_CONFIG_PAGE     0x5001
#define HASH_FREE_PAGE       0x6001
#define HASH_DELETED_PAGE    0x8001
#define HASH_RAW_PAGE        0xffee

/**
 * User data information page format
 *
 * ---- ------------------      ------  --------        --------------
 * 0    current page number     4       pgno_t          CURR_PGNO(p)
 * 4    previous page number    4       pgno_t          PREV_PGNO(P)
 * 8    next page number        4       pgno_t          NEXT_PGNO(P)
 * 12   page signature          4       pgno_t          PAGE_SIGN(P)
 * 16   # pairs on page         2       indx_t          NUM_ENT(P)
 * 18   page type               2       indx_t          TYPE(P)
 * 20   check sum (crc)         4       pgno_t          CHKSUM(P)
 * 24   highest free byte       4       pgno_t          OFFSET(P)
 * 28   user data len           4       pgno_t          USER_DATA_LEN(P)
 */
#define I_USER_DATA_LEN         28
#define USER_DATA_LEN(P) (FIND_VALUE((P), pgno_t, I_USER_DATA_LEN))

/**
 * Configuration information page format
 *
 * ---- ------------------      ------  --------        --------------
 * 0    current page number     4       pgno_t          CURR_PGNO(p)
 * 4    previous page number    4       pgno_t          PREV_PGNO(P)
 * 8    next page number        4       pgno_t          NEXT_PGNO(P)
 * 12   page signature          4       pgno_t          PAGE_SIGN(P)
 * 16   # pairs on page         2       indx_t          NUM_ENT(P)
 * 18   page type               2       indx_t          TYPE(P)
 * 20   check sum (crc)         4       pgno_t          CHKSUM(P)
 * 24   highest free byte       4       pgno_t          OFFSET(P)
 * 28   configuration data 
 */
#define CONFIG_INFO(P, N) ((ffdb_config_info_t *)((unsigned char *)(P) + PAGE_OVERHEAD + (N) * sizeof(ffdb_config_info_t)))

/**
 * Page in and out routines
 */
extern void ffdb_pgin_routine (void* arg, pgno_t pgno, void* page);
extern void ffdb_pgout_routine (void* arg, pgno_t pgno, void* page);



#endif
