#ifndef LIME_BINARY_HEADER_H
#define LIME_BINARY_HEADER_H

#include "lime_fixed_types.h"

/** LIME Binary Header Format **/

/**< Header length in 64-bit words */
#define MAX_HDR64 18

static union { 
  n_uint64_t int64[MAX_HDR64];
  n_uint32_t int32[2*MAX_HDR64];
  n_uint16_t int16[4*MAX_HDR64];
  unsigned char uchr[8*MAX_HDR64];
} lime_header;

/* Total LIME header length in bytes */
#define LIME_HDR_LENGTH 8*MAX_HDR64

/* First 64-bit word 
 * int64:    0
 * int32:    0                1
 * Shorts:   0           1    2         3
 * Bytes:    0           3    4     5   6               7
 * Bits:     0          31    32   47   48  49    50 - 63
 * Meaning  |  magic number | version | MB | ME | Reserved |
 */

static n_uint32_t *lime_hdr_magic_no = &lime_header.int32[0];
static n_uint16_t *lime_hdr_version  = &lime_header.int16[2];
static unsigned char *lime_hdr_mbme  = &lime_header.uchr[6];

/**< Message Begin Mask (Internal) */
#define MB_MASK                   ((unsigned char)0x80)

/**< Message End Mask (Internal) */
#define ME_MASK                   ((unsigned char)0x40)

/* Second 64-bit word
 * int64:      1
 * Meaning   | data length in bytes |
 */

static n_uint64_t *lime_hdr_data_len = &lime_header.int64[1];

/* Next 16 64-bit words 
 * int64:    2                     17
 * Bytes:    16                        143
 * Meaning | string specifying record type |
 */

static unsigned char *lime_hdr_rec_type     = &lime_header.uchr[16];

/* Note: max 127 chars in string to permit null termination */
#define MAX_LIME_HDR_REC_TYPE 128  

#endif
