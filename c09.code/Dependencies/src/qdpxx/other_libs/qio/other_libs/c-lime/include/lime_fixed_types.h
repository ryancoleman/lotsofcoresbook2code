#ifndef LIME_FIXED_TYPES_H
#define LIME_FIXED_TYPES_H

#include "lime_config.h"



#if ( HAVE_STDINT_H )
/* We have either got <inttypes.h> or <stdint.h> or both */
#include <stdint.h>

#ifdef HAVE_UINT16_T
typedef uint16_t n_uint16_t;
#else
#error "stdint.h found but uint16_t undefined"
#endif

#ifdef HAVE_UINT32_T
typedef uint32_t n_uint32_t;
#else
#error "stdint.h found but uint32_t undefined"
#endif

#ifdef HAVE_UINT64_T
typedef uint64_t n_uint64_t;
#else 
#error "stdint.h found but uint64_t undefined"
#endif


#else 
  
/* stdint not found. Try other ways */

/* n_uint_16_t */
#if ( SIZEOF_UNSIGNED_CHAR == 2 )
typedef unsigned char n_uint16_t;
#elif ( SIZEOF_UNSIGNED_SHORT == 2 )
typedef unsigned short n_uint16_t;
#elif ( SIZEOF_UNSIGNED_INT == 2 )
typedef unsigned int   n_uint16_t;
#elif ( SIZEOF_UNSIGNED_LONG == 2 ) 
typedef unsigned long  n_uint16_t;
#elif ( SIZEOF_UNSIGNED_LONG_LONG == 2 ) 
typedef unsigned long long n_uint16_t;
#else
#error "Could not find unsigned type of length 2 bytes = 16 bits"
#endif

/* n_uint_32_t */
#if ( SIZEOF_UNSIGNED_CHAR == 4 )
typedef unsigned char n_uint32_t;
#elif ( SIZEOF_UNSIGNED_SHORT == 4 )
typedef unsigned short n_uint32_t;
#elif ( SIZEOF_UNSIGNED_INT == 4 )
typedef unsigned int   n_uint32_t;
#elif ( SIZEOF_UNSIGNED_LONG == 4 ) 
typedef unsigned long  n_uint32_t;
#elif ( SIZEOF_UNSIGNED_LONG_LONG == 4 ) 
typedef unsigned long long n_uint32_t;
#else
#error "Could not find unsigned type of length 4 bytes = 32 bits"
#endif

/* n_uint_64_t */
#if ( SIZEOF_UNSIGNED_CHAR == 8 )
typedef unsigned char n_uint64_t;
#elif ( SIZEOF_UNSIGNED_SHORT == 8 )
typedef unsigned short n_uint64_t;
#elif ( SIZEOF_UNSIGNED_INT == 8 )
typedef unsigned int   n_uint64_t;
#elif ( SIZEOF_UNSIGNED_LONG == 8 ) 
typedef unsigned long  n_uint64_t;
#elif ( SIZEOF_UNSIGNED_LONG_LONG == 8 ) 
typedef unsigned long long n_uint64_t;
#else
#error "Could not find unsigned type of length 4 bytes = 32 bits"
#endif

#endif /* ifdef ( HAVE_STDINT_H ) */


#endif /* initial trigger guard */
