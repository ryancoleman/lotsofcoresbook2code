#ifndef QIO_STDINT_H
#define QIO_STDINT_H

#include "qio_config.h"


#if ( HAVE_STDINT_H )
#include <stdint.h>

#ifndef HAVE_UINT16_T
#error "stdint.h found but uint16_t undefined"
#endif

#ifndef HAVE_UINT32_T
#error "stdint.h found but uint32_t undefined"
#endif

#ifndef HAVE_UINT64_T
#error "stdint.h found but uint64_t undefined"
#endif


#elif ( HAVE_INTTYPES_H )
#include <inttypes.h>
  
#ifndef HAVE_UINT16_T
#error "inttypes.h found but uint16_t undefined"
#endif

#ifndef HAVE_UINT32_T
#error "inttypes.h found but uint32_t undefined"
#endif

#ifndef HAVE_UINT64_T
#error "inttypes.h found but uint64_t undefined"
#endif

#else

/* stdint not found. Try other ways */

/* n_uint_16_t */
#if ( SIZEOF_UNSIGNED_CHAR == 2 )
typedef unsigned char uint16_t;
#elif ( SIZEOF_UNSIGNED_SHORT == 2 )
typedef unsigned short uint16_t;
#elif ( SIZEOF_UNSIGNED_INT == 2 )
typedef unsigned int   uint16_t;
#elif ( SIZEOF_UNSIGNED_LONG == 2 ) 
typedef unsigned long  uint16_t;
#elif ( SIZEOF_UNSIGNED_LONG_LONG == 2 ) 
typedef unsigned long long uint16_t;
#else
#error "Could not find unsigned type of length 2 bytes = 16 bits"
#endif

/* n_uint_32_t */
#if ( SIZEOF_UNSIGNED_CHAR == 4 )
typedef unsigned char uint32_t;
#elif ( SIZEOF_UNSIGNED_SHORT == 4 )
typedef unsigned short uint32_t;
#elif ( SIZEOF_UNSIGNED_INT == 4 )
typedef unsigned int   uint32_t;
#elif ( SIZEOF_UNSIGNED_LONG == 4 ) 
typedef unsigned long  uint32_t;
#elif ( SIZEOF_UNSIGNED_LONG_LONG == 4 ) 
typedef unsigned long long uint32_t;
#else
#error "Could not find unsigned type of length 4 bytes = 32 bits"
#endif

/* n_uint_64_t */
#if ( SIZEOF_UNSIGNED_CHAR == 8 )
typedef unsigned char uint64_t;
#elif ( SIZEOF_UNSIGNED_SHORT == 8 )
typedef unsigned short uint64_t;
#elif ( SIZEOF_UNSIGNED_INT == 8 )
typedef unsigned int   uint64_t;
#elif ( SIZEOF_UNSIGNED_LONG == 8 ) 
typedef unsigned long  uint64_t;
#elif ( SIZEOF_UNSIGNED_LONG_LONG == 8 ) 
typedef unsigned long long uint64_t;
#else
#error "Could not find unsigned type of length 4 bytes = 32 bits"
#endif

#endif /* ifdef ( HAVE_STDINT_H || HAVE_INTTYPES_H ) */


#endif
