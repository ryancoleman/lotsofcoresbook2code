#ifndef _DSLASH_SCALAR_COMPLETE_SPECIALIZATIONS_H_
#define _DSLASH_SCALAR_COMPLETE_SPECIALIZATIONS_H_

/* Disgusting hack to get rid of _mm_prefetches left in the generated scalar code */
#define _mm_prefetch(a,b) {}


#include "qphix/geometry.h"
#define QUOTEME(M)       #M
#define INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#define INCLUDE_FILE_VAR(PRE,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)


// No SOALEN, COMPRESS12, COMPRESS_SUFFIX -> generic template
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX
#include "qphix/scalar/dslash_scalar_complete_specialization_form.h"

/* SINGLE PRECISION */
#define FPTYPE float
#define VEC 1

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 1
#include "qphix/scalar/dslash_scalar_complete_specialization_form.h"
#undef SOA


#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 1
#include "qphix/scalar/dslash_scalar_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX
#undef VEC
#undef FPTYPE
/* --- END OF SINGLE PRECISION */


/* DOUBLE PRECISION */
#define FPTYPE double
#define VEC 1

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 1
#include "qphix/scalar/dslash_scalar_complete_specialization_form.h"
#undef SOA


#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 1
#include "qphix/scalar/dslash_scalar_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX
#undef VEC
#undef FPTYPE
/* --- END OF DOUBLE PRECISION */

#undef QUOTEME
#undef INCLUDE_FILE
#undef INCLUDE_FILE_VAR 

#undef _mm_prefetch
#endif
