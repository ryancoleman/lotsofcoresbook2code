#ifndef QPHIX_CLOV_DSLASH_AVX_COMPLETE_SPECIALIZATIONS_H
#define QPHIX_CLOV_DSLASH_AVX_COMPLETE_SPECIALIZATIONS_H

#include "qphix/geometry.h"
#include "qphix/avx/avx_utils.h"
#define QUOTEME(M)       #M
#define INCLUDE_FILE(PRE,PRE2,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## _ ## FPTYPE ## _ ## PRE2 ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#define INCLUDE_FILE_VAR(PRE,PRE2,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE(PRE,PRE2,FPTYPE,VLEN,SLEN,POST)


/* No SOALEN, COMPRESS12 COMPRESS_SUFFIX defined so this will include the generic template definitions */
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"

/* --------  SINGLE PRECISION  ----------- */
#define FPTYPE float
#define VEC 8

/* Uncompressed */
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */

#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------------- END OF SINGLE PRECISION ----------- */

/* --------  DOUBLE PRECISION  ----------- */
#define FPTYPE double
#define VEC 4

/* Uncompressed */
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 2
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA

#define SOA 4
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */

#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 2
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA

#define SOA 4
#include "qphix/avx/clov_dslash_avx_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------------- END OF DOUBLE PRECISION ----------- */

#undef QUOTEME
#undef INCLUDE_FILE
#undef INCLUDE_FILE_VAR 

#endif
