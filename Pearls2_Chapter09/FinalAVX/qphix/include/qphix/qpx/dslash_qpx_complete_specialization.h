#ifndef QPHIX_DSLASH_QPX_COMPLETE_SPECIALIZATIONS_H
#define QPHIX_DSLASH_QPX_COMPLETE_SPECIALIZATIONS_H

#include "qphix/geometry.h"
#include "qphix/qpx/qpx_utils.h"

#define QUOTEME(M)       #M
#define INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#define INCLUDE_FILE_VAR(PRE,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)


// No SOALEN, COMPRESS12, COMPRESS_SUFFIX -> generic template
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"

// No single precision for now
#if 0
/* SINGLE PRECISION */
#define FPTYPE float
#define VEC 8

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA


#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX
#undef VEC
#undef FPTYPE
/* --- END OF SINGLE PRECISION */
#endif

/* DOUBLE PRECISION */
#define FPTYPE double
#define VEC 4

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

// No SOA 2 for now
#if 0
#define SOA 2
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA
#endif
#define SOA 4
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA


#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

// No SOA 2 for now
#if 0
#define SOA 2
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA
#endif

#define SOA 4
#include "qphix/qpx/dslash_qpx_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX
#undef VEC
#undef FPTYPE
/* --- END OF DOUBLE PRECISION */

#undef QUOTEME
#undef INCLUDE_FILE
#undef INCLUDE_FILE_VAR 


#endif


