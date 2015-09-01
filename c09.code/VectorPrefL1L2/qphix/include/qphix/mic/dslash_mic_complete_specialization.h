#ifndef _DSLASH_MIC_COMPLETE_SPECIALIZATIONS_H_
#define _DSLASH_MIC_COMPLETE_SPECIALIZATIONS_H_

#include "qphix/geometry.h"


#define QUOTEME(M)       #M
#define INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#define INCLUDE_FILE_VAR(PRE,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)


#undef FPTYPE
#undef VEC
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX
// No SOALEN, COMPRESS12, COMPRESS_SUFFIX -> generic template
#include "qphix/mic/dslash_mic_complete_specialization_form.h"




/* --------------- SINGLE PRECISION ------------------- */
#define FPTYPE float
#define VEC 16

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------- END OF SINGLE PRECISION --------------- */


/* --------------- DOUBLE PRECISION ------------------- */
#define FPTYPE double
#define VEC 8

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------- END OF DOUBLE PRECISION --------------- */

/* --------------- HALF PRECISION ------------------- */
/* HALF PRECISION IS FUNNY. SAME VECLEN/SOALEN AS SINGLE 
 * BUT THE DATA IS KEPT AS a 16 bit type. I WILL USE
 * SHORT HERE
 */

#define FPTYPE half
#define VEC 16

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/dslash_mic_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------- END OF HALF PRECISION --------------- */

#undef QUOTEME
#undef INCLUDE_FILE
#undef INCLUDE_FILE_VAR 


#endif


