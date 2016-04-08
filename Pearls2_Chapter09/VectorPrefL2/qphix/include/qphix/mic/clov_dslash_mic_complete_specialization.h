#ifndef CLOV_DSLASH_MIC_COMPLETE_SPECIALIZATIONS_H
#define CLOV_DSLASH_MIC_COMPLETE_SPECIALIZATIONS_H

#include "qphix/geometry.h"
#define QUOTEME(M)       #M
#define INCLUDE_FILE(PRE,PRE2,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## _ ## FPTYPE ## _ ## PRE2 ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#define INCLUDE_FILE_VAR(PRE,PRE2,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE(PRE,PRE2,FPTYPE,VLEN,SLEN,POST)


/* No SOALEN, COMPRESS12 COMPRESS_SUFFIX defined so this will include the generic template definitions */
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"

/* --------  SINGLE PRECISION  ----------- */
#define FPTYPE float
#define VEC 16

/* Uncompressed */
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */

#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA


#define SOA 16
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------------- END OF SINGLE PRECISION ----------- */


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
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 16
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA
#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------- END OF HALF PRECISION --------------- */


/* --------------- DOUBLE PRECISION ------------------- */
#define FPTYPE double
#define VEC 8

// Uncompressed
#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

/* Compressed */
#define COMPRESS12 true
#define COMPRESS_SUFFIX _12

#define SOA 4
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/mic/clov_dslash_mic_complete_specialization_form.h"
#undef SOA

#undef COMPRESS12
#undef COMPRESS_SUFFIX

#undef VEC
#undef FPTYPE
/* -------------- END OF DOUBLE PRECISION --------------- */

#undef QUOTEME
#undef INCLUDE_FILE
#undef INCLUDE_FILE_VAR 

#endif
