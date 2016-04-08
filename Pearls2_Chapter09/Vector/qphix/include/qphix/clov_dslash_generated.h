#ifndef QPHIX_CLOV_DSLASH_C_GENERATED
#define QPHIX_CLOV_DSLASH_C_GENERATED

#if defined (QPHIX_MIC_SOURCE)
#include "qphix/mic/clov_dslash_mic_complete_specialization.h"

#elif defined(QPHIX_AVX_SOURCE)
#include "qphix/avx/clov_dslash_avx_complete_specialization.h"

#elif defined(QPHIX_SCALAR_SOURCE)
#include "qphix/scalar/clov_dslash_scalar_complete_specialization.h"

#elif defined(QPHIX_QPX_SOURCE)
#include "qphix/qpx/clov_dslash_qpx_complete_specialization.h"

#endif

#endif
