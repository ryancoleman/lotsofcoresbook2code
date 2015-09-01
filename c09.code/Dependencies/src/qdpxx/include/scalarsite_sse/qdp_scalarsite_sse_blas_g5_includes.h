#ifndef QDP_SCALARSITE_SSE_BLAS_G5_INCLDUES_H
#define QDP_SCALARSITE_SSE_BLAS_G5_INCLUDES_H

#include "qdp_config.h"

#ifndef QDP_USE_SSE
#include "scalarsite_generic/generic_blas_vaxpy3_g5.h"
#include "scalarsite_generic/generic_blas_vaypx3_g5.h"
#include "scalarsite_generic/generic_blas_vadd3_g5.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
#include "scalarsite_generic/generic_blas_vaxpby3_g5.h"
#include "scalarsite_generic/generig_blas_g5.h"
#else

#if BASE_PRECISION == 32

#if defined (__GNUC__)
#include "scalarsite_sse/sse_blas_vaxpy3_g5.h"
#include "scalarsite_sse/sse_blas_vaypx3_g5.h"
#include "scalarsite_sse/sse_blas_vadd3_g5.h"
#include "scalarsite_sse/sse_blas_vscal3_g5.h"
#include "scalarsite_sse/sse_blas_vaxpby3_g5.h"
#include "scalarsite_generic/generic_blas_g5.h"
#else
#include "scalarsite_generic/generic_blas_vaxpy3_g5.h"
#include "scalarsite_generic/generic_blas_vaypx3_g5.h"
#include "scalarsite_generic/generic_blas_vadd3_g5.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
#include "scalarsite_generic/generic_blas_vaxpby3_g5.h"
#include "scalarsite_generic/generic_blas_g5.h"
#endif // GNUC

#else
#include "scalarsite_generic/generic_blas_vaxpy3_g5.h"
#include "scalarsite_generic/generic_blas_vaypx3_g5.h"
#include "scalarsite_generic/generic_blas_vadd3_g5.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
#include "scalarsite_generic/generic_blas_vaxpby3_g5.h"
#include "scalarsite_generic/generic_blas_g5.h"
#endif // BASE_PRECISION

#endif // QDP_USE_SSE

#endif
