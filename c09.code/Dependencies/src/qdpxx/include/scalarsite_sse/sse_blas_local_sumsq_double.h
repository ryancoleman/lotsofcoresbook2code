// $Id: sse_blas_local_sumsq_double.h,v 1.2 2008-06-18 16:03:02 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_LOCAL_SUMSQ_DOUBLE
#define QDP_SSE_BLAS_LOCAL_SUMSQ_DOUBLE

#include "qdp_precision.h"

namespace QDP {

#include <xmmintrin.h>

  void local_sumsq4(REAL64 *sum, REAL64 *vecptr, int n_4spin);


} // namespace QDP;

#endif // guard
