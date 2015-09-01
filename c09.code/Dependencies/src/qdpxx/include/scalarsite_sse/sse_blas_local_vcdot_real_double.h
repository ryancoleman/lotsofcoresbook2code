// $Id: sse_blas_local_vcdot_real_double.h,v 1.2 2008-06-18 16:03:02 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_LOCAL_VCDOT_REAL_DOUBLE
#define QDP_SSE_BLAS_LOCAL_VCDOT_REAL_DOUBLE

#include "qdp_precision.h"

namespace QDP {

  void local_vcdot_real4(REAL64 *sum, REAL64 *y, REAL64* x,int n_4spin);



} // namespace QDP;

#endif // guard
