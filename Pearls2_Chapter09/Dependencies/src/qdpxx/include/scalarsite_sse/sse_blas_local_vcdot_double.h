// $Id: sse_blas_local_vcdot_double.h,v 1.2 2008-06-18 16:03:02 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_LOCAL_VCDOT_DOUBLE
#define QDP_SSE_BLAS_LOCAL_VCDOT_DOUBLE

#include "qdp_precision.h"

namespace QDP {


  void local_vcdot4(REAL64 *sum, REAL64 *y, REAL64* x,int n_4spin);



} // namespace QDP;

#endif // guard
