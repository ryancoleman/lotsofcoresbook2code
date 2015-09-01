// $Id: sse_blas_vscal4_double.h,v 1.3 2008-06-18 16:03:02 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VSCAL_DOUBLE
#define QDP_SSE_BLAS_VSCAL_DOUBLE

#include "qdp_precision.h"

namespace QDP {

  void vscal4(REAL64 *z,REAL64 *a,REAL64 *x, int n_4spin);


} // namespace QDP;

#endif // guard
