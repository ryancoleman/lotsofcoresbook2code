// $Id: sse_blas_vaypx4_double.h,v 1.1 2008-06-20 20:03:01 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VAYPX4_DOUBLE
#define QDP_SSE_BLAS_VAYPX4_DOUBLE

#include "qdp_precision.h"
namespace QDP {

  void vaypx4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, int n_4spin);

} // namespace QDP;

#endif // guard
