// $Id: sse_blas_vaxpy4_double.h,v 1.2 2008-06-18 16:03:02 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VAXPY4_DOUBLE
#define QDP_SSE_BLAS_VAXPY4_DOUBLE

#include "qdp_precision.h"
namespace QDP {

  void vaxpy4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, int n_4spin);
  void vaxpyz4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, REAL64 *Add,int n_4vec);


} // namespace QDP;

#endif // guard
