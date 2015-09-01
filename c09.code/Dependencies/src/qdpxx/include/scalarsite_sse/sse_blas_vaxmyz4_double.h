// $Id: sse_blas_vaxmyz4_double.h,v 1.3 2008-06-23 14:19:43 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VAXMYZ4_DOUBLE
#define QDP_SSE_BLAS_VAXMYZ4_DOUBLE

#include "qdp_precision.h"
namespace QDP {

  void vaxmyz4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, REAL64 *Add,int n_4vec);
  void vaxmy4(REAL64 *Out,REAL64 *scalep,REAL64 *InScale, int n_4spin);

} // namespace QDP;

#endif // guard
