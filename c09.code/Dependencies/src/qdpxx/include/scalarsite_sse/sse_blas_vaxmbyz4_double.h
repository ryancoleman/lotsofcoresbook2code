// $Id: sse_blas_vaxmbyz4_double.h,v 1.3 2008-06-23 14:19:43 bjoo Exp $

/*! @file
 *  @brief Generic Scalar VAXPY routine
 *
 */

#ifndef QDP_SSE_BLAS_VAXMBYZ4_DOUBLE
#define QDP_SSE_BLAS_VAXMBYZ4_DOUBLE

#include "qdp_precision.h"
namespace QDP {

  void vaxmbyz4(REAL64 *z,REAL64 *a,REAL64 *x, REAL64 *b, REAL64 *y, int n_4vec);

  void vaxmby4(REAL64 *z,REAL64 *a,REAL64 *x, REAL64 *b, int n_4vec);


} // namespace QDP;

#endif // guard
