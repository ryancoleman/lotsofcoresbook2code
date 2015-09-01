// -*- C++ -*-

/*! @file
 * @brief Intel SSE optimizations
 *
 * SSE optimizations of basic operations
 */

#ifndef QDP_SCALARSITE_SSE_H
#define QDP_SCALARSITE_SSE_H


#warning "Using SSE BLAS. If your compiler cant handle intrinsics your build will break"
#include "scalarsite_sse/qdp_scalarsite_sse_linalg.h"
#include "scalarsite_sse/qdp_scalarsite_sse_blas.h"
#include "scalarsite_sse/qdp_scalarsite_sse_blas_double.h"
#include "scalarsite_sse/qdp_scalarsite_sse_linalg_double.h"
#include "scalarsite_generic/qdp_scalarsite_generic_cblas.h"

//#if BASE_PRECISION == 32
//#include "scalarsite_sse/sse_spin_aggregate.h"
//#else
//#include "scalarsite_sse/sse_spin_aggregate.h"
//#include "scalarsite_generic/generic_spin_aggregate.h"
//#endif

#endif  // guard

