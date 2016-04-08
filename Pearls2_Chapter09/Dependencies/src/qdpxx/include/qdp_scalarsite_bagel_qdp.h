// -*- C++ -*-

/*! @file
 * @brief Qcdoc optimizations
 *
 * Qcdoc version of optimized basic operations
 */

#ifndef QDP_SCALARSITE_BAGEL_QDP_H
#define QDP_SCALARSITE_BAGEL_QDP_H

#include "qdp_config.h"

#include "scalarsite_bagel_qdp/qdp_scalarsite_bagel_qdp_linalg.h"

// Use QCDOC specific BLAS for now -- use Pete's assembler
#include "scalarsite_bagel_qdp/qdp_scalarsite_bagel_qdp_blas.h"
// Use GENERIC Chiral Projector BLAS for now, BAGEL should generate
// this eventually

#include "scalarsite_bagel_qdp/qdp_scalarsite_bagel_qdp_blas_g5.h"

// Use GENERIC Complex BLAS for now as there is no other yet
#include "scalarsite_generic/qdp_scalarsite_generic_cblas.h"

// #include "scalarsite_generic/generic_spin_aggregate.h"
#endif
