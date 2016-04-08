#ifndef QPHIX_BLAS_H
#define QPHIX_BLAS_H

#include "qphix/qphix_config.h"
#include "qphix/comm.h"

// Generic OpenMP templated
#include "qphix/blas_c.h"
#include "qphix/blas_new_c.h"

// MIC Specializations
#ifdef QPHIX_MIC_SOURCE
#include "qphix/blas_mic.h"
#endif

// SSE Specializations
#ifdef QPHIX_SSE_SOURCE
#endif

// AVX Specialization
#ifdef QPHIX_AVX_SOURCE
#endif

#endif
