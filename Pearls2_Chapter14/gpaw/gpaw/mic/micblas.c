/******************************************************************************
 * Copyright 2014 Intel Corporation All Rights Reserved.                      *
 *                                                                            *
 * The source code, information and material ("Material") contained herein    *
 * is owned by Intel Corporation or its suppliers or licensors, and title to  *
 * such Material remains with Intel Corporation or its suppliers or           *
 * licensors.                                                                 *
 * The Material contains proprietary information of Intel or its suppliers    *
 * and licensors. The Material is protected by worldwide copyright laws and   *
 * treaty provisions. No part of the Material may be used, copied,            *
 * reproduced, modified, published, uploaded, posted, transmitted,            *
 * distributed or disclosed in any way without Intel's prior express written  *
 * permission. No license under any patent, copyright or other intellectual   *
 * property rights in the material is granted to or conferred upon you,       *
 * either expressly, by implication, inducement, estoppel or otherwise. Any   *
 * license under such intellectual property rights must be express and        *
 * approved by Intel in writing.                                              *
 ******************************************************************************/

#include <pymic_kernel.h>

#include <stdio.h>
#include <complex.h>

#include <mkl.h>

#define DTYPE_INT       0
#define DTYPE_FLOAT     1
#define DTYPE_COMPLEX   2
 

 PYMIC_KERNEL
void mic_gemm(const int64_t *dtype, const void *A_, const void *B_,
              void *C_, const int64_t *m_, const int64_t *n_, const int64_t *k_,
              const int64_t *lda_, const int64_t *ldb_, const int64_t *ldc_,
              const void *alpha_, const void *beta_, const int64_t *trans) {
    int i;
#if CHECK_KERNEL        
    double sum;
#endif       

    int m = (int) *m_;
    int n = (int) *n_;
    int k = (int) *k_;
    int lda = (int) *lda_;
    int ldb = (int) *ldb_;
    int ldc = (int) *ldc_;


    switch(*dtype) {
    case DTYPE_FLOAT:
        {
            const double *A = (const double *) A_;
            const double *B = (const double *) B_;
            double *C = (double *) C_;
            const double *alpha = (const double *) alpha_;
            const double *beta = (const double *) beta_;

            CBLAS_TRANSPOSE ctrans = (*trans == 0 ? CblasNoTrans : CblasTrans);
            
            cblas_dgemm(CblasColMajor, ctrans, CblasNoTrans,
                        m, n, k, *alpha, A, lda, B, ldb, *beta, C, ldc);
        }
        break;
    case DTYPE_COMPLEX:
        {
            const double complex *A = (const double complex *) A_;
            const double complex *B = (const double complex *) B_;
            double complex *C = (double complex *) C_;
            const double complex *alpha = (const double complex *) alpha_;
            const double complex *beta = (const double complex *) beta_;

            CBLAS_TRANSPOSE ctrans = (*trans == 0 ? CblasNoTrans : CblasConjTrans);
            
            cblas_zgemm(CblasColMajor, ctrans, CblasNoTrans,
                        m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        }
        break;
    }
}


PYMIC_KERNEL
void mic_syrk(const int64_t *dtype, const void *A_, void *C_, 
              const int64_t *n, const int64_t *k, const int64_t *ldc,
              const void *alpha_, const void *beta_) {
    switch(*dtype) {
    case DTYPE_FLOAT:
        {
            const double *A = (const double *) A_;
            double *C = (double *) C_;
            const double *alpha = (const double *) alpha_;
            const double *beta = (const double *) beta_;

            cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 
                        *n, *k, *alpha, A, *k, *beta, C, *ldc);
        }
        break;
    case DTYPE_COMPLEX:
        {
            const double complex *A = (const double complex *) A_;
            double complex *C = (double complex *) C_;
            const double complex *alpha = (const double complex *) alpha_;
            const double complex *beta = (const double complex *) beta_;

            cblas_zsyrk(CblasColMajor, CblasUpper, CblasTrans, 
                        *n, *k, alpha, A, *k, beta, C, *ldc);
        }
        break;
    }
}


void mic_dsyr2k(const double *A, const double *B, double *C,
                const int64_t *n, const int64_t *k, const int64_t *ldc,
                double *beta, double *alpha) {
        cblas_dsyr2k(CblasColMajor, CblasUpper, CblasTrans, 
                     *n, *k, *alpha, A, *k, B, *k, *beta, C, *ldc);
}
