/* FILE LICENSE TAG: SAMPLE */

#include <mkl.h>
#include <hStreams_common.h>


HSTR_RESULT
hStreams_custom_dtrsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                      const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                      const CBLAS_DIAG Diag, const MKL_INT M, const MKL_INT N,
                      const double alpha, const double *A, const MKL_INT lda,
                      double *B, const MKL_INT ldb, const uint32_t in_LogicalStr,
                      HSTR_EVENT *out_pEvent);

HSTR_RESULT
hStreams_custom_dpotrf(int matrix_order, char uplo, lapack_int n,
                       double *a, lapack_int lda, const uint32_t in_LogicalStr,
                       HSTR_EVENT *out_pEvent);

HSTR_RESULT
hStreams_custom_dsyrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                      const CBLAS_TRANSPOSE Trans, const MKL_INT N, const MKL_INT K,
                      const double alpha, const double *A, const MKL_INT lda,
                      const double beta, double *C, const MKL_INT ldc,
                      const uint32_t in_LogicalStr, HSTR_EVENT *out_pEvent);

HSTR_RESULT
hStreams_custom_dgetrf_square_nopiv_rowMaj(int size, double *mat,
        const uint32_t in_LogicalStr, HSTR_EVENT *out_pEvent);

HSTR_RESULT
hStreams_custom_dgetrf_square_nopiv_colMaj(int size, double *mat,
        const uint32_t in_LogicalStr, HSTR_EVENT *out_pEvent);
