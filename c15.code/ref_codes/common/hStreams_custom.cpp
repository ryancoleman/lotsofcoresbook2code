/* FILE LICENSE TAG: SAMPLE */

#include "hStreams_custom.h"
#include <hStreams_source.h>
#include <hStreams_common.h>
#include <intel-coi/common/COIMacros_common.h>
#include <mkl.h>

#if COI_LIBRARY_VERSION >= 2
#define COI_LIB_FLAGS_HACK NULL,
#else
#define COI_LIB_FLAGS_HACK /* nothing */
#endif


HSTR_RESULT
hStreams_custom_dtrsm(const CBLAS_ORDER Order, const CBLAS_SIDE Side,
                      const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                      const CBLAS_DIAG Diag, const MKL_INT M, const MKL_INT N,
                      const double alpha, const double *A, const MKL_INT lda,
                      double *B, const MKL_INT ldb, const uint32_t in_LogicalStr,
                      HSTR_EVENT *out_pEvent)
{
    doubleToUint64_t uAlpha;

    uAlpha.Set(alpha);

    uint64_t args[12];
    // Set up scalar args, then heap args
    args[ 0] = (uint64_t)(Order);
    args[ 1] = (uint64_t)(Side);
    args[ 2] = (uint64_t)(Uplo);
    args[ 3] = (uint64_t)(TransA);
    args[ 4] = (uint64_t)(Diag);
    args[ 5] = (uint64_t)(M);
    args[ 6] = (uint64_t)(N);
    args[ 7] = (uint64_t)(uAlpha.Get_uint64_t());
    args[ 8] = (uint64_t)(lda);
    args[ 9] = (uint64_t)(ldb);
    args[10] = (uint64_t)(A);
    args[11] = (uint64_t)(B);


    CHECK_HSTR_RESULT(
        hStreams_EnqueueCompute(
            in_LogicalStr,
            "hStreams_custom_dtrsm_sink",
            10,            // scalar args
            2,             // heap args
            args,          // arg array
            out_pEvent,    // event
            NULL,          // return value pointer
            0));           // return value size

    return HSTR_RESULT_SUCCESS;
}

HSTR_RESULT
hStreams_custom_dpotrf(int matrix_order, char uplo, lapack_int n,
                       double *a, lapack_int lda, const uint32_t in_LogicalStr,
                       HSTR_EVENT *out_pEvent)
{
    charToUint64_t uUplo;

    uUplo.Set(uplo);

    uint64_t args[5];
    // Set up scalar args, then heap args
    args[ 0] = (uint64_t)(matrix_order);
    args[ 1] = (uint64_t)(uUplo.Get_uint64_t());
    args[ 2] = (uint64_t)(n);
    args[ 3] = (uint64_t)(lda);
    args[ 4] = (uint64_t)(a);

    CHECK_HSTR_RESULT(
        hStreams_EnqueueCompute(
            in_LogicalStr,
            "hStreams_custom_dpotrf_sink",
            4,            // scalar args
            1,             // heap args
            args,          // arg array
            out_pEvent,    // event
            NULL,          // return value pointer
            0));           // return value size

    return HSTR_RESULT_SUCCESS;
}

HSTR_RESULT
hStreams_custom_dsyrk(const CBLAS_ORDER Order, const CBLAS_UPLO Uplo,
                      const CBLAS_TRANSPOSE Trans, const MKL_INT N, const MKL_INT K,
                      const double alpha, const double *A, const MKL_INT lda,
                      const double beta, double *C, const MKL_INT ldc,
                      const uint32_t in_LogicalStr, HSTR_EVENT *out_pEvent)
{
    doubleToUint64_t uAlpha, uBeta;

    uAlpha.Set(alpha);
    uBeta.Set(beta);

    uint64_t args[11];
    // Set up scalar args, then heap args
    args[ 0] = (uint64_t)(Order);
    args[ 1] = (uint64_t)(Uplo);
    args[ 2] = (uint64_t)(Trans);
    args[ 3] = (uint64_t)(N);
    args[ 4] = (uint64_t)(K);
    args[ 5] = (uint64_t)(uAlpha.Get_uint64_t());
    args[ 6] = (uint64_t)(lda);
    args[ 7] = (uint64_t)(uBeta.Get_uint64_t());
    args[ 8] = (uint64_t)(ldc);
    args[ 9] = (uint64_t)(A);
    args[10] = (uint64_t)(C);


    CHECK_HSTR_RESULT(
        hStreams_EnqueueCompute(
            in_LogicalStr,
            "hStreams_custom_dsyrk_sink",
            9,            // scalar args
            2,             // heap args
            args,          // arg array
            out_pEvent,    // event
            NULL,          // return value pointer
            0));           // return value size
    return HSTR_RESULT_SUCCESS;

}

HSTR_RESULT
hStreams_custom_dgetrf_square_nopiv_rowMaj(int size, double *mat,
        const uint32_t in_LogicalStr, HSTR_EVENT *out_pEvent)
{
    uint64_t args[2];
    // Set up scalar args, then heap args
    args[ 0] = (uint64_t)(size);
    args[ 1] = (uint64_t)(mat);

    CHECK_HSTR_RESULT(
        hStreams_EnqueueCompute(
            in_LogicalStr,
            "hStreams_custom_dgetrf_square_nopiv_rowMaj_sink",
            1,            // scalar args
            1,             // heap args
            args,          // arg array
            out_pEvent,    // event
            NULL,          // return value pointer
            0));           // return value size

    return HSTR_RESULT_SUCCESS;
}

HSTR_RESULT
hStreams_custom_dgetrf_square_nopiv_colMaj(int size, double *mat,
        const uint32_t in_LogicalStr, HSTR_EVENT *out_pEvent)
{
    uint64_t args[2];
    // Set up scalar args, then heap args
    args[ 0] = (uint64_t)(size);
    args[ 1] = (uint64_t)(mat);

    CHECK_HSTR_RESULT(
        hStreams_EnqueueCompute(
            in_LogicalStr,
            "hStreams_custom_dgetrf_square_nopiv_colMaj_sink",
            1,            // scalar args
            1,             // heap args
            args,          // arg array
            out_pEvent,    // event
            NULL,          // return value pointer
            0));           // return value size

    return HSTR_RESULT_SUCCESS;
}
