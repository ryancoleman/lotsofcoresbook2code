/* FILE LICENSE TAG: SAMPLE */
#include <stdio.h>
#include <mkl.h>
#include <math.h>
#ifndef WIN32
#include <sys/time.h>
#endif
#include <omp.h>
#include <stdlib.h>
#include <cmath>
#include <string.h>
#include <hStreams_source.h>
#include <hStreams_app_api.h>
#include <intel-coi/common/COIMacros_common.h>
#include <omp.h>
#include "helper_functions.h"
#include "hStreams_custom.h"
#include "dtime.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define SWITCH_CHAR  '-'

extern double *dpo_generate(size_t);

HSTR_CPU_MASK avail_cpu_mask, avoid_cpu_mask;

static int loc_verbose = 0;

//custom dgetrf nopiv for square matrices
//This function does the lu factorization without pivoting for
//a col-major storage of input matrix
int custom_dgetrf_square_nopiv_colMaj(int size, double *mat)
{
    double invDiag;

    #pragma omp parallel
    {
        int i, j, k, ivec;
        double factor;
        int off1, off2;

        for (k = 0; k < size; ++k) {
            off1 = k * size;
            invDiag = 1 / mat[off1 + k];

            #pragma omp for simd
            for (i = k + 1; i < size; ++i)
                //get L entries below the diagonal
            {
                mat[off1 + i] = invDiag * mat[off1 + i];
            }

            //get Schur Complement
            #pragma omp for
            for (j = k + 1; j < size; ++j) {
                off2 = j * size;
                factor = mat[off2 + k];
#pragma simd
                for (i = k + 1; i < size; ++i) {
                    mat[off2 + i] -= factor * mat[off1 + i];
                }
            }
        }
    }
    return 0;
}

//custom dgetrf nopiv for square matrices
//This function does the lu factorization without pivoting for
//a row-major storage of input matrix
int custom_dgetrf_square_nopiv_rowMaj(int size, double *mat)
{
    double invDiag;

    #pragma omp parallel
    {
        int i, j, k, ivec;
        double factor;
        int off1, off2;

        for (k = 0; k < size; ++k) {
            off1 = k * size;
            invDiag = 1 / mat[off1 + k];

            #pragma omp for
            for (i = k + 1; i < size; ++i) {
                off2 = i * size;
                //get L entries below the diagonal
                mat[off2 + k] = invDiag * mat[off2 + k];

                //get Schur Complement
                factor =  mat[off2 + k];
#pragma simd
                for (j = k + 1; j < size; ++j) {
                    mat[off2 + j] -= factor * mat[off1 + j];
                }
            }
        }
    }
    return 0;
}


void lu_tiled(double *mat, int tile_size, int num_tiles, int mat_size,
              int niter, int max_log_str, bool layRow, int verify)
{
    //total number of tiles
    int tot_tiles = num_tiles * num_tiles;

    //memory allocation for matrix for tiled-lu
    double *A_my = (double *)malloc(mat_size * mat_size * sizeof(double));

    //memory allocation for matrix for MKL lu (for comparison)
    double *A_MKL = (double *)malloc(mat_size * mat_size * sizeof(double));

    //used for verification of results
    double *Cout = (double *)malloc(mat_size * mat_size * sizeof(double));
    double *eye = (double *)malloc(mat_size * mat_size * sizeof(double));

    //memory allocation for tiled matrix
    double **Asplit = new double* [tot_tiles];
    int mem_size_tile = tile_size * tile_size * sizeof(double);

    for (int i = 0; i < tot_tiles; ++i) {
        //Buffer per tile, host allocation
        Asplit[i] = (double *)_mm_malloc(mem_size_tile, 64);

        //Buffer creation and allocation on the card
        hStreams_app_create_buf((void *)Asplit[i], mem_size_tile);
    }

    double tbegin, tend;

    //initializing the identity matrix
    //used for verification of results
    eyeInit(eye, mat_size);

    int iter;
    int info;

    //Events are needed for various synchronizations to enforce
    //data dependence between and among data-transfers/computes
    HSTR_EVENT *eventcpyto = new HSTR_EVENT[tot_tiles];
    HSTR_EVENT *eventcpyfr = new HSTR_EVENT[tot_tiles];
    HSTR_EVENT *eventgetrf = new HSTR_EVENT[tot_tiles];
    HSTR_EVENT *eventtrsm = new HSTR_EVENT[tot_tiles];
    HSTR_EVENT *eventgemm = new HSTR_EVENT[tot_tiles];

    //for timing tiled lu
    double *totTimeMsec = new double [niter];

    //for timing MKL lu
    double *totTimeMsecMKL = new double [niter];

    HSTR_RESULT res;

    int *ipiv_MKL = new int [mat_size];
    int *ipiv_my = new int [mat_size];

    //these queues are used for queining up compute on the card and
    //data transfers to/from the card.
    //q_trsm for dtrsm, q_getrf for custom_getrf, q_gemm for dgemm.
    //The queues q_getrf and q_gemm are incremented by one for
    //every compute queued, q_trsm is incremented by one for every
    //two computes. The queues wrap around the max_log_str available.
    //This ensures good load-balancing.
    int q_trsm, q_getrf, q_gemm;

    CBLAS_ORDER blasLay;
    int lapackLay;

    if (layRow) {
        blasLay = CblasRowMajor;
        lapackLay = LAPACK_ROW_MAJOR;
    } else {
        blasLay = CblasColMajor;
        lapackLay = LAPACK_COL_MAJOR;
    }

    for (iter = 0; iter < niter; ++iter) {

        //copying matrices into separate variables for tiled lu (A_my)
        //and MKL lu (A_MKL)
        //The output overwrites the matrices and hence the need to copy
        //for each iteration
        copy_mat(mat, A_my, mat_size);
        copy_mat(mat, A_MKL, mat_size);

        unsigned int m, n, k;

        printf("\nIteration = %d, Tbegin\n", iter);

        //beginning of timing
        tbegin = dtimeGet();

        //splitting time included in the timing
        //This splits the input matrix into tiles (or blocks)
        split_into_blocks(A_my, Asplit, num_tiles, tile_size,
                          mat_size, layRow);

        q_getrf = 0;
        for (k = 0; k < num_tiles; ++k) {
            //GETRF
            //custom dgetrf (no pivoting) is executed on the host on the
            //diagonal tile. The results are then sent to the card
            if (k > 0) {
                hStreams_app_event_wait(1, &eventgemm[k * num_tiles + k]);
                if (loc_verbose > 0)
                    printf("Sending tile[%d][%d] to host in queue %d\n",
                           k, k, (int)(q_getrf % max_log_str)) ;

                hStreams_app_xfer_memory(Asplit[k * num_tiles + k],
                                         Asplit[k * num_tiles + k], mem_size_tile,
                                         (int)(q_getrf % max_log_str), HSTR_SINK_TO_SRC,
                                         &eventcpyfr[k * num_tiles + k]);

                hStreams_app_event_wait(1, &eventcpyfr[k * num_tiles + k]);
            }

            if (loc_verbose > 0) {
                printf("Executing getrf on host for tile[%d][%d]\n", k, k);
            }

            if (layRow)
                info = custom_dgetrf_square_nopiv_rowMaj(tile_size,
                        Asplit[k * num_tiles + k]);
            else
                info = custom_dgetrf_square_nopiv_colMaj(tile_size,
                        Asplit[k * num_tiles + k]);

            if (k < num_tiles - 1) {
                if (loc_verbose > 0)
                    printf("Sending tile[%d][%d] to card in queue %d\n",
                           k, k, (int)(q_getrf % max_log_str));

                hStreams_app_xfer_memory(Asplit[k * num_tiles + k],
                                         Asplit[k * num_tiles + k], mem_size_tile,
                                         (int)(q_getrf % max_log_str), HSTR_SRC_TO_SINK,
                                         &eventcpyto[k * num_tiles + k]);
            }
            q_getrf++;

            q_trsm = 0;
            for (m = k + 1; m < num_tiles; ++m) {
                if (k == 0) {
                    if (loc_verbose > 0)
                        printf("Sending tile[%d][%d] to card in queue %d\n",
                               m, k, (int)(q_trsm % max_log_str));

                    hStreams_app_xfer_memory(Asplit[m * num_tiles + k],
                                             Asplit[m * num_tiles + k], mem_size_tile,
                                             (int)(q_trsm % max_log_str), HSTR_SRC_TO_SINK,
                                             &eventcpyto[m * num_tiles + k]);

                    if (loc_verbose > 0)
                        printf("Sending tile[%d][%d] to card in queue %d\n",
                               k, m, (int)(q_trsm % max_log_str));

                    hStreams_app_xfer_memory(Asplit[k * num_tiles + m],
                                             Asplit[k * num_tiles + m], mem_size_tile,
                                             (int)(q_trsm % max_log_str), HSTR_SRC_TO_SINK,
                                             &eventcpyto[k * num_tiles + m]);
                }

                //DTRSM
                hStreams_app_event_wait(1, &eventcpyto[k * num_tiles + k]);

                if (k > 0) {
                    hStreams_app_event_wait(1, &eventgemm[m * num_tiles + k]);
                    hStreams_app_event_wait(1, &eventgemm[k * num_tiles + m]);
                }

                if (loc_verbose > 0)
                    printf("Executing trsm for tile[%d][%d] on card in queue %d\n",
                           m, k, (int)(q_trsm % max_log_str));

                //dtrsm is executed on the card
                res = hStreams_custom_dtrsm(blasLay, CblasRight,
                                            CblasUpper, CblasNoTrans, CblasNonUnit,
                                            tile_size, tile_size, 1.0, Asplit[k * num_tiles + k],
                                            tile_size, Asplit[m * num_tiles + k], tile_size,
                                            (int)(q_trsm % max_log_str),
                                            &eventtrsm[m * num_tiles + k]);

                if (loc_verbose > 0)
                    printf("Executing trsm for tile[%d][%d] on card in queue %d\n",
                           k, m, (int)(q_trsm % max_log_str));

                res = hStreams_custom_dtrsm(blasLay, CblasLeft,
                                            CblasLower, CblasNoTrans, CblasUnit, tile_size,
                                            tile_size, 1.0, Asplit[k * num_tiles + k],
                                            tile_size, Asplit[k * num_tiles + m], tile_size,
                                            (int)(q_trsm % max_log_str),
                                            &eventtrsm[m * num_tiles + k]);

                if (loc_verbose > 0)
                    printf("Sending tile[%d][%d] back to host in queue %d\n",
                           m, k, (int)(q_trsm % max_log_str));

                hStreams_app_xfer_memory(Asplit[m * num_tiles + k],
                                         Asplit[m * num_tiles + k], mem_size_tile,
                                         (int)(q_trsm % max_log_str), HSTR_SINK_TO_SRC,
                                         &eventcpyfr[m * num_tiles + k]);

                if (loc_verbose > 0)
                    printf("Sending tile[%d][%d] back to host in queue %d\n",
                           k, m, (int)(q_trsm % max_log_str));

                hStreams_app_xfer_memory(Asplit[k * num_tiles + m],
                                         Asplit[k * num_tiles + m], mem_size_tile,
                                         (int)(q_trsm % max_log_str), HSTR_SINK_TO_SRC,
                                         &eventcpyfr[k * num_tiles + m]);

                q_trsm++;
            }

            q_gemm = 0;
            for (m = k + 1; m < num_tiles; ++m) {
                for (n = k + 1; n < num_tiles; ++n) {
                    if (k == 0) {
                        if (loc_verbose > 0)
                            printf("Sending tile[%d][%d] to card in queue %d\n",
                                   m, n, (int)(q_gemm % max_log_str));

                        hStreams_app_xfer_memory(Asplit[m * num_tiles + n],
                                                 Asplit[m * num_tiles + n], mem_size_tile,
                                                 (int)(q_gemm % max_log_str), HSTR_SRC_TO_SINK,
                                                 &eventcpyto[m * num_tiles + n]);

                    }

                    //DGEMM
                    hStreams_app_event_wait(1, &eventtrsm[m * num_tiles + k]);
                    hStreams_app_event_wait(1, &eventtrsm[k * num_tiles + n]);

                    if (k > 0) {
                        hStreams_app_event_wait(1, &eventgemm[m * num_tiles + n]);
                    }
                    //dgemm is executed on the card
                    if (loc_verbose > 0)
                        printf("Executing gemm for tile[%d][%d] on card in queue %d\n",
                               m, n, (int)(q_gemm % max_log_str));

                    res = hStreams_app_dgemm(blasLay, CblasNoTrans,
                                             CblasNoTrans, tile_size, tile_size, tile_size,
                                             -1.0, Asplit[m * num_tiles + k], tile_size,
                                             Asplit[k * num_tiles + n], tile_size, 1.0,
                                             Asplit[m * num_tiles + n], tile_size,
                                             (int)(q_gemm % max_log_str),
                                             &eventgemm[m * num_tiles + n]);

                    q_gemm++;
                }
            }
        }

        //syncrhonizing all the streams
        hStreams_app_thread_sync();

        //end of timing
        tend = dtimeGet();

        totTimeMsec[iter] = 1e3 * (tend - tbegin);
        printf("time for Tiled hstreams LU for iteration %d = %.2f msec\n",
               iter, totTimeMsec[iter]);

        //assembling of tiles back to form the full matrix
        assemble(Asplit, A_my, num_tiles, tile_size, mat_size, layRow);

        //calling mkl lu on host - AUTO OFFLOAD
        //this is to time the mkl lu for timing comparison
#ifndef _WIN32
        //FIXME: Calling this function on Windows causes a crash
        mkl_mic_enable();
#endif
        tbegin = dtimeGet();

        info = LAPACKE_dgetrf(lapackLay, mat_size, mat_size, A_MKL, mat_size, ipiv_MKL);

        tend = dtimeGet();
        totTimeMsecMKL[iter] = 1e3 * (tend - tbegin);
        printf("time for MKL LU (dgetrf) (AO) for iteration %d = %.2f msec\n",
               iter, totTimeMsecMKL[iter]);
        if (info != 0) {
            printf("error with dgetrf\n");
        }

        mkl_mic_disable();


        //verification of tiled-lu
        //for verification, we
        //1) obtain inverse of the input matrix using the lu factored matrix
        //2) multiply this inverse with the input matrix
        //3) compare the output of #2 with identity matrix
        //
        //Note: we cannot use the output from MKL dgetrf for verification
        //because MKL dgetrf does pivoting, and tiled-LU does not do pivoting.
        if (verify == 1) {
            //1. Generate ipiv for no pivoting
            for (int i = 0; i < mat_size; ++i) {
                ipiv_my[i] = i + 1;
            }

            //2. call getri for tiled LU to get matrix inverse
            info = LAPACKE_dgetri(lapackLay, mat_size, A_my,
                                  mat_size, ipiv_my);
            if (info != 0) {
                printf("error with dgetri, info = %d\n", info);
            }

            //3. multiply Ainv with Original A
            cblas_dgemm(blasLay, CblasNoTrans, CblasNoTrans, mat_size,
                        mat_size, mat_size, 1.0, A_my, mat_size, mat,
                        mat_size, 0.0, Cout, mat_size);

            //3. call verify_results (compare Cout with identity matrix)
            bool result = verify_results(Cout, eye, mat_size * mat_size);
            if (result == true) {
                printf("Tiled LU successful\n");
            } else {
                printf("Tiled LU failed\n");
            }
        }
    }

    double meanTimeMsec, stdDevMsec;
    double meanTimeMsecMKL, stdDevMsecMKL;

    double gflops = 2 * pow(mat_size, 3.0) / 3.0 * 1e-9;

    mean_and_stdev(totTimeMsec, meanTimeMsec, stdDevMsec, niter);
    mean_and_stdev(totTimeMsecMKL, meanTimeMsecMKL, stdDevMsecMKL, niter);

    printf("\nMatrix size = %d\n", mat_size);

    printf("Tiled hStreams LU: for %d iterations (ignoring first),\n"
           "mean Time = %.2f msec, stdDev Time = %.2f msec,\n"
           "Mean Gflops (using mean Time) = %.2f\n",
           niter - 1, meanTimeMsec, stdDevMsec,
           gflops / (meanTimeMsec * 1e-3));

    printf("\nMKL AO LU (dgetrf): For %d iterations (ignoring first),\n"
           "mean Time = %.2f msec, stdDev Time = %.2f msec,\n"
           "Mean Gflops (using mean Time) = %.2f\n\n",
           niter - 1, meanTimeMsecMKL, stdDevMsecMKL,
           gflops / (meanTimeMsecMKL * 1e-3));

    //Free
    free(A_my);
    free(A_MKL);
    free(Cout);
    free(eye);

    for (int i = 0; i < tot_tiles; ++i) {
        _mm_free(Asplit[i]);
    }

    delete [] Asplit;
    delete [] eventcpyto;
    delete [] eventcpyfr;
    delete [] eventgetrf;
    delete [] eventtrsm;
    delete [] eventgemm;
    delete [] totTimeMsec;
    delete [] totTimeMsecMKL;
    delete [] ipiv_MKL;
    delete [] ipiv_my;

}

int main(int argc, char **argv)
{
    HSTR_OPTIONS hstreams_options;
    hStreams_GetCurrentOptions(&hstreams_options, sizeof(HSTR_OPTIONS));

    hstreams_options.verbose = 0;
    char *libNames[200] = {NULL, NULL};

    //Library to be loaded for sink-side code
    libNames[0] = "lu_sink_1.so";
    hstreams_options.libNameCnt = 1;
    hstreams_options.libNames = libNames;
    hstreams_options.libFlags = NULL;

    int mat_size_m, num_tiles, niter, tile_size;
    niter = 5;
    num_tiles = 1;
    mat_size_m = 0; //must be an input
    bool layRow = true;

    //max_log_str defines the no. of physical partitions on the card
    int max_log_str = 5;

    int verify = 1;
    int use_num_threads;

    hStreams_SetOptions(&hstreams_options);
    for (int i = 1; i < argc; i++) {
        if (*argv[i] == SWITCH_CHAR) {
            switch (*(argv[i] + 1)) {
            case 'm':
                mat_size_m = (int)atol(argv[i] + 3);
                break;

            case 't':
                num_tiles = (int)atol(argv[i] + 3);
                break;

            case 's':
                max_log_str = (int)atol(argv[i] + 3);
                break;

            case 'l':
                if ((strcmp("row", argv[i] + 3) == 0) ||
                        (strcmp("ROW", argv[i] + 3) == 0)) {
                    layRow = true;
                    printf("matrix is in Row major format\n");
                } else {
                    layRow = false;
                    printf("matrix is in Col major format\n");
                }
                break;

            case 'i':
                niter = (int)atol(argv[i] + 3);
                if (niter < 3) {
                    niter = 3;
                }
                break;

            case 'v':
                verify = (int)atol(argv[i] + 3);
                break;

            default:
                break;
            }
        }
    }
    dtimeInit();

    printf("no. of streams (partitions) = %d, mat_size = %d, num_tiles = %d, "
           "niter = %d\n\n", max_log_str, mat_size_m, num_tiles, niter);

    //Check that mat_size is divisible by num_tiles
    if (mat_size_m % num_tiles != 0) {
        printf("matrix size MUST be divisible by num_tiles.. aborting\n");
        exit(0);
    }

    if (mat_size_m == 0) {
        printf("mat_size_m is not defined\n");
        exit(0);
    }

    tile_size = mat_size_m / num_tiles;

    //This allocates memory for the full input matrix
    double *A = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //Generate a symmetric positve-definite matrix
    A = dpo_generate(mat_size_m);

    //No. of PlacesPerDomain is same as no. of logical streams since LogStreamsPerPlace is 1.
    uint32_t PlacesPerDomain = max_log_str;
    uint32_t LogStreamsPerPlace = 1;
    hStreams_app_init(PlacesPerDomain, LogStreamsPerPlace);

    //Calling the tiled LU function. This does the factorization of the full matrix using a tiled implementation.
    lu_tiled(A, tile_size, num_tiles, mat_size_m, niter,
             max_log_str, layRow, verify);

    hStreams_app_fini();

    free(A);

    return 0;
}
