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
#include "dtime.h"
#include <string.h>
#include "helper_functions.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define SWITCH_CHAR  '-'

extern double *dpo_generate(size_t);

//Implements the tiled version of Cholesky factorization
void tiled_cholesky(double *mat_sp[], int tile_size, int num_tiles,
                    CBLAS_ORDER blasLay, int lapackLay)
{
    unsigned int m, n, k;
    int info;

    for (k = 0; k < num_tiles; ++k) {
        //POTRF - MKL call
        info = LAPACKE_dpotrf(lapackLay, 'L', tile_size, mat_sp[k * num_tiles + k],
                              tile_size);

        for (m = k + 1; m < num_tiles; ++m) {
            //DTRSM - MKL call
            cblas_dtrsm(blasLay, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
                        tile_size, tile_size, 1.0, mat_sp[k * num_tiles + k], tile_size,
                        mat_sp[m * num_tiles + k], tile_size);
        }

        for (n = k + 1; n < num_tiles; ++n) {
            //DSYRK - MKL call
            cblas_dsyrk(blasLay, CblasLower, CblasNoTrans, tile_size, tile_size,
                        -1.0, mat_sp[n * num_tiles + k], tile_size, 1.0,
                        mat_sp[n * num_tiles + n], tile_size);

            for (m = n + 1; m < num_tiles; ++m) {
                //DGEMM - MKL call
                cblas_dgemm(blasLay, CblasNoTrans, CblasTrans, tile_size,
                            tile_size, tile_size, -1.0, mat_sp[m * num_tiles + k],
                            tile_size, mat_sp[n * num_tiles + k], tile_size, 1.0,
                            mat_sp[m * num_tiles + n], tile_size);
            }
        }
    }
}


int main(int argc, char **argv)
{
    int info;
    double tbegin, tend;

    int mat_size_m, num_tiles, niter, tile_size, tot_tiles;
    niter = 5;
    num_tiles = 1;
    mat_size_m = 0; //must be an input

    bool layRow = true;
    int verify = 1;

    for (int i = 1; i < argc; i++) {
        if (*argv[i] == SWITCH_CHAR) {
            switch (*(argv[i] + 1)) {
            case 'm':
                mat_size_m = (int)atol(argv[i] + 3);
                break;

            case 't':
                num_tiles = (int)atol(argv[i] + 3);
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
    printf("mat_size = %d, num_tiles = %d, niter = %d\n\n",
           mat_size_m, num_tiles, niter);

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
    tot_tiles = num_tiles * num_tiles;

    //allocating memory for input matrix (full matrix)
    double *A = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //allocating memory for tiled_cholesky for the full matrix
    double *A_my = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //allocating memory for MKL cholesky for the full matrix
    double *A_MKL = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //memory allocation for tiled matrix
    double **Asplit = new double* [tot_tiles];

    for (int i = 0; i < tot_tiles; ++i) {
        //Buffer per tile
        Asplit[i] = (double *)malloc(tile_size * tile_size * sizeof(double));
    }

    //Generate a symmetric positve-definite matrix
    A = dpo_generate(mat_size_m);

    double *totTimeMsecTile = new double [niter];
    double *totTimeMsecMKL = new double [niter];

    CBLAS_ORDER blasLay;
    int lapackLay;

    if (layRow) {
        blasLay = CblasRowMajor;
        lapackLay = LAPACK_ROW_MAJOR;
    } else {
        blasLay = CblasColMajor;
        lapackLay = LAPACK_COL_MAJOR;
    }

    for (int iter = 0; iter < niter; ++iter) {
        printf("\niter %d\n", iter);

        //copying matrices into separate variables for tiled cholesky (A_my)
        //and MKL cholesky (A_MKL)
        //The output overwrites the matrices and hence the need to copy
        //for each iteration
        copy_mat(A, A_my, mat_size_m);
        copy_mat(A, A_MKL, mat_size_m);

        mkl_mic_disable(); //disabling auto offload feature of MKL

        //beginning of timing
        tbegin = dtimeGet();

        //splitting time included in the timing
        //This splits the input matrix into tiles (or blocks)
        split_into_blocks(A_my, Asplit, num_tiles, tile_size, mat_size_m, layRow);

        //Calling the tiled Cholesky function. This does the factorization of the full matrix using a tiled implementation.
        tiled_cholesky(Asplit, tile_size, num_tiles, blasLay, lapackLay);

        //end of timing
        tend = dtimeGet();

        totTimeMsecTile[iter] = 1e3 * (tend - tbegin);
        printf("time for Tiled Cholesky = %.2f msec\n",
               totTimeMsecTile[iter]);

        //assembling of tiles back into full matrix
        assemble(Asplit, A_my, num_tiles, tile_size, mat_size_m, layRow);

        tbegin = dtimeGet();

        //calling mkl cholesky for verification and timing comparison
        info = LAPACKE_dpotrf(lapackLay, 'L', mat_size_m, A_MKL, mat_size_m);

        tend = dtimeGet();
        totTimeMsecMKL[iter] = 1e3 * (tend - tbegin);
        printf("time for MKL Cholesky (host, NO AO) = %.2f msec\n",
               totTimeMsecMKL[iter]);
        if (info != 0) {
            printf("error with dpotrf, info = %d\n", info);
        }

        if (verify == 1) {
            bool res = verify_results(A_my, A_MKL, mat_size_m * mat_size_m);
            if (res == true) {
                printf("Tiled Cholesky successful\n");
            } else {
                printf("Tiled Chloesky failed\n");
            }
        }
    }

    double meanTimeMsecTile, stdDevMsecTile;
    double meanTimeMsecMKL, stdDevMsecMKL;

    mean_and_stdev(totTimeMsecTile, meanTimeMsecTile, stdDevMsecTile, niter);
    mean_and_stdev(totTimeMsecMKL, meanTimeMsecMKL, stdDevMsecMKL, niter);

    double gflops = pow(mat_size_m, 3.0) / 3.0 * 1e-9;

    printf("\nTiled Cholesky, for %d iterations (ignoring first),\n"
           "mean Time = %.2f msec, stdDev Time = %.2f msec,\n"
           "Mean Gflops (using mean Time) = %.2f\n",
           niter - 1, meanTimeMsecTile, stdDevMsecTile,
           gflops / (meanTimeMsecTile * 1e-3));

    printf("\nMKL Cholesky (host, NO AO), for %d iterations (ignoring first),\n"
           "mean Time = %.2f msec, stdDev Time = %.2f msec,\n"
           "Mean Gflops (using mean Time) = %.2f\n\n",
           niter - 1, meanTimeMsecMKL, stdDevMsecMKL,
           gflops / (meanTimeMsecMKL * 1e-3));

    //free
    free(A);
    free(A_my);
    free(A_MKL);
    for (int i = 0; i < tot_tiles; ++i) {
        free(Asplit[i]);
    }

    delete [] Asplit;
    delete [] totTimeMsecTile;
    delete [] totTimeMsecMKL;

    return 0;
}
