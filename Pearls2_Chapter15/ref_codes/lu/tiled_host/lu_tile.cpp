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
#include <iostream>
#include <string.h>
#include "helper_functions.h"
#include "dtime.h"

using namespace std;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define SWITCH_CHAR  '-'

extern double *dpo_generate(size_t);

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

//Implements the tiled version of the LU factorization
void tiled_lu(double *mat_sp[], int tile_size, int num_tiles,
              CBLAS_ORDER blasLay, int lapackLay, bool layRow)
{
    unsigned int m, n, k;
    int info;

    for (k = 0; k < num_tiles; ++k) {
        //GETRF - MKL call
        if (layRow)
            info = custom_dgetrf_square_nopiv_rowMaj(tile_size,
                    mat_sp[k * num_tiles + k]);
        else
            info = custom_dgetrf_square_nopiv_colMaj(tile_size,
                    mat_sp[k * num_tiles + k]);

        //DTRSM - MKL call
        for (m = k + 1; m < num_tiles; ++m) {
            //Get L column
            cblas_dtrsm(blasLay, CblasRight, CblasUpper, CblasNoTrans,
                        CblasNonUnit, tile_size, tile_size, 1.0,
                        mat_sp[k * num_tiles + k], tile_size,
                        mat_sp[m * num_tiles + k], tile_size);

            //Get U row
            cblas_dtrsm(blasLay, CblasLeft, CblasLower, CblasNoTrans,
                        CblasUnit, tile_size, tile_size, 1.0,
                        mat_sp[k * num_tiles + k], tile_size,
                        mat_sp[k * num_tiles + m], tile_size);
        }

        //Get Schur Complement
        for (m = k + 1; m < num_tiles; ++m) {
            for (n = k + 1; n < num_tiles; ++n) {
                //DGEMM - MKL call
                cblas_dgemm(blasLay, CblasNoTrans, CblasNoTrans,
                            tile_size, tile_size, tile_size, -1.0,
                            mat_sp[m * num_tiles + k], tile_size,
                            mat_sp[k * num_tiles + n], tile_size, 1.0,
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

    //allocating memory for tiled_lu for the full matrix
    double *A_my = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //allocating memory for MKL lu for the full matrix
    double *A_MKL = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //used for verification of results
    double *eye = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));
    double *Cout = (double *)malloc(mat_size_m * mat_size_m * sizeof(double));

    //memory allocation for tiled matrix
    double **Asplit = new double* [tot_tiles];

    for (int i = 0; i < tot_tiles; ++i) {
        //Buffer per tile
        Asplit[i] = (double *)malloc(tile_size * tile_size * sizeof(double));
    }

    //Generate a symmetric positve-definite matrix
    A = dpo_generate(mat_size_m);

    double totTimeMsecTile[niter];
    double totTimeMsecMKL[niter];

    int *ipiv_MKL = new int [mat_size_m];
    int *ipiv_my = new int [mat_size_m];

    eyeInit(eye, mat_size_m);

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

        //copying matrices into separate variables for tiled lu (A_my)
        //and MKL lu (A_MKL)
        //The output overwrites the matrices and hence the need to copy
        //for each iteration
        copy_mat(A, A_my, mat_size_m);
        copy_mat(A, A_MKL, mat_size_m);

        mkl_mic_disable(); //disabling auto-offload feature of MKL

        //beginning of timing
        tbegin = dtimeGet();

        //splitting time included in the timing
        //This splits the input matrix into tiles (or blocks)
        split_into_blocks(A_my, Asplit, num_tiles, tile_size, mat_size_m, layRow);

        //Calling the tiled lu function. This does the factorization of the full matrix using a tiled implementation.
        tiled_lu(Asplit, tile_size, num_tiles, blasLay, lapackLay, layRow);

        //end of timing
        tend = dtimeGet();

        totTimeMsecTile[iter] = 1e3 * (tend - tbegin);
        printf("time for Tiled LU = %f msec\n",
               totTimeMsecTile[iter]);

        //assembling of tiles back into full matrix
        assemble(Asplit, A_my,  num_tiles, tile_size, mat_size_m, layRow);

        tbegin = dtimeGet();

        //calling mkl LU for timing comparison
        info = LAPACKE_dgetrf(lapackLay, mat_size_m, mat_size_m, A_MKL, mat_size_m, ipiv_MKL);

        tend = dtimeGet();
        totTimeMsecMKL[iter] = 1e3 * (tend - tbegin);
        printf("time for MKL LU (dgetrf) (host, NO AO) = %f msec\n",
               totTimeMsecMKL[iter]);
        if (info != 0) {
            printf("error with getrf\n");
        }

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
            for (int i = 0; i < mat_size_m; ++i) {
                ipiv_my[i] = i + 1;
            }

            //2. call getri for tiled LU to get matrix inverse
            info = LAPACKE_dgetri(lapackLay, mat_size_m, A_my, mat_size_m, ipiv_my);
            if (info != 0) {
                printf("error with dgetri, info = %d\n", info);
            }

            //3. multiply Ainv with Original A
            cblas_dgemm(blasLay, CblasNoTrans, CblasNoTrans, mat_size_m, mat_size_m, mat_size_m, 1.0, A_my, mat_size_m, A, mat_size_m, 0.0, Cout, mat_size_m);

            //3. call verify_results (compare Cout with identity matrix)
            bool result = verify_results(Cout, eye, mat_size_m * mat_size_m);
            if (result == true) {
                printf("Tiled LU successful\n");
            } else {
                printf("Tiled-LU failed\n");
            }
        }
    }

    double meanTimeMsecTile, stdDevMsecTile;
    double meanTimeMsecMKL, stdDevMsecMKL;

    double gflops = 2 * pow(mat_size_m, 3.0) / 3.0 * 1e-9;

    mean_and_stdev(totTimeMsecTile, meanTimeMsecTile, stdDevMsecTile, niter);
    mean_and_stdev(totTimeMsecMKL, meanTimeMsecMKL, stdDevMsecMKL, niter);

    printf("\nTiled LU, for %d iterations (ignoring first),\n"
           "mean Time = %.2f msec, stdDev Time = %.2f msec,\n"
           "Mean Gflops (using mean Time) = %.2f\n",
           niter - 1, meanTimeMsecTile, stdDevMsecTile,
           gflops / (meanTimeMsecTile * 1e-3));

    printf("\nMKL LU, For %d iterations (ignoring first),\n"
           "mean Time = %.2f msec, stdDev Time = %.2f msec,\n"
           "Mean Gflops (using mean Time) = %.2f\n\n",
           niter - 1, meanTimeMsecMKL, stdDevMsecMKL,
           gflops / (meanTimeMsecMKL * 1e-3));

    //Free
    free(A);
    free(A_my);
    free(A_MKL);
    free(eye);
    free(Cout);

    for (int i = 0; i < tot_tiles; ++i) {
        free(Asplit[i]);
    }

    delete [] Asplit;
    delete [] ipiv_MKL;
    delete [] ipiv_my;

    return 0;
}
