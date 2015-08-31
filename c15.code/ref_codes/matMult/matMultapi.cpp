/* FILE LICENSE TAG: SAMPLE */

//******************************************************************************************
// Author:   Gaurav Bansal
// Purpose:  Matrix multiply using hStreams
// Created:  5/2/14
//********************************************************************************************

#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>

#include "dtime.h"
#include <math.h>

//********************************************************************************************
#include <mkl.h>
#include <hStreams_source.h>
#include <hStreams_app_api.h>

#define SWITCH_CHAR  '-'
typedef FP_DATA_TYPE MATRIX_TYPE;

//********************************************************************************************
//
typedef struct _mtxSize { // Optional Command-line multiplier for matrix sizes
    unsigned int numColsA, numRowsA, numColsB, numRowsB, numColsC, numRowsC;
}   sMtxSize;

static int partitions  = 4;  // computing partitions per card
static int max_log_str = 4;  // number of logical streams: one can make this bigger
static int pass        = 0;  // number of passes in former_main

// Calculate mean gigaFlops
double mean_and_stddev(double *gflops, int size, double &mean, double &stddev, double &max, double &min)
{
    mean = stddev = 0.0;
    if (size < 3) {
        return (0.0);
    }
    min = 1.0e38;
    max = -min;
    for (int i = 1; i < size; ++i) {
        mean = mean + gflops[i];
        if (gflops[i] > max) {
            max = gflops[i];
        }
        if (gflops[i] < min) {
            min = gflops[i];
        }
    }
    mean = mean / (size - 1);

    for (int i = 1; i < size; ++i) {
        stddev = stddev + pow(gflops[i] - mean, 2);
    }

    stddev = stddev / (size - 2);
    stddev = sqrt(stddev);
    return (stddev);
}

// Initializes a matrix with random MATRIX_TYPE entries.
void randomInit(MATRIX_TYPE *data, int size)
{
    for (int i = 0; i < size; ++i) {
        data[i] = rand() / (MATRIX_TYPE)RAND_MAX;
    }
}

void initializeMatrices(sMtxSize &allSizes, int m, int k, int n)
{
    allSizes.numRowsA =  m;
    allSizes.numColsA =  k;
    allSizes.numRowsB =  k;
    allSizes.numColsB =  n;
    allSizes.numRowsC =  m;
    allSizes.numColsC =  n;
}

//
// Do Block Matrix Multiply
//
int matrixMultiply(sMtxSize &allSizes,
                   int mblocks, int kblocks, int nblocks,
                   int nIter)
{
    int totalblocks_A = mblocks * kblocks;
    int totalblocks_B = kblocks * nblocks;
    int totalblocks_C = mblocks * nblocks;

    unsigned int size_A = allSizes.numColsA * allSizes.numRowsA;
    unsigned int size_B = allSizes.numColsB * allSizes.numRowsB;
    unsigned int size_C = allSizes.numColsC * allSizes.numRowsC;

    // allocate host memory for matrices A and B and C
    MATRIX_TYPE *h_A = new MATRIX_TYPE[size_A];
    MATRIX_TYPE *h_B = new MATRIX_TYPE[size_B];
    MATRIX_TYPE *h_C = new MATRIX_TYPE[size_C];

    // host memory for matrix blocks for A and B and C
    MATRIX_TYPE **h_Ablock = new MATRIX_TYPE*[totalblocks_A];
    MATRIX_TYPE **h_Bblock = new MATRIX_TYPE*[totalblocks_B];
    MATRIX_TYPE **h_Cblock = new MATRIX_TYPE*[totalblocks_C];

    unsigned int mem_size_Ablock = sizeof(MATRIX_TYPE) * size_A / totalblocks_A;
    unsigned int mem_size_Bblock = sizeof(MATRIX_TYPE) * size_B / totalblocks_B;
    unsigned int mem_size_Cblock = sizeof(MATRIX_TYPE) * size_C / totalblocks_C;

    // malloc & initialize blocks of A/B/C on host and device
    for (int i = 0; i < totalblocks_A; ++i) {
        h_Ablock[i] = new MATRIX_TYPE[size_A / totalblocks_A];
        CHECK_HSTR_RESULT(hStreams_app_create_buf((void *)h_Ablock[i], mem_size_Ablock);
                          _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                "Block A[%d]: ",
                                                HSTR_INFO_TYPE_TRACE,
                                                i));

        // set seed for rand()
        srand(2006 + i);

        // initialize host memory
        randomInit(h_Ablock[i], size_A / totalblocks_A);
    }

    for (int i = 0; i < totalblocks_B; ++i) {
        h_Bblock[i] = new MATRIX_TYPE[size_B / totalblocks_B];
        CHECK_HSTR_RESULT(hStreams_app_create_buf((void *) h_Bblock[i], mem_size_Bblock);
                          _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                "Block B[%d]: ",
                                                HSTR_INFO_TYPE_TRACE,
                                                i));

        // set seed for rand()
        srand(201 + i);

        // initialize host memory
        randomInit(h_Bblock[i], size_B / totalblocks_B);
    }

    // Create buffers on the card
    for (int i = 0; i < totalblocks_C; ++i) {
        h_Cblock[i] = new MATRIX_TYPE[size_C / totalblocks_C];
        CHECK_HSTR_RESULT(hStreams_app_create_buf((void *) h_Cblock[i], mem_size_Cblock);
                          _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                "Block C[%d]: ",
                                                HSTR_INFO_TYPE_TRACE,
                                                i));
    }

    // Assemble matrices A and B for later host computation check
    int Bindex;
    for (int i = 0; i < kblocks; i++) {
        for (unsigned int k = 0; k < allSizes.numRowsB / kblocks; k++) {
            for (int j = i * nblocks; j < (i + 1)*nblocks; j++) {
                for (unsigned int l = 0; l < allSizes.numColsB / nblocks; l++) {
                    Bindex = allSizes.numColsB * (i * allSizes.numRowsB / kblocks + k) +
                             ((j - i * nblocks) * allSizes.numColsB / nblocks + l);
                    h_B[Bindex] = h_Bblock[j][l + k * allSizes.numColsB / nblocks];
                }
            }
        }
    }
    int Aindex;
    for (int i = 0; i < mblocks; i++) {
        for (unsigned int k = 0; k < allSizes.numRowsA / mblocks; k++) {
            for (int j = i; j <= i + (mblocks) * (kblocks - 1); j = j + mblocks) {
                for (unsigned int l = 0; l < allSizes.numColsA / kblocks; l++) {
                    Aindex = allSizes.numColsA * (i * allSizes.numRowsA / mblocks + k) +
                             (((j - i) / mblocks) * allSizes.numColsA / kblocks + l);
                    h_A[Aindex] = h_Ablock[j][l + k * allSizes.numColsA / kblocks];
                }
            }
        }
    }

    double *gigaFlops = new double[nIter];
    const MATRIX_TYPE alpha = 1.0;
    const MATRIX_TYPE beta  = 1.0;

    HSTR_EVENT *eventAcpy = new HSTR_EVENT[totalblocks_A];
    HSTR_EVENT *eventBcpy = new HSTR_EVENT[totalblocks_B];
    HSTR_EVENT *eventCcpy = new HSTR_EVENT[totalblocks_C];
    HSTR_EVENT *eventCset = new HSTR_EVENT[totalblocks_C];

    for (int iter = 0; iter < nIter; iter++) {
        double CdataRate = 0.0;
        double BdataRate = 0.0;
        double AdataRate = 0.0;

        double iterStart = dtimeGet();

        // Memset h_Cblock to zero initially
        for (int i = 0; i < totalblocks_C; ++i) {
            CHECK_HSTR_RESULT(hStreams_app_memset(h_Cblock[i], 0, mem_size_Cblock, (i % max_log_str), &eventCset[i]));
        }

        int qindex;
        int cblock_index = 0;
        int countA = -1, countB = -1;
        div_t divresult1, divresult2;

        while (countA < totalblocks_A - 1 || countB < totalblocks_B - 1) {
            if ((countA <= countB || countB == totalblocks_B - 1) && (countA < totalblocks_A - 1)) {
                countA++;
                double tbegin = dtimeGet();
                CHECK_HSTR_RESULT(hStreams_app_xfer_memory(
                                      h_Ablock[countA], h_Ablock[countA], mem_size_Ablock,
                                      (int)countA % max_log_str, HSTR_SRC_TO_SINK,
                                      &eventAcpy[countA]));
                double tend = dtimeGet();
                _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                      "Memcpy A Mic time= %.3f milliseconds\n",
                                      HSTR_INFO_TYPE_TRACE,
                                      1.0e3 * (tend - tbegin));
                if (tend != tbegin) {
                    AdataRate = 1.0e-9 * mem_size_Ablock / (tend - tbegin);
                    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                          "A Data Rate= %.3f GB/s\n",
                                          HSTR_INFO_TYPE_TRACE,
                                          AdataRate);
                }
            } else if ((countB < countA || countA == totalblocks_A - 1) && (countB < totalblocks_B - 1)) {
                countB++;
                double tbegin = dtimeGet();
                CHECK_HSTR_RESULT(hStreams_app_xfer_memory(
                                      h_Bblock[countB], h_Bblock[countB], mem_size_Bblock,
                                      (int)countB % max_log_str, HSTR_SRC_TO_SINK,
                                      &eventBcpy[countB]));
                double tend = dtimeGet();
                _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                      "Memcpy B Mic time= %.3f milliseconds\n",
                                      HSTR_INFO_TYPE_TRACE,
                                      1.0e3 * (tend - tbegin));

                if (tend != tbegin) {
                    BdataRate = 1.0e-9 * mem_size_Bblock / (tend - tbegin);
                    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                          "B Data Rate= %.3f GB/s\n",
                                          HSTR_INFO_TYPE_TRACE,
                                          BdataRate);
                }
            }

            // start multiplications
            if ((countA + countB) >= 0) {
                if (countA > countB) {
                    //received a new A block
                    for (int j = 0; j <= countB; j++) {
                        divresult1 = div(countA, mblocks);
                        divresult2 = div(j, nblocks);
                        if (divresult1.quot == divresult2.quot) {
                            //validity
                            cblock_index = divresult1.rem * nblocks + divresult2.rem;
                            qindex = (int)cblock_index % max_log_str;

                            _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                  "Multiplying A[%d]xB[%d] = C[%d].. Using logical stream %d...",
                                                  HSTR_INFO_TYPE_TRACE,
                                                  countA, j, cblock_index, qindex);

                            CHECK_HSTR_RESULT(hStreams_app_event_wait(1, &eventBcpy[j]));
                            CHECK_HSTR_RESULT(hStreams_app_event_wait(1, &eventAcpy[countA]));
                            CHECK_HSTR_RESULT(hStreams_app_event_wait(1, &eventCset[cblock_index]));

                            double tbegin = dtimeGet();
                            //Async compute on the device
                            CHECK_HSTR_RESULT(hStreams_app_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                                                 allSizes.numColsB / nblocks, allSizes.numRowsA / mblocks,
                                                                 allSizes.numRowsB / kblocks,
                                                                 alpha, h_Bblock[j], allSizes.numColsB / nblocks,
                                                                 h_Ablock[countA], allSizes.numColsA / kblocks,
                                                                 beta, h_Cblock[cblock_index], allSizes.numColsC / nblocks,
                                                                 qindex, &eventCcpy[cblock_index]));
                            double tend = dtimeGet();

                            _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                  "GEMM Time= %.3f microseconds\n",
                                                  HSTR_INFO_TYPE_TRACE,
                                                  1.0e6 * (tend - tbegin));
                        }
                    }
                } else { //received a new B block
                    for (int i = 0; i <= countA; i++) {
                        divresult1 = div(i, mblocks);
                        divresult2 = div(countB, nblocks);

                        if (divresult1.quot == divresult2.quot) {
                            //validity
                            cblock_index = divresult1.rem * nblocks + divresult2.rem;
                            qindex = (int)cblock_index % max_log_str;

                            _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                  "Multiplying A[%d]xB[%d] = C[%d].. Using logical stream %d...",
                                                  HSTR_INFO_TYPE_TRACE,
                                                  i, countB, cblock_index, qindex);

                            CHECK_HSTR_RESULT(hStreams_app_event_wait(1, &eventBcpy[countB]));
                            CHECK_HSTR_RESULT(hStreams_app_event_wait(1, &eventAcpy[i]));
                            CHECK_HSTR_RESULT(hStreams_app_event_wait(1, &eventCset[cblock_index]));

                            double tbegin = dtimeGet();
                            //Async compute
                            CHECK_HSTR_RESULT(hStreams_app_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                                                 allSizes.numColsB / nblocks, allSizes.numRowsA / mblocks,
                                                                 allSizes.numRowsB / kblocks,
                                                                 alpha, h_Bblock[countB], allSizes.numColsB / nblocks,
                                                                 h_Ablock[i], allSizes.numColsA / kblocks,
                                                                 beta, h_Cblock[cblock_index], allSizes.numColsC / nblocks,
                                                                 qindex, &eventCcpy[cblock_index]));
                            double tend = dtimeGet();

                            _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                                  "GEMM Time= %.3f microseconds\n",
                                                  HSTR_INFO_TYPE_TRACE,
                                                  1.0e6 * (tend - tbegin));
                        }
                    }
                }
            }
        }

        // copy result from device to host
        for (int i = 0; i < totalblocks_C; ++i) {
            _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                  "transferring matrix h_Cblock[%d] from device to host..asynchronously..using logical stream %d.",
                                  HSTR_INFO_TYPE_TRACE,
                                  i, i % max_log_str);

            //Synchronize before transferring the result back
            CHECK_HSTR_RESULT(hStreams_app_stream_sync(i % max_log_str));

            double tbegin = dtimeGet();
            CHECK_HSTR_RESULT(hStreams_app_xfer_memory(
                                  h_Cblock[i], h_Cblock[i], mem_size_Cblock, (int)i % max_log_str,
                                  HSTR_SINK_TO_SRC, &eventCcpy[i]));
            double tend = dtimeGet();
            _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                  "Memcpy C Time= %7.3f milliseconds --> ",
                                  HSTR_INFO_TYPE_TRACE,
                                  1.0e3 * (tend - tbegin));
            if (tend != tbegin) {
                CdataRate = 1.0e-9 * mem_size_Cblock / (tend - tbegin);
                _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                                      "C Data Rate= %7.3f GB/s\n",
                                      HSTR_INFO_TYPE_TRACE,
                                      CdataRate);
            }
        }

        //Synchronize after transferring the result back
        CHECK_HSTR_RESULT(hStreams_app_thread_sync());

        double iterStop = dtimeGet();

        //get Gflops
        double flopsPerMatrixMul = 2.0 * (double)allSizes.numColsA *
                                   (double)allSizes.numRowsA *
                                   (double)allSizes.numColsB;
        if (iterStop != iterStart) {
            gigaFlops[iter] = (flopsPerMatrixMul * 1.0e-9f) / (iterStop - iterStart);
        } else {
            gigaFlops[iter] = 0.0;
        }

        _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                              "Perf: %.3f Gflops/sec, Time= %.3f msec\r\n",
                              HSTR_INFO_TYPE_ALWAYSEMIT,
                              gigaFlops[iter], 1000.0 * (iterStop - iterStart));

    }

    double mean_gflops, stddev_gflops, max_gflops, min_gflops;
    mean_and_stddev(gigaFlops, nIter, mean_gflops, stddev_gflops, max_gflops, min_gflops);

    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                          "Size=, %d, %d, %d, %d, Pass=, %d, Iters=, %d, Max=, %.2f, GF/s, Avg_DGEMM=, %.2f, GFlop/s, StdDev=, %.2f, GFlop/s, %.2f, percent (Ignoring first iteration)\n",
                          HSTR_INFO_TYPE_ALWAYSEMIT,
                          allSizes.numRowsA, allSizes.numColsA, allSizes.numRowsB, allSizes.numColsB,
                          pass, nIter, max_gflops, mean_gflops, stddev_gflops, 100.0 * stddev_gflops / mean_gflops);

    int Cindex;
    for (int i = 0; i < mblocks; i++) {
        for (unsigned int k = 0; k < allSizes.numRowsC / mblocks; k++) {
            for (int j = i * nblocks; j < (i + 1)*nblocks; j++) {
                for (unsigned int l = 0; l < allSizes.numColsC / nblocks; l++) {
                    Cindex = allSizes.numColsC * (i * allSizes.numRowsC / mblocks + k) +
                             ((j - i * nblocks) * allSizes.numColsC / nblocks + l);
                    h_C[Cindex] = h_Cblock[j][l + k * allSizes.numColsC / nblocks];
                }
            }
        }
    }

    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                          "Computing result using host CPU...",
                          HSTR_INFO_TYPE_ALWAYSEMIT);

    MATRIX_TYPE *reference = new MATRIX_TYPE[size_C];
    bool resACBLAS = true;

    // faster MKL host compute reference solution
    memset(reference, 0, size_C * sizeof(MATRIX_TYPE));
    double refFlops = (2.0 * allSizes.numRowsC * allSizes.numColsC * allSizes.numColsA);
    double trefStart = dtimeGet();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                allSizes.numRowsA, allSizes.numColsB, allSizes.numColsA,
                1.0, h_A, allSizes.numColsA, h_B, allSizes.numColsB,
                0.0, reference, allSizes.numColsC);
    double trefStop = dtimeGet();
    resACBLAS = true;
    uint32_t num_mismatch = 0;
    for (unsigned int i = 0; i < size_C; ++i) {
        FP_DATA_TYPE diff = fabs(reference[i] - h_C[i]);
        if (reference[i] != 0) {
            diff /= reference[i];
        }
        if (diff > 1.0e-7) {
            num_mismatch++;
            if (resACBLAS) {
                printf("MKL dgemm: Detected problem at i=%d: ref %g actual %g\n",
                       i, reference[i], h_C[i]);
            }
            resACBLAS = false;
        }
    }

    if (resACBLAS) {
        printf("[If no MKLdgemm failures, then] Block Multiplication was successful.\r\n");
    } else {
        printf("ERROR: Mismatch between MKL's result and hStreams' result.\n"
               "ERROR: Number of mismatching elements: %u, number of matching elements: %u\n",
               num_mismatch, size_C - num_mismatch);
    }

    double mklHostRefFlops = (refFlops * 1.0e-9) / (trefStop - trefStart);
    printf("MKL Host DGEMM Perf, %.3f, GFlops/sec, Time= %.3f msec\r\n",
           mklHostRefFlops, 1000.0 * (trefStop - trefStart));
    //
    // Free Memory
    //
    delete [] reference;
    delete [] h_A;
    delete [] h_B;
    delete [] h_C;
    delete [] gigaFlops;
    delete [] eventAcpy;
    delete [] eventBcpy;
    delete [] eventCcpy;
    delete [] eventCset;

    for (int i = 0; i < totalblocks_A; ++i) {
        delete [] h_Ablock[i];
    }
    for (int i = 0; i < totalblocks_B; ++i) {
        delete [] h_Bblock[i];
    }
    for (int i = 0; i < totalblocks_C; ++i) {
        delete [] h_Cblock[i];
    }

    return 0;
}

static HSTR_OPTIONS hstreams_options;

static int formerMain(int argc, char **argv)
{
    hStreams_GetCurrentOptions(&hstreams_options, sizeof(hstreams_options));
    hstreams_options.verbose = 0;
    hstreams_options.phys_domains_limit = 256;  // Limit to 256 domains
    int nIter = 5, blocksize = 128;
    int mat_size_m = 512, mat_size_k = 512, mat_size_n = 5120;
    sMtxSize allSizes;

    ++pass;
    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                          "*** C(mxn) = A(mxk) x B(kxn)\n", HSTR_INFO_TYPE_ALWAYSEMIT);

    while (argc) {
        if (*argv[0] == SWITCH_CHAR) {
            switch (*(argv[0] + 1)) {
            case 'b':
                blocksize = (int)atol(argv[0] + 2);
                break;

            case 'm':
                mat_size_m = (int)atol(argv[0] + 2);
                break;

            case 'k':
                mat_size_k = (int)atol(argv[0] + 2);
                break;

            case 'n':
                mat_size_n = (int)atol(argv[0] + 2);
                break;

            case 'i':
                nIter = (int)atol(argv[0] + 2);
                if (nIter < 3) {
                    nIter = 3;
                }
                break;

            case 'v':
                hstreams_options.verbose = (int)atol(argv[0] + 2);
                break;

            default:
                break;
            }
        }
        argc -= 1;
        argv++;
    }

    hStreams_SetOptions(&hstreams_options);
    int mblocks, nblocks, kblocks;
    mblocks = (int)mat_size_m / blocksize;
    kblocks = (int)mat_size_k / blocksize;
    nblocks = (int)mat_size_n / blocksize;

    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                          "blocksize = %d x %d, m = %d, k = %d, n = %d\n",
                          HSTR_INFO_TYPE_ALWAYSEMIT,
                          blocksize, blocksize, mat_size_m, mat_size_k, mat_size_n);

    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                          "mblocks = %d, kblocks = %d, nblocks = %d\n",
                          HSTR_INFO_TYPE_ALWAYSEMIT,
                          mblocks, kblocks, nblocks);

    _hStreams_EmitMessage(HSTR_SEVERITY_INFO, __FUNCTION__,
                          "Matrix in blocks A(%d,%d), B(%d,%d), C(%d,%d)\n",
                          HSTR_INFO_TYPE_ALWAYSEMIT,
                          mblocks, kblocks, kblocks, nblocks, mblocks, nblocks);

    //the second arg is logStrPerParition! not perCard
    int iret = hStreams_app_init(partitions, 1);
    if (iret != 0) {
        printf("hstreams_app_init failed!\r\n");
        exit(-1);
    }

    initializeMatrices(allSizes, mat_size_m, mat_size_k, mat_size_n);

    matrixMultiply(allSizes, mblocks, kblocks, nblocks, nIter);

    CHECK_HSTR_RESULT(hStreams_app_thread_sync());

    CHECK_HSTR_RESULT(hStreams_app_fini());

    return 0;

}

int main(int argc, char **argv)
{
    dtimeInit();
    return formerMain(argc, argv);
}
