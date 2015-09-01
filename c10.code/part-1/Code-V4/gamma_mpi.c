#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "global.h"

_OFFLOADABLE void gamma_3d_offload(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare)
{
    int i, j;
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    int xsize_spl = get_b_xsize();
    double* xvec_spl = (double*) malloc(xsize_spl * sizeof(double));
    get_b_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) malloc(xsize_spl * sizeof(double));
    double* ixdiff_spl = (double*) malloc(xsize_spl * sizeof(double));
    for (i = 0; i < xsize_spl-1; i++)
    {
        xdiff_spl[i] = xvec_spl[i+1] - xvec_spl[i];
        ixdiff_spl[i] = 1.0 / xdiff_spl[i];
    }

    // jb zero the gamma array
    int m, n;
    for (m=0; m<gamma_size; m++)
        for (n=0; n<gamma_size; n++)
            gamma[m][n] = 0.0;

    // JB thread data for CO algorithm
    thread_g3d_data_t* thread_data;
    //corebarrier_t* core_barriers;

    // JB start omp parallel region here
    #pragma omp parallel default(none) private(n,m,i,j) \
    shared(rank, numranks, gamma, gamma_size, xsize_spl, xvec_spl, xdiff_spl, ixdiff_spl, thread_data, workshare)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        // JB allocate the thread data for all threads
        #pragma omp barrier
        #pragma omp master
        {
            thread_data = (thread_g3d_data_t*) malloc(nthreads * sizeof(thread_g3d_data_t));
            //core_barriers = (corebarrier_t*) malloc(nthreads * sizeof(corebarrier_t));
        }
        #pragma omp barrier

        // allocate private thread data
        thread_data[tid].mvec = (double*) malloc(gamma_size * sizeof(double));
        thread_data[tid].intgrlvec = (double*) malloc(gamma_size * sizeof(double));
        thread_data[tid].yvec = (double*) malloc(xsize_spl * sizeof(double));
        thread_data[tid].xsize = xsize_spl;
        thread_data[tid].xvec = xvec_spl;
        thread_data[tid].xdiff = xdiff_spl;
        thread_data[tid].ixdiff = ixdiff_spl;

        // jb v3 mkl spline objects 
        double* scoeff = (double*) malloc(4*(xsize_spl-1) * sizeof *scoeff);

        int status = dfdNewTask1D(&(thread_data[tid].task_spl), xsize_spl, thread_data[tid].xvec, DF_NO_HINT, 1, thread_data[tid].yvec, DF_NO_HINT);
        status = dfdEditPPSpline1D(thread_data[tid].task_spl, DF_PP_CUBIC, DF_PP_NATURAL, DF_BC_FREE_END, \
                0, DF_NO_IC, 0, scoeff, DF_NO_HINT);

        // zero the memory
        for (m=0; m<gamma_size; m++)
            thread_data[tid].mvec[m]=0.0;

        #pragma omp master
        {
            decompose_gamma_3d_mpi(thread_data, rank, numranks, workshare);
        }
        #pragma omp barrier

        double t0, t1;

        // SJP: For timing accuracy, ensure that all threads start and end sync'd
        #pragma omp barrier
        t0 = omp_get_wtime();

        // main calculation
        for (n=0; n < gamma_size; n++)
        {
            calculate_gamma_3d(n, thread_data);
            #pragma omp barrier

            // SJP: Just do this as a critical for now.  Swap for an array reduction later.
            #pragma omp critical
            {
                for (m = 0; m < gamma_size; m++)
                {
                    gamma[m][n] += thread_data[tid].mvec[m];
                }
            }
            for (m = 0; m < gamma_size; m++)
            {
                thread_data[tid].mvec[m]=0.0;
            }
        }
        #pragma omp barrier
        t1 = omp_get_wtime();

        /* JB free spline data */
        free(thread_data[tid].mvec);
        free(thread_data[tid].intgrlvec);
        //free(thread_data[tid].plijkz);
        free(thread_data[tid].yvec);
        //free(thread_data[tid].bar);
    }
    free(xvec_spl);
    free(xdiff_spl);
    free(ixdiff_spl);
    free(thread_data);
    printf("[%d] finished\n", rank);
}
