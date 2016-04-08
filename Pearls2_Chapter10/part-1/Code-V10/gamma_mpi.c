#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#include "global.h"
#include "local_barrier.h"

// blocksize used in gamma_3d
#define MAXBLOCKSIZE 64

_OFFLOADABLE void gamma_3d_offload(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare)
{
    int i, j;
    int gamma_size_pad = (gamma_size+7) & ~7;
    double (*restrict gamma)[gamma_size_pad] = (double (*restrict)[gamma_size_pad]) gamma_flat;

    int xsize_spl = get_b_xsize();
    double* xvec_spl = (double*) _mm_malloc(xsize_spl * sizeof(double),64);
    get_b_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double),64);
    double* ixdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double),64);
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
    corebarrier_t* core_barriers;

    // JB start omp parallel region here
    #pragma omp parallel default(none) private(n,m,i,j) \
    shared(rank, numranks, gamma, gamma_size, gamma_size_pad, xsize_spl, xvec_spl, xdiff_spl, ixdiff_spl, thread_data, core_barriers, workshare)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int maxblocksize = MAXBLOCKSIZE;

        #ifdef __MIC__
        int ntasks = N_MIC_CORES;
        int nworkers = nthreads/ntasks;
        int taskid = tid / nworkers;
        int workerid = tid % nworkers;
        #else
        int ntasks = nthreads;
        int nworkers = 1;
        int taskid = tid;
        int workerid = 0;
        #endif

        // JB allocate the thread data for all threads
        #pragma omp barrier
        #pragma omp master
        {
            thread_data = (thread_g3d_data_t*) _mm_malloc(nthreads * sizeof(thread_g3d_data_t),64);
            core_barriers = (corebarrier_t*) _mm_malloc(ntasks * sizeof(corebarrier_t),64);
        }
        #pragma omp barrier

        int task_master_tid = taskid*nworkers;

        // worker master threads to allocate task-shared work arrays
        if (workerid == 0)
        {
            thread_data[tid].mmat = (double*) _mm_malloc(gamma_size * gamma_size_pad * sizeof(double),64);
            thread_data[tid].intgrlvec = (double*) _mm_malloc(maxblocksize*gamma_size_pad * sizeof(double),64);
            thread_data[tid].plijkz = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double),64);
            core_barriers[taskid].userbarrier_arrive=0;
            core_barriers[taskid].userbarrier_depart=0;
        }
        #pragma omp barrier

        // other worker threads copy their masters pointer to shared work arrays
        if (workerid != 0)
        {
            thread_data[tid].mmat = thread_data[task_master_tid].mmat;
            thread_data[tid].intgrlvec = thread_data[task_master_tid].intgrlvec;
            thread_data[tid].plijkz = thread_data[task_master_tid].plijkz;
        }

        // allocate private thread data
        thread_data[tid].yvec = (double*) _mm_malloc(xsize_spl * sizeof(double),64);
        thread_data[tid].xsize = xsize_spl;
        thread_data[tid].xvec = xvec_spl;
        thread_data[tid].xdiff = xdiff_spl;
        thread_data[tid].ixdiff = ixdiff_spl;

        // core barriers (though not useful on the host)
        // every worker has a private local barrier struct 
        thread_data[tid].bar = (barrier_t*) _mm_malloc(sizeof(barrier_t),64);
        thread_data[tid].bar->me = &(core_barriers[taskid]);
        thread_data[tid].bar->usersense = 1;
        thread_data[tid].bar->mycoretid = workerid;
        #ifdef __MIC__
        thread_data[tid].bar->coreval = 0x01010101;    // assuming 4 threads per core
        #else
        thread_data[tid].bar->coreval = 0x00000001;    // assuming 1 threads per core
        #endif

        // zero the memory
        for (m=0; m<gamma_size; m++)
            for (n=0; n<gamma_size; n++)
                thread_data[tid].mmat[m*gamma_size_pad+n]=0.0;

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
        calculate_gamma_3d(thread_data, maxblocksize);

        // sum up all the results into gamma
        #pragma omp barrier
        t0 = omp_get_wtime();
        #pragma omp critical
        if (workerid == 0)
        {
            for (m = 0; m < gamma_size; m++)
            {
                for (n = 0; n < gamma_size; n++)
                {
                    gamma[m][n] += thread_data[tid].mmat[m*gamma_size_pad + n];
                }
            }
        }
        #pragma omp barrier
        #pragma omp barrier
        t1 = omp_get_wtime();
        //#pragma omp master
        //printf("[%d]\treduction\t%d\t%e\n", rank, n, t1-t0);

        #pragma omp master
        {
            _mm_free(core_barriers);
        }
        /* JB free spline data */
        if (workerid == 0)
        {
            _mm_free(thread_data[tid].mmat);
            _mm_free(thread_data[tid].intgrlvec);
            _mm_free(thread_data[tid].plijkz);
        }
        _mm_free(thread_data[tid].yvec);
        _mm_free(thread_data[tid].bar);
    }
    _mm_free(xvec_spl);
    _mm_free(xdiff_spl);
    _mm_free(ixdiff_spl);
    _mm_free(thread_data);
    printf("[%d] finished\n", rank);
}
