#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#pragma offload_attribute(push, target(mic))
#include <omp.h>
#pragma offload_attribute(pop)
#include <mkl.h>

#include "global.h"

#define MAXBLOCKSIZE 64

_OFFLOADABLE void gamma_3d_host(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare)
{
    int i, j;
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    int xsize_spl = get_b_xsize();
    int xsize_pad = (xsize_spl + 7) & ~7;	
    double* xvec_spl = (double*) _mm_malloc(xsize_pad * sizeof(double) + 64, 64);
    get_b_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
    double* ixdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
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
    shared(rank, numranks, gamma, gamma_size, xsize_spl, xsize_pad, xvec_spl, xdiff_spl, ixdiff_spl, thread_data, core_barriers, workshare)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int maxblocksize = MAXBLOCKSIZE;

        // JB allocate the thread data for all threads
        #pragma omp barrier
        #pragma omp master
        {
            thread_data = (thread_g3d_data_t*) _mm_malloc(nthreads * sizeof(thread_g3d_data_t),64);
            core_barriers = (corebarrier_t*) _mm_malloc(nthreads * sizeof(corebarrier_t), 64);
        }
        #pragma omp barrier

        // allocate private thread data
        int gamma_size_pad = (gamma_size+7) & ~7;
        thread_data[tid].mvec = (double*) _mm_malloc(gamma_size_pad * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].intgrlvec = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].plijkz = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].yvec[0] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
        thread_data[tid].xsize = xsize_spl;
        thread_data[tid].xvec = xvec_spl;
        thread_data[tid].xdiff = xdiff_spl;
        thread_data[tid].ixdiff = ixdiff_spl;

        // core barriers (though not useful on the host)
        // TODO JB get rid of core_barriers on the host side?
        core_barriers[tid].userbarrier_arrive=0;
        core_barriers[tid].userbarrier_depart=0;
        thread_data[tid].bar = (barrier_t*) _mm_malloc(sizeof(barrier_t), 64);
        thread_data[tid].bar->me = &(core_barriers[tid]);
        thread_data[tid].bar->usersense = 1;
        thread_data[tid].bar->mycoretid = 0;
        thread_data[tid].bar->coreval = 0x00000001;    // assuming 1 threads per core

        // zero the memory
        for (m=0; m<gamma_size; m++)
            for (n=0; n<gamma_size; n++)
                thread_data[tid].mvec[m*gamma_size_pad+n]=0.0;

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
        calculate_gamma_3d_nested(thread_data,maxblocksize);

        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\tgamma_3D \t%d\t%e\n", rank, n, t1-t0);

        // SJP: Just do this as a critical for now.  Swap for an array reduction later.
        #pragma omp barrier
        t0 = omp_get_wtime();
        #pragma omp critical
        for (m = 0; m < gamma_size; m++)
        {
            for (n = 0; n < gamma_size; n++)
            {
                gamma[m][n] += thread_data[tid].mvec[m*gamma_size_pad + n];
            }
        }
        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\treduction\t%d\t%e\n", rank, n, t1-t0);

        /* JB free spline data */

        #pragma omp master
        {
            _mm_free(core_barriers);
        }
        _mm_free(thread_data[tid].mvec);
        _mm_free(thread_data[tid].intgrlvec);
        _mm_free(thread_data[tid].plijkz);
        _mm_free(thread_data[tid].yvec[0]);
        _mm_free(thread_data[tid].bar);
    }
    _mm_free(xvec_spl);
    _mm_free(xdiff_spl);
    _mm_free(ixdiff_spl);
    _mm_free(thread_data);
}

void gamma_3d_offload(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare)
{
    int i, j;
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    int xsize_spl = get_b_xsize();
    int xsize_pad = (xsize_spl + 7) & ~7;
    double* xvec_spl = (double*) _mm_malloc(xsize_pad * sizeof(double) + 64, 64);
    get_b_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
    double* ixdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
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
    shared(rank, numranks, gamma, gamma_size, xsize_spl, xsize_pad, xvec_spl, xdiff_spl, ixdiff_spl, thread_data, core_barriers, workshare)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int nworkers = nthreads/59;
        int taskid = tid / nworkers;
        int workerid = tid % nworkers;
        int maxblocksize = MAXBLOCKSIZE;

        // JB allocate the thread data for all threads
        #pragma omp barrier
        #pragma omp master
        {
            thread_data = (thread_g3d_data_t*) _mm_malloc(nthreads * sizeof(thread_g3d_data_t),64);
            core_barriers = (corebarrier_t*) _mm_malloc(59 * sizeof(corebarrier_t), 64);
        }
        #pragma omp barrier

        int task_master_tid = taskid*nworkers;

        // worker master threads to allocate task-shared work arrays
        if (workerid == 0)
        {
            int gamma_size_pad = (gamma_size+7) & ~7;
            thread_data[tid].mvec = (double*) _mm_malloc(gamma_size_pad * gamma_size_pad * sizeof(double), 64);
            thread_data[tid].intgrlvec = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
            thread_data[tid].plijkz = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
            core_barriers[taskid].userbarrier_arrive=0;
            core_barriers[taskid].userbarrier_depart=0;
        }
        #pragma omp barrier

        // other worker threads copy their masters pointer to shared work arrays
        if (workerid != 0)
        {
            thread_data[tid].mvec = thread_data[task_master_tid].mvec;
            thread_data[tid].intgrlvec = thread_data[task_master_tid].intgrlvec;
            thread_data[tid].plijkz = thread_data[task_master_tid].plijkz;
        }

        // every worker malloc a private yvec array
        thread_data[tid].yvec[0] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);

        // every worker has a private local barrier struct 
        thread_data[tid].bar = (barrier_t*) _mm_malloc(sizeof(barrier_t), 64);
        thread_data[tid].bar->me = &(core_barriers[taskid]);
        thread_data[tid].bar->usersense = 1;
        thread_data[tid].bar->mycoretid = workerid;
        thread_data[tid].bar->coreval = 0x01010101;    // assuming 4 threads per core

        // shared ready-only data for all workers
        thread_data[tid].xsize = xsize_spl;
        thread_data[tid].xvec = xvec_spl;
        thread_data[tid].xdiff = xdiff_spl;
        thread_data[tid].ixdiff = ixdiff_spl;

        if (workerid == 0)
        {
            // zero the memory
            int gamma_size_pad = (gamma_size+7) & ~7;
            for (m=0; m<gamma_size; m++)
                for (n=0; n<gamma_size; n++)
                    thread_data[tid].mvec[m*gamma_size_pad+n]=0.0;
        }
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
        calculate_gamma_3d_nested(thread_data,maxblocksize);

        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\tgamma_3D \t%d\t%e\n", rank, n, t1-t0);

        // SJP: Just do this as a critical for now.  Swap for an array reduction later.
        #pragma omp barrier
        t0 = omp_get_wtime();
        #pragma omp critical
        if (workerid == 0)
        {
            int gamma_size_pad = (gamma_size+7) & ~7;
            for (m = 0; m < gamma_size; m++)
            {
                for (n = 0; n < gamma_size; n++)
                {
                    gamma[m][n] += thread_data[tid].mvec[m*gamma_size_pad + n];
                }
            }
        }
        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\treduction\t%d\t%e\n", rank, n, t1-t0);

        /* JB free spline data */

        #pragma omp master
        {
            _mm_free(core_barriers);
        }
        if(workerid==0)
        {
            _mm_free(thread_data[tid].mvec);
            _mm_free(thread_data[tid].intgrlvec);
            _mm_free(thread_data[tid].plijkz);
        }
        _mm_free(thread_data[tid].yvec[0]);
        _mm_free(thread_data[tid].bar);
    }
    _mm_free(xvec_spl);
    _mm_free(xdiff_spl);
    _mm_free(ixdiff_spl);
    _mm_free(thread_data);
}

#ifdef __INTEL_OFFLOAD
void gamma_3d_offload_teams(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare)
{
    int i, j;
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    int xsize_spl = get_b_xsize();
    int xsize_pad = (xsize_spl + 7) & ~7;
    double* xvec_spl = (double*) _mm_malloc(xsize_pad * sizeof(double) + 64, 64);
    get_b_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
    double* ixdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
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

    // JB start omp parallel region here
    int tid;
    const int nteams = omp_get_max_teams();
    printf("nteams = %d\n", nteams);
    int nthreads = nteams;  // lazy james
    const int nworkers = 1;
    int maxblocksize = MAXBLOCKSIZE;

    // JB allocate the thread data for all threads
    thread_data = (thread_g3d_data_t*) _mm_malloc(nthreads * sizeof(thread_g3d_data_t),64);
    printf("sizeof(thread_data) = %d bytes\n", sizeof(thread_g3d_data_t));

    #pragma omp target data device(0) map(to: xvec_spl[0:xsize_pad], xdiff_spl[0:xsize_spl], ixdiff_spl[0:xsize_spl], thread_data) map(tofrom: gamma_flat[0:gamma_size*gamma_size])
    {

        // SJP: Two target regions here are necessary because:
        // 1) decompose_gamma_3d_mpi doesn't need to be called for all teams (although it could be, thinking about it);
        // 2) MKL initialization needs to take place on the core where we want MKL to run.  I think.
        #pragma omp target device(0)
        {
            decompose_gamma_3d_mpi(thread_data, rank, numranks, workshare);
        }

        #pragma omp target device(0)
        #pragma omp teams default(none) shared(xvec_spl, xdiff_spl, ixdiff_spl, gamma_flat, gamma_size, xsize_spl, xsize_pad, thread_data, maxblocksize) num_teams(59) thread_limit(4)
        {

            int tid = omp_get_team_num();

            // worker master threads to allocate task-shared work arrays
            int gamma_size_pad = (gamma_size+7) & ~7;
            thread_data[tid].mvec = (double*) _mm_malloc(gamma_size_pad * gamma_size_pad * sizeof(double), 64);
            thread_data[tid].intgrlvec = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
            thread_data[tid].plijkz = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);

            // every worker malloc a private yvec array
            thread_data[tid].yvec[0] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
            thread_data[tid].yvec[1] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
            thread_data[tid].yvec[2] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
            thread_data[tid].yvec[3] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);

            // shared ready-only data for all workers
            thread_data[tid].xsize = xsize_spl;
            thread_data[tid].xvec = xvec_spl;
            thread_data[tid].xdiff = xdiff_spl;
            thread_data[tid].ixdiff = ixdiff_spl;

            // zero the memory
            int m, n;
            for (m=0; m<gamma_size; m++)
                for (n=0; n<gamma_size; n++)
                    thread_data[tid].mvec[m*gamma_size_pad+n]=0.0;
            int i;
            for (i = 0; i < maxblocksize * gamma_size_pad; i++)
            {
                thread_data[tid].intgrlvec[i] = 0.0;
                thread_data[tid].plijkz[i] = 0.0;
            }
            int j;
            for (j = 0; j < 4; j++)
            {
                for (i = 0; i < xsize_pad; i++)
                {
                    thread_data[tid].yvec[j][i] = 0.0;
                }
            }
        }

        // SJP: For timing accuracy, ensure that all threads start and end sync'd
        double t0, t1;

        t0 = omp_get_wtime();

        // main calculation
        #pragma omp target device(0)
        #pragma omp teams default(none) shared(thread_data, maxblocksize) num_teams(59) thread_limit(4)
        {
            calculate_gamma_3d_nested(thread_data, maxblocksize);
        }

        t1 = omp_get_wtime();

        printf("[%d]\tgamma_3D \t%d\t%e\n", rank, n, t1-t0);

        // SJP: Just do this as a critical for now.  Swap for an array reduction later.
        t0 = omp_get_wtime();

        #pragma omp target device(0)
        {
            int gamma_size_pad = (gamma_size+7) & ~7;
            for (tid = 0; tid < nteams; tid++)
            {
                double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat; // I think this is necessary because host ptr != device ptr.
                for (m = 0; m < gamma_size; m++)
                {
                    for (n = 0; n < gamma_size; n++)
                    {
                        gamma[m][n] += thread_data[tid].mvec[m*gamma_size_pad+ n];
                    }
                }
            }
        }

        t1 = omp_get_wtime();

        printf("[%d]\treduction\t%d\t%e\n", rank, n, t1-t0);

        /* JB free spline data */
        #pragma omp target device(0)
        {
            for (tid = 0; tid < nteams; tid++)
            {
                _mm_free(thread_data[tid].mvec);
                _mm_free(thread_data[tid].intgrlvec);
                _mm_free(thread_data[tid].plijkz);
                _mm_free(thread_data[tid].yvec[0]);
                _mm_free(thread_data[tid].yvec[1]);
                _mm_free(thread_data[tid].yvec[2]);
                _mm_free(thread_data[tid].yvec[3]);
            }
            // These are managed by the runtime now.
            /*_mm_free(xvec_spl);
            _mm_free(xdiff_spl);
            _mm_free(ixdiff_spl);*/
        }

    }

}
#endif

void gamma_3d_offload_nested(double * restrict gamma_flat, int gamma_size, int rank, int numranks, double workshare)
{
    int i, j;
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    int xsize_spl = get_b_xsize();
    int xsize_pad = (xsize_spl + 7) & ~7;
    double* xvec_spl = (double*) _mm_malloc(xsize_pad * sizeof(double) + 64, 64);
    get_b_xvec(xvec_spl);

    // SJP: Pre-compute the xdiff and inverse xdiff for the inline spline calculation.
    double* xdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
    double* ixdiff_spl = (double*) _mm_malloc(xsize_spl * sizeof(double), 64);
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
    shared(rank, numranks, gamma, gamma_size, xsize_spl, xsize_pad, xvec_spl, xdiff_spl, ixdiff_spl, thread_data, core_barriers, workshare) num_threads(59)
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        int nworkers = 4;
        int maxblocksize = MAXBLOCKSIZE;

        // JB allocate the thread data for all threads
        #pragma omp barrier
        #pragma omp master
        {
            thread_data = (thread_g3d_data_t*) _mm_malloc(nthreads * sizeof(thread_g3d_data_t),64);
        }
        #pragma omp barrier

        // worker master threads to allocate task-shared work arrays
        int gamma_size_pad = (gamma_size+7) & ~7;
        thread_data[tid].mvec = (double*) _mm_malloc(gamma_size_pad * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].intgrlvec = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);
        thread_data[tid].plijkz = (double*) _mm_malloc(maxblocksize * gamma_size_pad * sizeof(double), 64);

        // every worker malloc a private yvec array
        thread_data[tid].yvec[0] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
        thread_data[tid].yvec[1] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
        thread_data[tid].yvec[2] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);
        thread_data[tid].yvec[3] = (double*) _mm_malloc(xsize_pad * sizeof(double), 64);

        // shared ready-only data for all workers
        thread_data[tid].xsize = xsize_spl;
        thread_data[tid].xvec = xvec_spl;
        thread_data[tid].xdiff = xdiff_spl;
        thread_data[tid].ixdiff = ixdiff_spl;

        // zero the memory
        for (m=0; m<gamma_size; m++)
            for (n=0; n<gamma_size; n++)
                thread_data[tid].mvec[m*gamma_size_pad+n]=0.0;

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
        calculate_gamma_3d_nested(thread_data,maxblocksize);

        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\tgamma_3D \t%d\t%e\n", rank, n, t1-t0);

        // SJP: Just do this as a critical for now.  Swap for an array reduction later.
        #pragma omp barrier
        t0 = omp_get_wtime();
        #pragma omp critical
        {
            int gamma_size_pad = (gamma_size+7) & ~7;
            for (m = 0; m < gamma_size; m++)
            {
                for (n = 0; n < gamma_size; n++)
                {
                    gamma[m][n] += thread_data[tid].mvec[m*gamma_size_pad + n];
                }
            }
        }
        #pragma omp barrier
        t1 = omp_get_wtime();
        #pragma omp master
        printf("[%d]\treduction\t%d\t%e\n", rank, n, t1-t0);

        _mm_free(thread_data[tid].mvec);
        _mm_free(thread_data[tid].intgrlvec);
        _mm_free(thread_data[tid].plijkz);
        _mm_free(thread_data[tid].yvec[0]);
        _mm_free(thread_data[tid].yvec[1]);
        _mm_free(thread_data[tid].yvec[2]);
        _mm_free(thread_data[tid].yvec[3]);
    }
    _mm_free(xvec_spl);
    _mm_free(xdiff_spl);
    _mm_free(ixdiff_spl);
    _mm_free(thread_data);
    printf("offload finished...\n");
}
