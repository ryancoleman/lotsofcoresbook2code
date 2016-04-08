#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#include "global.h"

// original l1 upper bound
#ifndef L1_CUTOFF
#define L1_CUTOFF (lmax+1)
#endif

void gamma_3d_host(double * restrict gamma_flat, int gamma_size, int lmax)
{
    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);


    // decompose loop amongst mpi tasks
    //int gamma_3D = (lmax-1)*gamma_size;
    int gamma_3D = (L1_CUTOFF-2)*gamma_size;
    int gamma_pairs_3D[gamma_3D][2];

    int m=0;
    int i,j,n,r;
    int l1;
    for(n=0;n<gamma_size;n++){
        for(l1=2;l1<L1_CUTOFF;l1++){
            gamma_pairs_3D[m][0] = n;
            gamma_pairs_3D[m][1] = l1;
            m++;
        }   
    }   
    int loops=gamma_3D/nproc;
    int auxloop = fmod(gamma_3D,nproc);
    int start_loop = rank*loops;
    int end_loop = (rank+1)*loops-1;

    if (auxloop != 0){ 
        if (rank < auxloop){
            start_loop = start_loop + rank;
            end_loop = end_loop + rank + 1;
        }else{
            start_loop = start_loop + auxloop;
            end_loop = end_loop + auxloop;
        }   
    }   

    // collect partial result in mvec
    double *mvec = (double *) malloc(gamma_size * sizeof(double));

    // main loop
    for(n=start_loop;n<end_loop+1;n++)
    {
        m = gamma_pairs_3D[n][0];
        l1 = gamma_pairs_3D[n][1];

        calculate_gamma_3D(m,l1,mvec);

        // write partial result to gamma matrix
        for (r=0;r<gamma_size;r++) {
            gamma_flat[r*gamma_size+m] += mvec[r];
        }
    }
    free(mvec);
    printf("[%d] finished\n", rank);
}

