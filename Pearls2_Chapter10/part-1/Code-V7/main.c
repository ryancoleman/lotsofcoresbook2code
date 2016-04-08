#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
//#include <mkl.h>

#include "global.h"

int main( int argc, char *argv[] ){

    // set up timing
    double time1, time2, time3, time4, time5, duration;
    double time_total;

    if (argc != 2) {
        printf ( "**** Incorrect number of arguments	****\n" );
        printf( "Usage is >:%s dir\n", argv[0] );
        printf ( "**** Program terminated ****\n" );
        exit (1);
    }

    #ifdef _OPENMP
    time1 = omp_get_wtime();
    #else
    time1 = csecond();
    #endif
    time_total = time1;

    // mpi vars
    int rank, nproc, lnext;

    // mpi init
    int required=MPI_THREAD_SERIALIZED;
    int provided;
    MPI_Init_thread(&argc, &argv, required, &provided); // jb
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // jb check the threading support level
    if (provided < required)
    {
        // insufficient support, degrade to 1 hread and warn the user
        if (rank==0)
        {
            printf("Warning: This MPI implementation provides insufficient threading support.\n");
        }
        omp_set_num_threads(1);
    }

    // read in ini file
    char inifile[MAXLEN];
    strcpy(inifile, argv[1]);
    initilise(inifile);

    // sets the dimension of the matrix and loads the ordering files	
    set_terms_prim();
    set_terms_late();

    int i,j,k,r,s,t,l,m,n;

    // reads in the projected primordial basis functions \tilde{q}
    read_beta();

    int terms = get_terms_prim();
    int lsize = get_b_lsize();
    int *lvec = create_ivector(lsize);
    get_b_lvec(lvec);
    int lmax = lvec[lsize-1];
    int xsize = lmax+1;

    // loads in the cls and other bits and bobs
    init_lmax(lmax);
    create_cl(lsize);
    create_beam(lsize);
    create_noise(lsize);
    create_t_wgt(lsize);
    create_lens(lsize);
    load_cl(lsize, lvec);
    load_BN(lsize, lvec);
    load_TL(lsize, lvec);
    load_lens(lsize, lvec);

    double sum1,sum2,sum3,sum4;
    double x1,x2,x3,x4;

    int ortho_size = terms;

    if(rank==0){
        printf("lmax: %d\n", (int)lmax);
        printf("prim pmax: %d\n", get_pmax_prim());
        printf("late pmax: %d\n", get_pmax_late());
        printf("prim terms: %d\n", get_terms_prim());
        printf("late terms: %d\n", get_terms_late());
        printf("xsize: %d\n",get_b_xsize());
        for(n=0;n<ortho_size;n++){
            find_perm_prim(n,&i,&j,&k);
            find_perm_late(n,&r,&s,&t);
            //printf("%d\t(%d,%d,%d)\t(%d,%d,%d)\n",n,i,j,k,r,s,t);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    double xvec[lmax+1];
    for(i=0;i<xsize;i++){
        xvec[i] = (double)i;
    }
    // creates the late time basis Q
    double xmax = (double)lmax;
    create_basis_late(xsize, xmax, xvec);

    // paralelisation variables
    int loops;
    int auxloop ;
    int start_loop;
    int end_loop;

    MPI_Barrier(MPI_COMM_WORLD);

    // reads in orthogonalisation arrays 
    read_orthol();
    read_lambdal();
    double **orthoinv = (double **)create_array(ortho_size,ortho_size);

    for(i=0;i<ortho_size;i++){
        for(j=0;j<ortho_size;j++){
            orthoinv[i][j]=0.0;
            for(k=0;k<ortho_size;k++){
                orthoinv[i][j] += get_lambdal(k,i)*get_lambdal(k,j);
            }
        }
    }

    create_gamma();
    int gamma_size = terms;
    int gamma_total = terms*terms;
    double *gamma_flat = (double *) malloc(gamma_total*sizeof(double));
    double (*restrict gamma)[gamma_size] = (double (*restrict)[gamma_size]) gamma_flat;

    /*
     * GAMMA3D CALCULATION
     */

    int numranks = nproc;
    MPI_Barrier(MPI_COMM_WORLD);

    #ifdef __INTEL_OFFLOAD
    double workshare_offload = 1.0;
    #else
    double workshare_offload = 0.0;
    #endif

    // Precomputation
    precompute_gamma_3d();

    // Main calculation.
    if (rank==0)
    {
        #ifdef __INTEL_OFFLOAD
        printf("Running on device.\n");
        #else
        printf("Running on host.\n");
        #endif
    }
    double t2 = omp_get_wtime();
    int offload_target =0;
    #pragma offload_transfer target(mic:offload_target) in(gamma_flat[0:gamma_total] : ALLOC RETAIN)
    #pragma offload target(mic:offload_target) in(rank, numranks, gamma_size, workshare_offload) \
    out(gamma_flat[0:gamma_total] : REUSE FREE)
    {
        gamma_3d_offload(gamma_flat, gamma_size, rank, numranks, workshare_offload);
    }
    double thost = omp_get_wtime() - t2;
    printf("[%d] finished\n", rank);

    // now we're done and collecting results into one big array
    double* gamma_send = (double *)malloc( gamma_total * sizeof(double) );
    double* gamma_recv = (double *)malloc( gamma_total * sizeof(double) );

    MPI_Barrier(MPI_COMM_WORLD);

    n=0;
    for (i=0;i<gamma_size;i++) {
        for (j=0;j<gamma_size;j++) {
            gamma_send[n] = gamma[i][j];
            n++;
        }
    }

    MPI_Reduce(gamma_send,gamma_recv,gamma_total,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);

    #ifdef _OPENMP
    duration = omp_get_wtime() - time_total;
    #else
    duration = csecond() - time_total;
    #endif

    printf("[%d] Total Time = %f\n", rank, duration);
    //return 0;

    if(rank==0)
    {
        puts("Writing gamma output to ./Output/rmpi_V7.output...");
        FILE *foutput = fopen("./Output/rmpi_V7.output", "w");

        if (foutput == NULL)
        {
            perror("Failed to open file for output. Exiting.");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        n=0;
        for (i=0;i<gamma_size;i++) {
            for (j=0;j<gamma_size;j++) {
                gamma[i][j] = gamma_recv[n];
                n++;
                fprintf(foutput, "%.15e\t",gamma[i][j]);
            }
            fprintf(foutput, "\n");
        }
        // doing last matrix multiplication then outputing the result.

        double **gamma_temp = (double **)create_array(gamma_size,gamma_size);

        for (i=0;i<gamma_size;i++) {
            for (j=0;j<gamma_size;j++) {
                gamma_temp[i][j] = 0.0;
                for (n=0;n<gamma_size;n++) {
                    gamma_temp[i][j] +=orthoinv[i][n]*gamma[n][j];
                }
            }
        }

        fprintf(foutput, "\n");

        double *results_g =  malloc( sizeof(double)*3);
        results_g[0] = 0.0;
        results_g[1] = 0.0;
        results_g[2] = 0.0;

        for (i=0;i<gamma_size;i++) {
            for (j=0;j<gamma_size;j++) {
                results_g[0] = (double)i;
                results_g[1] = (double)j;
                results_g[2] = gamma_temp[i][j];
                update_gamma(results_g);
                // 				printf("Gamma:\t%d\t%d\t%e\n", i, j, results_g[2]);
                fprintf(foutput, "%e\t",results_g[2]);
            }
            fprintf(foutput, "\n");
        }
        output_gamma();

        for (n=0;n<gamma_size;n++) {
            for (i=0;i<gamma_size;i++) {
                gamma[i][n] = 0.0;
                for (j=i;j<gamma_size;j++) {
                    gamma[i][n] += gamma_temp[j][n]*get_orthol(j,i);
                }
                fprintf(foutput, "%e\t",gamma[i][n]);
            }
            fprintf(foutput, "\n");
        }

        fprintf(foutput, "\n");
        for (i=0;i<gamma_size;i++) {
            for (j=0;j<gamma_size;j++) {
                x1=x2=x3=0.0;
                for (n=0;n<gamma_size;n++) {
                    x1 += gamma[n][i]*gamma[n][i];
                    x2 += gamma[n][j]*gamma[n][j];
                    x3 += gamma[n][i]*gamma[n][j];
                }
                x4 = x3 / sqrt(x1*x2);
                fprintf(foutput, "%e\t",x4);
            }
            fprintf(foutput, "\n");
        }
        fclose(foutput);
    }

    MPI_Finalize();
    return 0;
}
