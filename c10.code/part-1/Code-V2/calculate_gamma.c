#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

static double* svec;
void precompute_gamma_3d() {

	int lsize = get_lmax()+1;
	double* restrict cl = get_cl_array();
	double* restrict noise = get_noise_array();
	double* restrict beam = get_beam_array();

	svec = (double*) malloc(lsize*sizeof(double));
	if (!svec) printf("ERROR: Couldn't malloc svec.\n");
	int l;
	const double third = 1.0/3.0;
	for (l = 0; l < lsize; l++) {
		svec[l] = pow(2.0*l+1.0,third)*(cl[l]+noise[l]/(beam[l]*beam[l]));
	}
}

void calculate_gamma_3d(int n, int i, double *mvec) {
	
	int j,k,m,t1,t2,t3;
	double x,y,z;
	int terms = get_terms_prim();
	int lsize = get_lmax()+1;
	int lmax = get_lmax();
	double s1,s2,s3;
	
    int xsize = get_b_xsize();
    double *xvec = create_vector(xsize);
    get_b_xvec(xvec);
    double xmax = xvec[xsize-1];
	
    for(m=0;m<terms;m++) mvec[m] = 0.0;

    s1 = svec[i];
    for(j=i; j<lsize; j++)
    {   
        s2 = svec[j];
        for(k=j+i%2; k<min(i+j,lmax)+1; k+=2)
        {
            x = calculate_xint(i,j,k,n,xsize,xvec);
            s3 = svec[k];
            z = permsix(i,j,k)*calculate_geometric(i,j,k)/sqrt(s1 * s2 * s3);

            for(m=0;m<terms;m++)
            {
                y = plijk(m,i,j,k);
                mvec[m] += x*y*z;
            }
        }
    }
    for(m=0;m<terms;m++) 
        mvec[m] *= 6.0*deltaphi*deltaphi;
	return;
}
