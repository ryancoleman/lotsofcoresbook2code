#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include "global.h"

double calculate_xint(int l1, int l2, int l3, int n, int xsize, double *xvec) {

	double result;
	int i,j,k;
	int a1,a2,a3;
	find_perm_prim(n,&a1,&a2,&a3);
	
	int size = xsize;
	double x[size];
	double y[size];
	
	for(i=0;i<size;i++){
		x[i] = xvec[i];
	}
	
	double min = x[0];
	double max = x[size-1];
	
    for (i=0;i<size;i++)
    {
        double p1,p2,p3,p4,p5,p6,p7,p8,p9;
        double t1,t2,t3,t4,t5,t6;

        p1 = get_beta(l1,a1,i);
        p2 = get_beta(l1,a2,i);
        p3 = get_beta(l1,a3,i);
        p4 = get_beta(l2,a1,i);
        p5 = get_beta(l2,a2,i);
        p6 = get_beta(l2,a3,i);
        p7 = get_beta(l3,a1,i);
        p8 = get_beta(l3,a2,i);
        p9 = get_beta(l3,a3,i);

        t1 = p1*p5*p9;
		t2 = p2*p6*p7;
		t3 = p3*p4*p8;
		t4 = p3*p5*p7;
		t5 = p2*p4*p9;
		t6 = p1*p6*p8;
		
		y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6)/6.0;
		if( i<10 && fabs(y[i])>1e-18 )y[i]=0;
// 		if(l1==2&&l2==2000)printf("%d\t%d\t%e\t%e\n",n,i,x[i],y[i]);
	}
	
	gsl_spline* sp =  gsl_spline_alloc (gsl_interp_cspline, size);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	
	gsl_spline_init(sp,x,y,size);
	result = gsl_spline_eval_integ(sp,min,max,acc);
	
	gsl_spline_free(sp);
	gsl_interp_accel_free(acc);
	
// 	return result*pow(kmax,-3);
	return result;
}
