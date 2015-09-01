#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <gsl/gsl_spline.h>
#include <mkl.h>
#include "global.h"

static const double sixth = (1.0/6.0);
//static gsl_interp* xint_interp = NULL;
//static gsl_interp_accel* xint_acc = NULL;
double calculate_xint(int l1, int l2, int l3, int n, int xsize, double* x,    double* y, DFTaskPtr task) {

    double result;
    int i,j,k;
    int a1,a2,a3;
    find_perm_prim(n,&a1,&a2,&a3);

    int size = xsize;

    double min = x[0];
    double max = x[size-1];

    int pmax = get_pmax_prim();
    double* restrict beta_flat = get_beta_array();
    double (*restrict beta)[pmax+1][xsize] = (double (*restrict)[pmax+1][xsize])  beta_flat;  
    for (i=0; i < size; i++) {
                
        double p1 = beta[l1][a1][i];
        double p2 = beta[l1][a2][i];
        double p3 = beta[l1][a3][i];
        double p4 = beta[l2][a1][i];
        double p5 = beta[l2][a2][i];
        double p6 = beta[l2][a3][i];
        double p7 = beta[l3][a1][i];
        double p8 = beta[l3][a2][i];
        double p9 = beta[l3][a3][i];

        double t1 = p1*p5*p9;
        double t2 = p2*p6*p7;
        double t3 = p3*p4*p8;
        double t4 = p3*p5*p7;
        double t5 = p2*p4*p9;
        double t6 = p1*p6*p8;

        double _y = x[i] * x[i] * (t1+t2+t3+t4+t5+t6) * sixth;
        if (i < 10 && fabs(_y) > 1e-18) _y = 0;
        y[i] = _y;

    }

    double *llim = &(x[0]);
    double *rlim = &(x[size-1]);
    int status;
    status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
    status = dfdIntegrate1D(task, DF_METHOD_PP, 1, llim, DF_NO_HINT, rlim, \
                            DF_NO_HINT, 0, 0, &result, DF_NO_HINT);

    return result;
}

