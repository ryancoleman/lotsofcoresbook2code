#include <math.h>
//#include <mkl.h>
#include "global.h"

_OFFLOADABLE
inline double Delta(double* xdiff, int i)
{
	return xdiff[i];
}
_OFFLOADABLE
inline double delta(double* restrict xdiff, double* restrict ixdiff, double* restrict y, int i)
{
	return (y[i+1] - y[i]) * ixdiff[i];
}
_OFFLOADABLE
inline double deriv(double* restrict xdiff, double* restrict ixdiff, double* restrict y, int i)
{
        //return (Delta(xdiff,i-1)*delta(xdiff, ixdiff,y,i) + Delta(xdiff,i)*delta(xdiff,ixdiff,y,i-1)) / (Delta(xdiff,i-1) + Delta(xdiff,i));
	return (y[i+1] - y[i]) * ixdiff[i];
}

#define TRAPEZIUM
//#define HERMITE
_OFFLOADABLE static const double sixth = (1.0/6.0);
_OFFLOADABLE double integrate(int size, double* restrict xdiff, double* restrict ixdiff, double* restrict y)
{

#ifdef HERMITE
	__assume_aligned(xdiff, 64);
	__assume_aligned(y, 64);
	__assume_aligned(ixdiff, 64);

	int i;
	double r = 0;
	#pragma loop count (213)
	for (i = 1; i < size-2; i++) {
		r += 0.5 * xdiff[i] * (y[i] + y[i+1] + xdiff[i] * (deriv(xdiff,ixdiff,y,i) - deriv(xdiff,ixdiff,y,i+1)) * sixth);
	}
	return r;
#endif

#ifdef TRAPEZIUM
	int i;
	double r = 0;
	for (i = 0; i < size-1; i++)
    {
		r += (xdiff[i]) * (y[i] + y[i+1]) * 0.5;
	}
	return r;
#endif
}

_OFFLOADABLE double calculate_xint(int n, int xsize, double* restrict x, double* restrict y, int l1, int l2, int l3, double* restrict xdiff, double* restrict ixdiff)
{
	int i;
	int a1, a2, a3;
	int xsize_pad = (xsize + 7) & ~7;
	int pmax = get_pmax_prim();

	double* restrict beta_flat = get_beta_array();
	double (*restrict beta)[pmax+1][xsize_pad] = (double (*restrict)[pmax+1][xsize_pad]) beta_flat;

    find_perm_prim(n,&a1,&a2,&a3);

    #pragma loop count (216)
    #pragma prefetch beta:0:1
    #pragma prefetch y:0:1
    #pragma prefetch x:0:1
    #pragma vector aligned
    for (i = 0; i < xsize; i++) {

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

        y[i] = x[i] * x[i] * (t1+t2+t3+t4+t5+t6) * sixth;
    }

    for (i = 0; i < 8; i++)
        if (fabs(y[i]) > 1e-18) y[i] = 0;

    return integrate(xsize, xdiff, ixdiff, y);
}
