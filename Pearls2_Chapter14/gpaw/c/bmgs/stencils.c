/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2005       CSC - IT Center for Science Ltd.
 *  Please see the accompanying LICENSE file for further information. */

#include <stdlib.h>
#include "bmgs.h"


// Expansion coefficients for finite difference Laplacian.  The numbers are
// from J. R. Chelikowsky et al., Phys. Rev. B 50, 11355 (1994):

bmgsstencil bmgs_stencil(int ncoefs, const double* coefs, const long* offsets,
			 int r, const long n[3])
{
  bmgsstencil stencil = 
    {ncoefs,
     (double*)malloc(ncoefs * sizeof(double)),
     (long*)malloc(ncoefs * sizeof(long)),
     {n[0], n[1], n[2]},
     {2 * r * (n[2] + 2 * r) * (n[1] + 2 * r),
     2 * r * (n[2] + 2 * r),
     2 * r}};
  assert((stencil.coefs != NULL) && (stencil.offsets != NULL));
  memcpy(stencil.coefs, coefs, ncoefs * sizeof(double));
  memcpy(stencil.offsets, offsets, ncoefs * sizeof(long));
  return stencil;
}


static const double laplace[4][5] = 
  {{-2.0,        1.0,      0.0,      0.0,        0.0},
   {-5.0/2.0,    4.0/3.0, -1.0/12.0, 0.0,        0.0},
   {-49.0/18.0,  3.0/2.0, -3.0/20.0, 1.0/90.0,   0.0},
   {-205.0/72.0, 8.0/5.0, -1.0/5.0,  8.0/315.0, -1.0/560.0}};


bmgsstencil bmgs_laplace(int k, double scale, 
				const double h[3],
				const long n[3])
{
  int ncoefs = 3 * k - 2;
  double* coefs = (double*)malloc(ncoefs * sizeof(double));
  long* offsets = (long*)malloc(ncoefs * sizeof(long));
  assert((coefs != NULL) && (offsets != NULL));
  double f1 = 1.0 / (h[0] * h[0]);
  double f2 = 1.0 / (h[1] * h[1]);
  double f3 = 1.0 / (h[2] * h[2]);
  int r = (k - 1) / 2;   // range
  double s[3] = {(n[2] + 2 * r) * (n[1] + 2 * r), n[2] + 2 * r, 1};
  int m = 0;
  for (int j = 1; j <= r; j++)
    {
      double c = scale * laplace[r - 1][j];
      coefs[m] = c * f1; offsets[m++] = -j * s[0];
      coefs[m] = c * f1; offsets[m++] = +j * s[0];
      coefs[m] = c * f2; offsets[m++] = -j * s[1];
      coefs[m] = c * f2; offsets[m++] = +j * s[1];
      coefs[m] = c * f3; offsets[m++] = -j;
      coefs[m] = c * f3; offsets[m++] = +j;
    }
  double c = scale * laplace[r - 1][0];
  coefs[m] = c * (f1 + f2 + f3); offsets[m] = 0;
  bmgsstencil stencil = 
    {ncoefs, coefs, offsets,
     {n[0], n[1], n[2]},
     {2 * r * (n[2] + 2 * r) * (n[1] + 2 * r),
     2 * r * (n[2] + 2 * r),
     2 * r}};
  return stencil;
}


bmgsstencil bmgs_mslaplaceA(double scale, 
				   const double h[3],
				   const long n[3])
{
  int ncoefs = 19;
  double* coefs = (double*)malloc(ncoefs * sizeof(double));
  long* offsets = (long*)malloc(ncoefs * sizeof(long));
  assert((coefs != NULL) && (offsets != NULL));
  double e[3]  = {-scale / (12.0 * h[0] * h[0]),
		  -scale / (12.0 * h[1] * h[1]),
		  -scale / (12.0 * h[2] * h[2])};
  double f = -16.0 * (e[0] + e[1] + e[2]);
  double g[3] = {10.0 * e[0] + 0.125 * f,
		 10.0 * e[1] + 0.125 * f,
		 10.0 * e[2] + 0.125 * f};
  double s[3] = {(n[2] + 2) * (n[1] + 2), n[2] + 2, 1};
  int m = 0;
  coefs[m] = f;
  offsets[m++] = 0;
  for (int j = -1; j <= 1; j += 2)
    {
      coefs[m] = g[0];
      offsets[m++] = j * s[0];
      coefs[m] = g[1];
      offsets[m++] = j * s[1];
      coefs[m] = g[2];
      offsets[m++] = j * s[2];
    }
  for (int j = -1; j <= 1; j += 2)
    for (int k = -1; j <= 1; j += 2)
      {
	coefs[m] = e[1] + e[2];
	offsets[m++] = -j * s[1] - k * s[2];
	coefs[m] = e[0] + e[2];
	offsets[m++] = -j * s[0] - k * s[2];
	coefs[m] = e[0] + e[1];
	offsets[m++] = -j * s[0] - k * s[1];
      }
  bmgsstencil stencil = 
    {ncoefs, coefs, offsets,
     {n[0], n[1], n[2]},
     {2 * s[0], 2 * s[1], 2}};
  return stencil;
}


bmgsstencil bmgs_mslaplaceB(const long n[3])
{
  int ncoefs = 7;
  double* coefs = (double*)malloc(ncoefs * sizeof(double));
  long* offsets = (long*)malloc(ncoefs * sizeof(long));
  assert((coefs != NULL) && (offsets != NULL));
  double s[3] = {(n[2] + 2) * (n[1] + 2), n[2] + 2, 1};
  int k = 0;
  coefs[k] = 0.5;
  offsets[k++] = 0;
  for (int j = -1; j <= 1; j += 2)
    {
      coefs[k] = 1.0 / 12.0;
      offsets[k++] = j * s[0];
      coefs[k] = 1.0 / 12.0;
      offsets[k++] = j * s[1];
      coefs[k] = 1.0 / 12.0;
      offsets[k++] = j * s[2];
    }
  bmgsstencil stencil = 
    {ncoefs, coefs, offsets,
     {n[0], n[1], n[2]},
     {2 * s[0], 2 * s[1], 2}};
  return stencil;
}


bmgsstencil bmgs_gradient(int k, int i, double h, 
				 const long n[3])
{
  int ncoefs = k - 1;
  double* coefs = (double*)malloc(ncoefs * sizeof(double));
  long* offsets = (long*)malloc(ncoefs * sizeof(long));
  assert((coefs != NULL) && (offsets != NULL));
  int r = 1;
  double s[3] = {(n[2] + 2 * r) * (n[1] + 2 * r), n[2] + 2 * r, 1};
  double c = 0.5 / h;
  coefs[0] = +c; offsets[0] = +s[i];
  coefs[1] = -c; offsets[1] = -s[i];
  bmgsstencil stencil = 
    {ncoefs, coefs, offsets,
     {n[0], n[1], n[2]},
     {2 * r * s[0], 2 * r * s[1], 2 * r}};
  return stencil;
}

void bmgs_deletestencil(bmgsstencil* stencil)
{
  free(stencil->coefs);
  free(stencil->offsets);
}
