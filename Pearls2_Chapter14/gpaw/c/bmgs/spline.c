/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "bmgs.h"


bmgsspline bmgs_spline(int l, double dr, int nbins, double* f)
{
  double c = 3.0 / (dr * dr);
  double* f2 = (double*)malloc((nbins + 1) * sizeof(double));
  assert(f2 != NULL);
  double* u = (double*)malloc(nbins * sizeof(double));
  assert(u != NULL);
  f2[0] = -0.5;
  u[0] = (f[1] - f[0]) * c;
  for (int b = 1; b < nbins; b++)
    {
      double p = 0.5 * f2[b - 1] + 2.0;
      f2[b] = -0.5 / p;
      u[b] = ((f[b + 1] - 2.0 * f[b] + f[b - 1]) * c - 0.5 * u[b - 1]) / p;
    }
  f2[nbins] = ((f[nbins - 1] * c - 0.5 * u[nbins - 1]) /
               (0.5 * f2[nbins - 1] + 1.0));
  for (int b = nbins - 1; b >= 0; b--)
    f2[b] = f2[b] * f2[b + 1] + u[b];
  double* data = (double*)malloc(4 * (nbins + 1) * sizeof(double));
  assert(data != NULL);
  bmgsspline spline = {l, dr, nbins, data};
  for (int b = 0; b < nbins; b++)
    {
      *data++ = f[b];
      *data++ = (f[b + 1] - f[b]) / dr - (f2[b] / 3 + f2[b + 1] / 6) * dr;
      *data++ = 0.5 * f2[b];
      *data++ = (f2[b + 1] - f2[b]) / (6 * dr);
    }
  data[0] = 0.0;
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 0.0;
  free(u);
  free(f2);
  return spline;
}


double bmgs_splinevalue(const bmgsspline* spline, double r)
{
  int b = r / spline->dr;
  if (b >= spline->nbins)
    return 0.0;
  double u = r - b * spline->dr;
  double* s = spline->data + 4 * b;
  return  s[0] + u * (s[1] + u * (s[2] + u * s[3]));
}


void bmgs_get_value_and_derivative(const bmgsspline* spline, double r,
				   double *f, double *dfdr)
{
  int b = r / spline->dr;
  if (b >= spline->nbins)
    {
      *f = 0.0;
      *dfdr = 0.0;
      return;
    }
  double u = r - b * spline->dr;
  double* s = spline->data + 4 * b;
  *f = s[0] + u * (s[1] + u * (s[2] + u * s[3]));
  *dfdr = s[1] + u * (2.0 * s[2] + u * 3.0 * s[3]);
}


void bmgs_deletespline(bmgsspline* spline)
{
  free(spline->data);
}


void bmgs_radial1(const bmgsspline* spline, 
		  const int n[3], const double C[3],
		  const double h[3],
		  int* b, double* d)
{
  int nbins = spline->nbins;
  double dr = spline->dr;
  double x = C[0];
  for (int i0 = 0; i0 < n[0]; i0++)
    {
      double xx = x * x;
      double y = C[1];
      for (int i1 = 0; i1 < n[1]; i1++)
	{
	  double xxpyy = xx + y * y;
	  double z = C[2];
	  for (int i2 = 0; i2 < n[2]; i2++)
	    {
	      double r = sqrt(xxpyy + z * z);
	      int j = r / dr;
	      if (j < nbins)
		{
		  *b++ = j;
		  *d++ = r - j * dr;
		}
	      else
		{
		  *b++ = nbins;
		  *d++ = 0.0;
		}
	      z += h[2];
	    }
	  y += h[1];
	}
      x += h[0];
    }
}


void bmgs_radial2(const bmgsspline* spline, const int n[3],
		  const int* b, const double* d, 
		  double* f, double* g)
{
  double dr = spline->dr;
  for (int q = 0; q < n[0] * n[1] * n[2]; q++)
    {
      int j = b[q];
      const double* s = spline->data + 4 * j;
      double u = d[q];
      f[q] = s[0] + u * (s[1] + u * (s[2] + u * s[3]));
      if (g != 0)
	{
	  if (j == 0)
	    g[q] = 2.0 * s[2] + u * 3.0 * s[3];
	  else
	    g[q] = (s[1] + u * (2.0 * s[2] + u * 3.0 * s[3])) / (j * dr + u);
	}
    }
}

//Computer generated code! Hands off!
    
// inserts values of f(r) r^l Y_lm(theta, phi) in elements of input array 'a'
void bmgs_radial3(const bmgsspline* spline, int m, 
		  const int n[3], 
		  const double C[3],
		  const double h[3],
		  const double* f, double* a)
{
  int l = spline->l;
  if (l == 0)
    for (int q = 0; q < n[0] * n[1] * n[2]; q++)
      a[q] = 0.28209479177387814 * f[q];
  else if (l == 1)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == -1)
                    a[q] = f[q] * 0.48860251190291992 * y;
                  else if (m == 0)
                    a[q] = f[q] * 0.48860251190291992 * z;
                  else
                    a[q] = f[q] * 0.48860251190291992 * x;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (l == 2)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -2)
                    a[q] = f[q] * 1.0925484305920792 * x*y;
                  else if (m == -1)
                    a[q] = f[q] * 1.0925484305920792 * y*z;
                  else if (m == 0)
                    a[q] = f[q] * 0.31539156525252005 * (3*z*z-r2);
                  else if (m == 1)
                    a[q] = f[q] * 1.0925484305920792 * x*z;
                  else
                    a[q] = f[q] * 0.54627421529603959 * (x*x-y*y);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (l == 3)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -3)
                    a[q] = f[q] * 0.59004358992664352 * (-y*y*y+3*x*x*y);
                  else if (m == -2)
                    a[q] = f[q] * 2.8906114426405538 * x*y*z;
                  else if (m == -1)
                    a[q] = f[q] * 0.45704579946446577 * (-y*r2+5*y*z*z);
                  else if (m == 0)
                    a[q] = f[q] * 0.3731763325901154 * (5*z*z*z-3*z*r2);
                  else if (m == 1)
                    a[q] = f[q] * 0.45704579946446577 * (5*x*z*z-x*r2);
                  else if (m == 2)
                    a[q] = f[q] * 1.4453057213202769 * (x*x*z-y*y*z);
                  else
                    a[q] = f[q] * 0.59004358992664352 * (x*x*x-3*x*y*y);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (l == 4)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -4)
                    a[q] = f[q] * 2.5033429417967046 * (x*x*x*y-x*y*y*y);
                  else if (m == -3)
                    a[q] = f[q] * 1.7701307697799307 * (-y*y*y*z+3*x*x*y*z);
                  else if (m == -2)
                    a[q] = f[q] * 0.94617469575756008 * (-x*y*r2+7*x*y*z*z);
                  else if (m == -1)
                    a[q] = f[q] * 0.66904654355728921 * (-3*y*z*r2+7*y*z*z*z);
                  else if (m == 0)
                    a[q] = f[q] * 0.10578554691520431 * (-30*z*z*r2+3*r2*r2+35*z*z*z*z);
                  else if (m == 1)
                    a[q] = f[q] * 0.66904654355728921 * (7*x*z*z*z-3*x*z*r2);
                  else if (m == 2)
                    a[q] = f[q] * 0.47308734787878004 * (-x*x*r2+7*x*x*z*z+y*y*r2-7*y*y*z*z);
                  else if (m == 3)
                    a[q] = f[q] * 1.7701307697799307 * (x*x*x*z-3*x*y*y*z);
                  else
                    a[q] = f[q] * 0.62583573544917614 * (-6*x*x*y*y+x*x*x*x+y*y*y*y);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else
    assert(0 == 1);
}


// insert values of
// d( f(r) * r^l Y_l^m )                           d( r^l Y_l^m )
// --------------------- = g(r) q r^l Y_l^m + f(r) --------------
//        dq                                             dq
// where q={x, y, z} and g(r) = 1/r*(df/dr)
void bmgs_radiald3(const bmgsspline* spline, int m, int c, 
		  const int n[3], 
		  const double C[3],
		  const double h[3],
		  const double* f, const double* g, double* a)
{
  int l = spline->l;
  // x
  if (c == 0 && l == 0)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == 0)
                    a[q] = g[q] * 0.28209479177387814 * x;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 0 && l == 1)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == -1)
                    a[q] = g[q] * 0.48860251190291992 * x*y;
                  else if (m == 0)
                    a[q] = g[q] * 0.48860251190291992 * x*z;
                  else
                    a[q] = g[q] * 0.48860251190291992 * x*x + f[q] * 0.48860251190291992;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 0 && l == 2)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -2)
                    a[q] = g[q] * 1.0925484305920792 * x*x*y + f[q] * 1.0925484305920792 * y;
                  else if (m == -1)
                    a[q] = g[q] * 1.0925484305920792 * x*y*z;
                  else if (m == 0)
                    a[q] = g[q] * 0.31539156525252005 * (3*x*z*z-x*r2) + f[q] * 0.63078313050504009 * -x;
                  else if (m == 1)
                    a[q] = g[q] * 1.0925484305920792 * x*x*z + f[q] * 1.0925484305920792 * z;
                  else
                    a[q] = g[q] * 0.54627421529603959 * (x*x*x-x*y*y) + f[q] * 1.0925484305920792 * x;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 0 && l == 3)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -3)
                    a[q] = g[q] * 0.59004358992664352 * (3*x*x*x*y-x*y*y*y) + f[q] * 3.5402615395598613 * x*y;
                  else if (m == -2)
                    a[q] = g[q] * 2.8906114426405538 * x*x*y*z + f[q] * 2.8906114426405538 * y*z;
                  else if (m == -1)
                    a[q] = g[q] * 0.45704579946446577 * (-x*y*r2+5*x*y*z*z) + f[q] * 0.91409159892893155 * -x*y;
                  else if (m == 0)
                    a[q] = g[q] * 0.3731763325901154 * (5*x*z*z*z-3*x*z*r2) + f[q] * 2.2390579955406924 * -x*z;
                  else if (m == 1)
                    a[q] = g[q] * 0.45704579946446577 * (-x*x*r2+5*x*x*z*z) + f[q] * 0.45704579946446577 * (5*z*z-r2-2*x*x);
                  else if (m == 2)
                    a[q] = g[q] * 1.4453057213202769 * (x*x*x*z-x*y*y*z) + f[q] * 2.8906114426405538 * x*z;
                  else
                    a[q] = g[q] * 0.59004358992664352 * (-3*x*x*y*y+x*x*x*x) + f[q] * 1.7701307697799307 * (x*x-y*y);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 0 && l == 4)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -4)
                    a[q] = g[q] * 2.5033429417967046 * (-x*x*y*y*y+x*x*x*x*y) + f[q] * 2.5033429417967046 * (-y*y*y+3*x*x*y);
                  else if (m == -3)
                    a[q] = g[q] * 1.7701307697799307 * (-x*y*y*y*z+3*x*x*x*y*z) + f[q] * 10.620784618679583 * x*y*z;
                  else if (m == -2)
                    a[q] = g[q] * 0.94617469575756008 * (7*x*x*y*z*z-x*x*y*r2) + f[q] * 0.94617469575756008 * (-y*r2+7*y*z*z-2*x*x*y);
                  else if (m == -1)
                    a[q] = g[q] * 0.66904654355728921 * (-3*x*y*z*r2+7*x*y*z*z*z) + f[q] * 4.0142792613437353 * -x*y*z;
                  else if (m == 0)
                    a[q] = g[q] * 0.10578554691520431 * (-30*x*z*z*r2+3*x*r2*r2+35*x*z*z*z*z) + f[q] * 1.2694265629824517 * (-5*x*z*z+x*r2);
                  else if (m == 1)
                    a[q] = g[q] * 0.66904654355728921 * (-3*x*x*z*r2+7*x*x*z*z*z) + f[q] * 0.66904654355728921 * (7*z*z*z-6*x*x*z-3*z*r2);
                  else if (m == 2)
                    a[q] = g[q] * 0.47308734787878004 * (-x*x*x*r2+x*y*y*r2+7*x*x*x*z*z-7*x*y*y*z*z) + f[q] * 0.94617469575756008 * (-x*x*x+7*x*z*z-x*r2+x*y*y);
                  else if (m == 3)
                    a[q] = g[q] * 1.7701307697799307 * (-3*x*x*y*y*z+x*x*x*x*z) + f[q] * 5.3103923093397913 * (x*x*z-y*y*z);
                  else
                    a[q] = g[q] * 0.62583573544917614 * (-6*x*x*x*y*y+x*y*y*y*y+x*x*x*x*x) + f[q] * 2.5033429417967046 * (-3*x*y*y+x*x*x);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  // y
  else if (c == 1 && l == 0)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == 0)
                    a[q] = g[q] * 0.28209479177387814 * y;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 1 && l == 1)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == -1)
                    a[q] = g[q] * 0.48860251190291992 * y*y + f[q] * 0.48860251190291992;
                  else if (m == 0)
                    a[q] = g[q] * 0.48860251190291992 * y*z;
                  else
                    a[q] = g[q] * 0.48860251190291992 * x*y;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 1 && l == 2)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -2)
                    a[q] = g[q] * 1.0925484305920792 * x*y*y + f[q] * 1.0925484305920792 * x;
                  else if (m == -1)
                    a[q] = g[q] * 1.0925484305920792 * y*y*z + f[q] * 1.0925484305920792 * z;
                  else if (m == 0)
                    a[q] = g[q] * 0.31539156525252005 * (-y*r2+3*y*z*z) + f[q] * 0.63078313050504009 * -y;
                  else if (m == 1)
                    a[q] = g[q] * 1.0925484305920792 * x*y*z;
                  else
                    a[q] = g[q] * 0.54627421529603959 * (-y*y*y+x*x*y) + f[q] * 1.0925484305920792 * -y;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 1 && l == 3)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -3)
                    a[q] = g[q] * 0.59004358992664352 * (3*x*x*y*y-y*y*y*y) + f[q] * 1.7701307697799307 * (x*x-y*y);
                  else if (m == -2)
                    a[q] = g[q] * 2.8906114426405538 * x*y*y*z + f[q] * 2.8906114426405538 * x*z;
                  else if (m == -1)
                    a[q] = g[q] * 0.45704579946446577 * (-y*y*r2+5*y*y*z*z) + f[q] * 0.45704579946446577 * (5*z*z-r2-2*y*y);
                  else if (m == 0)
                    a[q] = g[q] * 0.3731763325901154 * (-3*y*z*r2+5*y*z*z*z) + f[q] * 2.2390579955406924 * -y*z;
                  else if (m == 1)
                    a[q] = g[q] * 0.45704579946446577 * (-x*y*r2+5*x*y*z*z) + f[q] * 0.91409159892893155 * -x*y;
                  else if (m == 2)
                    a[q] = g[q] * 1.4453057213202769 * (-y*y*y*z+x*x*y*z) + f[q] * 2.8906114426405538 * -y*z;
                  else
                    a[q] = g[q] * 0.59004358992664352 * (x*x*x*y-3*x*y*y*y) + f[q] * 3.5402615395598613 * -x*y;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 1 && l == 4)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -4)
                    a[q] = g[q] * 2.5033429417967046 * (x*x*x*y*y-x*y*y*y*y) + f[q] * 2.5033429417967046 * (x*x*x-3*x*y*y);
                  else if (m == -3)
                    a[q] = g[q] * 1.7701307697799307 * (3*x*x*y*y*z-y*y*y*y*z) + f[q] * 5.3103923093397913 * (x*x*z-y*y*z);
                  else if (m == -2)
                    a[q] = g[q] * 0.94617469575756008 * (-x*y*y*r2+7*x*y*y*z*z) + f[q] * 0.94617469575756008 * (-2*x*y*y+7*x*z*z-x*r2);
                  else if (m == -1)
                    a[q] = g[q] * 0.66904654355728921 * (-3*y*y*z*r2+7*y*y*z*z*z) + f[q] * 0.66904654355728921 * (7*z*z*z-3*z*r2-6*y*y*z);
                  else if (m == 0)
                    a[q] = g[q] * 0.10578554691520431 * (3*y*r2*r2-30*y*z*z*r2+35*y*z*z*z*z) + f[q] * 1.2694265629824517 * (y*r2-5*y*z*z);
                  else if (m == 1)
                    a[q] = g[q] * 0.66904654355728921 * (-3*x*y*z*r2+7*x*y*z*z*z) + f[q] * 4.0142792613437353 * -x*y*z;
                  else if (m == 2)
                    a[q] = g[q] * 0.47308734787878004 * (-7*y*y*y*z*z+y*y*y*r2+7*x*x*y*z*z-x*x*y*r2) + f[q] * 0.94617469575756008 * (y*r2+y*y*y-7*y*z*z-x*x*y);
                  else if (m == 3)
                    a[q] = g[q] * 1.7701307697799307 * (x*x*x*y*z-3*x*y*y*y*z) + f[q] * 10.620784618679583 * -x*y*z;
                  else
                    a[q] = g[q] * 0.62583573544917614 * (x*x*x*x*y-6*x*x*y*y*y+y*y*y*y*y) + f[q] * 2.5033429417967046 * (y*y*y-3*x*x*y);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  // z
  else if (c == 2 && l == 0)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == 0)
                    a[q] = g[q] * 0.28209479177387814 * z;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 2 && l == 1)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  if (m == -1)
                    a[q] = g[q] * 0.48860251190291992 * y*z;
                  else if (m == 0)
                    a[q] = g[q] * 0.48860251190291992 * z*z + f[q] * 0.48860251190291992;
                  else
                    a[q] = g[q] * 0.48860251190291992 * x*z;
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 2 && l == 2)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -2)
                    a[q] = g[q] * 1.0925484305920792 * x*y*z;
                  else if (m == -1)
                    a[q] = g[q] * 1.0925484305920792 * y*z*z + f[q] * 1.0925484305920792 * y;
                  else if (m == 0)
                    a[q] = g[q] * 0.31539156525252005 * (3*z*z*z-z*r2) + f[q] * 1.2615662610100802 * z;
                  else if (m == 1)
                    a[q] = g[q] * 1.0925484305920792 * x*z*z + f[q] * 1.0925484305920792 * x;
                  else
                    a[q] = g[q] * 0.54627421529603959 * (x*x*z-y*y*z);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 2 && l == 3)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -3)
                    a[q] = g[q] * 0.59004358992664352 * (-y*y*y*z+3*x*x*y*z);
                  else if (m == -2)
                    a[q] = g[q] * 2.8906114426405538 * x*y*z*z + f[q] * 2.8906114426405538 * x*y;
                  else if (m == -1)
                    a[q] = g[q] * 0.45704579946446577 * (-y*z*r2+5*y*z*z*z) + f[q] * 3.6563663957157262 * y*z;
                  else if (m == 0)
                    a[q] = g[q] * 0.3731763325901154 * (-3*z*z*r2+5*z*z*z*z) + f[q] * 1.1195289977703462 * (3*z*z-r2);
                  else if (m == 1)
                    a[q] = g[q] * 0.45704579946446577 * (5*x*z*z*z-x*z*r2) + f[q] * 3.6563663957157262 * x*z;
                  else if (m == 2)
                    a[q] = g[q] * 1.4453057213202769 * (x*x*z*z-y*y*z*z) + f[q] * 1.4453057213202769 * (x*x-y*y);
                  else
                    a[q] = g[q] * 0.59004358992664352 * (x*x*x*z-3*x*y*y*z);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else if (c == 2 && l == 4)
    {
      int q = 0;
      double x = C[0];
      for (int i0 = 0; i0 < n[0]; i0++)
        {
          double y = C[1];
          for (int i1 = 0; i1 < n[1]; i1++)
            {
              double z = C[2];
	      for (int i2 = 0; i2 < n[2]; i2++, q++)
		{
                  double r2 = x*x+y*y+z*z;
                  if (m == -4)
                    a[q] = g[q] * 2.5033429417967046 * (x*x*x*y*z-x*y*y*y*z);
                  else if (m == -3)
                    a[q] = g[q] * 1.7701307697799307 * (-y*y*y*z*z+3*x*x*y*z*z) + f[q] * 1.7701307697799307 * (-y*y*y+3*x*x*y);
                  else if (m == -2)
                    a[q] = g[q] * 0.94617469575756008 * (-x*y*z*r2+7*x*y*z*z*z) + f[q] * 11.354096349090721 * x*y*z;
                  else if (m == -1)
                    a[q] = g[q] * 0.66904654355728921 * (-3*y*z*z*r2+7*y*z*z*z*z) + f[q] * 2.0071396306718676 * (-y*r2+5*y*z*z);
                  else if (m == 0)
                    a[q] = g[q] * 0.10578554691520431 * (-30*z*z*z*r2+3*z*r2*r2+35*z*z*z*z*z) + f[q] * 1.6925687506432689 * (5*z*z*z-3*z*r2);
                  else if (m == 1)
                    a[q] = g[q] * 0.66904654355728921 * (-3*x*z*z*r2+7*x*z*z*z*z) + f[q] * 2.0071396306718676 * (5*x*z*z-x*r2);
                  else if (m == 2)
                    a[q] = g[q] * 0.47308734787878004 * (-x*x*z*r2+7*x*x*z*z*z+y*y*z*r2-7*y*y*z*z*z) + f[q] * 5.6770481745453605 * (x*x*z-y*y*z);
                  else if (m == 3)
                    a[q] = g[q] * 1.7701307697799307 * (x*x*x*z*z-3*x*y*y*z*z) + f[q] * 1.7701307697799307 * (x*x*x-3*x*y*y);
                  else
                    a[q] = g[q] * 0.62583573544917614 * (x*x*x*x*z-6*x*x*y*y*z+y*y*y*y*z);
                  z += h[2];
		}
	      y += h[1];
	    }
	  x += h[0];
	}
    }
  else
      assert(0 == 1);
}

