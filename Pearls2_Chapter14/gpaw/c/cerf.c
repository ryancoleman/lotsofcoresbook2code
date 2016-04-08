#include <Python.h>
#include "extensions.h"
#include <float.h>
#include <math.h>

#define eps 1.e-15

double_complex itpp_erf(double_complex z);

PyObject* cerf(PyObject *self, PyObject *args)
{
  double complex z, res;
  if (!PyArg_ParseTuple(args, "D", &z)) 
    return NULL;

  res = itpp_erf(z);
  return Py_BuildValue("D", &res);
}

/* taken from 
   http://prdownloads.sourceforge.net/itpp/itpp-3.10.7.tar.bz2
   and transformed to C */ 

 /*!
 * \file
 * \brief Implementation of scalar functions
 * \author Tony Ottosson, Pal Frenger and Adam Piatyszek
 *
 * $Date: 2006-08-19 10:53:33 +0200 (sob, 19 sie 2006) $
 * $Revision: 643 $
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

double_complex itpp_cerf_series(double_complex z);
double_complex itpp_cerfc_continued_fraction(double_complex z);
double_complex itpp_cerf_continued_fraction(double_complex z);
double_complex itpp_cerf_rybicki(double_complex z);

double cabs(double_complex z) {
  return sqrt(creal(z) * creal(z) + cimag(z) * cimag(z));
}

/*
 * This function calculates a well known error function erf(z) for
 * complex z. Three methods are implemented. Which one is used
 * depends on z.
 */
double_complex itpp_erf(double_complex z)
{
    // Use the method appropriate to size of z -
    // there probably ought to be an extra option for NaN z, or infinite z
  if (cabs(z) < 2.0)
    return itpp_cerf_series(z);
  else {
    if (fabs(creal(z)) < 0.5)
      // XXX neither rybicki nor continued_fraction seem to work here
      //     return itpp_cerf_rybicki(z);
      //return itpp_cerf_continued_fraction(z);
      return itpp_cerf_series(z);
    else
      return itpp_cerf_continued_fraction(z);
  }
}

/*
 * Abramawitz and Stegun: Eq. (7.1.5) gives a series for erf(z) good
 * for all z, but converges faster for smallish abs(z), say abs(z) < 2.
 */
double_complex itpp_cerf_series(double_complex z)
{
  double_complex sum, term, z2, oldsum;
  double error;

  sum = 0.0;
  term = z;
  z2 = z * z;

  oldsum = 1.e32;
  for (int n = 0; 1; n++) {
    sum += term / (2. * n + 1);
    term *= -z2 / (1. * n + 1);
    error = cabs(sum / oldsum - 1.);
    if (error < eps) {
      return sum * (2.0 / sqrt(M_PI)); }
    oldsum = sum;
  }
}

/*
 * Abramowitz and Stegun: Eq. (7.1.14) gives this continued fraction
 * for erfc(z)
 *
 * erfc(z) = sqrt(pi).exp(-z^2).  1   1/2   1   3/2   2   5/2
 *                               ---  ---  ---  ---  ---  --- ...
 *                               z +  z +  z +  z +  z +  z +
 *
 * This is evaluated using Lentz's method, as described in the
 * narative of Numerical Recipes in C.
 *
 * The continued fraction is true providing real(z) > 0. In practice
 * we like real(z) to be significantly greater than 0, say greater
 * than 0.5.
 */
double_complex itpp_cerfc_continued_fraction(double_complex z)
{
  // first calculate z+ 1/2   1
  //                    ---  --- ...
  //                    z +  z +
  double_complex f, C, D, delta;
  double a;
  //  printf("itpp_cerfc_continued_fraction\n");
    
  f = z;
  C = f;
  D = 0.0;
  
  a = 0.0;
  do {
    a += 0.5;
    D = z + a * D;
    C = z + a / C;
    if ((creal(D) == 0.0) && (cimag(D) == 0.0))
      D = DBL_MIN;
    D = 1.0 / D;
    delta = C * D;
    f = f * delta;
  } while (cabs(1.0 - delta) > eps);

  // Do the first term of the continued fraction
  f = 1.0 / f;

  // and do the final scaling
  f = f * exp(-z * z) / sqrt(M_PI);
  
  return f;
}

double_complex itpp_cerf_continued_fraction(double_complex z)
{
  if (creal(z) > 0)
    return 1.0 - itpp_cerfc_continued_fraction(z);
  else
    return -1.0 + itpp_cerfc_continued_fraction(-z);
}

/*
 * Numerical Recipes quotes a formula due to Rybicki for evaluating
 * Dawson's Integral:
 *
 * exp(-x^2) integral exp(t^2).dt = 1/sqrt(pi) lim  sum  exp(-(z-n.h)^2) / n
 *            0 to x                           h->0 n odd
 *
 * This can be adapted to erf(z).
 */
double_complex itpp_cerf_rybicki(double_complex z)
{
  double h = 0.2; // numerical experiment suggests this is small enough
  printf("itpp_cerf_rybicki");

  // choose an even n0, and then shift z->z-n0.h and n->n-h.
  // n0 is chosen so that real((z-n0.h)^2) is as small as possible.
  int n0 = 2 * ((int)(cimag(z) / (2 * h) + 0.5));

  double_complex z0 = I * n0 * h;
  double_complex zp = z - z0;
  double_complex sum = 0.0;

  // limits of sum chosen so that the end sums of the sum are
  // fairly small. In this case exp(-(35.h)^2)=5e-22
  for (int np = -35; np <= 35; np += 2) {
    double_complex t = creal(zp) + I * (cimag(zp) - np * h);
    double_complex b = (exp(t * t) / ((double)(np + n0)));
    sum += b;
  }

  sum *= 2.0 * exp(-z * z) / M_PI;
  sum = - cimag(sum) + creal(sum) * I;
  return sum;
}

