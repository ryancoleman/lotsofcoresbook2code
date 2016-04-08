/* Definitions needed by code transfered with f2c */
#include <stdio.h>
#include <math.h>

typedef int integer;
typedef double doublereal;
typedef struct { doublereal r, i; } doublecomplex;

#ifndef STATIC_NUMERIC
inline double pow_dd(double *x, double *y) {
  return pow(*x,*y);
}
#endif
