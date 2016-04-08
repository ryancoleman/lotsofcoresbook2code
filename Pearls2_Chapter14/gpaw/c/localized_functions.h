/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2008  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <Python.h>

typedef struct 
{
  PyObject_HEAD
  double dv;     // volume per grid point
  int size[3];   // dimensions of big box
  int start[3];  // corner of small box
  int size0[3];  // dimensions of small box
  int ng;        // number of grid points in big box
  int ng0;       // number of grid points in small box
  int nf;        // number of localized functions
  int nfd;       // number of derivatives: zero or 3*nf 
                 // pointers to size0 arrays:
  double* f;     // localized functions
  double* fd;    // xyz-derivatives of localized functions
  double* w;     // work array for one double or double complex array
} LocalizedFunctionsObject;

