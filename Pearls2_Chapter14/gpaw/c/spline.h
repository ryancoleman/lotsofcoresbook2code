/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include "extensions.h"
#include "bmgs/bmgs.h"

typedef struct 
{
  PyObject_HEAD
  bmgsspline spline;
} SplineObject;
