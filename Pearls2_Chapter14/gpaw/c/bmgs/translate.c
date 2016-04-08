/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <string.h>
#include "bmgs.h"

void bmgs_translate(double* a, const int sizea[3], const int size[3],
		    const int start1[3], const int start2[3])
{
  const double* restrict s = 
    a + start1[2] + (start1[1] + start1[0] * sizea[1]) * sizea[2];
  double* restrict d = 
    a + start2[2] + (start2[1] + start2[0] * sizea[1]) * sizea[2];
  for (int i0 = 0; i0 < size[0]; i0++)
    {
      for (int i1 = 0; i1 < size[1]; i1++)
        {
          memcpy(d, s, size[2] * sizeof(double));
          s += sizea[2];
          d += sizea[2];
        }
      s += sizea[2] * (sizea[1] - size[1]);
      d += sizea[2] * (sizea[1] - size[1]);
    }
}

void bmgs_translatemz(double_complex* a, const int sizea[3], const int size[3],
		      const int start1[3], const int start2[3],
		      double_complex phase)
{
  const double_complex* restrict s = 
    a + start1[2] + (start1[1] + start1[0] * sizea[1]) * sizea[2];
  double_complex* restrict d = 
    a + start2[2] + (start2[1] + start2[0] * sizea[1]) * sizea[2];
  for (int i0 = 0; i0 < size[0]; i0++)
    {
      for (int i1 = 0; i1 < size[1]; i1++)
        {
	  for (int i2 = 0; i2 < size[2]; i2++)
	    d[i2] = phase * s[i2];
          s += sizea[2];
          d += sizea[2];
        }
      s += sizea[2] * (sizea[1] - size[1]);
      d += sizea[2] * (sizea[1] - size[1]);
    }
}
