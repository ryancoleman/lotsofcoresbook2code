/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include "bmgs.h"

void Z(bmgs_paste)(const T* a, const int sizea[3],
		   T* b, const int sizeb[3], const int startb[3])
{
  b += startb[2] + (startb[1] + startb[0] * sizeb[1]) * sizeb[2];
  for (int i0 = 0; i0 < sizea[0]; i0++)
    {
      for (int i1 = 0; i1 < sizea[1]; i1++)
	{
	  memcpy(b, a, sizea[2] * sizeof(T));
	  a += sizea[2];
	  b += sizeb[2];
	}
      b += sizeb[2] * (sizeb[1] - sizea[1]);
    }
}

void Z(bmgs_pastep)(const T* a, const int sizea[3],
		    T* b, const int sizeb[3], const int startb[3])
{
  b += startb[2] + (startb[1] + startb[0] * sizeb[1]) * sizeb[2];
  for (int i0 = 0; i0 < sizea[0]; i0++)
    {
      for (int i1 = 0; i1 < sizea[1]; i1++)
	{
	  for (int i2 = 0; i2 < sizea[2]; i2++)
	    b[i2] += *a++;
	  b += sizeb[2];
	}
      b += sizeb[2] * (sizeb[1] - sizea[1]);
    }
}
