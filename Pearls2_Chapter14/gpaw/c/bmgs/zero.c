/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <string.h>
#include "bmgs.h"

void Z(bmgs_zero)(T* a, const int n[3], const int c[3],
		  const int s[3])
{
  a += c[2] + (c[1] + c[0] * n[1]) * n[2];
  for (int i0 = 0; i0 < s[0]; i0++)
    {
      for (int i1 = 0; i1 < s[1]; i1++)
	{
	  memset(a, 0, s[2] * sizeof(T));
	  a += n[2];
	}
      a += n[2] * (n[1] - s[1]);
    }
}
