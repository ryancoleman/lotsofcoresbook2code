/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <string.h>
#include "bmgs.h"

void Z(bmgs_cut)(const T* a, const int n[3], const int c[3],
                 T* b, const int m[3])
{
  a += c[2] + (c[1] + c[0] * n[1]) * n[2];
  for (int i0 = 0; i0 < m[0]; i0++)
    {
      for (int i1 = 0; i1 < m[1]; i1++)
        {
          memcpy(b, a, m[2] * sizeof(T));
          a += n[2];
          b += m[2];
        }
      a += n[2] * (n[1] - m[1]);
    }
}

#ifdef BMGSCOMPLEX
void bmgs_cutmz(const double_complex* a, const int sizea[3],
		const int start[3],
		double_complex* b, const int sizeb[3], double_complex p)
{
  a += start[2] + (start[1] + start[0] * sizea[1]) * sizea[2];
  for (int i0 = 0; i0 < sizeb[0]; i0++)
    {
      for (int i1 = 0; i1 < sizeb[1]; i1++)
	{
	  for (int i2 = 0; i2 < sizeb[2]; i2++)
	    b[i2] = p * a[i2];
	  a += sizea[2];
	  b += sizeb[2];
	}
      a += sizea[2] * (sizea[1] - sizeb[1]);
    }
}
#endif
