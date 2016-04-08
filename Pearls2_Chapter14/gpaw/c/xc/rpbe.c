/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <math.h>
#include "xc_gpaw.h"

double rpbe_exchange(const xc_parameters* par,
		     double n, double rs, double a2,
                     double* dedrs, double* deda2)
{
  double e = C1 / rs;
  *dedrs = -e / rs;
  if (par->gga)  // not really needed? XXX
    {
      double c = C2 * rs / n;
      c *= c;
      double s2 = a2 * c;
      double x = exp(-MU * s2 / 0.804);
      double Fx = 1.0 + 0.804 * (1 - x);
      double dFxds2 = MU * x;
      double ds2drs = 8.0 * c * a2 / rs;
      *dedrs = *dedrs * Fx + e * dFxds2 * ds2drs;
      *deda2 = e * dFxds2 * c;
      e *= Fx;
    }
  return e;
}
