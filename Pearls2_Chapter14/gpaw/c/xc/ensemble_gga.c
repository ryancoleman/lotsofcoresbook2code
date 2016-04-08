/*  Copyright (C) 2003-2007  CAMP
 *  Copyright (C) 2007-2009  CAMd
 *  Please see the accompanying LICENSE file for further information. */

#include <math.h>
#include "xc_gpaw.h"

double beefvdw_exchange(const xc_parameters* par,
                       double n, double rs, double a2,
                       double* dedrs, double* deda2)
{
  double e = C1 / rs;
  *dedrs = -e / rs;
  double c = C2 * rs / n;
  c *= c;
  double s2 = a2 * c;
  
  /* Legendre polynomial basis expansion */
  int parlen = par->nparameters-1;
  double p = par->parameters[0];
  double tmp = p + s2;
  double x = 2.0 * s2 / tmp - 1.0;
  double dxds2 = 2.0 * p / pow(tmp,2);
  double Fx = 0.0;
  double dFxds2 = 0.0;

  int max_order = par->parameters[parlen+1];
  double L[max_order+1];
  double dL[max_order+1];
  double coef;
  int m;
  int order;

  /* initializing */
  L[0] = 1.0;
  L[1] = x;
  dL[0] = 0.0;
  dL[1] = 1.0;

  /* recursively building polynomia and their derivatives */
  for(int i = 2; i < max_order+1; i++)
    {
    L[i] = 2.0 * x * L[i-1] - L[i-2] - (x * L[i-1] - L[i-2])/i;
    dL[i] = i * L[i-1] + x * dL[i-1];
    }

  /* building enhancement factor Fx and derivative dFxds2 */
  m = 0;
  for(int i = 0; i < max_order+1; i++)
    {
    order = par->parameters[2+m];
    if(order == i)
      {
      coef = par->parameters[2+parlen+m];
      Fx += coef * L[i];
      dFxds2 += coef * dL[i] * dxds2;
      m += 1;
      }
    }

  double ds2drs = 8.0 * c * a2 / rs;
  *dedrs = *dedrs * Fx + e * dFxds2 * ds2drs;
  *deda2 = e * dFxds2 * c;
  e *= Fx;
  return e;
}
