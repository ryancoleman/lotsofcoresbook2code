/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <math.h>
#include "xc_gpaw.h"

double pbe_exchange(const xc_parameters* par,
		    double n, double rs, double a2,
		    double* dedrs, double* deda2)
{
  double e = C1 / rs;
  *dedrs = -e / rs;
  if (par->gga)
    {
      double c = C2 * rs / n;
      c *= c;
      double s2 = a2 * c;
      double x = 1.0 + MU * s2 / par->kappa;
      double Fx = 1.0 + par->kappa - par->kappa / x;
      double dFxds2 = MU / (x * x);
      double ds2drs = 8.0 * c * a2 / rs;
      //double ds2drs = 8.0 * s2 / rs;
      *dedrs = *dedrs * Fx + e * dFxds2 * ds2drs;
      *deda2 = e * dFxds2 * c;
      e *= Fx;
    }
  return e;
}

/* inline */ double G(double rtrs, double A, double alpha1,
                double beta1, double beta2, double beta3, double beta4,
                double* dGdrs)
{
  double Q0 = -2.0 * A * (1.0 + alpha1 * rtrs * rtrs);
  double Q1 = 2.0 * A * rtrs * (beta1 + 
                                rtrs * (beta2 + 
                                        rtrs * (beta3 + 
                                                rtrs * beta4)));
  double G1 = Q0 * log(1.0 + 1.0 / Q1);
  double dQ1drs = A * (beta1 / rtrs + 2.0 * beta2 +
                       rtrs * (3.0 * beta3 + 4.0 * beta4 * rtrs));
  *dGdrs = -2.0 * A * alpha1 * G1 / Q0 - Q0 * dQ1drs / (Q1 * (Q1 + 1.0));
  return G1;
}

/*
 * In[1]:= H=g Log[1+b/g t^2(1+a t^2)/(1+a t^2 + a^2 t^4)]    
 * 
 *                        2         2
 *                     b t  (1 + a t )
 * Out[1]= g Log[1 + --------------------]
 *                             2    2  4
 *                   g (1 + a t  + a  t )  
 *
 * In[4]:= Simplify[D[H,t]] 
 *
 *                                             2
 *                           2 b g t (1 + 2 a t )
 * Out[4]= ---------------------------------------------------------
 *                 2    2  4          2        2        4    2    4
 *         (1 + a t  + a  t ) (g + b t  + a g t  + a b t  + a  g t )
 *
 */

double pbe_correlation(double n, double rs, double zeta, double a2, 
		       bool gga, bool spinpol,
		       double* dedrs, double* dedzeta, double* deda2)
{
  double rtrs = sqrt(rs);
  double de0drs;
  double e0 = G(rtrs, GAMMA, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294,
		&de0drs);
  double e;
  double xp = 117.0;
  double xm = 117.0;
  if (spinpol)
    {
      double de1drs;
      double e1 = G(rtrs, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662,
                    0.62517, &de1drs);
      double dalphadrs;
      double alpha = -G(rtrs, 0.016887, 0.11125, 10.357, 3.6231, 0.88026,
                        0.49671, &dalphadrs);
      dalphadrs = -dalphadrs;
      double zp = 1.0 + zeta;
      double zm = 1.0 - zeta;
      xp = pow(zp, THIRD);
      xm = pow(zm, THIRD);
      double f = CC1 * (zp * xp + zm * xm - 2.0);
      double f1 = CC2 * (xp - xm);
      double zeta2 = zeta * zeta;
      double zeta3 = zeta2 * zeta;
      double zeta4 = zeta2 * zeta2;
      double x = 1.0 - zeta4;
      *dedrs = (de0drs * (1.0 - f * zeta4) + 
               de1drs * f * zeta4 +
               dalphadrs * f * x * IF2);
      *dedzeta = (4.0 * zeta3 * f * (e1 - e0 - alpha * IF2) +
                 f1 * (zeta4 * e1 - zeta4 * e0 + x * alpha * IF2));
      e = e0 + alpha * IF2 * f * x + (e1 - e0) * f * zeta4;
    }
  else
    {
      *dedrs = de0drs;
      e = e0;
    }
  if (gga)
    {
      double n2 = n * n;
      double t2;
      double y;
      double phi = 117.0;
      double phi2 = 117.0;
      double phi3 = 117.0;
      if (spinpol)
        {
          phi = 0.5 * (xp * xp + xm * xm);
          phi2 = phi * phi;
          phi3 = phi * phi2;
          t2 = C3 * a2 * rs / (n2 * phi2);
          y = -e / (GAMMA * phi3);
        }
      else
        {
          t2 = C3 * a2 * rs / n2;
          y = -e / GAMMA;
        }
      double x = exp(y);
      double A;
      if (x != 1.0)
        A = BETA / (GAMMA * (x - 1.0)); 
      else
        A = BETA / (GAMMA * y);
      double At2 = A * t2;
      double nom = 1.0 + At2;
      double denom = nom + At2 * At2;
      double H = GAMMA * log( 1.0 + BETA * t2 * nom / (denom * GAMMA));
      double tmp = (GAMMA * BETA /
                    (denom * (BETA * t2 * nom + GAMMA * denom)));
      double tmp2 = A * A * x / BETA;
      double dAdrs = tmp2 * *dedrs;
      if (spinpol)
        {
          H *= phi3;
          tmp *= phi3;
          dAdrs /= phi3;
        }
      double dHdt2 = (1.0 + 2.0 * At2) * tmp;
      double dHdA = -At2 * t2 * t2 * (2.0 + At2) * tmp;
      *dedrs += dHdt2 * 7 * t2 / rs + dHdA * dAdrs;
      *deda2 = dHdt2 * C3 * rs / n2;
      if (spinpol)
        {
          double dphidzeta = (1.0 / xp - 1.0 / xm) / 3.0;
          double dAdzeta = tmp2 * (*dedzeta - 
				   3.0 * e * dphidzeta / phi) / phi3;
          *dedzeta += ((3.0 * H / phi - dHdt2 * 2.0 * t2 / phi ) * dphidzeta +
                      dHdA * dAdzeta);
          *deda2 /= phi2;
        }
      e += H;
    }
  return e;
}
