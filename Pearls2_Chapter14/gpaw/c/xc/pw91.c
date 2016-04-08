/*  Copyright (C) 2003-2007  CAMP
 *  Please see the accompanying LICENSE file for further information. */

#include <math.h>
#include "xc_gpaw.h"

double G(double rtrs, double A, double alpha1,
	 double beta1, double beta2, double beta3, double beta4,
	 double* dGdrs);

double pw91_exchange(const xc_parameters* par,
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
      double s = sqrt(s2);
      double f1 = 7.7956 * s;
      double f2 = 0.19645 * asinh(f1);
      double f3 = 0.1508 * exp(-100.0 * s2);
      double f4 = 0.004 * s2 * s2;
      double f5 = 1.0 + s * f2;
      double f6 = f5 + f4;
      double f7 = 0.2743 - f3;
      double f8 = f5 + f7 * s2;
      double Fx = f8 / f6;
      double f9 = 0.5 * 7.7956 * 0.19645 / sqrt(1.0 + f1 * f1);
      if (s < 0.00001)
	f9 += 0.5 * 7.7956 * 0.19645;
      else
	f9 += 0.5 * f2 / s;
      double dFxds2 = ((f9 + f7 + 100.0 * f3 * s2) * f6 -
		       f8 * (f9 + 0.008 * s2)) / (f6 * f6);
      double ds2drs = 8.0 * s2 / rs;
      *dedrs = *dedrs * Fx + e * dFxds2 * ds2drs;
      *deda2 = e * dFxds2 * c;
      e *= Fx;
    }
  return e;
}

double pw91_correlation(double n, double rs, double zeta, double a2, 
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
      double phi;
      double phi2;
      double phi3;
      double phi4;
      double GAMMAPW91 = BETA * BETA / 0.18;
      if (spinpol)
        {
          phi = 0.5 * (xp * xp + xm * xm);
          phi2 = phi * phi;
          phi3 = phi * phi2;
          phi4 = phi * phi3;

        }
      else
        {
          phi = 1.0;
          phi2 = 1.0;
          phi3 = 1.0;
          phi4 = 1.0;
        }
      t2 = C3 * a2 * rs / (n2 * phi2);
      y = -e / (GAMMAPW91 * phi3);
      double x = exp(y);
      double A = BETA / (GAMMAPW91 * (x - 1.0)); 
      double At2 = A * t2;
      double nom = 1.0 + At2;
      double denom = nom + At2 * At2;
      double H0 = (phi3 * GAMMAPW91 * 
		   log(1.0 + BETA * t2 * nom / (denom * GAMMAPW91)));
      double tmp = (phi3 * GAMMAPW91 * BETA /
                    (denom * (BETA * t2 * nom + GAMMAPW91 * denom)));
      double tmp2 = A * A * x / BETA;
      double dAdrs = tmp2 * *dedrs / phi3;

      const double KK = 66.343643960645011; // 100*4/pi*(4/pi/9)**(1/3.)
      const double XNU = 15.75592;
      const double Cc0 = 0.004235;
      const double Cx = -0.001667212;
      const double K1 = 0.002568;
      const double K2 = 0.023266;
      const double K3 = 7.389e-6;
      const double K4 = 8.723;
      const double K5 = 0.472;
      const double K6 = 7.389e-2;

      double f0 = XNU * exp(-KK * rs * phi4 * t2);
      double rs2 = rs * rs;
      double f1 = K1 + K2 * rs + K3 * rs2;
      double f2 = 1.0 + K4 * rs + K5 * rs2 + K6 * rs2 * rs;
      double f3 = -10.0 * Cx / 7.0 - Cc0 + f1 / f2;
      double H1 = f0 * phi3 * f3 * t2;
      double dH1drs = (-KK * phi4 * t2 * H1 + f0 * phi3 * t2 *
		       ((K2 + 2.0 * K3 * rs) * f2 -
			(K4 + 2.0 * K5 * rs + 3.0 * K6 * rs2) * f1) / (f2 * f2));
      double dH1dt2 = -KK * rs * phi4 * H1 + f0 * phi3 * f3;
      double dH1dphi = (-4.0 * KK * rs * phi3 * H1 + 3.0 * f0 * phi2 * f3) * t2;
      double dH0dt2 = (1.0 + 2.0 * At2) * tmp;
      double dH0dA = -At2 * t2 * t2 * (2.0 + At2) * tmp;
      *dedrs += (dH0dt2 + dH1dt2) * 7 * t2 / rs + dH0dA * dAdrs + dH1drs;
      *deda2 = (dH0dt2 + dH1dt2) * C3 * rs / n2;
      if (spinpol)
        {
          double dphidzeta = (1.0 / xp - 1.0 / xm) / 3.0;
          double dAdzeta = tmp2 * (*dedzeta - 
				   3.0 * e * dphidzeta / phi) / phi3;
          *dedzeta += ((3.0 * H0 / phi - dH0dt2 * 2.0 * t2 / phi ) * dphidzeta +
                      dH0dA * dAdzeta);
          *dedzeta += (dH1dphi - dH1dt2 * 2.0 * t2 / phi ) * dphidzeta;
          *deda2 /= phi2;
        }
      e += H0 + H1;
    }
  return e;
}
