/***********************************************************************
 *   DL_MESO       Version 2.6                                         *
 *   Authors   :   R. S. Qin, M. A. Seaton                             *
 *   Copyright :   STFC Daresbury Laboratory                           *
 *             :   07/08/2014                                          *
 ***********************************************************************
 *   This extract used with permission of the copyright owners         *
 ***********************************************************************
 *   The full DL_MESO package is supplied to individuals under an      *
 *   academic licence, which is free of cost to academic scientists    *
 *   pursuing scientific research of a non-commercial nature           *
 *   To register please visit www.ccp5.ac.uk/DL_MESO                   *
 *   Commercial organisations interested in acquiring the package      *
 *   should approach Dr. M. A. Seaton at Daresbury Laboratory in the   *
 *   first instance. Daresbury Laboratory is the sole centre for       * 
 *   distribution of the package                                       *
 ***********************************************************************/

long fNextStep(int q, int xpos, int ypos, int zpos);
long fNextStep(int dx, int dy, int dz, long tpos);
int fBounceBackF(long tpos);
int fBounceBackC(long tpos);
int fBounceBackT(long tpos);
int fMidBounceBackF(long tpos);
int fMidBounceBackC(long tpos);
int fMidBounceBackT(long tpos);
int fSiteBlankF(long tpos);
int fSiteBlankC(long tpos);
int fSiteBlankT(long tpos);

int fD2Q9VCE(double v0, double v1, double &f0, double &f1, double &f2,
	     double &f3, double &f4, double &f5, double &f6,double &f7, double &f8);
int fD2Q9VCEIncom(double v0, double v1, double rho0, double &f0, double &f1, double &f2,
	          double &f3, double &f4, double &f5, double &f6,double &f7, double &f8);
int fD2Q9VF(int tpos, int prop, double *uwall);
int fD2Q9VCC(double p, double *v, double * startpos);
int fD2Q9VCCIncom(double p, double p0, double *v, double * startpos);
int fD2Q9PCE(double p, double &f0, double &f1, double &f2, double &f3, 
	     double &f4, double &f5, double &f6, double &f7, double &f8, double& vel);
int fD2Q9PCC(double p, double *v, double * startpos);
int fD2Q9PCCIncom(double p, double p0, double *v, double * startpos);
int fD2Q9PF(int tpos, int prop, double *uwall);
int fD2Q9CCE(double p, double v0, double v1, double &f0, double &f1, double &f2,
             double &f3, double &f4, double &f5, double &f6, double &f7, double &f8);
int fD2Q9CCC(double p, double *v, double * startpos);
int fD2Q9PC(int tpos, int prop, double *uwall);
int fD2Q9TCC(double p, double *v, double * startpos);
int fD2Q9PT(int tpos, int prop, double *uwall);

int fD3Q15VPS(double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14);
int fD3Q15VPSIncom(double v0, double v1, double v2, double rho0, double &f0,
                   double &f1, double &f2, double &f3, double &f4, double &f5, 
	           double &f6, double &f7, double &f8, double &f9, double &f10,
                   double &f11, double &f12, double &f13, double &f14);
int fD3Q15VCE(double p, double *v, double * startpos);
int fD3Q15VCEIncom(double p, double p0, double *v, double * startpos);
int fD3Q15VF(int tpos, int prop, double *uwall);
int fD3Q15PPS(double p, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double& vel);
int fD3Q15PF(int tpos, int prop, double *uwall);
int fD3Q15CPS(double p, double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, double &f6, 
              double &f7, double &f8, double &f9, double &f10, double &f11,
              double &f12, double &f13, double &f14);
int fD3Q15CCE(double p, double *v, double * startpos);
int fD3Q15PC(int tpos, int prop, double *uwall);
int fD3Q15TCE(double p, double *v, double * startpos);
int fD3Q15PT(int tpos, int prop, double *uwall);

int fD3Q19VPS(double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double &f15, double &f16, double &f17, 
	      double &f18);
int fD3Q19VPSIncom(double v0, double v1, double v2, double rho0, 
                   double &f0, double &f1, double &f2, double &f3, double &f4,
                   double &f5, double &f6, double &f7, double &f8, double &f9,
	           double &f10, double &f11, double &f12, double &f13, double &f14,
                   double &f15, double &f16, double &f17, double &f18);
int fD3Q19VCE(double p, double *v, double * startpos);
int fD3Q19VCEIncom(double p, double p0, double *v, double * startpos);
int fD3Q19VF(int tpos, int prop, double *uwall);
int fD3Q19PPS(double p, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double &f15, double &f16, double &f17, 
	      double &f18, double& vel);
int fD3Q19PF(int tpos, int prop, double *uwall);
int fD3Q19CPS(double p, double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, double &f6,
              double &f7, double &f8, double &f9, double &f10, double &f11,
              double &f12, double &f13, double &f14, double &f15, double &f16,
              double &f17, double &f18);
int fD3Q19CCE(double p, double *v, double * startpos);
int fD3Q19PC(int tpos, int prop, double *uwall);
int fD3Q19TCE(double p, double *v, double * startpos);
int fD3Q19PT(int tpos, int prop, double *uwall);

int fD3Q27VPS(double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double &f15, double &f16, double &f17,
	      double &f18, double &f19, double &f20, double &f21,
	      double &f22, double &f23, double &f24, double &f25,
	      double &f26);
int fD3Q27VPSIncom(double v0, double v1, double v2, double rho0, 
                   double &f0, double &f1, double &f2, double &f3, double &f4,
                   double &f5, double &f6, double &f7, double &f8, double &f9,
	           double &f10, double &f11, double &f12, double &f13, double &f14,
                   double &f15, double &f16, double &f17, double &f18, double &f19,
                   double &f20, double &f21, double &f22, double &f23, double &f24,
                   double &f25, double &f26);
int fD3Q27VCE(double p, double *v, double * startpos);
int fD3Q27VCEIncom(double p, double *v, double * startpos);
int fD3Q27VF(int tpos, int prop, double *uwall);
int fD3Q27PPS(double p, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double &f15, double &f16, double &f17,
	      double &f18, double &f19, double &f20, double &f21,
	      double &f22, double &f23, double &f24, double &f25,
	      double &f26, double& vel);
int fD3Q27PF(int tpos, int prop, double *uwall);
int fD3Q27CPS(double p, double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, double &f6,
              double &f7, double &f8, double &f9, double &f10, double &f11,
              double &f12, double &f13, double &f14, double &f15, double &f16,
              double &f17, double &f18, double &f19, double &f20, double &f21,
	      double &f22, double &f23, double &f24, double &f25, double &f26);
int fD3Q27CCE(double p, double *v, double * startpos);
int fD3Q27PC(int tpos, int prop, double *uwall);
int fD3Q27TCE(double p, double *v, double * startpos);
int fD3Q27PT(int tpos, int prop, double *uwall);

int fFixedSpeedFluid(int tpos, int prop, double *uwall);
int fFixedDensityFluid(int tpos, int prop, double *uwall);
int fFixedSoluteConcen(int tpos, int prop, double *uwall);
int fFixedTemperature(int tpos, int prop, double *uwall);
int fPostCollBoundary();
int fPostPropBoundary();
int fNeighbourBoundary();

