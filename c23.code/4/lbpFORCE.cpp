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


int fInteractionForceZero()
{
  // reset all interaction forces to zero

  for(long k=0; k<3*lbsy.nf*lbdm.touter; k++)
    lbinterforce[k]=0;

  return 0;
}

int fCalcPotential_ShanChen()
{
  // calculate interaction potentials
  double *pt2 = &lbf[0];
  double *pt3 = &lbft[0];

  for(int j=0; j<lbsy.nf; j++) {
    
    int pottyp = lbscpot[j];
    switch (pottyp) {
      case 0:
      // Shan-Chen model (1993): psi = rho0 * (1 - exp(-rho/rho0))
        #pragma omp parallel
        {
          double rho0 = lbincp[j];
          long counter_pt2 = 0;
          long counter_pt3 = 0;
      
          #pragma omp for
          for(long i=0; i<lbdm.touter; i++) {
            counter_pt3 = i*lbsy.nf;
            counter_pt2 = i*lbsitelength;
            pt3[counter_pt3+j] = rho0 * (1.0 - exp(- fGetOneMassSite(&pt2[counter_pt2+j*lbsy.nq])/rho0));
          }
        }
        break;
      case 1:
      // Shan-Chen model (1994): psi = psi0 * exp (-rho0/rho)
        #pragma omp parallel
        {
          double rho0 = lbincp[j];
          double psi0 = lbpsi0[j];
          long counter_pt2 = 0;
          long counter_pt3 = 0;
      
          #pragma omp for
          for(long i=0; i<lbdm.touter; i++) {
            counter_pt3 = i*lbsy.nf;
            counter_pt2 = i*lbsitelength;
            pt3[counter_pt3+j] = psi0 * exp(- rho0/fGetOneMassSite(&pt2[counter_pt2+j*lbsy.nq]));
          }
        }
        break;
      case 2:
      // psi = rho
        #pragma omp parallel
        {
          long counter_pt2 = 0;
          long counter_pt3 = 0;
      
          #pragma omp for
          for(long i=0; i<lbdm.touter; i++) {
            counter_pt3 = i*lbsy.nf;
            counter_pt2 = i*lbsitelength;
            pt3[counter_pt3+j] = fGetOneMassSite(&pt2[counter_pt2+j*lbsy.nq]);
          }
        }
        break;
      // cases 3+: various Equations of State - TO DO
    }
  }

  pt2 = NULL;
  pt3 = NULL;
  return 0;
}

int fCalcInteraction_ShanChen(int xpos, int ypos, int zpos)
{
  long spos, tpos;
  int q;
  double factor0[lbsy.nf], factor1[lbsy.nf], factor2[lbsy.nf];
  double wfactor0 = 0.0;
  double wfactor1 = 0.0;
  double wfactor2 = 0.0;

  tpos = (xpos * lbdm.youter + ypos) * lbdm.zouter + zpos;
  for(int k=0; k<lbsy.nf; k++){
    factor0[k] =0.0; factor1[k] =0.0; factor2[k] =0.0;
  }
  for(int m=1; m<lbsy.nq; m++) {
    q = 3 * m;
    spos = ((xpos+lbv[q])*lbdm.youter + ypos+lbv[q+1])*lbdm.zouter + zpos+lbv[q+2];

    for(int k=0; k<lbsy.nf; k++) {
      double psi = lbft[spos * lbsy.nf + k];
      factor0[k] += lbvw[q]   * psi;
      factor1[k] += lbvw[q+1] * psi;
      factor2[k] += lbvw[q+2] * psi;
    }
    double spsi = (lbneigh[tpos]>0 && lbphi[spos]>10);
    wfactor0 += lbvw[q]   * spsi;
    wfactor1 += lbvw[q+1] * spsi;
    wfactor2 += lbvw[q+2] * spsi;
  }

  for(int k=0; k<lbsy.nf; k++) {
    double psik = fCppAbs(lbft[tpos*lbsy.nf+k]);
    for(int j=0; j<lbsy.nf; j++) {
      double gfluid = lbg[k*lbsy.nf+j]*psik;
      lbinterforce[tpos*3*lbsy.nf+3*k]  -=gfluid*factor0[j];
      lbinterforce[tpos*3*lbsy.nf+3*k+1]-=gfluid*factor1[j];
      lbinterforce[tpos*3*lbsy.nf+3*k+2]-=gfluid*factor2[j];
    }
    double gwall = lbgwall[k]*psik;
    lbinterforce[tpos*3*lbsy.nf+3*k]  -=gwall*wfactor0;
    lbinterforce[tpos*3*lbsy.nf+3*k+1]-=gwall*wfactor1;
    lbinterforce[tpos*3*lbsy.nf+3*k+2]-=gwall*wfactor2;
  }
  return 0;
}

int fCalcInteraction_ShanChen_Boundary(int xpos, int ypos, int zpos)
{
  long spos, tpos;
  int q;
  double factor0[lbsy.nf], factor1[lbsy.nf], factor2[lbsy.nf];
  double wfactor0 = 0.0;
  double wfactor1 = 0.0;
  double wfactor2 = 0.0;

  tpos = (xpos * lbdm.youter + ypos) * lbdm.zouter + zpos;
  for(int k=0; k<lbsy.nf; k++){
    factor0[k] =0.0; factor1[k] =0.0; factor2[k] =0.0;
  }
  for(int m=1; m<lbsy.nq; m++) {
    q = 3 * m;
    spos = (fCppMod(xpos+lbv[q]  , lbdm.xouter)  * lbdm.youter
          + fCppMod(ypos+lbv[q+1], lbdm.youter)) * lbdm.zouter
          + fCppMod(zpos+lbv[q+2], lbdm.zouter);

    for(int k=0; k<lbsy.nf; k++) {
      double psi = lbft[spos * lbsy.nf + k];
      factor0[k] += lbvw[q]   * psi;
      factor1[k] += lbvw[q+1] * psi;
      factor2[k] += lbvw[q+2] * psi;
    }    
    double spsi = (lbneigh[tpos]>0 && lbphi[spos]>10);
    wfactor0 += lbvw[q]   * spsi;
    wfactor1 += lbvw[q+1] * spsi;
    wfactor2 += lbvw[q+2] * spsi;
  }

  for(int k=0; k<lbsy.nf; k++) {
    double psik = fCppAbs(lbft[tpos*lbsy.nf+k]);
    for(int j=0; j<lbsy.nf; j++) {
      double gfluid = lbg[k*lbsy.nf+j]*psik;
      lbinterforce[tpos*3*lbsy.nf+3*k]  -=gfluid*factor0[j];
      lbinterforce[tpos*3*lbsy.nf+3*k+1]-=gfluid*factor1[j];
      lbinterforce[tpos*3*lbsy.nf+3*k+2]-=gfluid*factor2[j];
    }
    double gwall = lbgwall[k]*psik;
    lbinterforce[tpos*3*lbsy.nf+3*k]  -=gwall*wfactor0;
    lbinterforce[tpos*3*lbsy.nf+3*k+1]-=gwall*wfactor1;
    lbinterforce[tpos*3*lbsy.nf+3*k+2]-=gwall*wfactor2;
  }
  return 0;
}




int fsInteractionForceShanChen()
{
  unsigned long il, i, j, k;
  unsigned long Xmax = (unsigned long) lbdm.xouter,
                Ymax = (unsigned long) lbdm.youter,
                Zmax = (unsigned long) lbdm.zouter;

#pragma omp parallel
  {
    #pragma omp for
    for(i=0; i<lboutersize; i++) {
      if (lbphi[lbouter[4*i]] == 0)
        fCalcInteraction_ShanChen_Boundary((int) lbouter[4*i+1], (int) lbouter[4*i+2], (int) lbouter[4*i+3]);
    }
    #pragma omp for private(i,j,k,il) collapse(2)
    for(i=lbdm.owidx; i<(Xmax-lbdm.owidx); i++)
      for(j=lbdm.owidy; j<(Ymax-lbdm.owidy); j++)
        for(k=lbdm.owidz; k<(Zmax-lbdm.owidz); k++) {
          il = (i * Ymax + j) * Zmax + k;
          if(lbphi[il] == 0)
            fCalcInteraction_ShanChen((int) i, (int) j, (int) k);
      }

  }
  return 0;
}











int fCalcForce_Boussinesq(long tpos, double temph, double templ)
{
  double temp, rhotemp;
  double temp0 = 0.5 * (temph + templ);
  temp = fGetTemperatureSite(tpos);
  double *pt2 = &lbf[tpos*lbsitelength];
  if(!incompress) {
    for(int j=0; j<lbsy.nf; j++) {
//      rho = fGetOneMassSite(pt2);
      rhotemp = fGetOneMassSite(pt2) * (temp-temp0) * fReciprocal(temph-templ);
      lbinterforce[tpos*3*lbsy.nf+3*j]   -= rhotemp * lbbousforce[j*3];
      lbinterforce[tpos*3*lbsy.nf+3*j+1] -= rhotemp * lbbousforce[j*3+1];
      lbinterforce[tpos*3*lbsy.nf+3*j+2] -= rhotemp * lbbousforce[j*3+2];
      pt2 += lbsy.nq;
    }
  }
  else {
    for(int j=0; j<lbsy.nf; j++) {
//      rho = lbincp[j];
      rhotemp = lbincp[j] * (temp-temp0) * fReciprocal(temph-templ);
      lbinterforce[tpos*3*lbsy.nf+3*j]   -= rhotemp * lbbousforce[j*3];
      lbinterforce[tpos*3*lbsy.nf+3*j+1] -= rhotemp * lbbousforce[j*3+1];
      lbinterforce[tpos*3*lbsy.nf+3*j+2] -= rhotemp * lbbousforce[j*3+2];
      pt2 += lbsy.nq;
    }
  }
  pt2 = NULL;
  return 0;
}


int fConvectionForceBoussinesq(double temph, double templ)
{
  long il;
  int Xouter = lbdm.xouter;
  int Youter = lbdm.youter;
  int Zouter = lbdm.zouter;

  if(lbsy.nt==1) {
#pragma omp parallel private(il)
    {
      #pragma omp for collapse(3)
      for(int i=0; i<Xouter; i++) {
        for(int j=0; j<Youter; j++) {
          for(int k=0; k<Zouter; k++) {
            il = (i * Youter + j) * Zouter + k;
            if(lbphi[il]<11)
              fCalcForce_Boussinesq(il, temph, templ);
          }
        }
      }
    }
  }
  return 0;
}



