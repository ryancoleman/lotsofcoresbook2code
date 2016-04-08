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


int fSiteFluidCollisionBGK(double* startpos, double *sitespeed, double* bodyforce)
{
  
  // calculate fluid collisions at grid point: uses BGK single relaxation time

  double speed[3], feq[lbsy.nq];
  double onemass, invmassrelax, relax;
  double *pt2 = &startpos[0];
  double *pt3 = &bodyforce[0];
  int counter_pt2 = 0;
  int counter_pt3 = 0;

  for(int j=0; j<lbsy.nf; j++) {
    onemass = fGetOneMassSite(&pt2[counter_pt2]);
    invmassrelax = fReciprocal(onemass)/lbtf[j];
    relax = 1.0 - lbtf[j];
    speed[0] = sitespeed[0] + pt3[counter_pt3]   * invmassrelax;
    speed[1] = sitespeed[1] + pt3[counter_pt3+1] * invmassrelax;
    speed[2] = sitespeed[2] + pt3[counter_pt3+2] * invmassrelax;
    counter_pt3 += 3;
    fGetEquilibriumF(&feq[0], speed, onemass);
    for(int i=0; i<lbsy.nq; i++){
      pt2[counter_pt2] = relax * pt2[counter_pt2] + feq[i]*lbtf[j];
      counter_pt2 ++;
    }
  }

  pt2 = NULL;
  pt3 = NULL;  
  return 0;
}


int fSiteFluidIncomCollisionBGK(double* startpos, double *sitespeed, double* bodyforce)
{
  
  // calculate fluid collisions at grid point: uses BGK single relaxation time
  // with incompressible fluids

  double speed[3], feq[lbsy.nq];
  double onemass, density, invmassrelax, relax;
  double *pt2 = &startpos[0];
  double *pt3 = &bodyforce[0];
  int counter_pt2 = 0;
  int counter_pt3 = 0;

  for(int j=0; j<lbsy.nf; j++) {
    onemass = fGetOneMassSite(&pt2[counter_pt2]);
    density = lbincp[j];
    invmassrelax = fReciprocal(density * lbtf[j]);
    relax = 1.0 - lbtf[j];
    speed[0] = sitespeed[0] + pt3[counter_pt3]   * invmassrelax;
    speed[1] = sitespeed[1] + pt3[counter_pt3+1] * invmassrelax;
    speed[2] = sitespeed[2] + pt3[counter_pt3+2] * invmassrelax;
    counter_pt3 += 3;
    fGetEquilibriumFIncom(&feq[0], speed, onemass, density);
    for(int i=0; i<lbsy.nq; i++){
      pt2[counter_pt2] = relax * pt2[counter_pt2] + feq[i]*lbtf[j];
      counter_pt2 ++;
    }
  }

  pt2 = NULL;
  pt3 = NULL;  
  return 0;
}



int fSiteFluidCollisionBGKEDM(double* startpos, double *sitespeed, double* bodyforce)
{
  
  // calculate fluid collisions at grid point: uses BGK single relaxation time with
  // Exact Difference Method forcing term

  double dv[3], feq[lbsy.nq];
  double onemass, invmass, relax, source;
  double modv,eixux,eiyuy,eizuz,eixdux,eiyduy,eizduz;
  double *pt2 = &startpos[0];
  double *pt3 = &bodyforce[0];
  int counter_pt2 = 0;
  int counter_pt3 = 0;

  for(int j=0; j<lbsy.nf; j++) {
    onemass = fGetOneMassSite(&pt2[counter_pt2]);
    invmass = fReciprocal(onemass);
    dv[0] = pt3[counter_pt3]   * invmass;
    dv[1] = pt3[counter_pt3+1] * invmass;
    dv[2] = pt3[counter_pt3+2] * invmass;
    counter_pt3 += 3;
    fGetEquilibriumF(&feq[0], sitespeed, onemass);
    relax = 1.0 - lbtf[j];
    modv = dv[0]*(2.0*sitespeed[0]+dv[0]) + dv[1]*(2.0*sitespeed[1]+dv[1]) + dv[2]*(2.0*sitespeed[2]+dv[2]);
    for(int i=0; i<lbsy.nq; i++){
      eixux =  lbv[i*3  ]*sitespeed[0];
      eiyuy =  lbv[i*3+1]*sitespeed[1];
      eizuz =  lbv[i*3+2]*sitespeed[2];
      eixdux = lbv[i*3  ]*dv[0];
      eiyduy = lbv[i*3+1]*dv[1];
      eizduz = lbv[i*3+2]*dv[2];
      source = onemass * lbw[i] * (3.0 * (eixdux+eiyduy+eizduz) +
                                   4.5 * (eixdux*(eixdux+2.0*(eixux+eiyuy+eizuz+eiyduy+eizduz)) +
                                          eiyduy*(eiyduy+2.0*(eixux+eiyuy+eizuz+eizduz)) +
                                          eizduz*(eizduz+2.0*(eixux+eiyuy+eizuz))) -
                                   1.5 * modv);
      pt2[counter_pt2] = relax * pt2[counter_pt2] + source + feq[i]*lbtf[j];
      counter_pt2++;
    }
  }

  pt2 = NULL;
  pt3 = NULL;  
  return 0;
}


int fSiteFluidIncomCollisionBGKEDM(double* startpos, double *sitespeed, double* bodyforce)
{
  
  // calculate fluid collisions at grid point: uses BGK single relaxation time
  // with incompressible fluids and Exact Difference Method forcing term

  double dv[3], feq[lbsy.nq];
  double onemass, density, invmass, relax, source;
  double modv,eixux,eiyuy,eizuz,eixdux,eiyduy,eizduz;
  double *pt2 = &startpos[0];
  double *pt3 = &bodyforce[0];
  int counter_pt2 = 0;
  int counter_pt3 = 0;

  for(int j=0; j<lbsy.nf; j++) {
    onemass = fGetOneMassSite(&pt2[counter_pt2]);
    density = lbincp[j];
    invmass = fReciprocal(density);
    relax = 1.0 - lbtf[j];
    dv[0] = pt3[counter_pt3]   * invmass;
    dv[1] = pt3[counter_pt3+1] * invmass;
    dv[2] = pt3[counter_pt3+2] * invmass;
    counter_pt3 += 3;
    fGetEquilibriumFIncom(&feq[0], sitespeed, onemass, density);
    modv = dv[0]*(2.0*sitespeed[0]+dv[0]) + dv[1]*(2.0*sitespeed[1]+dv[1]) + dv[2]*(2.0*sitespeed[2]+dv[2]);
    for(int i=0; i<lbsy.nq; i++){
      eixux =  lbv[i*3  ]*sitespeed[0];
      eiyuy =  lbv[i*3+1]*sitespeed[1];
      eizuz =  lbv[i*3+2]*sitespeed[2];
      eixdux = lbv[i*3  ]*dv[0];
      eiyduy = lbv[i*3+1]*dv[1];
      eizduz = lbv[i*3+2]*dv[2];
      source = density * lbw[i] * (3.0 * (eixdux+eiyduy+eizduz) +
                                   4.5 * (eixdux*(eixdux+2.0*(eixux+eiyuy+eizuz+eiyduy+eizduz)) +
                                          eiyduy*(eiyduy+2.0*(eixux+eiyuy+eizuz+eizduz)) +
                                          eizduz*(eizduz+2.0*(eixux+eiyuy+eizuz))) -
                                   1.5 * modv);
      pt2[counter_pt2] = relax * pt2[counter_pt2] + source + feq[i]*lbtf[j];
      counter_pt2 ++;
    }
  }

  pt2 = NULL;
  pt3 = NULL;  
  return 0;
}



int fSiteFluidCollisionBGKGuo(double* startpos, double *sitespeed, double* bodyforce)
{
  
  // calculate fluid collisions at grid point: uses BGK single relaxation time with
  // Guo forcing term

  double speed[3], force[3], feq[lbsy.nq];
  double onemass, udot, source, invmass, relax, halfrelax, ex, ey, ez;
  double *pt2 = &startpos[0];
  double *pt3 = &bodyforce[0];
  int counter_pt2 = 0;
  int counter_pt3 = 0;

  for(int j=0; j<lbsy.nf; j++) {
    onemass = fGetOneMassSite(&pt2[counter_pt2]);
    invmass = fReciprocal(onemass);
    force[0] = pt3[counter_pt3];
    force[1] = pt3[counter_pt3+1];
    force[2] = pt3[counter_pt3+2];
    speed[0] = sitespeed[0] + 0.5 * force[0] * invmass;
    speed[1] = sitespeed[1] + 0.5 * force[1] * invmass;
    speed[2] = sitespeed[2] + 0.5 * force[2] * invmass;
    counter_pt3 += 3;
    fGetEquilibriumF(&feq[0], speed, onemass);
    relax = 1.0 - lbtf[j];
    halfrelax = 1.0 - 0.5*lbtf[j];
    for(int i=0; i<lbsy.nq; i++){
      ex = lbv[3*i];
      ey = lbv[3*i+1];
      ez = lbv[3*i+2];
      udot = ex * speed[0] + ey * speed[1] + ez * speed[2];
      source = (3.0*(ex-speed[0]) + 9.0*udot*ex) * force[0]
             + (3.0*(ey-speed[1]) + 9.0*udot*ey) * force[1]
             + (3.0*(ez-speed[2]) + 9.0*udot*ez) * force[2];
      pt2[counter_pt2] = relax * pt2[counter_pt2] + halfrelax * lbw[i] * source + feq[i]*lbtf[j];
      counter_pt2++;
    }
  }

  pt2 = NULL;
  pt3 = NULL;  
  return 0;
}


int fSiteFluidIncomCollisionBGKGuo(double* startpos, double *sitespeed, double* bodyforce)
{
  
  // calculate fluid collisions at grid point: uses BGK single relaxation time with
  // Guo forcing term for incompressible fluids

  double speed[3], force[3], feq[lbsy.nq];
  double onemass, density, udot, source, invmass, relax, halfrelax, ex, ey, ez;
  double *pt2 = &startpos[0];
  double *pt3 = &bodyforce[0];
  int counter_pt2 = 0;
  int counter_pt3 = 0;

  for(int j=0; j<lbsy.nf; j++) {
    onemass = fGetOneMassSite(pt2);    
    density = lbincp[j];
    invmass = fReciprocal(density);
    force[0] = pt3[counter_pt3];
    force[1] = pt3[counter_pt3+1];
    force[2] = pt3[counter_pt3+2];
    speed[0] = sitespeed[0] + 0.5 * force[0] * invmass;
    speed[1] = sitespeed[1] + 0.5 * force[1] * invmass;
    speed[2] = sitespeed[2] + 0.5 * force[2] * invmass;
    counter_pt3 += 3;
    fGetEquilibriumFIncom(&feq[0], speed, onemass, density);
    relax = 1.0 - lbtf[j];
    halfrelax = 1.0 - 0.5*lbtf[j];
    for(int i=0; i<lbsy.nq; i++){
      ex = lbv[3*i];
      ey = lbv[3*i+1];
      ez = lbv[3*i+2];
      udot = ex * speed[0] + ey * speed[1] + ez * speed[2];
      source = (3.0*(ex-speed[0]) + 9.0*udot*ex) * force[0]
             + (3.0*(ey-speed[1]) + 9.0*udot*ey) * force[1]
             + (3.0*(ez-speed[2]) + 9.0*udot*ez) * force[2];
      pt2[counter_pt2] = relax * pt2[counter_pt2] + halfrelax * lbw[i] * source + feq[i]*lbtf[j];
      counter_pt2++;
    }
  }

  pt2 = NULL;
  pt3 = NULL;  
  return 0;
}



int fSiteSoluteCollisionBGK(double* startpos, double *sitespeed)
{
  
  // calculate solute collisions at grid point: uses BGK single relaxation time

  double onemass, relax;
  double feq[lbsy.nq];
  double *pt2 = &startpos[0];
  int counter_pt2 = 0;

  for(int j=0; j<lbsy.nc; j++) {
    onemass = fGetOneMassSite(&pt2[counter_pt2]);
    fGetEquilibriumC(&feq[0], sitespeed, onemass);
    relax = 1.0 - lbtc[j];
    for(int i=0; i<lbsy.nq; i++) {
      pt2[counter_pt2] = relax * pt2[counter_pt2] + feq[i]*lbtc[j];
      counter_pt2++;
    }
  }

  pt2 = NULL;
  return 0;
}


int fSiteThermalCollisionBGK(double* startpos, double *sitespeed)
{
  
  // calculate thermal collisions at grid point: uses BGK single relaxation time

  double onemass, relax;
  double feq[lbsy.nq];
  double *pt2 = &startpos[0];
  int counter_pt2 = 0;

  if(lbsy.nt == 1) {
    onemass = fGetOneMassSite(pt2);    
    fGetEquilibriumT(&feq[0], sitespeed, onemass);
    relax = 1.0 - lbtt[0];
    for(int i=0; i<lbsy.nq; i++) {
      pt2[counter_pt2] = relax * pt2[counter_pt2] + feq[i]*lbtt[0];
      counter_pt2++;
    }
  }

  pt2 = NULL;
  return 0;
}


int fCollisionBGK()
{

  long il;
  int Xmax = (long) lbdm.xouter,
      Ymax = (long) lbdm.youter,
      Zmax = (long) lbdm.zouter; 
  int i,j,k;
  int ll;
  if(!incompress) {

   #pragma omp parallel private(i,j,k,il,ll)
   {
      double interforce_local_t[3*lbsy.nf], sitespeed_local_t[3];
      double *interforce = &interforce_local_t[0];
      double *sitespeed = &sitespeed_local_t[0];

      #pragma omp for collapse(2)
      for(i=0; i<Xmax; i++) {
        for(j=0; j<Ymax; j++) {
          for(k=0; k<Zmax; k++) {
            il = (i * Ymax + j) * Zmax + k;
            fGetSpeedSite(sitespeed, &lbf[il*lbsitelength]);
            if(lbphi[il] != 11) {
              for(ll=0; ll<3*lbsy.nf; ll++)
	            interforce [ll] = lbinterforce[il*3*lbsy.nf+ll] + postequil*lbbdforce[ll];
              if(lbphi[il] != 12 && lbphi[il] != 13)
                fSiteFluidCollisionBGK(&lbf[il*lbsitelength], sitespeed, interforce);
              fSiteSoluteCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*lbsy.nf], sitespeed);
              fSiteThermalCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*(lbsy.nf+lbsy.nc)], sitespeed);
            }
          }
        }
      }
    }
  }
  else {
      
   #pragma omp parallel private(i,j,k,il,ll)
   {
      double interforce_local_t[3*lbsy.nf], sitespeed_local_t[3];
      double *interforce = &interforce_local_t[0];
      double *sitespeed = &sitespeed_local_t[0];

      #pragma omp for collapse(2)
      for(int i=0; i<lbdm.xouter; i++) {
        for(int j=0; j<lbdm.youter; j++) {
          for(int k=0; k<lbdm.zouter; k++) {
            il = (i * lbdm.youter + j) * lbdm.zouter + k;
            fGetSpeedIncomSite(sitespeed, &lbf[il*lbsitelength]);
            if(lbphi[il] != 11) {
              for(ll=0; ll<3*lbsy.nf; ll++)
  	            interforce[ll] = lbinterforce[il*3*lbsy.nf+ll] + postequil*lbbdforce[ll];
  	          if(lbphi[il] != 12 && lbphi[il] != 13)
	            fSiteFluidIncomCollisionBGK(&lbf[il*lbsitelength], sitespeed, interforce);
              fSiteSoluteCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*lbsy.nf], sitespeed);
              fSiteThermalCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*(lbsy.nf+lbsy.nc)], sitespeed);
            }
          }
        }
      }
    }
  }
  return 0;
}



int fCollisionBGKEDM()
{
  long il;
  int Xmax = (long) lbdm.xouter, Ymax = (long) lbdm.youter, Zmax = (long) lbdm.zouter;
  int i, j, k, ll;

  if(!incompress) {

    #pragma omp parallel private(i,j,k,il,ll)
    {
      double interforce_t[3*lbsy.nf], sitespeed_t[3];
      double *interforce = &interforce_t[0];
      double *sitespeed = &sitespeed_t[0];
                                     
      #pragma omp for collapse(2)
      for(i=0; i<Xmax; i++) {
        for(j=0; j<Ymax; j++) {
          for(k=0; k<Zmax; k++) {
            il = (i * Ymax + j) * Zmax + k;
            fGetSpeedSite(sitespeed, &lbf[il*lbsitelength]);
            if(lbphi[il] != 11) {
              for(ll=0; ll<3*lbsy.nf; ll++)
	            interforce[ll] = lbinterforce[il*3*lbsy.nf+ll] + postequil*lbbdforce[ll];
              if(lbphi[il] != 12 && lbphi[il] != 13)
  	            fSiteFluidCollisionBGKEDM(&lbf[il*lbsitelength], sitespeed, interforce);
              fSiteSoluteCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*lbsy.nf], sitespeed);
              fSiteThermalCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*(lbsy.nf+lbsy.nc)], sitespeed);
            }
          }
        }
      }
    }
  }
  else {

    #pragma omp parallel private(i,j,k,il,ll)
    {
      double interforce_t[3*lbsy.nf], sitespeed_t[3];
      double *interforce = &interforce_t[0];
      double *sitespeed = &sitespeed_t[0];
                                     
      #pragma omp for collapse(2)
      for(i=0; i<Xmax; i++) {
        for(j=0; j<Ymax; j++) {
          for(k=0; k<Zmax; k++) {
            il = (i * Ymax + j) * Zmax + k;
            fGetSpeedIncomSite(sitespeed, &lbf[il*lbsitelength]);
            if(lbphi[il] != 11) {
              for(ll=0; ll<3*lbsy.nf; ll++)
	            interforce[ll] = lbinterforce[il*3*lbsy.nf+ll] + postequil*lbbdforce[ll];
              if(lbphi[il] != 12 && lbphi[il] != 13)
  	            fSiteFluidIncomCollisionBGKEDM(&lbf[il*lbsitelength], sitespeed, interforce);
              fSiteSoluteCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*lbsy.nf], sitespeed);
              fSiteThermalCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*(lbsy.nf+lbsy.nc)], sitespeed);
            }
          }
        }
      }
    }
  }
  return 0;
}



int fCollisionBGKGuo()
{
  long il;
  int Xmax = (long) lbdm.xouter, Ymax = (long) lbdm.youter, Zmax = (long) lbdm.zouter;
  int i, j, k, ll;

  if(!incompress) {

    #pragma omp parallel private(i,j,k,il,ll)
    {
      double interforce_t[3*lbsy.nf], sitespeed_t[3];
      double *interforce = &interforce_t[0];
      double *sitespeed = &sitespeed_t[0];
                                     
      #pragma omp for collapse(2)
      for(i=0; i<Xmax; i++) {
        for(j=0; j<Ymax; j++) {
          for(k=0; k<Zmax; k++) {
            il = (i * Ymax + j) * Zmax + k;
            fGetSpeedSite(sitespeed, &lbf[il*lbsitelength]);
            if(lbphi[il] != 11) {
              for(ll=0; ll<3*lbsy.nf; ll++)
	            interforce[ll] = lbinterforce[il*3*lbsy.nf+ll] + postequil*lbbdforce[ll];
              if(lbphi[il] != 12 && lbphi[il] != 13)
  	            fSiteFluidCollisionBGKGuo(&lbf[il*lbsitelength], sitespeed, interforce);
              fSiteSoluteCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*lbsy.nf], sitespeed);
              fSiteThermalCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*(lbsy.nf+lbsy.nc)], sitespeed);
            }
          }
        }
      }
    }
  }
  else {

    #pragma omp parallel private(i,j,k,il,ll)
    {
      double interforce_t[3*lbsy.nf], sitespeed_t[3];
      double *interforce = &interforce_t[0];
      double *sitespeed = &sitespeed_t[0];
                                     
      #pragma omp for collapse(2)
      for(i=0; i<Xmax; i++) {
        for(j=0; j<Ymax; j++) {
          for(k=0; k<Zmax; k++) {
            il = (i * Ymax + j) * Zmax + k;
            fGetSpeedIncomSite(sitespeed, &lbf[il*lbsitelength]);
            if(lbphi[il] != 11) {
              for(ll=0; ll<3*lbsy.nf; ll++)
	            interforce[ll] = lbinterforce[il*3*lbsy.nf+ll] + postequil*lbbdforce[ll];
              if(lbphi[il] != 12 && lbphi[il] != 13)
  	            fSiteFluidIncomCollisionBGKGuo(&lbf[il*lbsitelength], sitespeed, interforce);
              fSiteSoluteCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*lbsy.nf], sitespeed);
              fSiteThermalCollisionBGK(&lbf[il*lbsitelength + lbsy.nq*(lbsy.nf+lbsy.nc)], sitespeed);
            }
          }
        }
      }
    }
  }
  return 0;
}


