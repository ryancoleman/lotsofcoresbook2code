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


long fNextStep(int q, int xpos, int ypos, int zpos)
{

  // find position of particle for next time step when moving along q 
  // direction to nearest neighbour

  q *= 3;
  return (fCppMod(xpos+lbv[q], lbdm.xouter) * lbdm.youter 
    + fCppMod(ypos+lbv[q+1], lbdm.youter)) * lbdm.zouter 
    + fCppMod(zpos+lbv[q+2], lbdm.zouter);
}


long fNextStep(int dx, int dy, int dz, long tpos)
{

  // find position of particle for next time step when moving along q
  // direction to nearest neighbour
  int xpos, ypos, zpos;
  fGetCoord(tpos, xpos, ypos, zpos);
  return (fCppMod(xpos+dx, lbdm.xouter) * lbdm.youter 
    + fCppMod(ypos+dy, lbdm.youter)) * lbdm.zouter 
    + fCppMod(zpos+dz, lbdm.zouter);
}




int fBounceBackF(long tpos)
{

  // perform on-grid bounce-back for fluid distribution function
  
  long stpos=tpos * lbsitelength;
  int half = (lbsy.nq - 1) / 2;

  for(int j=0; j<lbsy.nf; j++)
    for(int m=1; m<=half; m++)
        fSwapPair (lbf[stpos + j*lbsy.nq + m], lbf[stpos + j*lbsy.nq + m + half]);
    
/*
  for(int j=0; j<lbsy.nf; j++) {
    for(int m=0; m<lbsy.nq; m++) 
      lbfeq[lbopv[m]] = lbf[stpos + j*lbsy.nq + m];
    for(int m=0; m<lbsy.nq; m++)
      lbf[stpos + j*lbsy.nq + m] = lbfeq[m];
  } 
*/
    
  return 0;
}


int fBounceBackC(long tpos)
{

  // perform on-grid bounce-back for solute distribution function
  
  long stpos=tpos * lbsitelength + lbsy.nf * lbsy.nq;
  int half = (lbsy.nq - 1) / 2;

  for(int j=0; j<lbsy.nc; j++)
    for(int m=1; m<=half; m++)
      fSwapPair (lbf[stpos + j*lbsy.nq + m], lbf[stpos + j*lbsy.nq + m + half]);
        
/*
  for(int j=0; j<lbsy.nc; j++) {
    for(int m=0; m<lbsy.nq; m++) 
      lbfeq[lbopv[m]] = lbf[stpos +j*lbsy.nq+ m];
    for(int m=0; m<lbsy.nq; m++)
      lbf[stpos +j*lbsy.nq+ m] = lbfeq[m];
  }    
*/

  return 0;
}


int fBounceBackT(long tpos)
{

  // perform on-grid bounce-back for temperature distribution function
  
  long stpos=tpos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq;
  int half = (lbsy.nq - 1) / 2;

  for(int j=0; j<lbsy.nt; j++)
    for(int m=1; m<=half; m++)
      fSwapPair (lbf[stpos + m], lbf[stpos + m + half]);
  
/*  for(int j=0; j<lbsy.nt; j++) {
      for(int m=0; m<lbsy.nq; m++)
        lbfeq[lbopv[m]] = lbf[stpos + m];
      for(int m=0; m<lbsy.nq; m++)
        lbf[stpos + m] = lbfeq[m];
      stpos += lbsy.nq;
    }
*/
    
  return 0;
}

int fMidBounceBackF(long tpos)
{

  // perform mid-link bounce-back for fluid distribution function
  
  long stpos= tpos * lbsitelength;
  long spos;
  int xpos, ypos, zpos;

  fGetCoord(tpos, xpos, ypos, zpos); 
  for(int j=0; j<lbsy.nf; j++) {

    for(int m=0; m<lbsy.nq; m++) {
      spos = fNextStep(m, xpos, ypos, zpos);
      lbf[stpos + j * lbsy.nq + m] = lbf[spos * lbsitelength + j * lbsy.nq + lbopv[m]];
    }

  }
  return 0;
}


int fMidBounceBackC(long tpos)
{

  // perform mid-link bounce-back for solute distribution function
  
  long stpos=tpos * lbsitelength + lbsy.nf * lbsy.nq;
  long spos;
  int xpos, ypos, zpos;

  fGetCoord(tpos, xpos, ypos, zpos);   
  for(int j=0; j<lbsy.nc; j++) {
 
    for(int m=0; m<lbsy.nq; m++) {
      spos = fNextStep(m, xpos, ypos, zpos);
      lbf[stpos + j * lbsy.nq + m] = lbf[spos * lbsitelength + (lbsy.nf + j) * lbsy.nq + lbopv[m]];
    }

  }    
  
  return 0;
}


int fMidBounceBackT(long tpos)
{

  // perform mid-link bounce-back for temperature distribution function
  
  long stpos=tpos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq;
  long spos;
  int xpos, ypos, zpos;

  fGetCoord(tpos, xpos, ypos, zpos);   
  for(int j=0; j<lbsy.nt; j++) {
     
    for(int m=0; m<lbsy.nq; m++) {
      spos = fNextStep(m, xpos, ypos, zpos);
      lbf[stpos + m] = lbf[spos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq + lbopv[m]];
    }

  }    
  
  return 0;
}

int fSiteBlankF(long tpos)
{

  // set fluid distribution function at grid point tpos to zero (solid block)

  double *pt2 = &lbf[tpos * lbsitelength];
  for(int j=0; j<lbsy.nf * lbsy.nq; j++){
    *pt2 =0;
    pt2 ++;
  }
  pt2 = NULL;
  return 0;
}


int fSiteBlankC(long tpos)
{

  // set solute distribution function at grid point tpos to zero (solid block)

  double *pt2 = &lbf[tpos * lbsitelength + lbsy.nf * lbsy.nq];
  for(int j=0; j<lbsy.nc * lbsy.nq; j++){
    *pt2 =0;
    pt2 ++;
  }
  pt2 = NULL;
  return 0;
}

int fSiteBlankT(long tpos)
{

  // set temperature distribution function at grid point tpos to zero (solid 
  // block)
  
  double *pt2 = &lbf[tpos * lbsitelength + (lbsy.nf+lbsy.nc)*lbsy.nq];
  for(int j=0; j<lbsy.nt * lbsy.nq; j++){
    *pt2 =0;
    pt2 ++;
  }
  pt2 = NULL;
  return 0;
}

// D2Q9

int fD2Q9VCE(double v0, double v1, double &f0, double &f1, double &f2,
	     double &f3, double &f4, double &f5, double &f6,
	     double &f7, double &f8)
{

  // produce fixed velocity at top boundary (PSDF) for compressible fluid

  double c1=2.0/3.0,c2=1.0/6.0;
  double rho = (f0+f2+f6+2*(f1+f7+f8))/(1.0+v1);
  f4=f8-c1*rho*v1;
  f3=f7+0.5*(f6-f2)-0.5*rho*v0-rho*v1*c2;
  f5=f1-0.5*(f6-f2)+0.5*rho*v0-rho*v1*c2;
  return 0;
}

int fD2Q9VCEIncom(double v0, double v1, double rho0, double &f0,
                 double &f1, double &f2, double &f3, double &f4,
                 double &f5, double &f6, double &f7, double &f8)
{

  // produce fixed velocity at top boundary (PSDF) for compressible fluid

  double c1=2.0/3.0,c2=1.0/6.0;
//  double rho = f0+f2+f6+2*(f1+f7+f8)-rho0*v1;
  f4=f8-c1*rho0*v1;
  f3=f7+0.5*(f6-f2)-0.5*rho0*v0-rho0*v1*c2;
  f5=f1-0.5*(f6-f2)+0.5*rho0*v0-rho0*v1*c2;
  return 0;
}

int fD2Q9VCC(double p, double *v, double * startpos)
{

  // produce fixed velocity at concave corner for compressible fluid

  double *pt1=startpos;
  fGetEquilibriumF(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD2Q9VCCIncom(double p, double p0, double *v, double * startpos)
{

  // produce fixed velocity at concave corner for incompressible fluid

  double *pt1=startpos;
  fGetEquilibriumFIncom(lbfeq, v, p, p0);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD2Q9VF(long tpos, int prop, double *uwall)
{
  double rho, rho0;
  long spos = tpos * lbsitelength, tpos1;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 49) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = 0.0;
        fD2Q9VCE(uwall[0], uwall[1], lbf[spos], lbf[spos+1], lbf[spos+2],
                 lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6], lbf[spos+7], lbf[spos+8]);
      }
      else if(prop == 47) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = 0.0;
        fD2Q9VCE(-uwall[0], -uwall[1], lbf[spos], lbf[spos+5], lbf[spos+6],
                 lbf[spos+7], lbf[spos+8], lbf[spos+1], lbf[spos+2], lbf[spos+3], lbf[spos+4]);
      }
      else if(prop == 50) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = 0.0;
        fD2Q9VCE(uwall[1], -uwall[0], lbf[spos], lbf[spos+3], lbf[spos+4],
                 lbf[spos+5], lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+1], lbf[spos+2]);
      }
      else if(prop == 48) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = 0.0;
        fD2Q9VCE(-uwall[1], uwall[0], lbf[spos], lbf[spos+7], lbf[spos+8],
                 lbf[spos+1], lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6]);
      }
      else if(prop == 31) {
        tpos1 = fNextStep(1, 1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = 0.0;
        fD2Q9VCC(rho, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        tpos1 = fNextStep(-1, 1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = 0.0;
        fD2Q9VCC(rho, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        tpos1 = fNextStep(-1, -1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = 0.0;
        fD2Q9VCC(rho, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        tpos1 = fNextStep(1, -1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = 0.0;
        fD2Q9VCC(rho, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 49) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = 0.0;
        fD2Q9VCEIncom(uwall[0], uwall[1], rho0, lbf[spos], lbf[spos+1], lbf[spos+2],
                      lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6], lbf[spos+7],
                      lbf[spos+8]);
      }
      else if(prop == 47) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = 0.0;
        fD2Q9VCEIncom(-uwall[0], -uwall[1], rho0, lbf[spos], lbf[spos+5], lbf[spos+6],
                      lbf[spos+7], lbf[spos+8], lbf[spos+1], lbf[spos+2], lbf[spos+3],
                      lbf[spos+4]);
      }
      else if(prop == 50) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = 0.0;
        fD2Q9VCEIncom(uwall[1], -uwall[0], rho0, lbf[spos], lbf[spos+3], lbf[spos+4],
                      lbf[spos+5], lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+1],
                      lbf[spos+2]);
      }
      else if(prop == 48) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = 0.0;
        fD2Q9VCEIncom(-uwall[1], uwall[0], rho0, lbf[spos], lbf[spos+7], lbf[spos+8],
                      lbf[spos+1], lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                      lbf[spos+6]);
      }
      else if(prop == 31) {
        tpos1 = fNextStep(1, 1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = 0.0;
        fD2Q9VCCIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        tpos1 = fNextStep(-1, 1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = 0.0;
        fD2Q9VCCIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        tpos1 = fNextStep(-1, -1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = 0.0;
        fD2Q9VCCIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        tpos1 = fNextStep(1, -1, 0, tpos);
        rho = fEvapLimit(fGetOneMassSite(&lbf[tpos1*lbsitelength+ktt*lbsy.nq]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = 0.0;
        fD2Q9VCCIncom(rho, rho0, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  return 0;  
}

int fD2Q9PCE(double p, double &f0, double &f1, double &f2, double &f3, 
	     double &f4, double &f5, double &f6, double &f7, double &f8, double &vel)
{

  // produce fixed density/pressure boundary for concave edge (PCED)

  double c1=2.0/3.0,c2=1.0/6.0;
  vel=f0+f2+f6+2*(f1+f7+f8)-p;  // vel = p*v1
  f4=f8-c1*vel;
  f3=f7+0.5*(f6-f2)-vel*c2;
  f5=f1-0.5*(f6-f2)-vel*c2;
  return 0;
}


int fD2Q9PCC(double p, double *v, double * startpos)
{

  // produce fixed velocity at concave corner for compressible fluid

  double *pt1=startpos;
  fGetEquilibriumF(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}


int fD2Q9PCCIncom(double p, double p0, double *v, double * startpos)
{

  // produce fixed velocity at concave corner for incompressible fluid

  double *pt1=startpos;
  fGetEquilibriumFIncom(lbfeq, v, p, p0);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}


int fD2Q9PF(long tpos, int prop, double *uwall)
{
  double rho0,moment,onemass,mass=0.0;
  long spos = tpos * lbsitelength;
  uwall[0]=0.0; uwall[1]=0.0; uwall[2]=0.0;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 49) {
        onemass=lbtopp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+1], lbf[spos+2],
	         lbf[spos+3], lbf[spos+4], lbf[spos+5], 
	         lbf[spos+6], lbf[spos+7], lbf[spos+8], moment);
        mass += onemass;
        uwall[1] += moment;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+5], lbf[spos+6],
	         lbf[spos+7], lbf[spos+8], lbf[spos+1],
                 lbf[spos+2], lbf[spos+3], lbf[spos+4], moment);
        mass += onemass;
        uwall[1] -= moment;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+3], lbf[spos+4],
	         lbf[spos+5], lbf[spos+6], lbf[spos+7], 
                 lbf[spos+8], lbf[spos+1], lbf[spos+2], moment);
        mass += onemass;
        uwall[0] -= moment;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+7], lbf[spos+8],
	         lbf[spos+1], lbf[spos+2], lbf[spos+3],
                 lbf[spos+4], lbf[spos+5], lbf[spos+6], moment);
        mass += onemass;
        uwall[0] += moment;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD2Q9PCC(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD2Q9PCC(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD2Q9PCC(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD2Q9PCC(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 49) {
        onemass=lbtopp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+1], lbf[spos+2],
	         lbf[spos+3], lbf[spos+4], lbf[spos+5],
  	         lbf[spos+6], lbf[spos+7], lbf[spos+8], moment);
        mass += rho0;
        uwall[1] += moment;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+5], lbf[spos+6],
	         lbf[spos+7], lbf[spos+8], lbf[spos+1],
                 lbf[spos+2], lbf[spos+3], lbf[spos+4], moment);
        mass += rho0;
        uwall[1] -= moment;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+3], lbf[spos+4],
	         lbf[spos+5], lbf[spos+6], lbf[spos+7], 
                 lbf[spos+8], lbf[spos+1], lbf[spos+2], moment);
        mass += rho0;
        uwall[0] -= moment;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD2Q9PCE(onemass, lbf[spos], lbf[spos+7], lbf[spos+8],
  	         lbf[spos+1], lbf[spos+2], lbf[spos+3],
                 lbf[spos+4], lbf[spos+5], lbf[spos+6], moment);
        mass += rho0;
        uwall[0] += moment;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD2Q9PCCIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD2Q9PCCIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD2Q9PCCIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD2Q9PCCIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      spos += lbsy.nq;
    }
  }
  uwall[0] *= fReciprocal(mass);
  uwall[1] *= fReciprocal(mass);
  return 0;  
}


int fD2Q9CCE(double p, double v0, double v1, double &f0, double &f1, double &f2,
             double &f3, double &f4, double &f5, double &f6, double &f7, double &f8)
{

  // produce fixed concentration boundary for concave edge (CCED)

  double c1=1.0/9.0,c2=1.0/36.0;
  double rho=6.0*(p-f0-f1-f2-f6-f7-f8)/(1.0-3.0*v1);
  f4=c1*rho*(1.0-3.0*v1);
  f3=c2*rho*(1.0-3.0*v0-3.0*v1);
  f5=c2*rho*(1.0+3.0*v0-3.0*v1);
  return 0;
}


int fD2Q9CCC(double p, double *v, double * startpos)
{

  // produce fixed solute concentration at concave corner

  double *pt1=startpos;
  fGetEquilibriumC(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD2Q9PC(long tpos, int prop, double *uwall)
{
  double oneconc;
  long spos = tpos * lbsitelength + lbsy.nf * lbsy.nq;
  for(int ktt=0; ktt<lbsy.nc; ktt++) {
    if(prop == 49) {
      oneconc=lbtopc[ktt];
      fD2Q9CCE(oneconc, uwall[0], uwall[1], lbf[spos], lbf[spos+1],
               lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
               lbf[spos+6], lbf[spos+7], lbf[spos+8]);
    }
    else if(prop == 47) {
      oneconc=lbbotc[ktt];
      fD2Q9CCE(oneconc, -uwall[0], -uwall[1], lbf[spos], lbf[spos+5],
               lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+1],
               lbf[spos+2], lbf[spos+3], lbf[spos+4]);
    }
    else if(prop == 50) {
      oneconc=lblefc[ktt];
      fD2Q9CCE(oneconc, uwall[1], -uwall[0], lbf[spos], lbf[spos+3],
               lbf[spos+4], lbf[spos+5], lbf[spos+6], lbf[spos+7],
               lbf[spos+8], lbf[spos+1], lbf[spos+2]);
    }
    else if(prop == 48) {
      oneconc=lbrigc[ktt];
      fD2Q9CCE(oneconc, -uwall[1], uwall[0], lbf[spos], lbf[spos+7],
               lbf[spos+8], lbf[spos+1], lbf[spos+2], lbf[spos+3],
               lbf[spos+4], lbf[spos+5], lbf[spos+6]);
    }
    else if(prop == 31) {
      oneconc=lblefc[ktt];
      fD2Q9CCC(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      oneconc=lbrigc[ktt];
      fD2Q9CCC(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      oneconc=lbrigc[ktt];
      fD2Q9CCC(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      oneconc=lblefc[ktt];
      fD2Q9CCC(oneconc, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}


int fD2Q9TCC(double p, double *v, double * startpos)
{

  // produce fixed temperature at concave corner

  double *pt1=startpos;
  fGetEquilibriumT(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}


int fD2Q9PT(long tpos, int prop, double *uwall)
{
  double onetemp;
  long spos = tpos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq;
  if(lbsy.nt == 1) {
    if(prop == 49) {
      onetemp=lbtopt;
      fD2Q9CCE(onetemp, uwall[0], uwall[1], lbf[spos], lbf[spos+1],
               lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
	       lbf[spos+6], lbf[spos+7], lbf[spos+8]);
    }
    else if(prop == 47) {
      onetemp=lbbott;
      fD2Q9CCE(onetemp, -uwall[0], -uwall[1], lbf[spos], lbf[spos+5],
               lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+1],
               lbf[spos+2], lbf[spos+3], lbf[spos+4]);
    }
    else if(prop == 50) {
      onetemp=lbleft;
      fD2Q9CCE(onetemp, uwall[1], -uwall[0], lbf[spos], lbf[spos+3],
               lbf[spos+4], lbf[spos+5], lbf[spos+6], lbf[spos+7],
               lbf[spos+8], lbf[spos+1], lbf[spos+2]);
    }
    else if(prop == 48) {
      onetemp=lbrigt;
      fD2Q9CCE(onetemp, -uwall[1], uwall[0], lbf[spos], lbf[spos+7],
               lbf[spos+8], lbf[spos+1], lbf[spos+2], lbf[spos+3],
               lbf[spos+4], lbf[spos+5], lbf[spos+6]);
    }
    else if(prop == 31) {
      onetemp=lbleft;
      fD2Q9TCC(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      onetemp=lbrigt;
      fD2Q9TCC(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      onetemp=lbrigt;
      fD2Q9TCC(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      onetemp=lbleft;
      fD2Q9TCC(onetemp, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}


// D3Q15

int fD3Q15VPS(double v0, double v1, double v2, double &f0, double &f1, double &f2,
              double &f3, double &f4, double &f5, double &f6, double &f7, double &f8,
              double &f9, double &f10, double &f11, double &f12, double &f13, double &f14)
{

  // produce fixed velocity at planar surface boundary for compressible fluid:
  // expressed for top wall
  
  double rho, n1, n2;
  double c1=2.0/3.0,c2=1.0/12.0,c3=1.0/6.0;
  rho=(f0+f1+f3+f8+f10+2.0*(f6+f7+f9+f11+f12))/(1+v1);
  n1=0.25*(f8-f1)-c3*rho*v0;
  n2=0.25*(f10-f3)-c3*rho*v2;
  f2=f9-c1*rho*v1;
  f4=f11-c2*rho*(v0+v1+v2)+n1+n2;
  f5=f12-c2*rho*(v0+v1-v2)+n1-n2;  
  f13=f6+c2*rho*(v0-v1+v2)-n1-n2;
  f14=f7+c2*rho*(v0-v1-v2)-n1+n2;
  return 0;
}

int fD3Q15VPSIncom(double v0, double v1, double v2, double rho0, double &f0, double &f1,
	           double &f2, double &f3, double &f4, double &f5, double &f6, double &f7,
                   double &f8, double &f9, double &f10, double &f11, double &f12, double &f13,
	           double &f14)
{

  // produce fixed velocity at planar surface boundary for incompressible fluid:
  // expressed for top wall
  
  double n1, n2;
  double c1=2.0/3.0,c2=1.0/12.0,c3=1.0/6.0;
//  double rho=f0+f1+f3+f8+f10+2.0*(f6+f7+f9+f11+f12)-rho0*v1;
  n1=0.25*(f8-f1)-c3*rho0*v0;
  n2=0.25*(f10-f3)-c3*rho0*v2;
  f2=f9-c1*rho0*v1;
  f4=f11-c2*rho0*(v0+v1+v2)+n1+n2;
  f5=f12-c2*rho0*(v0+v1-v2)+n1-n2;  
  f13=f6+c2*rho0*(v0-v1+v2)-n1-n2;
  f14=f7+c2*rho0*(v0-v1-v2)-n1+n2;
  return 0;
}

int fD3Q15VCE(double p, double *v, double * startpos)
{

  // produce fixed velocity at convex or concave corner for compressible fluid:
  // expressed for left-bottom-back corner (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumF(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q15VCEIncom(double p, double p0, double *v, double * startpos)
{

  // produce fixed velocity at convex or concave corner for incompressible fluid:
  // expressed for left-bottom-back corner (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumFIncom(lbfeq, v, p, p0);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q15VF(long tpos, int prop, double *uwall)
{
  long spos=tpos * lbsitelength;
  double rho, rho0;
  long rpos;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 22) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VPS(uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14]);
      }
      else if(prop == 21) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VPS(uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                  lbf[spos+9], lbf[spos+3], lbf[spos+6], lbf[spos+7],
                  lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+2],
                  lbf[spos+10], lbf[spos+13], lbf[spos+14], lbf[spos+11],
                  lbf[spos+12]);
      }
      else if(prop == 23) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q15VPS(uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                  lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+14], lbf[spos+13], lbf[spos+9], lbf[spos+8],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+7],
                  lbf[spos+6]);
      }
      else if(prop == 24) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q15VPS(uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                  lbf[spos+8], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+1],
                  lbf[spos+10], lbf[spos+6], lbf[spos+7], lbf[spos+11],
                  lbf[spos+12]);
      }
      else if(prop == 25) {
        uwall[0] = lbbacv[0]; uwall[1] = lbbacv[1]; uwall[2] = lbbacv[2];
        fD3Q15VPS(uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                  lbf[spos+10], lbf[spos+2], lbf[spos+5], lbf[spos+7],
                  lbf[spos+4], lbf[spos+6], lbf[spos+8], lbf[spos+3],
                  lbf[spos+9], lbf[spos+12], lbf[spos+14], lbf[spos+11],
                  lbf[spos+13]);
      }
      else if(prop == 26) {
        uwall[0] = lbfrov[0]; uwall[1] = lbfrov[1]; uwall[2] = lbfrov[2];
        fD3Q15VPS(uwall[0], uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+2], lbf[spos+4], lbf[spos+6],
                  lbf[spos+5], lbf[spos+7], lbf[spos+8], lbf[spos+10],
                  lbf[spos+9], lbf[spos+11], lbf[spos+13], lbf[spos+12],
                  lbf[spos+14]);
      }
      else if(prop == 27) {
        rpos = fNextStep(1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 28) {
        rpos = fNextStep(-1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 29) {
        rpos = fNextStep(-1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 30) {
        rpos = fNextStep(1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 31) {
        rpos = fNextStep(1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        rpos = fNextStep(-1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        rpos = fNextStep(-1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        rpos = fNextStep(1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 43) {
        rpos = fNextStep(1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 44) {
        rpos = fNextStep(-1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 45) {
        rpos = fNextStep(-1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 46) {
        rpos = fNextStep(1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 47) {
        rpos = fNextStep(0, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 48) {
        rpos = fNextStep(-1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 49) {
        rpos = fNextStep(0, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 50) {
        rpos = fNextStep(1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 51) {
        rpos = fNextStep(0, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 52) {
        rpos = fNextStep(-1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 53) {
        rpos = fNextStep(0, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 54) {
        rpos = fNextStep(1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q15VCE(rho, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 22) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VPSIncom(uwall[0], uwall[1], uwall[2], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                       lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                       lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                       lbf[spos+14]);
      }
      else if(prop == 21) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VPSIncom(uwall[0], -uwall[1], uwall[2], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+9], lbf[spos+3], lbf[spos+6], lbf[spos+7],
                       lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+2],
                       lbf[spos+10], lbf[spos+13], lbf[spos+14], lbf[spos+11],
                       lbf[spos+12]);
      }
      else if(prop == 23) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q15VPSIncom(uwall[1], uwall[0], uwall[2], rho0, lbf[spos], lbf[spos+2],
                       lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                       lbf[spos+14], lbf[spos+13], lbf[spos+9], lbf[spos+8],
                       lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+7],
                       lbf[spos+6]);
      }
      else if(prop == 24) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q15VPSIncom(uwall[1], -uwall[0], uwall[2], rho0, lbf[spos], lbf[spos+2],
                       lbf[spos+8], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                       lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+1],
                       lbf[spos+10], lbf[spos+6], lbf[spos+7], lbf[spos+11],
                       lbf[spos+12]);
      }
      else if(prop == 25) {
        uwall[0] = lbbacv[0]; uwall[1] = lbbacv[1]; uwall[2] = lbbacv[2];
        fD3Q15VPSIncom(uwall[0], -uwall[2], uwall[1], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+10], lbf[spos+2], lbf[spos+5], lbf[spos+7],
                       lbf[spos+4], lbf[spos+6], lbf[spos+8], lbf[spos+3],
                       lbf[spos+9], lbf[spos+12], lbf[spos+14], lbf[spos+11],
                       lbf[spos+13]);
      }
      else if(prop == 26) {
        uwall[0] = lbfrov[0]; uwall[1] = lbfrov[1]; uwall[2] = lbfrov[2];
        fD3Q15VPSIncom(uwall[0], uwall[2], uwall[1], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+3], lbf[spos+2], lbf[spos+4], lbf[spos+6],
                       lbf[spos+5], lbf[spos+7], lbf[spos+8], lbf[spos+10],
                       lbf[spos+9], lbf[spos+11], lbf[spos+13], lbf[spos+12],
                       lbf[spos+14]);
      }
      else if(prop == 27) {
        rpos = fNextStep(1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 28) {
        rpos = fNextStep(-1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 29) {
        rpos = fNextStep(-1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 30) {
        rpos = fNextStep(1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 31) {
        rpos = fNextStep(1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        rpos = fNextStep(-1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        rpos = fNextStep(-1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        rpos = fNextStep(1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 43) {
        rpos = fNextStep(1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 44) {
        rpos = fNextStep(-1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 45) {
        rpos = fNextStep(-1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 46) {
        rpos = fNextStep(1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 47) {
        rpos = fNextStep(0, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 48) {
        rpos = fNextStep(-1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 49) {
        rpos = fNextStep(0, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 50) {
        rpos = fNextStep(1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 51) {
        rpos = fNextStep(0, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 52) {
        rpos = fNextStep(-1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 53) {
        rpos = fNextStep(0, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 54) {
        rpos = fNextStep(1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q15VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  return 0;
  
}

int fD3Q15PPS(double p, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double& vel)
{

  // produce fixed pressure/density at planar surface: expressed for top wall

  double rx, rz;
  double c1=2.0/3.0,c2=1.0/12.0;
  rx=0.25*(f8-f1);
  vel=(f0+f1+f3+f8+f10+2.0*(f6+f7+f9+f11+f12)) -p;
  rz=0.25*(f10-f3);
  f2=f9-c1*vel;
  f4=f11-c2*vel+rx+rz; 
  f5=f12-c2*vel+rx-rz;
  f13=f6-c2*vel-rx-rz; 
  f14=f7-c2*vel-rx+rz; 
  
  return 0;
}


int fD3Q15PF(long tpos, int prop, double *uwall)
{
  double moment,rho0,onemass,mass=0.0;
  long spos=tpos * lbsitelength;
  uwall[0]=0.0; uwall[1]=0.0; uwall[2]=0.0;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 22) {
        onemass=lbtopp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], moment);
        mass += onemass;
        uwall[1] += moment;     
      }    
      else if(prop == 21) {
        onemass=lbbotp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+9], lbf[spos+3], lbf[spos+6], lbf[spos+7],
                  lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+2],
                  lbf[spos+10], lbf[spos+13], lbf[spos+14], lbf[spos+11],
                  lbf[spos+12], moment);
        mass += onemass;
        uwall[1] -= moment;     
      }       
      else if(prop == 23) {
        onemass=lbrigp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+14], lbf[spos+13], lbf[spos+9], lbf[spos+8],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+7],
                  lbf[spos+6], moment);
        mass += onemass;
        uwall[0] += moment;     
      }    
      else if(prop == 24) {
        onemass=lblefp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+8], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+1],
                  lbf[spos+10], lbf[spos+6], lbf[spos+7], lbf[spos+11],
                  lbf[spos+12], moment);
        mass += onemass;
        uwall[1] -= moment;     
      }      
      else if(prop == 25) {
        onemass=lbbacp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+10], lbf[spos+2], lbf[spos+5], lbf[spos+7],
                  lbf[spos+4], lbf[spos+6], lbf[spos+8], lbf[spos+3],
                  lbf[spos+9], lbf[spos+12], lbf[spos+14], lbf[spos+11],
                  lbf[spos+13], moment);
        mass += onemass;
        uwall[2] -= moment;     
      }    
      else if(prop == 26) {
        onemass=lbfrop[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+2], lbf[spos+4], lbf[spos+6],
                  lbf[spos+5], lbf[spos+7], lbf[spos+8], lbf[spos+10],
                  lbf[spos+9], lbf[spos+11], lbf[spos+13], lbf[spos+12],
                  lbf[spos+14], moment);
        mass += onemass;
        uwall[2] += moment;     
      }    
      else if(prop == 27) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 28) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 29) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 30) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 43) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 44) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 45) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 46) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 49) {
        onemass=lbtopp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 51) {
        onemass=lbbotp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 52) {
        onemass=lbrigp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 53) {
        onemass=lbtopp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 54) {
        onemass=lblefp[ktt];
        fD3Q15VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 22) {
        onemass=lbtopp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], moment);
        mass += rho0;
        uwall[1] += moment;     
      }    
      else if(prop == 21) {
        onemass=lbbotp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+9], lbf[spos+3], lbf[spos+6], lbf[spos+7],
                  lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+2],
                  lbf[spos+10], lbf[spos+13], lbf[spos+14], lbf[spos+11],
                  lbf[spos+12], moment);
        mass += rho0;
        uwall[1] -= moment;     
      }       
      else if(prop == 23) {
        onemass=lbrigp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+14], lbf[spos+13], lbf[spos+9], lbf[spos+8],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+7],
                  lbf[spos+6], moment);
        mass += rho0;
        uwall[0] += moment;     
      }    
      else if(prop == 24) {
        onemass=lblefp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+8], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+1],
                  lbf[spos+10], lbf[spos+6], lbf[spos+7], lbf[spos+11],
                  lbf[spos+12], moment);
        mass += rho0;
        uwall[0] -= moment;     
      }      
      else if(prop == 25) {
        onemass=lbbacp[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+10], lbf[spos+2], lbf[spos+5], lbf[spos+7],
                  lbf[spos+4], lbf[spos+6], lbf[spos+8], lbf[spos+3],
                  lbf[spos+9], lbf[spos+12], lbf[spos+14], lbf[spos+11],
                  lbf[spos+13], moment);
        mass += rho0;
        uwall[2] -= moment;     
      }    
      else if(prop == 26) {
        onemass=lbfrop[ktt];
        fD3Q15PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+2], lbf[spos+4], lbf[spos+6],
                  lbf[spos+5], lbf[spos+7], lbf[spos+8], lbf[spos+10],
                  lbf[spos+9], lbf[spos+11], lbf[spos+13], lbf[spos+12],
                  lbf[spos+14], moment);
        mass += rho0;
        uwall[2] += moment;     
      }    
      else if(prop == 27) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 28) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 29) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 30) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 43) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 44) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 45) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 46) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 47) {
        onemass=lbtopp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 49) {
        onemass=lbtopp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 51) {
        onemass=lbbotp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 52) {
        onemass=lbrigp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 53) {
          onemass=lbtopp[ktt];
          fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
          mass += rho0;
      }
      else if(prop == 54) {
        onemass=lblefp[ktt];
        fD3Q15VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      spos += lbsy.nq;
    }
  }
  uwall[0] *= fReciprocal(mass);
  uwall[1] *= fReciprocal(mass);
  uwall[2] *= fReciprocal(mass);
  return 0;  
}


int fD3Q15CPS(double p, double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, double &f6, 
              double &f7, double &f8, double &f9, double &f10, double &f11,
              double &f12, double &f13, double &f14)
{

  // produce fixed concentration at planar surface: expressed for top wall

  double rho;
  double c1=1.0/9.0,c2=1.0/72.0;
  rho=6.0*(p-f0-f1-f3-f6-f7-f8-f9-f10-f11-f12)/(1.0-3.0*v1);
  
  f2=c1*rho*(1.0-3.0*v1);
  f4=c2*rho*(1.0-3.0*v0-3.0*v1-3.0*v2); 
  f5=c2*rho*(1.0-3.0*v0-3.0*v1+3.0*v2);
  f13=c2*rho*(1.0+3.0*v0-3.0*v1+3.0*v2); 
  f14=c2*rho*(1.0+3.0*v0-3.0*v1-3.0*v2); 
  
  return 0;
}


int fD3Q15CCE(double p, double *v, double * startpos)
{

  // produce fixed solute concentration at convex or concave corner:
  // expressed for left-bottom-back corner (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumC(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q15PC(long tpos, int prop, double *uwall)
{
  double oneconc;
  long spos=tpos * lbsitelength + lbsy.nf * lbsy.nq;
  for(int ktt=0; ktt<lbsy.nc; ktt++) {
    if(prop == 22) {
      oneconc=lbtopc[ktt];
      fD3Q15CPS(oneconc, uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
		        lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6],
                lbf[spos+7], lbf[spos+8], lbf[spos+9], lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+13], lbf[spos+14]);
    }
    else if(prop == 21) {
      oneconc=lbbotc[ktt];
      fD3Q15CPS(oneconc, uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+1],
		        lbf[spos+9], lbf[spos+3], lbf[spos+6], lbf[spos+7], lbf[spos+4],
                lbf[spos+5], lbf[spos+8], lbf[spos+2], lbf[spos+10], lbf[spos+13],
                lbf[spos+14], lbf[spos+11], lbf[spos+12]);
    }
    else if(prop == 23) {
      oneconc=lbrigc[ktt];
      fD3Q15CPS(oneconc, uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+14],
                lbf[spos+13], lbf[spos+9], lbf[spos+8],	lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+7], lbf[spos+6]);
    }
    else if(prop == 24) {
      oneconc=lblefc[ktt];
      fD3Q15CPS(oneconc, uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+8], lbf[spos+3], lbf[spos+13], lbf[spos+14], lbf[spos+4],
                lbf[spos+5], lbf[spos+9], lbf[spos+1], lbf[spos+10], lbf[spos+6],
                lbf[spos+7], lbf[spos+11], lbf[spos+12]);
    }
    else if(prop == 25) {
      oneconc=lbbacc[ktt];
      fD3Q15CPS(oneconc, uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+10], lbf[spos+2], lbf[spos+5], lbf[spos+7], lbf[spos+4],
                lbf[spos+6], lbf[spos+8], lbf[spos+3], lbf[spos+9], lbf[spos+12],
                lbf[spos+14], lbf[spos+11], lbf[spos+13]);
    }
    else if(prop == 26) {
      oneconc=lbfroc[ktt];
      fD3Q15CPS(oneconc, uwall[0], uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+3], lbf[spos+2], lbf[spos+4], lbf[spos+6], lbf[spos+5],
                lbf[spos+7], lbf[spos+8], lbf[spos+10],	lbf[spos+9], lbf[spos+11],
                lbf[spos+13], lbf[spos+12], lbf[spos+14]);
    }
    else if(prop == 27) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 28) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 29) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 30) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 31) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 43) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 44) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 45) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 46) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 47) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 48) {
      oneconc=lbrigc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 49) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 50) {
      oneconc=lblefc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 51) {
      oneconc=lbbotc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 52) {
      oneconc=lbrigc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 53) {
      oneconc=lbtopc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 54) {
      oneconc=lblefc[ktt];
      fD3Q15CCE(oneconc, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}


int fD3Q15TCE(double p, double *v, double * startpos)
{

  // produce fixed temperature at convex or concave corner:
  // expressed for left-bottom-back corner (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumT(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q15PT(long tpos, int prop, double *uwall)
{
  double onetemp;
  long spos=tpos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq;
  if(lbsy.nt == 1) {
    if(prop == 22) {
      onetemp=lbtopt;
      fD3Q15CPS(onetemp, uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6],
                lbf[spos+7], lbf[spos+8], lbf[spos+9], lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+13], lbf[spos+14]);
    }
    else if(prop == 21) {
      onetemp=lbbott;
      fD3Q15CPS(onetemp, uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+9], lbf[spos+3], lbf[spos+6], lbf[spos+7], lbf[spos+4],
                lbf[spos+5], lbf[spos+8], lbf[spos+2], lbf[spos+10], lbf[spos+13],
                lbf[spos+14], lbf[spos+11], lbf[spos+12]);
    }
    else if(prop == 23) {
      onetemp=lbrigt;
      fD3Q15CPS(onetemp, uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+14],
                lbf[spos+13], lbf[spos+9], lbf[spos+8],	lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+7], lbf[spos+6]);
    }
    else if(prop == 24) {
      onetemp=lbleft;
      fD3Q15CPS(onetemp, uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+8], lbf[spos+3], lbf[spos+13], lbf[spos+14], lbf[spos+4],
                lbf[spos+5], lbf[spos+9], lbf[spos+1], lbf[spos+10], lbf[spos+6],
                lbf[spos+7], lbf[spos+11], lbf[spos+12]);
    }
    else if(prop == 25) {
      onetemp=lbbact;
      fD3Q15CPS(onetemp, uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+10], lbf[spos+2], lbf[spos+5], lbf[spos+7], lbf[spos+4],
                lbf[spos+6], lbf[spos+8], lbf[spos+3], lbf[spos+9], lbf[spos+12],
                lbf[spos+14], lbf[spos+11], lbf[spos+13]);
    }
    else if(prop == 26) {
      onetemp=lbfrot;
      fD3Q15CPS(lbfrot, uwall[0], uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+3], lbf[spos+2], lbf[spos+4], lbf[spos+6], lbf[spos+5],
                lbf[spos+7], lbf[spos+8], lbf[spos+10],	lbf[spos+9], lbf[spos+11],
                lbf[spos+13], lbf[spos+12], lbf[spos+14]);
    }
    else if(prop == 27) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 28) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 29) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 30) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 31) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 43) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 44) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 45) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 46) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 47) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 48) {
      onetemp=lbrigt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 49) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 50) {
      onetemp=lbleft;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 51) {
      onetemp=lbbott;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 52) {
      onetemp=lbrigt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 53) {
      onetemp=lbtopt;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 54) {
      onetemp=lbleft;
      fD3Q15TCE(onetemp, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}


// D3Q19

int fD3Q19VPS(double v0, double v1, double v2, double &f0, double &f1, double &f2, 
              double &f3, double &f4, double &f5, double &f6, double &f7, double &f8,
              double &f9, double &f10, double &f11, double &f12, double &f13, double &f14,
              double &f15, double &f16, double &f17, double &f18)
{

  // produce fixed velocity for compressible fluid at planar surface boundary:
  // expressed for top wall
  
  double rho, n1, n2;
  double c1=1.0/3.0, c2=1.0/6.0;
  rho=(f0+f1+f3+f6+f7+f10+f12+f15+f16+2.0*(f5+f11+f13+f17+f18))/(1+v1);
  n1=0.5*(f10-f1-f6-f7+f15+f16)-c1*rho*v0;
  n2=0.5*(f12-f3-f6+f7+f15-f16)-c1*rho*v2;
  f2=f11-c1*rho*v1;
  f4=f13-c2*rho*(v0+v1)+n1;
  f8=f17-c2*rho*(v1+v2)+n2;
  f9=f18-c2*rho*(v1-v2)-n2;  
  f14=f5+c2*rho*(v0-v1)-n1;
  return 0;
}

int fD3Q19VPSIncom(double v0, double v1, double v2, double rho0, double &f0, double &f1,
	           double &f2, double &f3, double &f4, double &f5, double &f6, double &f7,
                   double &f8, double &f9, double &f10, double &f11, double &f12, double &f13,
	           double &f14, double &f15, double &f16, double &f17, double &f18)
{

  // produce fixed velocity for incompressible fluid at planar surface boundary:
  // expressed for top wall
  
  double n1, n2;
  double c1=1.0/3.0, c2=1.0/6.0;
//  double rho=f0+f1+f3+f6+f7+f10+f12+f15+f16+2.0*(f5+f11+f13+f17+f18)-rho0*v1;
  n1=0.5*(f10-f1-f6-f7+f15+f16)-c1*rho0*v0;
  n2=0.5*(f12-f3-f6+f7+f15-f16)-c1*rho0*v2;
  f2=f11-c1*rho0*v1;
  f4=f13-c2*rho0*(v0+v1)+n1;
  f8=f17-c2*rho0*(v1+v2)+n2;
  f9=f18-c2*rho0*(v1-v2)-n2;  
  f14=f5+c2*rho0*(v0-v1)-n1;
  return 0;
}

int fD3Q19VCE(double p, double *v, double * startpos)
{

  // produce fixed velocity at convex or concave corner for compressible fluid:
  // expressed for left-bottom edge (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumF(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q19VCEIncom(double p, double p0, double *v, double * startpos)
{

  // produce fixed velocity at convex or concave corner for incompressible fluid: 
  // expressed for left-bottom edge (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumFIncom(lbfeq, v, p, p0);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q19VF(long tpos, int prop, double *uwall)
{
  long spos=tpos * lbsitelength;
  double rho,rho0;
  long rpos;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 22) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VPS(uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                  lbf[spos+18]);
      }
      else if(prop == 21) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VPS(-uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+10],
                  lbf[spos+11], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                  lbf[spos+16], lbf[spos+15], lbf[spos+18], lbf[spos+17],
                  lbf[spos+1], lbf[spos+2], lbf[spos+12], lbf[spos+4],
                  lbf[spos+5], lbf[spos+7], lbf[spos+6], lbf[spos+9],
                  lbf[spos+8]);
      }
      else if(prop == 23) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q19VPS(-uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+11],
                  lbf[spos+1], lbf[spos+3], lbf[spos+5], lbf[spos+13],
                  lbf[spos+18], lbf[spos+17], lbf[spos+6], lbf[spos+7],
                  lbf[spos+2], lbf[spos+10], lbf[spos+12], lbf[spos+14],
                  lbf[spos+4], lbf[spos+9], lbf[spos+8], lbf[spos+15],
                  lbf[spos+16]);
      }
      else if(prop == 24) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q19VPS(uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                  lbf[spos+10], lbf[spos+3], lbf[spos+14], lbf[spos+4],
                  lbf[spos+8], lbf[spos+9], lbf[spos+16], lbf[spos+15],
                  lbf[spos+11], lbf[spos+1], lbf[spos+12], lbf[spos+5],
                  lbf[spos+13], lbf[spos+17], lbf[spos+18], lbf[spos+7],
                  lbf[spos+6]);
      }
      else if(prop == 25) {
        uwall[0] = lbbacv[0]; uwall[1] = lbbacv[1]; uwall[2] = lbbacv[2];
        fD3Q19VPS(uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                  lbf[spos+12], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+17],
                  lbf[spos+10], lbf[spos+3], lbf[spos+11], lbf[spos+16],
                  lbf[spos+15], lbf[spos+13], lbf[spos+14], lbf[spos+18],
                  lbf[spos+8]);
      }
      else if(prop == 26) {
        uwall[0] = lbfrov[0]; uwall[1] = lbfrov[1]; uwall[2] = lbfrov[2];
        fD3Q19VPS(uwall[0], uwall[2], -uwall[1], lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+11], lbf[spos+6], lbf[spos+7],
                  lbf[spos+5], lbf[spos+4], lbf[spos+18], lbf[spos+8],
                  lbf[spos+10], lbf[spos+12], lbf[spos+2], lbf[spos+15],
                  lbf[spos+16], lbf[spos+14], lbf[spos+13], lbf[spos+9],
                  lbf[spos+17]);
      }
      else if(prop == 27) {
        rpos = fNextStep(1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 28) {
        rpos = fNextStep(-1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 29) {
        rpos = fNextStep(-1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 30) {
        rpos = fNextStep(1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 31) {
        rpos = fNextStep(1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        rpos = fNextStep(-1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        rpos = fNextStep(-1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        rpos = fNextStep(1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 43) {
        rpos = fNextStep(1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 44) {
        rpos = fNextStep(-1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 45) {
        rpos = fNextStep(-1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 46) {
        rpos = fNextStep(1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 47) {
        rpos = fNextStep(0, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 48) {
        rpos = fNextStep(-1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 49) {
        rpos = fNextStep(0, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 50) {
        rpos = fNextStep(1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 51) {
        rpos = fNextStep(0, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 52) {
        rpos = fNextStep(-1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 53) {
        rpos = fNextStep(0, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 54) {
        rpos = fNextStep(1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q19VCE(rho, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 22) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VPSIncom(uwall[0], uwall[1], uwall[2], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                       lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                       lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                       lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                       lbf[spos+18]);
      }
      else if(prop == 21) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VPSIncom(-uwall[0], -uwall[1], uwall[2], rho0, lbf[spos], lbf[spos+10],
                       lbf[spos+11], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                       lbf[spos+16], lbf[spos+15], lbf[spos+18], lbf[spos+17],
                       lbf[spos+1], lbf[spos+2], lbf[spos+12], lbf[spos+4],
                       lbf[spos+5], lbf[spos+7], lbf[spos+6], lbf[spos+9],
                       lbf[spos+8]);
      }
      else if(prop == 23) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q19VPSIncom(-uwall[1], uwall[0], uwall[2], rho0, lbf[spos], lbf[spos+11],
                       lbf[spos+1], lbf[spos+3], lbf[spos+5], lbf[spos+13],
                       lbf[spos+18], lbf[spos+17], lbf[spos+6], lbf[spos+7],
                       lbf[spos+2], lbf[spos+10], lbf[spos+12], lbf[spos+14],
                       lbf[spos+4], lbf[spos+9], lbf[spos+8], lbf[spos+15],
                       lbf[spos+16]);
      }
      else if(prop == 24) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q19VPSIncom(uwall[1], -uwall[0], uwall[2], rho0, lbf[spos], lbf[spos+2],
                       lbf[spos+10], lbf[spos+3], lbf[spos+14], lbf[spos+4],
                       lbf[spos+8], lbf[spos+9], lbf[spos+16], lbf[spos+15],
                       lbf[spos+11], lbf[spos+1], lbf[spos+12], lbf[spos+5],
                       lbf[spos+13], lbf[spos+17], lbf[spos+18], lbf[spos+7],
                       lbf[spos+6]);
      }
      else if(prop == 25) {
        uwall[0] = lbbacv[0]; uwall[1] = lbbacv[1]; uwall[2] = lbbacv[2];
        fD3Q19VPSIncom(uwall[0], -uwall[2], uwall[1], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+12], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                       lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+17],
                       lbf[spos+10], lbf[spos+3], lbf[spos+11], lbf[spos+16],
                       lbf[spos+15], lbf[spos+13], lbf[spos+14], lbf[spos+18],
                       lbf[spos+8]);
      }
      else if(prop == 26) {
        uwall[0] = lbfrov[0]; uwall[1] = lbfrov[1]; uwall[2] = lbfrov[2];
        fD3Q19VPSIncom(uwall[0], uwall[2], -uwall[1], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+3], lbf[spos+11], lbf[spos+6], lbf[spos+7],
                       lbf[spos+5], lbf[spos+4], lbf[spos+18], lbf[spos+8],
                       lbf[spos+10], lbf[spos+12], lbf[spos+2], lbf[spos+15],
                       lbf[spos+16], lbf[spos+14], lbf[spos+13], lbf[spos+9],
                       lbf[spos+17]);
      }
      else if(prop == 27) {
        rpos = fNextStep(1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 28) {
        rpos = fNextStep(-1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 29) {
        rpos = fNextStep(-1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 30) {
        rpos = fNextStep(1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 31) {
        rpos = fNextStep(1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        rpos = fNextStep(-1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        rpos = fNextStep(-1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        rpos = fNextStep(1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 43) {
        rpos = fNextStep(1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 44) {
        rpos = fNextStep(-1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 45) {
        rpos = fNextStep(-1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 46) {
        rpos = fNextStep(1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 47) {
        rpos = fNextStep(0, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 48) {
        rpos = fNextStep(-1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 49) {
        rpos = fNextStep(0, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 50) {
        rpos = fNextStep(1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 51) {
        rpos = fNextStep(0, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 52) {
        rpos = fNextStep(-1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 53) {
        rpos = fNextStep(0, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 54) {
        rpos = fNextStep(1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q19VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  
  return 0;
  
}


int fD3Q19PPS(double p, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, 
	      double &f6, double &f7, double &f8, double &f9,
	      double &f10, double &f11, double &f12, double &f13,
	      double &f14, double &f15, double &f16, double &f17, 
	      double &f18, double& vel)
{

  // produce fixed pressure/density at planar surface: expressed for top wall

  double rx, rz;
  double c1=1.0/3.0,c2=1.0/6.0;
  rx=0.5*(f10-f1-f6-f7+f15+f16);
  vel=(f0+f1+f3+f6+f7+f10+f12+f15+f16+2.0*(f5+f11+f13+f17+f18)) -p;
  rz=0.5*(f12-f3-f6+f7+f15-f16);
  f2=f11-c1*vel;
  f4=f13-c2*vel+rx;
  f8=f17-c2*vel+rz; 
  f9=f18-c2*vel-rz; 
  f14=f5-c2*vel-rx; 
  
  return 0;
}


int fD3Q19PF(long tpos, int prop, double *uwall)
{
  double moment,rho0,onemass,mass=0.0;
  long spos=tpos * lbsitelength;
  uwall[0]=0.0; uwall[1]=0.0; uwall[2]=0.0;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 22) {
        onemass=lbtopp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                  lbf[spos+18], moment);
        mass += onemass;
        uwall[1] += moment;
      }
      else if(prop == 21) {
        onemass=lbbotp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+10],
                  lbf[spos+11], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                  lbf[spos+16], lbf[spos+15], lbf[spos+18], lbf[spos+17],
                  lbf[spos+1], lbf[spos+2], lbf[spos+12], lbf[spos+4],
                  lbf[spos+5], lbf[spos+7], lbf[spos+6], lbf[spos+9],
                  lbf[spos+8], moment);
        mass += onemass;
        uwall[1] -= moment;
      }  
      else if(prop == 23) {
        onemass=lbrigp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+11],
                  lbf[spos+1], lbf[spos+3], lbf[spos+5], lbf[spos+13],
                  lbf[spos+18], lbf[spos+17], lbf[spos+6], lbf[spos+7],
                  lbf[spos+2], lbf[spos+10], lbf[spos+12], lbf[spos+14],
                  lbf[spos+4], lbf[spos+9], lbf[spos+8], lbf[spos+15],
                  lbf[spos+16], moment);
        mass += onemass;
        uwall[0] += moment;
      }  
      else if(prop == 24) {
        onemass=lblefp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+10], lbf[spos+3], lbf[spos+14], lbf[spos+4],
                  lbf[spos+8], lbf[spos+9], lbf[spos+16], lbf[spos+15],
                  lbf[spos+11], lbf[spos+1], lbf[spos+12], lbf[spos+5],
                  lbf[spos+13], lbf[spos+17], lbf[spos+18], lbf[spos+7],
                  lbf[spos+6], moment);
        mass += onemass;
        uwall[0] -= moment;
      }  
      else if(prop == 25) {
        onemass=lbbacp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+12], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+17],
                  lbf[spos+10], lbf[spos+3], lbf[spos+11], lbf[spos+16],
                  lbf[spos+15], lbf[spos+13], lbf[spos+14], lbf[spos+18],
                  lbf[spos+8], moment);
        mass += onemass;
        uwall[2] -= moment;
      }   
      else if(prop == 26) {
        onemass=lbfrop[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+11], lbf[spos+6], lbf[spos+7],
                  lbf[spos+5], lbf[spos+4], lbf[spos+18], lbf[spos+8],
                  lbf[spos+10], lbf[spos+12], lbf[spos+2], lbf[spos+15],
                  lbf[spos+16], lbf[spos+14], lbf[spos+13], lbf[spos+9],
                  lbf[spos+17], moment);
        mass += onemass;
        uwall[2] += moment;
      }
      else if(prop == 27) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 28) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 29) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 30) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 43) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 44) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 45) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 46) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 49) {
        onemass=lbtopp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 51) {
        onemass=lbbotp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 52) {
        onemass=lbrigp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 53) {
        onemass=lbtopp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 54) {
        onemass=lblefp[ktt];
        fD3Q19VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 22) {
        onemass=lbtopp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                  lbf[spos+18], moment);
        mass += rho0;
        uwall[1] += moment;
      }
      else if(prop == 21) {
        onemass=lbbotp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+10],
                  lbf[spos+11], lbf[spos+3], lbf[spos+13], lbf[spos+14],
                  lbf[spos+16], lbf[spos+15], lbf[spos+18], lbf[spos+17],
                  lbf[spos+1], lbf[spos+2], lbf[spos+12], lbf[spos+4],
                  lbf[spos+5], lbf[spos+7], lbf[spos+6], lbf[spos+9],
                  lbf[spos+8], moment);
        mass += rho0;
        uwall[1] -= moment;
      }  
      else if(prop == 23) {
        onemass=lbrigp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+11],
                  lbf[spos+1], lbf[spos+3], lbf[spos+5], lbf[spos+13],
                  lbf[spos+18], lbf[spos+17], lbf[spos+6], lbf[spos+7],
                  lbf[spos+2], lbf[spos+10], lbf[spos+12], lbf[spos+14],
                  lbf[spos+4], lbf[spos+9], lbf[spos+8], lbf[spos+15],
                  lbf[spos+16], moment);
        mass += rho0;
        uwall[0] += moment;
      }  
      else if(prop == 24) {
        onemass=lblefp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+10], lbf[spos+3], lbf[spos+14], lbf[spos+4],
                  lbf[spos+8], lbf[spos+9], lbf[spos+16], lbf[spos+15],
                  lbf[spos+11], lbf[spos+1], lbf[spos+12], lbf[spos+5],
                  lbf[spos+13], lbf[spos+17], lbf[spos+18], lbf[spos+7],
                  lbf[spos+6], moment);
        mass += rho0;
        uwall[0] -= moment;
      }  
      else if(prop == 25) {
        onemass=lbbacp[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+12], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+17],
                  lbf[spos+10], lbf[spos+3], lbf[spos+11], lbf[spos+16],
                  lbf[spos+15], lbf[spos+13], lbf[spos+14], lbf[spos+18],
                  lbf[spos+8], moment);
        mass += rho0;
        uwall[2] -= moment;
      }   
      else if(prop == 26) {
        onemass=lbfrop[ktt];
        fD3Q19PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+11], lbf[spos+6], lbf[spos+7],
                  lbf[spos+5], lbf[spos+4], lbf[spos+18], lbf[spos+8],
                  lbf[spos+10], lbf[spos+12], lbf[spos+2], lbf[spos+15],
                  lbf[spos+16], lbf[spos+14], lbf[spos+13], lbf[spos+9],
                  lbf[spos+17], moment);
        mass += rho0;
        uwall[2] += moment;
      }

      else if(prop == 27) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 28) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 29) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 30) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 43) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 44) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 45) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 46) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 49) {
        onemass=lbtopp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 51) {
        onemass=lbbotp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 52) {
        onemass=lbrigp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 53) {
        onemass=lbtopp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 54) {
        onemass=lblefp[ktt];
        fD3Q19VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      spos += lbsy.nq;
    }
  }
  uwall[0] *= fReciprocal(mass);
  uwall[1] *= fReciprocal(mass);
  uwall[2] *= fReciprocal(mass);
  return 0;  
}


int fD3Q19CPS(double p, double v0, double v1, double v2, double &f0, double &f1,
	      double &f2, double &f3, double &f4, double &f5, double &f6,
              double &f7, double &f8, double &f9, double &f10, double &f11,
              double &f12, double &f13, double &f14, double &f15, double &f16,
              double &f17, double &f18)
{

  // produce fixed concentration at planar surface: expressed for top wall

  double rho;
  double c1=1.0/18.0,c2=1.0/36.0;
  rho=6.0*(p-f0-f1-f3-f5-f6-f7-f10-f11-f12-f13-f15-f16-f17-f18)/(1.0-3.0*v1);
  f2=c1*rho*(1.0-3.0*v1);
  f4=c2*rho*(1.0-3.0*v0-3.0*v1);
  f8=c2*rho*(1.0-3.0*v1-3.0*v2);
  f9=c2*rho*(1.0-3.0*v1+3.0*v2);
  f14=c2*rho*(1.0+3.0*v0-3.0*v1);
  
  return 0;
}


int fD3Q19CCE(double p, double *v, double * startpos)
{

  // produce fixed solute concentration at convex or concave corner:
  // expressed for left-bottom edge (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumC(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q19PC(long tpos, int prop, double *uwall)
{
  double oneconc;
  long spos=tpos * lbsitelength + lbsy.nf * lbsy.nq;
  for(int ktt=0; ktt<lbsy.nc; ktt++) {
    if(prop == 22) {
      oneconc=lbtopc[ktt];
      fD3Q19CPS(oneconc, uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6],
                lbf[spos+7], lbf[spos+8], lbf[spos+9], lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+13], lbf[spos+14], lbf[spos+15], lbf[spos+16],
                lbf[spos+17], lbf[spos+18]);
    }
    else if(prop == 21) {
      oneconc=lbbotc[ktt];
      fD3Q19CPS(oneconc, -uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+10],
                lbf[spos+11], lbf[spos+3], lbf[spos+13], lbf[spos+14], lbf[spos+16],
                lbf[spos+15], lbf[spos+18], lbf[spos+17], lbf[spos+1], lbf[spos+2],
                lbf[spos+12], lbf[spos+4], lbf[spos+5], lbf[spos+7], lbf[spos+6],
                lbf[spos+9], lbf[spos+8]);
    }
    else if(prop == 23) {
      oneconc=lbrigc[ktt];
      fD3Q19CPS(oneconc, -uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+11],
                lbf[spos+1], lbf[spos+3], lbf[spos+5], lbf[spos+13], lbf[spos+18],
                lbf[spos+17], lbf[spos+6], lbf[spos+7],	lbf[spos+2], lbf[spos+10],
                lbf[spos+12], lbf[spos+14], lbf[spos+4], lbf[spos+9], lbf[spos+8],
                lbf[spos+15], lbf[spos+16]);
    }
    else if(prop == 24) {
      oneconc=lblefc[ktt];
      fD3Q19CPS(oneconc, uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+10], lbf[spos+3], lbf[spos+14], lbf[spos+4], lbf[spos+8],
                lbf[spos+9], lbf[spos+16], lbf[spos+15], lbf[spos+11], lbf[spos+1],
                lbf[spos+12], lbf[spos+5], lbf[spos+13], lbf[spos+17], lbf[spos+18],
                lbf[spos+7], lbf[spos+6]);
    }
    else if(prop == 25) {
      oneconc=lbbacc[ktt];
      fD3Q19CPS(oneconc, uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+12], lbf[spos+2], lbf[spos+7], lbf[spos+6], lbf[spos+4],
                lbf[spos+5], lbf[spos+9], lbf[spos+17], lbf[spos+10], lbf[spos+3],
                lbf[spos+11], lbf[spos+16], lbf[spos+15], lbf[spos+13], lbf[spos+14],
                lbf[spos+18], lbf[spos+8]);
    }
    else if(prop == 26) {
      oneconc=lbfroc[ktt];
      fD3Q19CPS(lbfroc[ktt], uwall[0], uwall[2], -uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+3], lbf[spos+11], lbf[spos+6], lbf[spos+7], lbf[spos+5],
                lbf[spos+4], lbf[spos+18], lbf[spos+8],	lbf[spos+10], lbf[spos+12],
                lbf[spos+2], lbf[spos+15], lbf[spos+16], lbf[spos+14], lbf[spos+13],
                lbf[spos+9], lbf[spos+17]);   
    }
    else if(prop == 27) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 28) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 29) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 30) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 31) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 43) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 44) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 45) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 46) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 47) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 48) {
      oneconc=lbrigc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 49) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 50) {
      oneconc=lblefc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 51) {
      oneconc=lbbotc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 52) {
      oneconc=lbrigc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 53) {
      oneconc=lbtopc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 54) {
      oneconc=lblefc[ktt];
      fD3Q19CCE(oneconc, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}

int fD3Q19TCE(double p, double *v, double * startpos)
{

  // produce fixed temperature at convex or concave corner:
  // expressed for left-bottom edge (against velocity 7)

  double *pt1=startpos;
  fGetEquilibriumT(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q19PT(long tpos, int prop, double *uwall)
{
  double onetemp;
  long spos=tpos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq;
  if(lbsy.nt == 1) {
    if(prop == 22) {
      onetemp=lbtopt;
      fD3Q19CPS(onetemp, uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6],
                lbf[spos+7], lbf[spos+8], lbf[spos+9], lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+13], lbf[spos+14], lbf[spos+15], lbf[spos+16],
                lbf[spos+17], lbf[spos+18]);
    }
    else if(prop == 21) {
      onetemp=lbbott;
      fD3Q19CPS(onetemp, -uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+10],
                lbf[spos+11], lbf[spos+3], lbf[spos+13], lbf[spos+14], lbf[spos+16],
                lbf[spos+15], lbf[spos+18], lbf[spos+17], lbf[spos+1], lbf[spos+2],
                lbf[spos+12], lbf[spos+4], lbf[spos+5], lbf[spos+7], lbf[spos+6],
                lbf[spos+9], lbf[spos+8]);
    }
    else if(prop == 23) {
      onetemp=lbrigt;
      fD3Q19CPS(onetemp, -uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+11],
                lbf[spos+1], lbf[spos+3], lbf[spos+5], lbf[spos+13], lbf[spos+18],
                lbf[spos+17], lbf[spos+6], lbf[spos+7],	lbf[spos+2], lbf[spos+10],
                lbf[spos+12], lbf[spos+14], lbf[spos+4], lbf[spos+9], lbf[spos+8],
                lbf[spos+15], lbf[spos+16]);
    }
    else if(prop == 24) {
      onetemp=lbleft;
      fD3Q19CPS(onetemp, uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+10], lbf[spos+3], lbf[spos+14], lbf[spos+4], lbf[spos+8],
                lbf[spos+9], lbf[spos+16], lbf[spos+15], lbf[spos+11], lbf[spos+1],
                lbf[spos+12], lbf[spos+5], lbf[spos+13], lbf[spos+17], lbf[spos+18],
                lbf[spos+7], lbf[spos+6]);
    }
    else if(prop == 25) {
      onetemp=lbbact;
      fD3Q19CPS(onetemp, uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+12], lbf[spos+2], lbf[spos+7], lbf[spos+6], lbf[spos+4],
                lbf[spos+5], lbf[spos+9], lbf[spos+17], lbf[spos+10], lbf[spos+3],
                lbf[spos+11], lbf[spos+16], lbf[spos+15], lbf[spos+13], lbf[spos+14],
                lbf[spos+18], lbf[spos+8]);
    }
    else if(prop == 26) {
      onetemp=lbfrot;
      fD3Q19CPS(onetemp, uwall[0], uwall[2], -uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+3], lbf[spos+11], lbf[spos+6], lbf[spos+7], lbf[spos+5],
                lbf[spos+4], lbf[spos+18], lbf[spos+8],	lbf[spos+10], lbf[spos+12],
                lbf[spos+2], lbf[spos+15], lbf[spos+16], lbf[spos+14], lbf[spos+13],
                lbf[spos+9], lbf[spos+17]);   
    }
    else if(prop == 27) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 28) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 29) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 30) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 31) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 43) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 44) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 45) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 46) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 47) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 48) {
      onetemp=lbrigt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 49) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 50) {
      onetemp=lbleft;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 51) {
      onetemp=lbbott;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 52) {
      onetemp=lbrigt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 53) {
      onetemp=lbtopt;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 54) {
      onetemp=lbleft;
      fD3Q19TCE(onetemp, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}


// D3Q27

int fD3Q27VPS(double v0, double v1, double v2, double &f0, double &f1, double &f2,
              double &f3, double &f4, double &f5, double &f6, double &f7, double &f8,
              double &f9, double &f10, double &f11, double &f12, double &f13, double &f14,
              double &f15, double &f16, double &f17, double &f18, double &f19, double &f20,
              double &f21, double &f22, double &f23, double &f24, double &f25, double &f26)
{

  // produce fixed velocity for compressible fluid at planar surface boundary:
  // expressed for top wall
  
  double rho, n1, n2;
  double c1=4.0/9.0,c2=1.0/9.0,c3=1.0/36.0,c4=1.0/6.0;
  rho=(f0+f1+f3+f6+f7+f14+f16+f19+f20+2.0*(f5+f12+f13+f15+f17+f21+f22+f23+f24))/(1+v1);
  n1=c4*(f14-f1-f6-f7+f19+f20)-c2*rho*v0;
  n2=c4*(f16-f3-f6+f7+f19-f20)-c2*rho*v2;
  f2=f15-c1*rho*v1;
  f4=f17-c2*rho*(v0+v1)+n1;
  f8=f21-c2*rho*(v1+v2)+n2;
  f9=f22-c2*rho*(v1-v2)-n2;
  f18=f5+c2*rho*(v0-v1)-n1;
  f10=f23-c3*rho*(v0+v1+v2)+n1+n2;
  f11=f24-c3*rho*(v0+v1-v2)+n1-n2;
  f25=f12+c3*rho*(v0-v1+v2)-n1-n2;
  f26=f13+c3*rho*(v0-v1-v2)-n1+n2;
  return 0;
}

int fD3Q27VPSIncom(double v0, double v1, double v2, double rho0, double &f0, double &f1,
                   double &f2, double &f3, double &f4, double &f5, double &f6, double &f7,
                   double &f8, double &f9, double &f10, double &f11, double &f12, double &f13,
                   double &f14, double &f15, double &f16, double &f17, double &f18, double &f19,
                   double &f20, double &f21, double &f22, double &f23, double &f24, double &f25,
                   double &f26)
{

  // produce fixed velocity for incompressible fluid at planar surface boundary:
  // expressed for top wall
  
  double n1, n2;
  double c1=4.0/9.0,c2=1.0/9.0,c3=1.0/36.0,c4=1.0/6.0;
//  double rho=f0+f1+f3+f6+f7+f14+f16+f19+f20+2.0*(f5+f12+f13+f15+f17+f21+f22+f23+f24)-rho0*v1;
  n1=c4*(f14-f1-f6-f7+f19+f20)-c2*rho0*v0;
  n2=c4*(f16-f3-f6+f7+f19-f20)-c2*rho0*v2;
  f2=f15-c1*rho0*v1;
  f4=f17-c2*rho0*(v0+v1)+n1;
  f8=f21-c2*rho0*(v1+v2)+n2;
  f9=f22-c2*rho0*(v1-v2)-n2;
  f18=f5+c2*rho0*(v0-v1)-n1;
  f10=f23-c3*rho0*(v0+v1+v2)+n1+n2;
  f11=f24-c3*rho0*(v0+v1-v2)+n1-n2;
  f25=f12+c3*rho0*(v0-v1+v2)-n1-n2;
  f26=f13+c3*rho0*(v0-v1-v2)-n1+n2;
  return 0;
}

int fD3Q27VCE(double p, double *v, double * startpos)
{

  // produce fixed velocity for compressible fluid at convex or concave corner:
  // expressed for left-bottom-back corner (against velocity 19) 

  double *pt1=startpos;
  fGetEquilibriumF(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q27VCEIncom(double p, double p0, double *v, double * startpos)
{

  // produce fixed velocity for incompressible fluid at convex or concave corner:
  // expressed for left-bottom-back corner (against velocity 19) 

  double *pt1=startpos;
  fGetEquilibriumFIncom(lbfeq, v, p, p0);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q27VF(long tpos, int prop, double *uwall)
{
  long spos=tpos * lbsitelength;
  double rho,rho0;
  long rpos;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 22) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VPS(uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                  lbf[spos+18], lbf[spos+19], lbf[spos+20], lbf[spos+21],
                  lbf[spos+22], lbf[spos+23], lbf[spos+24], lbf[spos+25],
                  lbf[spos+26]);
      }
      else if(prop == 21) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VPS(uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                  lbf[spos+15], lbf[spos+3], lbf[spos+5], lbf[spos+4],
                  lbf[spos+6], lbf[spos+7], lbf[spos+22], lbf[spos+21],
                  lbf[spos+12], lbf[spos+13], lbf[spos+10], lbf[spos+11],
                  lbf[spos+14], lbf[spos+2], lbf[spos+16], lbf[spos+18],
                  lbf[spos+17], lbf[spos+19], lbf[spos+20], lbf[spos+9],
                  lbf[spos+8], lbf[spos+25], lbf[spos+26], lbf[spos+23],
                  lbf[spos+24]);
      }
      else if(prop == 23) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q27VPS(uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                  lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+18],
                  lbf[spos+8], lbf[spos+9], lbf[spos+6], lbf[spos+7],
                  lbf[spos+10], lbf[spos+11], lbf[spos+26], lbf[spos+25],
                  lbf[spos+15], lbf[spos+14], lbf[spos+16], lbf[spos+17],
                  lbf[spos+5], lbf[spos+21], lbf[spos+22], lbf[spos+19],
                  lbf[spos+20], lbf[spos+23], lbf[spos+24], lbf[spos+13],
                  lbf[spos+12]);
      }
      else if(prop == 24) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q27VPS(uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                  lbf[spos+14], lbf[spos+3], lbf[spos+18], lbf[spos+4],
                  lbf[spos+8], lbf[spos+9], lbf[spos+20], lbf[spos+19],
                  lbf[spos+26], lbf[spos+25], lbf[spos+10], lbf[spos+11],
                  lbf[spos+15], lbf[spos+1], lbf[spos+16], lbf[spos+5],
                  lbf[spos+17], lbf[spos+21], lbf[spos+22], lbf[spos+7],
                  lbf[spos+6], lbf[spos+13], lbf[spos+12], lbf[spos+23],
                  lbf[spos+24]);
      }
      else if(prop == 25) {
        uwall[0] = lbbacv[0]; uwall[1] = lbbacv[1]; uwall[2] = lbbacv[2];
        fD3Q27VPS(uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                  lbf[spos+16], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+21],
                  lbf[spos+11], lbf[spos+13], lbf[spos+10], lbf[spos+12],
                  lbf[spos+14], lbf[spos+3], lbf[spos+15], lbf[spos+20],
                  lbf[spos+19], lbf[spos+17], lbf[spos+18], lbf[spos+22],
                  lbf[spos+8], lbf[spos+24], lbf[spos+26], lbf[spos+23],
                  lbf[spos+25]);
      }
      else if(prop == 26) {
        uwall[0] = lbfrov[0]; uwall[1] = lbfrov[1]; uwall[2] = lbfrov[2];
        fD3Q27VPS(uwall[0], uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+2], lbf[spos+6], lbf[spos+7],
                  lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+22],
                  lbf[spos+10], lbf[spos+12], lbf[spos+11], lbf[spos+13],
                  lbf[spos+14], lbf[spos+16], lbf[spos+15], lbf[spos+19],
                  lbf[spos+20], lbf[spos+17], lbf[spos+18], lbf[spos+21],
                  lbf[spos+9], lbf[spos+23], lbf[spos+25], lbf[spos+24],
                  lbf[spos+26]);
      }
      else if(prop == 27)  {
        rpos = fNextStep(1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 28) {
        rpos = fNextStep(-1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 29) {
        rpos = fNextStep(-1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 30) {
        rpos = fNextStep(1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 31)  {
        rpos = fNextStep(1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        rpos = fNextStep(-1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        rpos = fNextStep(-1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        rpos = fNextStep(1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 43) {
        rpos = fNextStep(1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 44) {
        rpos = fNextStep(-1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 45) {
        rpos = fNextStep(-1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 46) {
        rpos = fNextStep(1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 47) {
        rpos = fNextStep(0, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 48) {
        rpos = fNextStep(-1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 49) {
        rpos = fNextStep(0, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 50) {
        rpos = fNextStep(1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 51) {
        rpos = fNextStep(0, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 52) {
        rpos = fNextStep(-1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 53) {
        rpos = fNextStep(0, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      else if(prop == 54) {
        rpos = fNextStep(1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q27VCE(rho, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 22) {
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VPSIncom(uwall[0], uwall[1], uwall[2], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                       lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                       lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                       lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                       lbf[spos+18], lbf[spos+19], lbf[spos+20], lbf[spos+21],
                       lbf[spos+22], lbf[spos+23], lbf[spos+24], lbf[spos+25],
                       lbf[spos+26]);
      }
      else if(prop == 21) {
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VPSIncom(uwall[0], -uwall[1], uwall[2], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+15], lbf[spos+3], lbf[spos+5], lbf[spos+4],
                       lbf[spos+6], lbf[spos+7], lbf[spos+22], lbf[spos+21],
                       lbf[spos+12], lbf[spos+13], lbf[spos+10], lbf[spos+11],
                       lbf[spos+14], lbf[spos+2], lbf[spos+16], lbf[spos+18],
                       lbf[spos+17], lbf[spos+19], lbf[spos+20], lbf[spos+9],
                       lbf[spos+8], lbf[spos+25], lbf[spos+26], lbf[spos+23],
                       lbf[spos+24]);
      }
      else if(prop == 23) {
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q27VPSIncom(uwall[1], uwall[0], uwall[2], rho0, lbf[spos], lbf[spos+2],
                       lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+18],
                       lbf[spos+8], lbf[spos+9], lbf[spos+6], lbf[spos+7],
                       lbf[spos+10], lbf[spos+11], lbf[spos+26], lbf[spos+25],
                       lbf[spos+15], lbf[spos+14], lbf[spos+16], lbf[spos+17],
                       lbf[spos+5], lbf[spos+21], lbf[spos+22], lbf[spos+19],
                       lbf[spos+20], lbf[spos+23], lbf[spos+24], lbf[spos+13],
                       lbf[spos+12]);
      }
      else if(prop == 24) {
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q27VPSIncom(uwall[1], -uwall[0], uwall[2], rho0, lbf[spos], lbf[spos+2],
                       lbf[spos+14], lbf[spos+3], lbf[spos+18], lbf[spos+4],
                       lbf[spos+8], lbf[spos+9], lbf[spos+20], lbf[spos+19],
                       lbf[spos+26], lbf[spos+25], lbf[spos+10], lbf[spos+11],
                       lbf[spos+15], lbf[spos+1], lbf[spos+16], lbf[spos+5],
                       lbf[spos+17], lbf[spos+21], lbf[spos+22], lbf[spos+7],
                       lbf[spos+6], lbf[spos+13], lbf[spos+12], lbf[spos+23],
                       lbf[spos+24]);
      }
      else if(prop == 25) {
        uwall[0] = lbbacv[0]; uwall[1] = lbbacv[1]; uwall[2] = lbbacv[2];
        fD3Q27VPSIncom(uwall[0], -uwall[2], uwall[1], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+16], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                       lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+21],
                       lbf[spos+11], lbf[spos+13], lbf[spos+10], lbf[spos+12],
                       lbf[spos+14], lbf[spos+3], lbf[spos+15], lbf[spos+20],
                       lbf[spos+19], lbf[spos+17], lbf[spos+18], lbf[spos+22],
                       lbf[spos+8], lbf[spos+24], lbf[spos+26], lbf[spos+23],
                       lbf[spos+25]);
      }
      else if(prop == 26) {
        uwall[0] = lbfrov[0]; uwall[1] = lbfrov[1]; uwall[2] = lbfrov[2];
        fD3Q27VPSIncom(uwall[0], uwall[2], uwall[1], rho0, lbf[spos], lbf[spos+1],
                       lbf[spos+3], lbf[spos+2], lbf[spos+6], lbf[spos+7],
                       lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+22],
                       lbf[spos+10], lbf[spos+12], lbf[spos+11], lbf[spos+13],
                       lbf[spos+14], lbf[spos+16], lbf[spos+15], lbf[spos+19],
                       lbf[spos+20], lbf[spos+17], lbf[spos+18], lbf[spos+21],
                       lbf[spos+9], lbf[spos+23], lbf[spos+25], lbf[spos+24],
                       lbf[spos+26]);
      }
      else if(prop == 27) {
        rpos = fNextStep(1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 28) {
        rpos = fNextStep(-1, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 29) {
        rpos = fNextStep(-1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 30) {
        rpos = fNextStep(1, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 31) {
        rpos = fNextStep(1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 32) {
        rpos = fNextStep(-1, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 33) {
        rpos = fNextStep(-1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 34) {
        rpos = fNextStep(1, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 43) {
        rpos = fNextStep(1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 44) {
        rpos = fNextStep(-1, 1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 45) {
        rpos = fNextStep(-1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 46) {
        rpos = fNextStep(1, -1, 0, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 47) {
        rpos = fNextStep(0, 1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 48) {
        rpos = fNextStep(-1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 49) {
        rpos = fNextStep(0, -1, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 50) {
        rpos = fNextStep(1, 0, 1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 51) {
        rpos = fNextStep(0, 1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbbotv[0]; uwall[1] = lbbotv[1]; uwall[2] = lbbotv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 52) {
        rpos = fNextStep(-1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbrigv[0]; uwall[1] = lbrigv[1]; uwall[2] = lbrigv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 53) {
        rpos = fNextStep(0, -1, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lbtopv[0]; uwall[1] = lbtopv[1]; uwall[2] = lbtopv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      else if(prop == 54) {
        rpos = fNextStep(1, 0, -1, tpos) * lbsitelength + ktt * lbsy.nq;
        rho = fEvapLimit(fGetOneMassSite(&lbf[rpos]));       
        uwall[0] = lblefv[0]; uwall[1] = lblefv[1]; uwall[2] = lblefv[2];
        fD3Q27VCEIncom(rho, rho0, uwall, &lbf[spos]);
      }
      spos += lbsy.nq;
    }
  }

  return 0;
  
}


int fD3Q27PPS(double p, double &f0, double &f1,
              double &f2, double &f3, double &f4, double &f5,
              double &f6, double &f7, double &f8, double &f9,
              double &f10, double &f11, double &f12, double &f13,
              double &f14, double &f15, double &f16, double &f17,
              double &f18, double &f19, double &f20, double &f21,
              double &f22, double &f23, double &f24, double &f25,
              double &f26, double& vel)
{

  // produce fixed pressure/density at planar surface: expressed for top wall

  double rx, rz;
  double c1=4.0/9.0,c2=1.0/9.0,c3=1.0/36.0,c4=1.0/6.0;
  rx=c4*(f14-f1-f6-f7+f19+f20);
  vel=(f0+f1+f3+f6+f7+f14+f16+f19+f20+2.0*(f5+f12+f13+f15+f17+f21+f22+f23+f24)) -p;
  rz=c4*(f16-f3-f6+f7+f19-f20);
  f2=f15-c1*vel;
  f4=f17-c2*vel+rx;
  f8=f21-c2*vel+rz;
  f9=f22-c2*vel-rz;
  f18=f5-c2*vel-rx;
  f10=f23-c3*vel+rx+rz;
  f11=f24-c3*vel+rx-rz;
  f25=f12-c3*vel-rx-rz;
  f26=f13-c3*vel-rx+rz;
  return 0;
}


int fD3Q27PF(long tpos, int prop, double *uwall)
{
  double moment,rho0,onemass,mass=0.0;
  long spos=tpos * lbsitelength;
  uwall[0]=0.0; uwall[1]=0.0; uwall[2]=0.0;
  if(!incompress) {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      if(prop == 22) {
        onemass=lbtopp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                  lbf[spos+18], lbf[spos+19], lbf[spos+20], lbf[spos+21],
                  lbf[spos+22], lbf[spos+23], lbf[spos+24], lbf[spos+25],
                  lbf[spos+26], moment);
        mass += onemass;
        uwall[1] += moment;
      }
      else if(prop == 21) {
        onemass=lbbotp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+15], lbf[spos+3], lbf[spos+5], lbf[spos+4],
                  lbf[spos+6], lbf[spos+7], lbf[spos+22], lbf[spos+21],
                  lbf[spos+12], lbf[spos+13], lbf[spos+10], lbf[spos+11],
                  lbf[spos+14], lbf[spos+2], lbf[spos+16], lbf[spos+18],
                  lbf[spos+17], lbf[spos+19], lbf[spos+20], lbf[spos+9],
                  lbf[spos+8], lbf[spos+25], lbf[spos+26], lbf[spos+23],
                  lbf[spos+24], moment);
        mass += onemass;
        uwall[1] -= moment;
      }
      else if(prop == 23) {
        onemass=lbrigp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+18],
                  lbf[spos+8], lbf[spos+9], lbf[spos+6], lbf[spos+7],
                  lbf[spos+10], lbf[spos+11], lbf[spos+26], lbf[spos+25],
                  lbf[spos+15], lbf[spos+14], lbf[spos+16], lbf[spos+17],
                  lbf[spos+5], lbf[spos+21], lbf[spos+22], lbf[spos+19],
                  lbf[spos+20], lbf[spos+23], lbf[spos+24], lbf[spos+13],
                  lbf[spos+12], moment);
        mass += onemass;
        uwall[0] += moment;
      }
      else if(prop == 24) {
        onemass=lblefp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+14], lbf[spos+3], lbf[spos+18], lbf[spos+4],
                  lbf[spos+8], lbf[spos+9], lbf[spos+20], lbf[spos+19],
                  lbf[spos+26], lbf[spos+25], lbf[spos+10], lbf[spos+11],
                  lbf[spos+15], lbf[spos+1], lbf[spos+16], lbf[spos+5],
                  lbf[spos+17], lbf[spos+21], lbf[spos+22], lbf[spos+7],
                  lbf[spos+6], lbf[spos+13], lbf[spos+12], lbf[spos+23],
                  lbf[spos+24], moment);
        mass += onemass;
        uwall[0] -= moment;
      }
      else if(prop == 25) {
        onemass=lbbacp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+16], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+21],
                  lbf[spos+11], lbf[spos+13], lbf[spos+10], lbf[spos+12],
                  lbf[spos+14], lbf[spos+3], lbf[spos+15], lbf[spos+20],
                  lbf[spos+19], lbf[spos+17], lbf[spos+18], lbf[spos+22],
                  lbf[spos+8], lbf[spos+24], lbf[spos+26], lbf[spos+23],
                  lbf[spos+25], moment);
        mass += onemass;
        uwall[2] -= moment;
      }   
      else if(prop == 26) {
        onemass=lbfrop[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+2], lbf[spos+6], lbf[spos+7],
                  lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+22],
                  lbf[spos+10], lbf[spos+12], lbf[spos+11], lbf[spos+13],
                  lbf[spos+14], lbf[spos+16], lbf[spos+15], lbf[spos+19],
                  lbf[spos+20], lbf[spos+17], lbf[spos+18], lbf[spos+21],
                  lbf[spos+9], lbf[spos+23], lbf[spos+25], lbf[spos+24],
                  lbf[spos+26], moment);
        mass += onemass;
        uwall[2] += moment;
      }    
      else if(prop == 27) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 28) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 29) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 30) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 43) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 44) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 45) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 46) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 49) {
        onemass=lbtopp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 51) {
        onemass=lbbotp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 52) {
        onemass=lbrigp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 53) {
        onemass=lbtopp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      else if(prop == 54) {
        onemass=lblefp[ktt];
        fD3Q27VCE(onemass, uwall, &lbf[spos]);
        mass += onemass;
      }
      spos += lbsy.nq;
    }
  }
  else {
    for(int ktt=0; ktt<lbsy.nf; ktt++) {
      rho0 = lbincp[ktt];
      if(prop == 22) {
        onemass=lbtopp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5],
                  lbf[spos+6], lbf[spos+7], lbf[spos+8], lbf[spos+9],
                  lbf[spos+10], lbf[spos+11], lbf[spos+12], lbf[spos+13],
                  lbf[spos+14], lbf[spos+15], lbf[spos+16], lbf[spos+17],
                  lbf[spos+18], lbf[spos+19], lbf[spos+20], lbf[spos+21],
                  lbf[spos+22], lbf[spos+23], lbf[spos+24], lbf[spos+25],
                  lbf[spos+26], moment);
        mass += rho0;
        uwall[1] += moment;
      }
      else if(prop == 21) {
        onemass=lbbotp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+15], lbf[spos+3], lbf[spos+5], lbf[spos+4],
                  lbf[spos+6], lbf[spos+7], lbf[spos+22], lbf[spos+21],
                  lbf[spos+12], lbf[spos+13], lbf[spos+10], lbf[spos+11],
                  lbf[spos+14], lbf[spos+2], lbf[spos+16], lbf[spos+18],
                  lbf[spos+17], lbf[spos+19], lbf[spos+20], lbf[spos+9],
                  lbf[spos+8], lbf[spos+25], lbf[spos+26], lbf[spos+23],
                  lbf[spos+24], moment);
        mass += rho0;
        uwall[1] -= moment;
      }
      else if(prop == 23) {
        onemass=lbrigp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+18],
                  lbf[spos+8], lbf[spos+9], lbf[spos+6], lbf[spos+7],
                  lbf[spos+10], lbf[spos+11], lbf[spos+26], lbf[spos+25],
                  lbf[spos+15], lbf[spos+14], lbf[spos+16], lbf[spos+17],
                  lbf[spos+5], lbf[spos+21], lbf[spos+22], lbf[spos+19],
                  lbf[spos+20], lbf[spos+23], lbf[spos+24], lbf[spos+13],
                  lbf[spos+12], moment);
        mass += rho0;
        uwall[0] += moment;
      }
      else if(prop == 24) {
        onemass=lblefp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+2],
                  lbf[spos+14], lbf[spos+3], lbf[spos+18], lbf[spos+4],
                  lbf[spos+8], lbf[spos+9], lbf[spos+20], lbf[spos+19],
                  lbf[spos+26], lbf[spos+25], lbf[spos+10], lbf[spos+11],
                  lbf[spos+15], lbf[spos+1], lbf[spos+16], lbf[spos+5],
                  lbf[spos+17], lbf[spos+21], lbf[spos+22], lbf[spos+7],
                  lbf[spos+6], lbf[spos+13], lbf[spos+12], lbf[spos+23],
                  lbf[spos+24], moment);
        mass += rho0;
        uwall[0] -= moment;
      }
      else if(prop == 25) {
        onemass=lbbacp[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+16], lbf[spos+2], lbf[spos+7], lbf[spos+6],
                  lbf[spos+4], lbf[spos+5], lbf[spos+9], lbf[spos+21],
                  lbf[spos+11], lbf[spos+13], lbf[spos+10], lbf[spos+12],
                  lbf[spos+14], lbf[spos+3], lbf[spos+15], lbf[spos+20],
                  lbf[spos+19], lbf[spos+17], lbf[spos+18], lbf[spos+22],
                  lbf[spos+8], lbf[spos+24], lbf[spos+26], lbf[spos+23],
                  lbf[spos+25], moment);
        mass += rho0;
        uwall[2] -= moment;
      }
      else if(prop == 26) {
        onemass=lbfrop[ktt];
        fD3Q27PPS(onemass, lbf[spos], lbf[spos+1],
                  lbf[spos+3], lbf[spos+2], lbf[spos+6], lbf[spos+7],
                  lbf[spos+4], lbf[spos+5], lbf[spos+8], lbf[spos+22],
                  lbf[spos+10], lbf[spos+12], lbf[spos+11], lbf[spos+13],
                  lbf[spos+14], lbf[spos+16], lbf[spos+15], lbf[spos+19],
                  lbf[spos+20], lbf[spos+17], lbf[spos+18], lbf[spos+21],
                  lbf[spos+9], lbf[spos+23], lbf[spos+25], lbf[spos+24],
                  lbf[spos+26], moment);
        mass += rho0;
        uwall[2] += moment;
      }
      else if(prop == 27) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 28) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 29) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 30) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 31) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 32) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 33) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 34) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 43) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 44) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 45) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 46) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 47) {
        onemass=lbbotp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 48) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 49) {
        onemass=lbtopp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 50) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 51) {
        onemass=lbbotp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 52) {
        onemass=lbrigp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 53) {
        onemass=lbtopp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      else if(prop == 54) {
        onemass=lblefp[ktt];
        fD3Q27VCEIncom(onemass, rho0, uwall, &lbf[spos]);
        mass += rho0;
      }
      spos += lbsy.nq;
    }
  }
  uwall[0] *= fReciprocal(mass);
  uwall[1] *= fReciprocal(mass);
  uwall[2] *= fReciprocal(mass);
  return 0;  
}


int fD3Q27CPS(double p, double v0, double v1, double v2, double &f0, double &f1,
              double &f2, double &f3, double &f4, double &f5, double &f6,
              double &f7, double &f8, double &f9, double &f10, double &f11,
              double &f12, double &f13, double &f14, double &f15, double &f16,
              double &f17, double &f18, double &f19, double &f20, double &f21,
              double &f22, double &f23, double &f24, double &f25, double &f26)
{

  // produce fixed concentration at planar surface: expressed for top wall

  double rho;
  double c1=2.0/27.0,c2=1.0/54.0,c3=1.0/216.0;
  rho=6.0*(p-f0-f1-f3-f5-f6-f7-f12-f13-f14-f15-f16-f17-f19-f20-f21-f22-f23-f24)/(1.0-3.0*v1);
  f2=c1*rho*(1.0-3.0*v1);
  f4=c2*rho*(1.0-3.0*v0-3.0*v1);
  f8=c2*rho*(1.0-3.0*v1-3.0*v2);
  f9=c2*rho*(1.0-3.0*v1+3.0*v2);
  f10=c3*rho*(1.0-3.0*v0-3.0*v1-3.0*v2);
  f11=c3*rho*(1.0-3.0*v0-3.0*v1+3.0*v2);
  f18=c2*rho*(1.0+3.0*v0-3.0*v1);
  f25=c3*rho*(1.0+3.0*v0-3.0*v1+3.0*v2);
  f26=c3*rho*(1.0+3.0*v0-3.0*v1-3.0*v2);
  return 0;
}


int fD3Q27CCE(double p, double *v, double * startpos)
{

  // produce fixed solute concentration at convex or concave corner:
  // expressed for left-bottom-back corner (against velocity 19) 

  double *pt1=startpos;
  fGetEquilibriumC(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}

int fD3Q27PC(long tpos, int prop, double *uwall)
{
  double oneconc;
  long spos=tpos * lbsitelength + lbsy.nf * lbsy.nq;
  for(int ktt=0; ktt<lbsy.nc; ktt++) {
    if(prop == 22) {
      oneconc=lbtopc[ktt];
      fD3Q27CPS(oneconc, uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6],
                lbf[spos+7], lbf[spos+8], lbf[spos+9], lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+13], lbf[spos+14], lbf[spos+15], lbf[spos+16],
                lbf[spos+17], lbf[spos+18], lbf[spos+19], lbf[spos+20], lbf[spos+21],
                lbf[spos+22], lbf[spos+23], lbf[spos+24], lbf[spos+25],	lbf[spos+26]);
    }
    else if(prop == 21) {
      oneconc=lbbotc[ktt];
      fD3Q27CPS(oneconc, uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+15], lbf[spos+3], lbf[spos+5], lbf[spos+4], lbf[spos+6],
                lbf[spos+7], lbf[spos+22], lbf[spos+21], lbf[spos+12], lbf[spos+13],
                lbf[spos+10], lbf[spos+11], lbf[spos+14], lbf[spos+2], lbf[spos+16],
                lbf[spos+18], lbf[spos+17], lbf[spos+19], lbf[spos+20], lbf[spos+9],
                lbf[spos+8], lbf[spos+25], lbf[spos+26], lbf[spos+23], lbf[spos+24]);
    }
    else if(prop == 23) {
      oneconc=lbrigc[ktt];
      fD3Q27CPS(oneconc, uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+18], lbf[spos+8],
                lbf[spos+9], lbf[spos+6], lbf[spos+7], lbf[spos+10], lbf[spos+11],
                lbf[spos+26], lbf[spos+25], lbf[spos+15], lbf[spos+14], lbf[spos+16],
                lbf[spos+17], lbf[spos+5], lbf[spos+21], lbf[spos+22], lbf[spos+19],
                lbf[spos+20], lbf[spos+23], lbf[spos+24], lbf[spos+13], lbf[spos+12]);
    }
    else if(prop == 24) {
      oneconc=lblefc[ktt];
      fD3Q27CPS(oneconc, uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+14], lbf[spos+3], lbf[spos+18], lbf[spos+4], lbf[spos+8],
                lbf[spos+9], lbf[spos+20], lbf[spos+19], lbf[spos+26], lbf[spos+25],
                lbf[spos+10], lbf[spos+11], lbf[spos+15], lbf[spos+1], lbf[spos+16],
                lbf[spos+5], lbf[spos+17], lbf[spos+21], lbf[spos+22], lbf[spos+7],
                lbf[spos+6], lbf[spos+13], lbf[spos+12], lbf[spos+23], lbf[spos+24]);
    }
    else if(prop == 25) {
      oneconc=lbbacc[ktt];
      fD3Q27CPS(oneconc, uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+16], lbf[spos+2], lbf[spos+7], lbf[spos+6], lbf[spos+4],
                lbf[spos+5], lbf[spos+9], lbf[spos+21],	lbf[spos+11], lbf[spos+13],
                lbf[spos+10], lbf[spos+12], lbf[spos+14], lbf[spos+3], lbf[spos+15],
                lbf[spos+20], lbf[spos+19], lbf[spos+17], lbf[spos+18], lbf[spos+22],
                lbf[spos+8], lbf[spos+24], lbf[spos+26], lbf[spos+23], lbf[spos+25]);
    }
    else if(prop == 26) {
      oneconc=lbfroc[ktt];
      fD3Q27CPS(oneconc, uwall[0], uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+3], lbf[spos+2], lbf[spos+6], lbf[spos+7], lbf[spos+4],
                lbf[spos+5], lbf[spos+8], lbf[spos+22], lbf[spos+10], lbf[spos+12],
                lbf[spos+11], lbf[spos+13], lbf[spos+14], lbf[spos+16], lbf[spos+15],
                lbf[spos+19], lbf[spos+20], lbf[spos+17], lbf[spos+18], lbf[spos+21],
                lbf[spos+9], lbf[spos+23], lbf[spos+25], lbf[spos+24], lbf[spos+26]);
    }
    else if(prop == 27) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 28) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 29) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 30) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 31) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 43) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 44) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 45) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 46) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 47) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 48) {
      oneconc=lbrigc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 49) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 50) {
      oneconc=lblefc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 51) {
      oneconc=lbbotc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 52) {
      oneconc=lbrigc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 53) {
      oneconc=lbtopc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    else if(prop == 54) {
      oneconc=lblefc[ktt];
      fD3Q27CCE(oneconc, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}



int fD3Q27TCE(double p, double *v, double * startpos)
{

  // produce fixed temperature at convex or concave corner:
  // expressed for left-bottom-back corner (against velocity 19) 

  double *pt1=startpos;
  fGetEquilibriumT(lbfeq, v, p);
  for(int i=0; i<lbsy.nq; i++){
    *pt1 =lbfeq[i];
    pt1 ++;
  }
  pt1 = NULL;
  return 0;
}


int fD3Q27PT(long tpos, int prop, double *uwall)
{
  double onetemp;
  long spos=tpos * lbsitelength + (lbsy.nf + lbsy.nc) * lbsy.nq;
  if(lbsy.nt == 1) {
    if(prop == 22) {
      onetemp=lbtopt;
      fD3Q27CPS(onetemp, uwall[0], uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+2], lbf[spos+3], lbf[spos+4], lbf[spos+5], lbf[spos+6],
                lbf[spos+7], lbf[spos+8], lbf[spos+9], lbf[spos+10], lbf[spos+11],
                lbf[spos+12], lbf[spos+13], lbf[spos+14], lbf[spos+15], lbf[spos+16],
                lbf[spos+17], lbf[spos+18], lbf[spos+19], lbf[spos+20], lbf[spos+21],
                lbf[spos+22], lbf[spos+23], lbf[spos+24], lbf[spos+25],	lbf[spos+26]);
    }
    else if(prop == 21) {
      onetemp=lbbott;
      fD3Q27CPS(onetemp, uwall[0], -uwall[1], uwall[2], lbf[spos], lbf[spos+1],
                lbf[spos+15], lbf[spos+3], lbf[spos+5], lbf[spos+4], lbf[spos+6],
                lbf[spos+7], lbf[spos+22], lbf[spos+21], lbf[spos+12], lbf[spos+13],
                lbf[spos+10], lbf[spos+11], lbf[spos+14], lbf[spos+2], lbf[spos+16],
                lbf[spos+18], lbf[spos+17], lbf[spos+19], lbf[spos+20], lbf[spos+9],
                lbf[spos+8], lbf[spos+25], lbf[spos+26], lbf[spos+23], lbf[spos+24]);
    }
    else if(prop == 23) {
      onetemp=lbrigt;
      fD3Q27CPS(onetemp, uwall[1], uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+1], lbf[spos+3], lbf[spos+4], lbf[spos+18], lbf[spos+8],
                lbf[spos+9], lbf[spos+6], lbf[spos+7], lbf[spos+10], lbf[spos+11],
                lbf[spos+26], lbf[spos+25], lbf[spos+15], lbf[spos+14], lbf[spos+16],
                lbf[spos+17], lbf[spos+5], lbf[spos+21], lbf[spos+22], lbf[spos+19],
                lbf[spos+20], lbf[spos+23], lbf[spos+24], lbf[spos+13], lbf[spos+12]);
    }
    else if(prop == 24) {
      onetemp=lbleft;
      fD3Q27CPS(onetemp, uwall[1], -uwall[0], uwall[2], lbf[spos], lbf[spos+2],
                lbf[spos+14], lbf[spos+3], lbf[spos+18], lbf[spos+4], lbf[spos+8],
                lbf[spos+9], lbf[spos+20], lbf[spos+19], lbf[spos+26], lbf[spos+25],
                lbf[spos+10], lbf[spos+11], lbf[spos+15], lbf[spos+1], lbf[spos+16],
                lbf[spos+5], lbf[spos+17], lbf[spos+21], lbf[spos+22], lbf[spos+7],
                lbf[spos+6], lbf[spos+13], lbf[spos+12], lbf[spos+23], lbf[spos+24]);
    }
    else if(prop == 25) {
      onetemp=lbbact;
      fD3Q27CPS(onetemp, uwall[0], -uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+16], lbf[spos+2], lbf[spos+7], lbf[spos+6], lbf[spos+4],
                lbf[spos+5], lbf[spos+9], lbf[spos+21],	lbf[spos+11], lbf[spos+13],
                lbf[spos+10], lbf[spos+12], lbf[spos+14], lbf[spos+3], lbf[spos+15],
                lbf[spos+20], lbf[spos+19], lbf[spos+17], lbf[spos+18], lbf[spos+22],
                lbf[spos+8], lbf[spos+24], lbf[spos+26], lbf[spos+23], lbf[spos+25]);
    }
    else if(prop == 26) {
      onetemp=lbfrot;
      fD3Q27CPS(onetemp, uwall[0], uwall[2], uwall[1], lbf[spos], lbf[spos+1],
                lbf[spos+3], lbf[spos+2], lbf[spos+6], lbf[spos+7], lbf[spos+4],
                lbf[spos+5], lbf[spos+8], lbf[spos+22], lbf[spos+10], lbf[spos+12],
                lbf[spos+11], lbf[spos+13], lbf[spos+14], lbf[spos+16], lbf[spos+15],
                lbf[spos+19], lbf[spos+20], lbf[spos+17], lbf[spos+18], lbf[spos+21],
                lbf[spos+9], lbf[spos+23], lbf[spos+25], lbf[spos+24], lbf[spos+26]);
    }
    else if(prop == 27) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 28) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 29) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 30) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 31) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 32) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 33) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 34) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 43) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 44) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 45) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 46) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 47) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 48) {
      onetemp=lbrigt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 49) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 50) {
      onetemp=lbleft;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 51) {
      onetemp=lbbott;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 52) {
      onetemp=lbrigt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 53) {
      onetemp=lbtopt;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    else if(prop == 54) {
      onetemp=lbleft;
      fD3Q27TCE(onetemp, uwall, &lbf[spos]);
    }
    spos += lbsy.nq;
  }
  return 0;  
}


int fFixedSpeedFluid(long tpos, int prop, double *uwall)
{

  // produce boundary with fixed speed

  if(lbsy.nq == 9)
    fD2Q9VF(tpos, prop, uwall); 
  else if(lbsy.nq == 15)
    fD3Q15VF(tpos, prop, uwall);
  else if(lbsy.nq == 19)
    fD3Q19VF(tpos, prop, uwall);
  else if(lbsy.nq == 27)
    fD3Q27VF(tpos, prop, uwall);
  else if(lbdm.rank ==0) {
    cout << "the boundary condition for D" << lbsy.nd << "Q" << lbsy.nq
         << " model has not been defined" << endl;
    if(lbsy.nd == 2)
      cout << "the D2Q9 model can be used" << endl;
    else
      cout << "the D3Q15, D3Q19 or D3Q27 models can be used" << endl;
    exit(1);
  }
  
  return 0;
}

int fFixedDensityFluid(long tpos, int prop, double *uwall)
{

  // produce boundary with fixed density/pressure
  
  if(lbsy.nq == 9)
    fD2Q9PF(tpos, prop, uwall);
  else if(lbsy.nq == 15)
    fD3Q15PF(tpos, prop, uwall);
  else if(lbsy.nq == 19)
    fD3Q19PF(tpos, prop, uwall);
  else if(lbsy.nq == 27)
    fD3Q27PF(tpos, prop, uwall);
  else if(lbdm.rank ==0) {
    cout << "the boundary condition for D" << lbsy.nd << "Q" << lbsy.nq
         << " model has not been defined" << endl;
    if(lbsy.nd == 2)
      cout << "the D2Q9 model can be used" << endl;
    else
      cout << "the D3Q15, D3Q19 or D3Q27 models can be used" << endl;
    exit(1);
  }
  
  return 0;
}

int fFixedSoluteConcen(long tpos, int prop, double *uwall)
{

  // produce boundary with fixed solute concentration
  
  if(lbsy.nq == 9)
    fD2Q9PC(tpos, prop, uwall);
  else if(lbsy.nq == 15)
    fD3Q15PC(tpos, prop, uwall);
  else if(lbsy.nq == 19)
    fD3Q19PC(tpos, prop, uwall);
  else if(lbsy.nq == 27)
    fD3Q27PC(tpos, prop, uwall);
  else if(lbdm.rank ==0) {
    cout << "the boundary condition for D" << lbsy.nd << "Q" << lbsy.nq
         << " model has not been defined" << endl;
    if(lbsy.nd == 2)
      cout << "the D2Q9 model can be used" << endl;
    else
      cout << "the D3Q15, D3Q19 or D3Q27 models can be used" << endl;
    exit(1);
  }
  return 0;
}

int fFixedTemperature(long tpos, int prop, double *uwall)
{

  // produce boundary with fixed temperature
  
  if(lbsy.nq == 9)
    fD2Q9PT(tpos, prop, uwall);
  else if(lbsy.nq == 15)
    fD3Q15PT(tpos, prop, uwall);
  else if(lbsy.nq == 19)
    fD3Q19PT(tpos, prop, uwall);
  else if(lbsy.nq == 27)
    fD3Q27PT(tpos, prop, uwall);
  else if(lbdm.rank ==0) {
    cout << "the boundary condition for D" << lbsy.nd << "Q" << lbsy.nq
         << " model has not been defined" << endl;
    if(lbsy.nd == 2)
      cout << "the D2Q9 model can be used" << endl;
    else
      cout << "the D3Q15, D3Q19 or D3Q27 models can be used" << endl;
    exit(1);
  }
  return 0;
}

int fPostCollBoundary()
{
  for(long il=0; il<lbdm.touter; il++){
    if(lbphi[il] == 13) {
      fMidBounceBackF(il);
      fMidBounceBackC(il);
      fMidBounceBackT(il);
    }
  }    
  return 0;
}

int fPostPropBoundary()
{
  int a3, a12;
  double uw[3];
  if(lbsy.nt == 1) {
    lbinit += lbsysdt * lbdt;
    lbtopt += lbtopdt * lbdt;
    lbbott += lbbotdt * lbdt;
    lbfrot += lbfrodt * lbdt;
    lbbact += lbbacdt * lbdt;
    lbleft += lblefdt * lbdt;
    lbrigt += lbrigdt * lbdt;
  }
  for(long il=0; il<lbdm.touter; il++){
    if(lbphi[il] == 0) ;
    else if(lbphi[il] == 10);
    else if(lbphi[il] == 11) {
      fSiteBlankF(il);
      fSiteBlankC(il);
      fSiteBlankT(il); 
    }
    else if(lbphi[il] == 12) {
      fBounceBackF(il);
      fBounceBackC(il);
      fBounceBackT(il);
    }
    else if(lbphi[il] == 13);
    else if (lbphi[il]>110 && lbphi[il]<890) {
      a12 = lbphi[il] % 100;
      a3=lbphi[il]/100;
      if(a3 == 1) {
        fFixedSpeedFluid(il, a12, uw);
	    fFixedSoluteConcen(il, a12, uw);
	    fFixedTemperature(il, a12, uw);
      }
      else if(a3 == 2) {
	    fFixedSpeedFluid(il, a12, uw);
	    fBounceBackC(il);
	    fBounceBackT(il);
      }
      else if(a3 == 3) {
	    fFixedSpeedFluid(il, a12, uw);
	    fFixedSoluteConcen(il, a12, uw);
	    fBounceBackT(il);
      }
      else if(a3 == 4) {
	    fFixedSpeedFluid(il, a12, uw);
	    fBounceBackC(il);
	    fFixedTemperature(il, a12, uw);
      }
      else if(a3 == 5) {
	    fFixedDensityFluid(il, a12, uw);
	    fFixedSoluteConcen(il, a12, uw);
	    fFixedTemperature(il, a12, uw);
      } 
      else if(a3 == 6) {
	    fFixedDensityFluid(il, a12, uw);
	    fBounceBackC(il);
	    fFixedTemperature(il, a12, uw);
      }    
      else if(a3 == 7) {
	    fFixedDensityFluid(il, a12, uw);
	    fFixedSoluteConcen(il, a12, uw);
	    fBounceBackT(il);
      }
      else if(a3 == 8) {
	    fFixedDensityFluid(il, a12, uw);
	    fBounceBackC(il);
	    fBounceBackT(il);
      }
    }
    else {
      if(lbdm.rank == 0)
	cout << "error: boundary condition not recognised" << endl;
      exit(1);
    }
  }
  return 0;
}

int fNeighbourBoundary()
{
  int test, testphi, xpos, ypos, zpos;
  for(long il=0; il<lbdm.touter; il++) {
    fGetCoord(il, xpos, ypos, zpos);
    lbneigh[il] = 0;
    for(int m=1; m<lbsy.nq; m++) {
      testphi = lbphi[fNextStep(m, xpos, ypos, zpos)];
      test = (testphi==11 || testphi==12 || testphi==13 || (lbphi[il]>100 && testphi>100));
      lbneigh[il] += test;
    }
  }
  return 0;
}














