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


inline long fGetNodePosi(int xpos, int ypos, int zpos)
{
  return (xpos * lbdm.youter + ypos) * lbdm.zouter + zpos;
}

inline long fGetNodePosi(int xpos, int ypos)
{
  return xpos * lbdm.youter + ypos;
}

int fGetCoord(long tpos, int& xpos, int& ypos, int& zpos)
{
  zpos = tpos % lbdm.zouter;
  xpos = int (tpos / (lbdm.youter * lbdm.zouter));
  ypos = int(tpos % (lbdm.youter * lbdm.zouter)) / lbdm.zouter;
  return 0;
}


double fGetOneMassSite(double* startpos)
{

  // calculate mass density at grid position starting from startpos in lbf[] 
  // array

  double mass=0.0;
  double *pt1; 
  pt1= startpos;
  for(int i=0; i<lbsy.nq; i++) {
    mass += *pt1;
    pt1 ++;
  }
  pt1 = NULL;
  return mass;
}




double fGetTotMassSite(double* startpos)
{

  // calculate total mass at grid position starting from startpos in lbf[] 
  // array

  double mass=0.0;
  double *pt1 = startpos;
  for(int i=0; i<lbsy.nf * lbsy.nq; i++) {
    mass += *pt1;
    pt1 ++;
  }
  pt1 = NULL;
  return mass;
}




double fGetOneMassDomain(int fpos)
{

  // calculate total mass of one fluid in the domain

  double totmass =0;
  for(long il=0; il<lbdm.touter; il++)
    if(lbphi[il] == 0)
      totmass += fGetOneMassSite(&lbf[il * lbsitelength + fpos * lbsy.nq]); 
  return totmass;
}


double fGetTotMassDomain()
{

  // calculate total mass of all fluids in the domain

  double totmass =0;
  for(long il=0; il<lbdm.touter; il++)
    if(lbphi[il] == 0)
      totmass +=  fGetTotMassSite(&lbf[il * lbsitelength]);
  return totmass;
}


double fGetFracSite(int fpos, double* startpos)
{

  // calculate mass fraction of fpos phase at grid position starting from 
  // startpos in lbf[] array

  int fs = fpos * lbsy.nq;
  int fe = fs + lbsy.nq;
  double mass=0.0, fmass=0.0;
  double *pt1 = startpos;
    
  for(int i=0; i<lbsy.nf * lbsy.nq; i++) {
    if(i >= fs && i < fe)
      fmass += *pt1;
    mass += *pt1;
    pt1 ++;
  }
  pt1 = NULL;

  if(mass<lbevaplim) return 0.0;
    
  return fmass*fReciprocal(mass);
}






int fGetTotMomentSite(double *momentum, double* startpos)
{

  // calculate momentum of all fluids at grid position starting from startpos 
  // in lbf[] array
  
  double *pt1 = startpos;
  momentum[0]=0; momentum[1]=0; momentum[2]=0;
  for(int j=0; j<lbsy.nf; j++)
    for(int i=0; i<lbsy.nq; i++) {
      momentum[0] += (*pt1) * lbv[3*i];
      momentum[1] += (*pt1) * lbv[3*i+1];
      momentum[2] += (*pt1) * lbv[3*i+2];
      pt1 ++;
    }
  pt1 = NULL;
  return 0;
}


int fGetTotMomentDomain(double *momentum)
{

  // calculate total momentum of all fluids in the domain

  double speed[3];
  for(int i=0; i<3; i++)
    momentum[i] = 0;
  for(long il=0; il<lbdm.touter; il++)
    if(lbphi[il] == 0) {
      fGetTotMomentSite(speed, &lbf[il * lbsitelength]);
      for(int i=0; i<3; i++) 
	momentum[i] += speed[i];
    }
  return 0;  
}


int fGetSpeedSite(double *speed, double* startpos)
{

  // calculate macroscopic speed of all fluids at grid position starting from 
  // startpos in lbf[] array
  
  double mass=0.0,invmass;
  double *pt1 = startpos;
  int counter_pt1 = 0;
  speed[0]=0; speed[1]=0; speed[2]=0;

  for(int j=0; j<lbsy.nf; j++)
    for(int i=0; i<lbsy.nq; i++,counter_pt1++) {
      mass += pt1[counter_pt1];
      speed[0] += pt1[counter_pt1] * lbv[3*i];
      speed[1] += pt1[counter_pt1] * lbv[3*i+1];
      speed[2] += pt1[counter_pt1] * lbv[3*i+2];
    }

  invmass = fReciprocal(mass);
  speed[0] *= invmass;
  speed[1] *= invmass;
  speed[2] *= invmass; 
  return 0;
}


int fGetSpeedIncomSite(double *speed, double* startpos)
{

  // calculate macroscopic speed of all incompressible fluids at grid position
  // starting from startpos in lbf[] array
  
  double mass=0.0;
  double *pt1 = startpos;
  speed[0]=0; speed[1]=0; speed[2]=0;

  for(int j=0; j<lbsy.nf; j++) {
    mass += lbincp[j];
    for(int i=0; i<lbsy.nq; i++) {
      speed[0] += (*pt1) * lbv[3*i];
      speed[1] += (*pt1) * lbv[3*i+1];
      speed[2] += (*pt1) * lbv[3*i+2];
      pt1 ++;
    }
  }
  pt1 = NULL;
  
  speed[0] *= fReciprocal(mass);
  speed[1] *= fReciprocal(mass);
  speed[2] *= fReciprocal(mass); 
  return 0;
}


float fGetOneDirecSpeedSite(int dire, double* startpos)
{

  // calculate site speed along direction dire at grid position starting from 
  // startpos in lbf[] array

  double speed = 0;
  double mass=0.0;
  double *pt1 = startpos;
  for(int j=0; j<lbsy.nf; j++)
    for(int i=0; i<lbsy.nq; i++){
      mass += *pt1;
      speed += (*pt1) * lbv[3*i + dire];
      pt1 ++;
    }
  pt1 = NULL;
  if(mass<lbevaplim) return 0.0;
  return float(speed*fReciprocal(mass));
}


float fGetOneDirecSpeedIncomSite(int dire, double* startpos)
{

  // calculate site speed for incompressible fluid along direction dire
  // at grid position starting from startpos in lbf[] array

  double speed = 0;
  double mass=0.0;
  double *pt1 = startpos;
  for(int j=0; j<lbsy.nf; j++) {
    mass += lbincp[j];
    for(int i=0; i<lbsy.nq; i++){
      speed += (*pt1) * lbv[3*i + dire];
      pt1 ++;
    }
  }
  pt1 = NULL;
  return float(speed*fReciprocal(mass));
}










double fGetOneConcSite(int cpos, long tpos)
{
  return fGetOneMassSite(&lbf[tpos*lbsitelength+(lbsy.nf + cpos)* lbsy.nq]);
}



double fGetTemperatureSite(long tpos)
{
  return fGetOneMassSite(&lbf[tpos*lbsitelength+
			    (lbsy.nf + lbsy.nc)* lbsy.nq]);
}























