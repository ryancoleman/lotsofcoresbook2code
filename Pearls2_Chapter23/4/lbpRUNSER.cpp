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


// No mesoscopic interactions

// option 0: BGK with standard forcing


// option 1: BGK with Guo forcing


// option 2: BGK with EDM forcing


// option 3: MRT with standard forcing


// option 4: MRT with Guo-like forcing


// option 5: MRT with EDM forcing


// Shan-Chen pseudopotential interactions

// option 0: BGK with standard forcing

int fsBGKShanChen()
{
  do {
    if(lbcurstep%lbsave == 0 || lbcurstep==lbequstep) {
      if(lbequstep>0 && lbcurstep==lbequstep) {
        fPrintEndEquilibration();
        postequil = 1;
        fReadSpaceParameter();
        fNeighbourBoundary();
      }
      if(lbcurstep>=lbequstep)
        fsOutput();
      cout << lbcurstep << " \t ";
      fPrintDomainMass();
      fPrintDomainMomentum();
    }
    fInteractionForceZero();
    fConvectionForceBoussinesq(lbbousth, lbboustl);
    fCalcPotential_ShanChen();
    fsInteractionForceShanChen();
    fCollisionBGK();
    fPostCollBoundary();
    fPropagationSwap();
    fPostPropBoundary();
    lbcurstep++;
  }
  while(lbcurstep<=lbtotstep);
  return 0;
}

// option 1: BGK with Guo forcing

int fsBGKGuoShanChen()
{
  do {
    if(lbcurstep%lbsave == 0 || lbcurstep==lbequstep) {
      if(lbequstep>0 && lbcurstep==lbequstep) {
        fPrintEndEquilibration();
        postequil = 1;
        fReadSpaceParameter();
        fNeighbourBoundary();
      }
      if(lbcurstep>=lbequstep)
        fsOutput();
      cout << lbcurstep << " \t ";
      fPrintDomainMass();
      fPrintDomainMomentum();
    }
    fInteractionForceZero();
    fConvectionForceBoussinesq(lbbousth, lbboustl);
    fCalcPotential_ShanChen();
    fsInteractionForceShanChen();
    fCollisionBGKGuo();
    fPostCollBoundary();
    fPropagationSwap();
    fPostPropBoundary();
    lbcurstep++;
  }
  while(lbcurstep<=lbtotstep);
  return 0;
}

// option 2: BGK with EDM forcing

int fsBGKEDMShanChen()
{
  do {
    if(lbcurstep%lbsave == 0 || lbcurstep==lbequstep) {
      if(lbequstep>0 && lbcurstep==lbequstep) {
        fPrintEndEquilibration();
        postequil = 1;
        fReadSpaceParameter();
        fNeighbourBoundary();
      }
      if(lbcurstep>=lbequstep)
        fsOutput();
      cout << lbcurstep << " \t ";
      fPrintDomainMass();
      fPrintDomainMomentum();
    }
    fInteractionForceZero();
    fConvectionForceBoussinesq(lbbousth, lbboustl);
    fCalcPotential_ShanChen();
    fsInteractionForceShanChen();
    fCollisionBGKEDM();
    fPostCollBoundary();
    fPropagationSwap();
    fPostPropBoundary();
    lbcurstep++;
  }
  while(lbcurstep<=lbtotstep);
  return 0;
}

// option 3: MRT with standard forcing


// option 4: MRT with Guo-like forcing


// option 5: MRT with EDM forcing


// Shan/Chen pseudopotential interactions with quadratic pseudopotential term
// TO DO

// Lishchuk continuum-based interactions (with interfacial normals calculated non-locally)

// option 0: BGK with standard forcing


// option 1: BGK with Guo forcing


// option 2: BGK with EDM forcing


// option 3: MRT with standard forcing


// option 4: MRT with Guo-like forcing


// option 5: MRT with EDM forcing


// Lishchuk continuum-based interactions (with interfacial normals calculated locally)
// TO DO


