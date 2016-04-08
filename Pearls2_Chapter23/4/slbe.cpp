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

#include "slbe.hpp"

int main(int argc, char* argv[])
{
  bigend = fBigEndian();
  lbcurstep = 0;
  fDefineSystem();
  fSetSerialDomain();
  fStartDLMESO();
  fMemoryAllocation();
  fPrintSystemInfo();
  fsPrintDomainInfo();
  fInputParameters();
  fReadSpaceParameter();
  fGetModel();
  fNeighbourBoundary();
  fInitializeSystem();
  fReadInitialState();
  if(outformat==2)  // Grid for Plot3D output files
    fsOutputGrid(); 
  if(interact==1 || interact==2)
    fCalcPotential_ShanChen();
  timetotal=fCheckTimeSerial();
  switch (interact) {
    case 0:
    // no mesoscopic interactions
      switch (collide) {
        case 0:
          break;
        case 1:
          break;
        case 2:
          break;
        case 3:
          break;
        case 4:
          break;
        case 5:
          break;
      }
      break;
    case 1:
    // Shan/Chen pseudopotential interactions
      switch (collide) {
        case 0:
          fsBGKShanChen();
          break;
        case 1:
          fsBGKEDMShanChen();
          break;
        case 2:
          fsBGKGuoShanChen();
          break;
        case 3:
          break;
        case 4:
          break;
        case 5:
          break;
      }
      break;
    case 2:
    // Shan/Chen pseudopotential interactions with quadratic pseudopotential term
    // TO DO
      break;
    case 3:
    // Lishchuk continuum-based interactions (with interfacial normals determined non-locally)
      switch (collide) {
        case 0:
          break;
        case 1:
          break;
        case 2:
          break;
        case 3:
          break;
        case 4:
          break;
        case 5:
          break;
      }
      break;
    case 4:
    // Lishchuk continuum-based interactions (with interfacial normals determined locally)
    // TO DO
      break;
  }
  timetotal=fCheckTimeSerial();
  fFreeMemory();
  fFinishDLMESO();
  return 0;
}

