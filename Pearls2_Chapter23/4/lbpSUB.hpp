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

inline void fWeakMemory();
int fMemoryAllocation();
int fFreeMemory();
int fSetSerialDomain();
int fStartDLMESO();
int fFinishDLMESO();
int fsPrintDomainInfo();
int fGetModel();
int fInitializeSystem();
int fGetEquilibriumF(double *feq, double *v, double rho);
int fGetEquilibriumFIncom(double *feq, double *v, double rho, double rho0);
int fGetEquilibriumC(double *feq, double *v, double rho);
int fGetEquilibriumT(double *feq, double *v, double rho);
int fGetMomentEquilibriumF(double *feq, double *p, double rho);
int fGetMomentEquilibriumFIncom(double *feq, double *p, double rho, double rho0);
int fGetEDMForce(double *source, double *v, double *dv, double rho);
int fPropagationSwap();

