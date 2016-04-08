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

int fSiteFluidCollisionBGK(double* startpos, double *sitespeed, double* bodyforce);
int fSiteFluidIncomCollisionBGK(double* startpos, double *sitespeed, double* bodyforce);
int fSiteFluidCollisionBGKEDM(double* startpos, double *sitespeed, double* bodyforce);
int fSiteFluidIncomCollisionBGKEDM(double* startpos, double *sitespeed, double* bodyforce);
int fSiteFluidCollisionBGKGuo(double* startpos, double *sitespeed, double* bodyforce);
int fSiteFluidIncomCollisionBGKGuo(double* startpos, double *sitespeed, double* bodyforce);
int fSiteSoluteCollisionBGK(double* startpos, double *sitespeed);
int fSiteThermalCollisionBGK(double* startpos, double *sitespeed);
int fCollisionBGK();
int fCollisionBGKEDM();
int fCollisionBGKGuo();

