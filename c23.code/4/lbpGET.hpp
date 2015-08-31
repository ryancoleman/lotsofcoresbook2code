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

inline long fGetNodePosi(int xpos, int ypos, int zpos);
inline long fGetNodePosi(int xpos, int ypos);
int fGetCoord(long tpos, int& xpos, int& ypos, int& zpos);
double fGetOneMassSite(double* startpos);
double fGetTotMassSite(double* startpos);
double fGetOneMassDomain(int fpos);
double fGetTotMassDomain();
double fGetFracSite(int fpos, double* startpos);
int fGetTotMomentSite(double *speed, double* startpos);
int fGetTotMomentDomain(double *momentum);
int fGetSpeedSite(double *speed, double* startpos);
int fGetSpeedIncomSite(double *speed, double* startpos);
float fGetOneDirecSpeedSite(int dire, double* startpos);
float fGetOneDirecSpeedIncomSite(int dire, double* startpos);
double fGetOneConcSite(int cpos, long tpos);
double fGetTemperatureSite(long tpos);
