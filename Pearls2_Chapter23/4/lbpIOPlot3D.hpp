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

int fsOutputGrid(const char* filename);

int fsOutputQ(const char* filename);





int fsOutputGrid3D(const char* filename);
int fsOutputGrid2D(const char* filename);

int fsOutputQP3D(const char* filename, int iprop);
int fsOutputQP2D(const char* filename, int iprop);

int fsOutputQCA3D(const char* filename, int iprop);
int fsOutputQCA2D(const char* filename, int iprop);

int fsOutputQCB3D(const char* filename, int iprop);
int fsOutputQCB2D(const char* filename, int iprop);

int fsOutputQT3D(const char* filename);
int fsOutputQT2D(const char* filename);
