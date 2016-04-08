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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#ifdef __linux__
#include <sys/time.h>
#include <sys/stat.h>
#else
#include <windows.h>
#endif

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_num_thread() 1
  #define omp_get_thread_num() 0
#endif
using namespace std;

/*
  DL_MESO_LBE
  Authors   :   R.S. Qin, M. A. Seaton
  Copyright :   STFC Daresbury Laboratory
            :   05/11/2004, rev. 18/07/2014
  
*/


// system information

struct sSystem
{
  int nd;
  int nq;
  int nf;
  int nc;
  int nt;
  int pf;
  int nx;
  int ny;
  int nz;
};


// domain information

struct sDomain
{
  int rank;
  int size;
  int bwid;
  int owidx;
  int owidy;
  int owidz;
  int xcor;
  int ycor;
  int zcor;
  int xdim;
  int ydim;
  int zdim;
  int xs;
  int xe;
  int ys;
  int ye;
  int zs;
  int ze;
  int xinner;
  int yinner;
  int zinner;
  int xouter;
  int youter;
  int zouter;
  int touter;
};

// neighbour information


/*
  DL_MESO_LBE
  Author    :   R.S. Qin
  Copyright :   Daresbury Laboratory
            :   2004
  
*/


// structure objects

sSystem lbsy;
sDomain lbdm;

// parameters

double lbiniv[3];
double lbtopv[3];
double lbbotv[3];
double lbfrov[3];
double lbbacv[3];
double lblefv[3];
double lbrigv[3];

double* lbincp;
double* lbinip;
double* lbtopp;
double* lbbotp;
double* lbfrop;
double* lbbacp;
double* lblefp;
double* lbrigp;

double* lbinic;
double* lbtopc;
double* lbbotc;
double* lbfroc;
double* lbbacc;
double* lblefc;
double* lbrigc;

double lbinit;
double lbtopt;
double lbbott;
double lbfrot;
double lbbact;
double lbleft;
double lbrigt;

double lbsysdt;
double lbtopdt;
double lbbotdt;
double lbfrodt;
double lbbacdt;
double lblefdt;
double lbrigdt;


double lbxsize;
int    lbtotstep;
int    lbequstep;
int    lbsave;
int    lbcurstep;
int    lbsteer;
int    lbdump;
int    lbrestart;

// double lbsmw;
// double lbsmv;
// int    *lbanifold;
// double *lbanitens; 
// double *lbemradius;
// int    *lbemnumber;

double lbcalctime;
double lbendtime;
double lbnoise;
double *lbtf;
double *lbtfbulk;
double *lbtc;
double *lbtt;
int *lbscpot;
double *lbg;
double *lbgwall;
double *lbseg;
double *lbpsi0;
double *lbbdforce;
double *lbbousforce;
double *lbinterforce;

double *lbf;
double *lbft;
double *lbfeq;
int    *lbphi;
int    *lbneigh;

int    *lbv;
int    *lbopv;
double *lbw;
double *lbvw;
double *lbtr;
double *lbtrinv;
double  lbmrts[3];
double  lbmrtw[3];

int    lbsitelength;
double lbdx;
double lbdt;
double lbcs;
double lbcssq;
double lbrcssq;
double lbsoundv;
double lbreynolds;
double lbkinetic;
int    lbyz;
double lbxst;
double lbtcoe; 
double lbbousth;
double lbboustl;
double lbevaplim;

unsigned long *lbouter;
int lboutersize;

// calculation parameters

int    bigend;
double timetotal;
int    qVersion;
int    collide;
int    interact;
int    incompress;
int    outformat;
int    postequil;

/*
  DL_MESO_LBE
  Author    :   R.S. Qin, M. A. Seaton
  Copyright :   STFC Daresbury Laboratory
            :   05/11/2004, rev. 29/07/2014

*/  
//	    DL_MESO_LBE functions
//	    camel
//          p-package
//          f-function 
//          s-structure
//          c-class

#include "lbpBASIC.hpp"
#include "lbpMODEL.hpp"
#include "lbpIO.hpp"
#include "lbpIOPlot3D.hpp"
#include "lbpIOLegacyVTK.hpp"
#include "lbpIOVTK.hpp"
#include "lbpGET.hpp"
#include "lbpFORCE.hpp"
#include "lbpBGK.hpp"
#include "lbpSUB.hpp"
#include "lbpBOUND.hpp"
#include "lbpRUNSER.hpp"
#include "lbpBASIC.cpp"
#include "lbpMODEL.cpp"
#include "lbpIO.cpp"
#include "lbpIOPlot3D.cpp"
#include "lbpIOLegacyVTK.cpp"
#include "lbpIOVTK.cpp"
#include "lbpGET.cpp"
#include "lbpFORCE.cpp"
#include "lbpBGK.cpp"
#include "lbpSUB.cpp"
#include "lbpBOUND.cpp"
#include "lbpRUNSER.cpp"

