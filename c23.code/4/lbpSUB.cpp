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


inline void fWeakMemory() 
{ 
  cout << "error: cannot allocate more memory" << endl;
  exit(1);
  return;
}

int fMemoryAllocation()
{

  // allocate memory for lattice Boltzmann calculation

  unsigned long dim = lbsitelength * lbdm.touter;
  unsigned long dim1 = lbsy.nf + lbsy.nc + lbsy.nt;
  unsigned long dim2 = 3 * lbsy.nf * (lbsy.nf-1) / 2;

  long i;
  unsigned long dimmax, j;

  if(dim != 0) {
    lbf = new double[dim];  
    if(lbf == NULL) fWeakMemory();
    for(j=0; j<dim; j++)
      lbf[j] = 0.0;
  }

  dimmax = dim1;
  if (interact==3 || interact==4)
    dimmax = fCppMax (dim1, dim2);

  if(lbdm.touter * dimmax != 0) {
    lbft = new double[lbdm.touter * dimmax];
    if(lbft == NULL) fWeakMemory();
    for(i=0; i< lbdm.touter * dimmax; i++)
      lbft[i] = 0.0;
  }

  if(lbsy.nq != 0) {
    lbfeq = new double[lbsy.nq];
    if(lbfeq == NULL) fWeakMemory();
    for(i=0; i< lbsy.nq; i++)
      lbfeq[i] = 0.0;
  }

  if(lbdm.touter != 0) {
    lbphi = new int[lbdm.touter];
    if(lbphi == NULL) fWeakMemory();
    for(i=0; i< lbdm.touter; i++)
      lbphi[i] = 0;
  }

  if(lbdm.touter != 0) {
    lbneigh = new int[lbdm.touter];
    if(lbneigh == NULL) fWeakMemory();
    for(i=0; i< lbdm.touter; i++)
      lbneigh[i] = 0;
  }

  if(lbsy.nq != 0) {
    lbw = new double[lbsy.nq];
    if(lbw == NULL) fWeakMemory();
    for(i=0; i< lbsy.nq; i++)
      lbw[i] = 0;
  }

  if(lbsy.nq != 0) {
    lbv = new int[3 * lbsy.nq];
    if(lbv == NULL) fWeakMemory();
    for(i=0; i< lbsy.nq; i++)
      lbv[i] = 0;
  }

  if(lbsy.nq != 0) {
    lbvw = new double[3*lbsy.nq];
    if(lbvw == NULL) fWeakMemory();
    for(i=0; i< lbsy.nq; i++)
      lbvw[i] = 0;
  }

  if(lbsy.nq != 0) {
    lbopv = new int[lbsy.nq];
    if(lbopv == NULL) fWeakMemory();
    for(i=0; i< lbsy.nq; i++)
      lbopv[i] = 0;
  }

  if(lbsy.nq != 0) {
    lbtr = new double[lbsy.nq*lbsy.nq];
    if(lbtr == NULL) fWeakMemory();
    for(i=0; i<lbsy.nq*lbsy.nq; i++)
      lbtr[i] = 0;
  }

  if(lbsy.nq != 0) {
    lbtrinv = new double[lbsy.nq*lbsy.nq];
    if(lbtrinv == NULL) fWeakMemory();
    for(i=0; i<lbsy.nf*lbsy.nq; i++)
      lbtrinv[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbtf = new double[lbsy.nf];
    if(lbtf == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbtf[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbtfbulk = new double[lbsy.nf];
    if(lbtfbulk == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbtfbulk[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbbdforce = new double[3 * lbsy.nf];
    if(lbbdforce == NULL) fWeakMemory();
    for(i=0; i< 3 * lbsy.nf; i++)
      lbbdforce[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbbousforce = new double[3 * lbsy.nf];
    if(lbbousforce == NULL) fWeakMemory();
    for(i=0; i< 3 * lbsy.nf; i++)
      lbbousforce[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbinterforce = new double[3 * lbsy.nf * lbdm.touter];
    if( lbinterforce== NULL) fWeakMemory();
    for(i=0; i< 3 * lbsy.nf * lbdm.touter; i++)
      lbinterforce[i] = 0;
  }


  if(lbsy.nc != 0) {
    lbtc = new double[lbsy.nc];
    if(lbtc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbtc[i] = 0;
  }

  if(lbsy.nt != 0) {
    lbtt = new double[lbsy.nt];
    if(lbtt == NULL) fWeakMemory();
    for(i=0; i< lbsy.nt; i++)
      lbtt[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbscpot = new int[lbsy.nf];
    if(lbscpot == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbscpot[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbg = new double[lbsy.nf * lbsy.nf];
    if(lbg == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf * lbsy.nf; i++)
      lbg[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbgwall = new double[lbsy.nf];
    if(lbgwall == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbgwall[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbseg = new double[lbsy.nf * lbsy.nf];
    if(lbseg == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf * lbsy.nf; i++)
      lbseg[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbpsi0 = new double[lbsy.nf];
    if(lbpsi0 == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbpsi0[i] = 0;
  }

/*
  if(lbsy.nf != 0) {
    lbemradius = new double[lbsy.nf];
    if(lbemradius == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbemradius[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbemnumber = new int[lbsy.nf];
    if(lbemnumber == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbemnumber[i] = 0;
  }
*/

  if(lbsy.nf != 0) {
    lbincp = new double[lbsy.nf];
    if(lbincp == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbincp[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbinip = new double[lbsy.nf];
    if(lbinip == NULL)  fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbinip[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbtopp = new double[lbsy.nf];
    if(lbtopp == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbtopp[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbbotp = new double[lbsy.nf];
    if(lbbotp == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbbotp[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbfrop = new double[lbsy.nf];
    if(lbfrop == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbfrop[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbbacp = new double[lbsy.nf];
    if(lbbacp == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbbacp[i] = 0;
  }

  if(lbsy.nf != 0) {
    lblefp = new double[lbsy.nf];
    if(lblefp == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lblefp[i] = 0;
  }

  if(lbsy.nf != 0) {
    lbrigp = new double[lbsy.nf];
    if(lbrigp == NULL) fWeakMemory();
    for(i=0; i< lbsy.nf; i++)
      lbrigp[i] = 0;
  }

/*
  if((lbsy.nf-1) != 0) {
    lbanifold = new int[(lbsy.nf * (lbsy.nf-1)) / 2];
    if(lbanifold == NULL) fWeakMemory();
    for(i=0; i< (lbsy.nf * (lbsy.nf-1)) / 2; i++)
      lbanifold[i] = 0;
  }

  if((lbsy.nf-1) != 0) {
    lbanitens = new double[(lbsy.nf * (lbsy.nf-1)) / 2];
    if(lbanitens == NULL) fWeakMemory();
    for(i=0; i< (lbsy.nf * (lbsy.nf-1)) / 2; i++)
      lbanitens[i] = 0;
  }
*/

  if(lbsy.nc != 0) {
    lbinic = new double[lbsy.nc];
    if(lbinic == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbinic[i] = 0;
  }

  if(lbsy.nc != 0) {
    lbtopc = new double[lbsy.nc];
    if(lbtopc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbtopc[i] = 0;
  }

  if(lbsy.nc != 0) {
    lbbotc = new double[lbsy.nc];
    if(lbbotc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbbotc[i] = 0;
  }

  if(lbsy.nc != 0) {
    lbfroc = new double[lbsy.nc];
    if(lbfroc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbfroc[i] = 0;
  }

  if(lbsy.nc != 0) {
    lbbacc = new double[lbsy.nc];
    if(lbbacc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbbacc[i] = 0;
  }

  if(lbsy.nc != 0) {
    lblefc = new double[lbsy.nc];
    if(lblefc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lblefc[i] = 0;
  }

  if(lbsy.nc != 0) {
    lbrigc = new double[lbsy.nc];
    if(lbrigc == NULL) fWeakMemory();
    for(i=0; i< lbsy.nc; i++)
      lbrigc[i] = 0;
  }

  // zero velocity arrays

  for(i=0; i<3; i++) {
    lbiniv[i] = 0;
    lbtopv[i] = 0;
    lbbotv[i] = 0;
    lbfrov[i] = 0;
    lbbacv[i] = 0;
    lblefv[i] = 0;
    lbrigv[i] = 0;
  }

  return 0;

}


int fFreeMemory()
{
  
  // free all allocated memory
  
  delete [] lbf;  
  delete [] lbft;
  delete [] lbphi;
  delete [] lbw;
  delete [] lbtf;
  delete [] lbtfbulk;
  delete [] lbtc;
  delete [] lbg;
  delete [] lbgwall;
  delete [] lbv;
  delete [] lbvw;
  delete [] lbopv;
  delete [] lbtr;
  delete [] lbtrinv;
  delete [] lbfeq;
//  delete [] lbemnumber;
//  delete [] lbemradius;
//  delete [] lbanitens;
//  delete [] lbanifold;
  delete [] lbinic;
  delete [] lbtopc;
  delete [] lbbotc;
  delete [] lbfroc;
  delete [] lbbacc;
  delete [] lblefc;
  delete [] lbrigc;
  delete [] lbincp;
  delete [] lbinip;
  delete [] lbtopp;
  delete [] lbbotp;
  delete [] lbfrop;
  delete [] lbbacp;
  delete [] lblefp;
  delete [] lbrigp;
  delete [] lbbdforce;
  delete [] lbbousforce;
  delete [] lbouter;

  return 0;
  
}

int fSetSerialDomain()
{
  int ii,jj,kk;
  lbdm.rank = 0;
  lbdm.size = 0;
  lbdm.xs = 0;
  lbdm.ys = 0;
  lbdm.xe = lbsy.nx;
  lbdm.ye = lbsy.ny;
  lbdm.zs = 0;
  lbdm.ze = lbsy.nz;
  lbdm.xinner = lbsy.nx;
  lbdm.yinner = lbsy.ny;
  lbdm.zinner = lbsy.nz;
  lbdm.xouter = lbsy.nx;
  lbdm.youter = lbsy.ny;
  lbdm.zouter = lbsy.nz;
  lbdm.xcor=0;
  lbdm.ycor=0;
  lbdm.zcor=0;
  lbdm.xdim=1;
  lbdm.ydim=1;
  lbdm.zdim=1;
  lbdm.touter = lbdm.xouter * lbdm.youter * lbdm.zouter;
  lbdm.owidx = lbdm.bwid;
  lbdm.owidy = lbdm.bwid;
  lbdm.owidz = lbdm.bwid;
  if (lbdm.owidx<1)
    lbdm.owidx = 1;
  if (lbdm.owidy<1)
    lbdm.owidy = 1;
  if (lbdm.owidz<1)
    lbdm.owidz = 1;
  lbdm.bwid = 0;
  if(lbsy.nd == 2) {
    lbdm.zs = 0;
    lbdm.ze = 1;
    lbdm.zinner = 1;
    lbdm.zouter = 1;
    lbdm.touter = lbdm.xouter * lbdm.youter;
    lbdm.owidz = 0;
  }
  // set grid boundary regions
  lboutersize = 0;
  ii = 2*lbdm.xouter*lbdm.youter*lbdm.owidz +
       2*lbdm.xouter*(lbdm.zouter-2*lbdm.owidz)*lbdm.owidy +
       2*(lbdm.youter-2*lbdm.owidy)*(lbdm.zouter-2*lbdm.owidz)*lbdm.owidx;
  lbouter = new unsigned long[4*ii];
  for(ii=0; ii<lbdm.xouter; ii++) {
    for(jj=0; jj<lbdm.youter; jj++) {
      for(kk=0; kk<lbdm.zouter; kk++) {
        if ((ii<lbdm.owidx) || (jj<lbdm.owidy) || (kk<lbdm.owidz) || (ii>=(lbdm.xouter-lbdm.owidx)) ||
            (jj>=(lbdm.youter-lbdm.owidy)) || (kk>=(lbdm.zouter-lbdm.owidz))) {
          lbouter[4*lboutersize  ] = (unsigned long) ((ii*lbdm.youter+jj)*lbdm.zouter+kk);
          lbouter[4*lboutersize+1] = (unsigned long) ii;
          lbouter[4*lboutersize+2] = (unsigned long) jj;
          lbouter[4*lboutersize+3] = (unsigned long) kk;
          lboutersize++;
        }
      }
    }
  }
  return 0;
}


int fStartDLMESO()
{
  
  if(lbdm.rank == 0) {
    cout << endl;
    cout << "Welcome to DL_MESO_LBE" << endl << endl;
    cout << "STFC/CCP5 Program Library Package" << endl;
    cout << "STFC Daresbury Laboratory Lattice Boltzmann Equation program" << endl << endl;
    cout << "DL_MESO version 2.6, July 2014" << endl;
    cout << "Authors: R. S. Qin & M. A. Seaton" << endl;
    cout << "Copyright (c) 2014 STFC Daresbury Laboratory, United Kingdom" << endl << endl;
    cout << "Any comments or queries, please contact Dr. M. A. Seaton at:" << endl;
    cout << "   Email: michael.seaton@stfc.ac.uk" << endl;
    cout << "   Tel:   +44 (0)1925 603317" << endl << endl;
    time_t starttime;
    struct tm * timeinfo;
    time ( &starttime);
    timeinfo = localtime ( &starttime );
    cout << "Program started at: " << asctime (timeinfo) << endl;
  }
  
  return 0;
}


int fFinishDLMESO()
{
  time_t endtime;
  struct tm * timeinfo;
  time ( &endtime);
  timeinfo = localtime ( &endtime );
  if(lbdm.rank == 0) {
    cout << "Calculation time elapsed: " << timetotal << " seconds" << endl;
    cout << "Efficiency measure: " << (0.000001*lbtotstep*(lbsy.nf+lbsy.nc+lbsy.nt)*lbsy.nx*lbsy.ny*lbsy.nz/timetotal) 
         << " MLUPS" << endl;
    cout << "Program finished at: " << asctime (timeinfo) << endl << endl;
    cout <<  "----------------------------------------------------------------------------" << endl;
    cout <<  "   This cut down version used with permission from copyright owners  " << endl;
    cout <<  "----------------------------------------------------------------------------" << endl;
    cout <<  "   The full DL_MESO package is supplied to individuals under an      " << endl;
    cout <<  "   academic licence, which is free of cost to academic scientists    " << endl;
    cout <<  "   pursuing scientific research of a non-commercial nature           " << endl;
    cout <<  "   To register please visit www.ccp5.ac.uk/DL_MESO                   " << endl;
    cout <<  "   Commercial organisations interested in acquiring the package      " << endl;
    cout <<  "   should approach Dr. M. A. Seaton at Daresbury Laboratory in the   " << endl;
    cout <<  "   first instance. Daresbury Laboratory is the sole centre for       " << endl;
    cout <<  "   distribution of the package                                       " << endl;
    cout <<  "----------------------------------------------------------------------------" << endl;

    cout << "Many thanks for using DL_MESO_LBE in your work. Please acknowledge" << endl;
    cout << "our efforts by including one of the following references when"<<endl;
    cout << "publishing data obtained using DL_MESO_LBE:" << endl << endl;
    cout << "   M. A. Seaton, R. L. Anderson, S. Metz and W. Smith, \"DL_MESO:" << endl;
    cout << "   highly scalable mesoscale simulations\", Mol. Sim. 39 (10), 796-821" << endl;
    cout << "   (2013), doi:10.1080/08927022.2013.772297" << endl << endl;
    cout << "   M. A. Seaton, \"The DL_MESO Mesoscale Simulation Package\", STFC" << endl;
    cout << "   Scientific Computing Department (2012), www.ccp5.ac.uk/DL_MESO" << endl << endl << endl;
  }

  return 0;

}


int fsPrintDomainInfo()
{

  // print domain information (serial version: prints number of threads if using OpenMP)

  int nthreads=1;

#ifdef _OPENMP
  #pragma omp parallel
  {
    #pragma omp single
    {
      nthreads = omp_get_num_threads();
      cout << "Running with " << nthreads << " threads" << endl;
    }
  }
#endif

  return 0;
}


int fGetModel()
{

  // initialize vectors lbv, lbw, lbopv

  if(lbsy.nq == 9) 
    D2Q9();

  else if(lbsy.nq == 15)
    D3Q15();

  else if(lbsy.nq == 19)
    D3Q19(); 

  else if(lbsy.nq == 27)
    D3Q27();

  else 
    {
      if(lbdm.rank == 0) {
        cout << "error: the requested model is not available" << endl;
        }
      exit(1);
    }
    
  lbdx = lbsoundv * lbdt / lbcs;
  
  return 0;
  
}








int fGetEquilibriumF(double *feq, double *v, double rho)
{

  // calculate equilibrium distribution function for compressible fluid:
  // only suitable for square lattices, not suitable for incompressible fluid

  double modv = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  double uv;

  for(int i=0; i<lbsy.nq; i++)
    {
      uv = lbv[i*3] * v[0] + lbv[i*3+1] * v[1] + lbv[i*3+2] * v[2];
      feq[i] = rho * lbw[i] * (1 + 3.0 * uv + 4.5 * uv * uv - 1.5 * modv);
    }
  
  return 0;
}


int fGetEquilibriumFIncom(double *feq, double *v, double rho, double rho0)
{

  // calculate equilibrium distribution function for incompressible fluid:
  // only suitable for square lattices

  double modv = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  double uv;
  
  for(int i=0; i<lbsy.nq; i++)
    {
      uv = lbv[i*3] * v[0] + lbv[i*3+1] * v[1] + lbv[i*3+2] * v[2];
      feq[i] = lbw[i] * (rho + rho0 * (3.0 * uv + 4.5 * uv * uv - 1.5 * modv));
    }
  
  return 0;
}


int fGetEquilibriumC(double *feq, double *v, double rho)
{

  // calculate equilibrium distribution function for solute:
  // only suitable for square lattices

  double uv;
  
  for(int i=0; i<lbsy.nq; i++)
    {
      uv = lbv[i*3] * v[0] + lbv[i*3+1] * v[1] + lbv[i*3+2] * v[2];
      feq[i] = rho * lbw[i] * (1.0 + 3.0 * uv);
    }
  
  return 0;
}


int fGetEquilibriumT(double *feq, double *v, double rho)
{

  // calculate equilibrium distribution function for temperature:
  // only suitable for square lattices

  double uv;
  
  for(int i=0; i<lbsy.nq; i++)
    {
      uv = lbv[i*3] * v[0] + lbv[i*3+1] * v[1] + lbv[i*3+2] * v[2];
      feq[i] = rho * lbw[i] * (1.0 + 3.0 * uv);
    }
  
  return 0;
}










int fInitializeSystem()
{

  // initialize system

  double tdone, ttodo;
  double *pt2=&lbf[0];
  int *pt3=&lbphi[0];  
  for(long i=0; i<lbdm.touter; i++) {
    if(lbphi[i] == 10 || lbphi[i] == 11 || lbphi[i] == 12 || lbphi[i] == 13) { 
      for(int l=0; l<lbsitelength; l++){
	    *pt2 = 0;
	    pt2++;
      }
      pt3++;
    }
    else {
      ttodo =0.0;
      for(int l=0; l<lbsy.nf; l++)
	    ttodo += lbinip[l];
      for(int l=1; l<lbsy.nf; l++) {
	    tdone = lbinip[l-1]*(1+lbnoise*fRandom());
	    ttodo -= tdone;
        if(!incompress)
  	      fGetEquilibriumF(lbfeq, lbiniv, tdone);
        else
  	      fGetEquilibriumFIncom(lbfeq, lbiniv, tdone, lbincp[l-1]);
        for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      if(lbsy.nf > 1) {
        if(!incompress)
  	      fGetEquilibriumF(lbfeq, lbiniv, ttodo);
        else
  	      fGetEquilibriumFIncom(lbfeq, lbiniv, ttodo, lbincp[lbsy.nf-1]);
	    for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      else if(lbsy.nf == 1) {
        if(!incompress)
	      fGetEquilibriumF(lbfeq, lbiniv, lbinip[0]*(1+lbnoise*fRandom()));
        else
	      fGetEquilibriumFIncom(lbfeq, lbiniv, lbinip[0]*(1+lbnoise*fRandom()), lbincp[0]);
	    for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      ttodo =0.0;
      for(int l=0; l<lbsy.nc; l++)
	    ttodo += lbinic[l];
      for(int l=1; l<lbsy.nc; l++) {
	    tdone = lbinic[l-1]*(1+lbnoise*fRandom());
	    ttodo -= tdone;
	    fGetEquilibriumC(lbfeq, lbiniv, tdone);
	    for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      if(lbsy.nc > 1) {
	    fGetEquilibriumC(lbfeq, lbiniv, ttodo);
	    for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      else if(lbsy.nc == 1) {
	    fGetEquilibriumC(lbfeq, lbiniv, lbinic[0]*(1+lbnoise*fRandom()));
	    for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      if(lbsy.nt ==1) {
	    fGetEquilibriumT(lbfeq, lbiniv, lbinit);
	    for(int m=0; m<lbsy.nq; m++) {
	      *pt2 = lbfeq[m];
	      pt2 ++;
	    }
      }
      if(lbsy.pf ==1) {
	    *pt2 = 0;
	    pt2 ++;
      }
      pt3++;
    }
  }
  pt2=NULL;
  pt3=NULL;
  return 0;
}



int fPropagationSwap()
{
  
  // move particles to neighbouring grid points using swap algorithm
  // separate loops for two swaps: suitable for lbdm.bwid >= 0
                                     
#pragma omp parallel
  {
    long il, ilnext;
    int qdim = (lbsy.nf + lbsy.nc + lbsy.nt);
    int half = (lbsy.nq - 1) / 2;
    int nextx, nexty, nextz;
    int Xmax = lbdm.xouter;
    int Ymax = lbdm.youter;
    int Zmax = lbdm.zouter;

    #pragma omp for
    for (il=0; il<lbdm.touter; il++)
      for (int l=0; l<qdim; l++)
        for (int m=1; m<=half; m++)
          fSwapPair (lbf[il*lbsitelength + l*lbsy.nq + m], lbf[il*lbsitelength + l*lbsy.nq + m + half]);

    #pragma omp for collapse(3)
    for(int i=0; i<Xmax; i++)
      for(int j=0; j<Ymax; j++)
        for(int k=0; k<Zmax; k++) {
          il = (i * Ymax + j) * Zmax + k;
          for (int l=0; l<qdim; l++) {
            for (int m=1; m<=half; m++) {
              nextx = fCppMod(i + lbv[3*m],   Xmax);
              nexty = fCppMod(j + lbv[3*m+1], Ymax);
              nextz = fCppMod(k + lbv[3*m+2], Zmax);
              ilnext = (nextx * Ymax + nexty) * Zmax + nextz;
              fSwapPair (lbf[il*lbsitelength + l*lbsy.nq + m + half], lbf[ilnext*lbsitelength + l*lbsy.nq + m]);
            }
          }
        }
                                     
  }
  return 0;    
}



