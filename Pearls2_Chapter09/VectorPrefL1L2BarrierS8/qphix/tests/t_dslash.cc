// $Id: t_dslash.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $



#include "qdp.h"
#include "unittest.h"

#include "testDslashFull.h"

#include <omp.h>
#include <iostream>
#include <cstdio>

using namespace QDP;
using namespace std;

int nrow_in[4]={4,4,4,4};
int iters=1;
int By_user = -1;
int Bz_user = -1;
int NCores_user = -1;
int Sy_user = -1;
int Sz_user = -1;

// Hardwire these for now.
int PadXY_user = 1;
int PadXYZ_user = 0;
int MinCt_user = 1;

bool compress12=false;
Prec prec_user = FLOAT_PREC;

void printHelp() 
{ 
       cout << "t_dslash -x Lx -y Ly -z Lz -t Lt -i iters -by BY -bz BZ -c NCores  -sy SY -sz SZ  -pxy Pxy -pxyz Pxyz -minct MinCt -compress12 -prec Prec" << endl;
       cout << "   Lx is the lattice size in X" << endl;
       cout << "   Ly is the lattice size in Y" << endl;
       cout << "   Lz is the lattice size in Z" << endl;
       cout << "   iters is the number of iterations " << endl;
       cout << "   BY is the block size in Y " << endl;
       cout << "   BZ is the block size in Z " << endl;
       cout << "   NCores is the number of cores " << endl;
       cout << "   Sy is the no of simt threads in y" << endl;
       cout << "   Sz is the no of simt threads in z" << endl;
       cout << "   Pxy is the extra pad in the XY plane" << endl;
       cout << "   Pxyz is the extra pad in the XYZ plane" << endl;
       cout << "   MinCt is the MinCt in the blocking scheme" << endl;
       cout << "   Prec for precision " << endl;
}

void processArgs(int argc, char *argv[]) {
    int i=1;
    if (argc == 1) {
      printHelp();
    }

    while( i < argc)  {
        if( string(argv[i]).compare("-x") == 0 ) {
          nrow_in[0]=atoi(argv[i+1]);
          i+=2;
        }
        else if ( string(argv[i]).compare("-y") == 0 ) {
           nrow_in[1]=atoi(argv[i+1]);
           i+=2;
        }
        else if ( string(argv[i]).compare("-z") == 0) {
           nrow_in[2]=atoi(argv[i+1]);
           i+=2;
        }
        else if ( string(argv[i]).compare("-t") == 0) {
           nrow_in[3]=atoi(argv[i+1]);
           i+=2;
        }
        else if ( string(argv[i]).compare("-i") == 0) {
           iters=atoi(argv[i+1]);
           i+=2;
        }
        else if ( string(argv[i]).compare("-by") == 0 ) {
	  By_user=atoi(argv[i+1]);
            i+=2;
        }
        else if (string(argv[i]).compare("-bz") == 0 ) {
          Bz_user=atoi(argv[i+1]);
            i+=2;
        }
        else if ( string(argv[i]).compare("-c") == 0 ) {
	  NCores_user=atoi(argv[i+1]);
            i+=2;
        }
        else if ( string(argv[i]).compare("-sy") == 0 ) {
	  Sy_user=atoi(argv[i+1]);
            i+=2;
        }
        else if (string(argv[i]).compare("-sz") == 0 ) {
          Sz_user=atoi(argv[i+1]);
            i+=2;
        }
	else if ( string(argv[i]).compare("-pxy") == 0 ) {
	  PadXY_user=atoi(argv[i+1]);
            i+=2;
        }
        else if (string(argv[i]).compare("-pxyz") == 0 ) {
          PadXYZ_user=atoi(argv[i+1]);
            i+=2;
        }
	else if (string(argv[i]).compare("-minct") == 0 ) {
	  MinCt_user=atoi(argv[i+1]);
	  i+=2;
        }
	else if (string(argv[i]).compare("-compress12") == 0 ) {
	  compress12 =true;
	  i++;
        }
	else if (string(argv[i]).compare("-prec") == 0 ) { 
	  string user_arg(argv[i+1]);
	  if( user_arg.compare("h") == 0 ) { 
	    prec_user = HALF_PREC;
	  }
	  if( user_arg.compare("f") == 0 ) { 
	    prec_user = FLOAT_PREC;
	  }
	  if( user_arg.compare("d") == 0 ) { 
	    prec_user = DOUBLE_PREC;
	  }
	  i+=2 ;
	}
        else {
           i++;
        }

    }

    if( NCores_user < 0 ) { printHelp(); exit(1); }
    if( Sy_user < 0 ) { printHelp(); exit(1); }
    if( Sz_user < 0 ) { printHelp(); exit(1); }


}

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc,argv);

  omp_set_num_threads(NCores_user*Sy_user*Sz_user);

  // QDP++ gets set up here
  TestRunner  tests(&argc, &argv, nrow_in);

  const multi1d<int>& localLattSize = Layout::subgridLattSize();

  // If user doesn't specify block size, use local volume dimensions.
  if( By_user < 0 ) { 
    By_user = localLattSize[1];
  }
  if( Bz_user < 0 ) { 
    Bz_user = localLattSize[2];
  }

  for(int i=0; i < iters; i++) { 
#if defined(QPHIX_SCALAR_SOURCE)
    tests.addTest(new testDslashFull(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user, compress12, prec_user, 1), "testDslashFull_S1\n" );

#elif defined(QPHIX_QPX_SOURCE)
    tests.addTest(new testDslashFull(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user, compress12, prec_user, 4), "testDslashFull_S4\n" );
#else
    tests.addTest(new testDslashFull(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user, compress12, prec_user, 4), "testDslashFull_S4\n" );
    
    tests.addTest(new testDslashFull(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user, compress12, prec_user, 8), "testDslashFull_S8\n" ); 
    tests.addTest(new testDslashFull(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user, compress12, prec_user, 16), "testDslashFull_S16\n" );
#endif
  }

  tests.run();
  
  // Testjig is destroyed
  tests.summary();

}

