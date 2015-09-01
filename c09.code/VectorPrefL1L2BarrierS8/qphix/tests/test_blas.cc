#include <iostream>
#include <cstdio>

using namespace std;

#include <omp.h>
#include "testBlas.h"

int qmp_geom[4]={-1,-1,-1,-1};
int nrow_in[4]={4,4,4,4};
int iters=10000;
int NCores_user = -1;
int Sy_user = -1;
int Sz_user = -1;

// Hardwire these for now.
int PadXY_user = 1;
int PadXYZ_user = 0;


void printHelp() 
{ 
       cout << "test_blas -x Lx -y Ly -z Lz -t Lt -i iters -c C -sy SY -sz SZ -geom Px Py Pz Pt" << endl;
       cout << "   Lx is the lattice size in X" << endl;
       cout << "   Ly is the lattice size in Y" << endl;
       cout << "   Lz is the lattice size in Z" << endl;
       cout << "   iters is the number of iterations " << endl;
       cout << "   Cy is the no of cores in Y" << endl;
       cout << "   Cz is the no of cores in Z" << endl;
       cout << "   Ct is the no of cores in T" << endl;
       cout << "   N is the number of SIMT threads per core " << endl;
       cout << "   Px Py Pz Pt define a 4D grid of MPI tasks" << endl;
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
	else if (string(argv[i]).compare("-geom") == 0 ) {
	  qmp_geom[0] = atoi(argv[i+1]);
	  qmp_geom[1] = atoi(argv[i+2]);
	  qmp_geom[2] = atoi(argv[i+3]);
	  qmp_geom[3] = atoi(argv[i+4]);
	  i+=4;

        }
        else {
           i++;
        }

    }

    if( NCores_user < 0 ) { printHelp(); exit(1); }
    if( Sy_user < 0 ) { printHelp(); exit(1); }
    if( Sz_user < 0 ) { printHelp(); exit(1); }

    

}

using namespace QPhiX;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc,argv);
  // Max no of threads = NCores * Sy * Sz
  omp_set_num_threads(NCores_user*Sy_user*Sz_user);

#ifdef QPHIX_QMP_COMMS

  // Initialize QMP
  if( QMP_is_initialized() == QMP_FALSE ) { 
    QMP_thread_level_t prv;
    if( QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS ) { 
      QMP_error("Failed to initialize QMP\n");
      abort();
  
    }
  }
  
  // Declare the logical topology
  if ( QMP_declare_logical_topology(qmp_geom, 4)!= QMP_SUCCESS ) { 
    QMP_error("Failed to declare QMP Logical Topology\n");
    abort();
  }

#endif
 
  masterPrintf("Declared QMP Topology: %d %d %d %d\n", 
	       qmp_geom[0], qmp_geom[1], qmp_geom[2], qmp_geom[3]);
  masterPrintf("Launching TestCase\n");
  
  // Launch the test case. 
  testBlas test(NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, iters);
  
  test.run(nrow_in, qmp_geom);

  masterPrintf("Test Case Done\n");

#ifdef QPHIX_QMP_COMMS
  QMP_finalize_msg_passing();
#endif
}

