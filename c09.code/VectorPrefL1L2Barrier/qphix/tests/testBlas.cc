#include "unittest.h"
#include "testBlas.h"
#include <omp.h>
#include "qphix/qphix_config.h"
#include "qphix/dslash_def.h"
#include "qphix/wilson.h"
#include "qphix/blas.h"
#include "qphix/invcg.h"

#include <cstdlib>

using namespace std;
using namespace QPhiX;
using namespace Assertions;

#include "qphix/blas_new_c.h"


#ifdef  QPHIX_MIC_SOURCE
#define VECLEN 16
#endif

#ifdef QPHIX_AVX_SOURCE
#define VECLEN 8
#endif

#ifdef QPHIX_SCALAR_SOURCE
#define VECLEN 1
#endif

#ifdef QPHIX_QPX_SOURCE
#define VECLEN 4
#endif


void copy(  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* x1,
	    Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* x2, 
	    int N_blocks) 
{

#pragma omp parallel for
 for(int block=0; block < N_blocks; block++) {
    for(int spin=0; spin < 4; spin++) { 
      for(int col=0; col < 3; col++) { 
	for(int site=0; site < QPHIX_SOALEN; site++){ 
	  x2[block][col][spin][RE][site] = x1[block][col][spin][RE][site];
	  x2[block][col][spin][IM][site] = x1[block][col][spin][IM][site];
	}
      }
    }
 }
}

  
void resetSpinors( Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*  x1,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* x2,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* y1,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* y2,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* z1,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* z2,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* t1,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* t2,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* w1,
		   Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* w2,
		   const Geometry<float,VECLEN,QPHIX_SOALEN,true>& geom)
{
  double start, end;
  int Nt=geom.Nt();
  int Nz=geom.Nz();
  int Ny=geom.Ny();
  int N_blocks = geom.getPxyz()*Nt;

  masterPrintf("Filling spinors: ");
  start=omp_get_wtime();
  
  
  
  
  // Now we need to fill the arrays with drand48 numbers
  // We could speed this up with a parallel RNG
#pragma omp parallel for collapse(4)
  for(int t=0; t < Nt; t++) { 
    for(int z=0; z < Nz; z++) { 
      for(int y=0; y < Ny; y++) { 
	for(int vec=0; vec < geom.nVecs(); vec++) { 
	  for(int col=0; col < 3; col++) {
	    for( int spin=0; spin < 4; spin++) { 
	      int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
	      for( int reim=0; reim < 2; reim++) { 
		for(int site=0; site < QPHIX_SOALEN; site++) { 
		  x1[block][col][spin][reim][site]=drand48()-0.5;
		  y1[block][col][spin][reim][site]=drand48()-0.5;
		  z1[block][col][spin][reim][site]=drand48()-0.5;
		  t1[block][col][spin][reim][site]=drand48()-0.5;
		  w1[block][col][spin][reim][site]=drand48()-0.5;
		}	    
	      }
	    }
	  }
	}
      }
    }
  }
	
  // Copy fields including zeros

#pragma omp parallel for 
  for(int block=0; block < N_blocks; block++) {
    for(int col=0; col < 3; col++) { 
      for(int spin=0; spin < 4; spin++) { 
	for(int site=0; site < QPHIX_SOALEN; site++){ 

	  // x2, y2, z2, t2 are copies of x1,y1,z1,t1
	  x2[block][col][spin][RE][site] = x1[block][col][spin][RE][site];
	  x2[block][col][spin][IM][site] = x1[block][col][spin][IM][site];

	  y2[block][col][spin][RE][site] = y1[block][col][spin][RE][site];
	  y2[block][col][spin][IM][site] = y1[block][col][spin][IM][site];

	  z2[block][col][spin][RE][site] = z1[block][col][spin][RE][site];
	  z2[block][col][spin][IM][site] = z1[block][col][spin][IM][site];

	  t2[block][col][spin][RE][site] = t1[block][col][spin][RE][site];
	  t2[block][col][spin][IM][site] = t1[block][col][spin][IM][site]
;
	  w2[block][col][spin][RE][site] = w1[block][col][spin][RE][site];
	  w2[block][col][spin][IM][site] = w1[block][col][spin][IM][site];
	}
      }
    }
  }

  end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);
}


void
testBlas::run(const int lattSize[], const int qmp_geom[]) 
{
  // Work out local lattice size
  int subLattSize[4];
  for(int mu=0; mu < 4; mu++){ 
    subLattSize[mu]=lattSize[mu]/qmp_geom[mu];
  }

  // Work out the size of checkerboarded X-dimension
  int Nxh = subLattSize[0]/2;
  int Nx = subLattSize[0];
  int Ny = subLattSize[1];
  int Nz = subLattSize[2];
  int Nt = subLattSize[3];

  // Diagnostic information:

  masterPrintf("Global Lattice Size = ");
  for(int mu=0; mu < 4; mu++){ 
    masterPrintf(" %d", lattSize[mu]);
  }
  masterPrintf("\n");
  
  masterPrintf("Local Lattice Size = ");
  for(int mu=0; mu < 4; mu++){ 
    masterPrintf(" %d", subLattSize[mu]);
  }
  masterPrintf("\n");
  
  
  masterPrintf("Cores=%d  SMT Grid Sy=%d Sz=%d\n", NCores, Sy, Sz);
  masterPrintf("Padding Factors:  PadXY=%d PadXYZ=%d\n", PadXY, PadXYZ);
  
  // Dummy Block Sizes
  int By = Ny;
  int Bz = Nz;
  int N_simt = Sy*Sz;
  Geometry<float, VECLEN, QPHIX_SOALEN,true> geom(subLattSize,
				       Ny,Nz,
				       NCores,
				       Sy,Sz,
				       PadXY,PadXYZ,1);

  typedef Geometry<float,VECLEN, QPHIX_SOALEN,true>::FourSpinorBlock Spinor;
  
 

  masterPrintf("Initializing Geometry\n");

  // Allocate data for the spinors
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* x1=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* x2=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* y1=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* y2=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* z1=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* z2=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* t1=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* t2=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* w1=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();
  Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock* w2=(Geometry<float,VECLEN,QPHIX_SOALEN,true>::FourSpinorBlock*)geom.allocCBFourSpinor();

  int N_blocks =  geom.getPxyz()*geom.Nt();
  int len = N_blocks*QPHIX_SOALEN*4*3*2;
  float ar=0.6;
  resetSpinors(x1,x2,y1,y2,z1,z2,t1,t2, w1,w2, geom);


#if 1
 // =============== COPY ==================
  {
    
    // copy to x2,y2,z2,t2
    copy(z1,x1,N_blocks);
    copy(t1,y1, N_blocks);
    copy(x1,x2, N_blocks);
    copy(y1,y2, N_blocks);
    
    masterPrintf("Testing copySpinor: ");
    
    float *xf1 = (float *)x1;
    float *yf1 = (float *)y1;

    // BASIC COPY
    for(int s=0; s < len; s++) { 
      yf1[s] = xf1[s];
    }
    float *yf2 = (float *)y2;
    copySpinor<float,VECLEN>((float *)yf2, (float *)xf1,len);
    
    try { 
      for(int s=0; s < len; s++) { 
	assertion( toBool( fabs(yf2[s]-yf1[s]) < 1.0e-6 ) );
      }
      masterPrintf("OK\n");
    }
    catch(std::exception) { 
      masterPrintf("FAILED \n");
    }
    
  }
  
  // ADVANCED COPY
  {
     // Optimized (?) it
    for( int bt=1; bt <=N_simt; bt++) { 
      masterPrintf("Testing copySpinor with %d threads per core: ", bt);
      copy(z1,x2,N_blocks);
      copy(t1,y2, N_blocks);
      float *yf2 = (float *)y2;
      float *yf1 = (float *)y1;

#if 0
      copySpinor<float,VECLEN>((float *)yf2, (float *)x1, len, NCores, N_simt, bt);
#else
      copySpinor<float,VECLEN,QPHIX_SOALEN,true>(y2, x1, geom, bt);
#endif

      try { 

#pragma omp parallel for collapse(6)
	for(int t=0; t < Nt; t++) { 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++) { 
	      for(int vec=0; vec < geom.nVecs(); vec++) { 
		for(int col=0; col < 3; col++) {
		  for( int spin=0; spin < 4; spin++) { 
		    for( int reim=0; reim < 2; reim++) { 
		      for(int s=0; s < QPHIX_SOALEN; s++) { 
			int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
			assertion( toBool( fabs(y2[block][col][spin][reim][s]
						-y1[block][col][spin][reim][s]) < 1.0e-6 ) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      
        masterPrintf("OK\n");
      }
      catch(std::exception) { 
	masterPrintf("FAILED \n");
      }
    }
    masterPrintf("Timing copySpinor\n");
    int titers=1000;
   
    double tstart = omp_get_wtime();
    for(int i=0; i < titers; i++) { 
      copySpinor<float,VECLEN>((float *)y1, (float *)x1,len);
    }
    double tstop = omp_get_wtime();
    double ttime=tstop-tstart;
    
    double BW = (2.0*len*sizeof(float)*(double)titers)/ttime;
    
    masterPrintf("copySpinor: len=%d time=%g iters=%d BW=%g\n", 
		 len, ttime, titers, BW*1.0e-9);
    
    for(int bt = 1; bt <=N_simt ; bt++) { 
      masterPrintf("Timing copySpinor: blasThreads = %d\n", bt);
      
      int titers=1000;
      double tstart = omp_get_wtime();
      
      for(int i=0; i < titers; i++) { 
#if 0
	copySpinor<float,VECLEN>((float *)y2, (float *)x1,len,NCores, N_simt, bt);
#else 
	copySpinor<float,VECLEN, QPHIX_SOALEN,true>(y2, x1, geom, bt);
#endif
      }
      double tstop = omp_get_wtime();
      double ttime=tstop-tstart;
#if 0
      double real_len = len*sizeof(float);
#else
      double real_len=(double)(Nxh*Ny*Nz*Nt*4*3*2*sizeof(float));
#endif

      double BW = (2.0*real_len*(double)titers)/ttime;
      masterPrintf("copySpinor: blasThreads=%d len=%g time=%g iters=%d BW=%g\n", 
		   bt, real_len, ttime, titers, BW*1.0e-9);
    }
  }
#endif
  //  ================== END COPY =======================

#if 1
 //   ====================   AYPX ========================
  {
    
    // copy to x2,y2,z2,t2
    copy(z1,x1,N_blocks);
    copy(t1,y1, N_blocks);
    copy(x1,x2, N_blocks);
    copy(y1,y2, N_blocks);
    masterPrintf("Testing aypx: ");
    
    float *xf1 = (float *)x1;
    float *yf1 = (float *)y1;
    
    // BASIC AYPX
    for(int s=0; s < len; s++) { 
      yf1[s] = ar*yf1[s] + xf1[s];
    }
    float *xf2 = (float *)x2;
    float *yf2 = (float *)y2;
    
    aypx<float,VECLEN>(ar, (float *)xf2, (float *)yf2, len);
    
    try { 
      for(int s=0; s < len; s++) { 
	assertion( toBool( fabs(yf2[s]-yf1[s]) < 1.0e-6 ) );
      }
      masterPrintf("OK\n");
    }
    catch(std::exception) { 
      masterPrintf("FAILED \n");
    }
    
    
   
    // Optimized (?) it
    for( int bt=1; bt <=N_simt; bt++) { 
      masterPrintf("Testing aypx with %d threads per core: ", bt);
      copy(z1,x2,N_blocks);
      copy(t1,y2, N_blocks);
      float *xf2 = (float *)x2;
      float *yf2 = (float *)y2;
      aypx<float,VECLEN,QPHIX_SOALEN,true>(ar, x2, y2, geom, bt);
      try { 
#pragma omp parallel for collapse(6)
	for(int t=0; t < Nt; t++) { 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++) { 
	      for(int vec=0; vec < geom.nVecs(); vec++) { 
		for(int col=0; col < 3; col++) {
		  for( int spin=0; spin < 4; spin++) { 
		    for( int reim=0; reim < 2; reim++) { 
		      for(int s=0; s < QPHIX_SOALEN; s++) { 
			int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
			assertion( toBool( fabs(y2[block][col][spin][reim][s]
						-y1[block][col][spin][reim][s]) < 1.0e-6 ) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	masterPrintf("OK\n");
      }
      catch(std::exception) { 
	masterPrintf("FAILED \n");
      }
    }
    
    masterPrintf("Timing aypx\n");
    int titers=1000;
    
    double tstart = omp_get_wtime();
    for(int i=0; i < titers; i++) { 
      aypx<float,VECLEN>(ar, (float *)xf1, (float *)yf1,len);
    }
    double tstop = omp_get_wtime();
    double ttime=tstop-tstart;
    double BW = (3.0*len*sizeof(float)*(double)titers)/ttime;
    
    masterPrintf("aypx: len=%d time=%g iters=%d BW=%g\n", 
		 len, ttime, titers, BW*1.0e-9);

    for(int bt = 1; bt <=N_simt ; bt++) { 
      masterPrintf("Timing aypx: blasThreads = %d\n", bt);
      
      int titers=1000;
      double tstart = omp_get_wtime();
      
      for(int i=0; i < titers; i++) { 
	aypx<float,VECLEN,QPHIX_SOALEN,true>(ar, x2, y2, geom, bt);
      }
      double tstop = omp_get_wtime();
      double ttime=tstop-tstart;
      double real_len=(double)(Nxh*Ny*Nz*Nt*4*3*2*sizeof(float));
      double BW = (3.0*real_len*(double)titers)/ttime;
      masterPrintf("aypx: blasThreads=%d len=%g time=%g iters=%d BW=%g\n", 
	       bt, real_len, ttime, titers, BW*1.0e-9);
    }
  }
    
#endif

//   ====================  END  AYPX ========================


#if 0

 //  ============ XMYNORM2SPINOR =========================
  {
    masterPrintf("Testing xmyNorm2Spinor: ");
    copy(z1,x1, N_blocks);
    copy(t1,y1, N_blocks);

    double norm=0;
    float *x1f = (float *)x1;
    float *y1f = (float *)y1;
    float *w1f = (float *)w1;
    
    for(int t=0; t < Nt; t++) { 
      for(int z=0; z < Nz; z++) { 
	for(int y=0; y < Ny; y++) { 
	  for(int vec=0; vec < geom.nVecs(); vec++) { 
	    for(int col=0; col < 3; col++) {
	      for( int spin=0; spin < 4; spin++) { 
		for( int reim=0; reim < 2; reim++) { 
		  for(int s=0; s < QPHIX_SOALEN; s++) { 
		    int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
		    w1[block][col][spin][reim][s] = x1[block][col][spin][reim][s] - y1[block][col][spin][reim][s];
		    double w1sd =  w1[block][col][spin][reim][s];
		    norm += w1sd*w1sd;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
   
    CommsUtils::sumDouble(&norm);

    copy(z1, x2, N_blocks);
    copy(t1, y2, N_blocks);

    double norm2=0;
    xmyNorm2Spinor<float,VECLEN>((float *)w2,(float *)x2, (float *)y2, norm2, len);


    float *w2f = (float *)w2;
    
    try { 
      for(int s=0; s < len; s++) { 
	assertion( toBool( fabs( w2f[s]-w1f[s] ) < 1.0e-6 ) );
      }
      assertion( toBool( fabs(norm2 - norm)/ norm < 1.0e-12 ));
      masterPrintf("OK\n"); 
      masterPrintf(" norm = %4.20e \n norm2= %4.20e\n relative diff =%4.16e\n", norm, norm2, fabs(norm2-norm)/norm);
      
    }
    catch(std::exception) { 
      masterPrintf("FAILED \n");
    }
    
    for(int bt=1; bt <=N_simt; bt++) {
      
      int vec_successes = 0;
      int sum_successes = 0;
      int vec_failures = 0;
      int sum_failures = 0;
      
      masterPrintf("Testing xmyNorm2Spinor with %d threads: ", bt);
      
      
      copy(z1, x2, N_blocks);
      copy(t1, y2, N_blocks);
      norm2=0;
      xmyNorm2Spinor<float,VECLEN,QPHIX_SOALEN,true>(w2,x2,y2, norm2, geom, bt);

      try { 
#pragma omp parallel for collapse(6)
	for(int t=0; t < Nt; t++) { 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++) { 
	      for(int vec=0; vec < geom.nVecs(); vec++) { 
		for(int col=0; col < 3; col++) {
		  for( int spin=0; spin < 4; spin++) { 
		    for( int reim=0; reim < 2; reim++) { 
		      for(int s=0; s < QPHIX_SOALEN; s++) { 
			int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
			assertion( toBool( fabs(w2[block][col][spin][reim][s]
						-w1[block][col][spin][reim][s]) < 1.0e-6 ) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	masterPrintf("OK\n");
      }
      catch( std::exception ) { 
	masterPrintf("FAILED\n");
      }
      
      try {  
	assertion( toBool( fabs(norm2 - norm)/ norm < 1.0e-12 ));
	masterPrintf ("Sum is OK\n");
      }
      catch(std::exception) { 
	masterPrintf("FAIL:  norm = %4.20e \n norm2= %4.20e\n relative diff =%4.16e\n", norm, norm2, fabs(norm2-norm)/norm); 
      }
    }
    
    masterPrintf("Timing xmyNorm2Spinor\n");
    
    int titers=1000;
    double tstart = omp_get_wtime();
    for(int i=0; i < titers; i++) { 
      xmyNorm2Spinor<float,VECLEN>((float *)w2, 
				   (float *)x2, 
				   (float *)y2, 
				   norm, len);
    }
    double tstop = omp_get_wtime();
    double ttime=tstop-tstart;
    double BW = (3.0*len*sizeof(float)*(double)titers)/ttime;
    
    masterPrintf("xmyNorm2Spinor: len=%d time=%g iters=%d BW=%g\n", 
	     len, ttime, titers, BW*1.0e-9);
        
    for( int bt =1 ; bt <= N_simt; bt++) { 
      masterPrintf("Timing xmyNorm2Spinor: blasThreads=%d \n", bt);
      
      int titers=1000;
      double tstart = omp_get_wtime();
      
      for(int i=0; i < titers; i++) { 
	xmyNorm2Spinor<float,VECLEN,QPHIX_SOALEN,true>(w2, 
					    x2, 
					    y2, 
					    norm, geom, bt );
      }
      double tstop = omp_get_wtime();
      double ttime=tstop-tstart;
      double real_len=(double)(Nxh*Ny*Nz*Nt*4*3*2*sizeof(float));
      double BW = (3.0*real_len*(double)titers)/ttime;
      
	masterPrintf("xmyNorm2Spinor: blasThreads = %d len=%d time=%g iters=%d BW=%g\n", 
	       bt, len, ttime, titers, BW*1.0e-9);
    }
  }
  
#endif
 
#if 1
// RMAMMPNORM2RXPAP
  {
    double norm;
    double norm2;
    resetSpinors(x1,x2,y1,y2,z1,z2,t1,t2, w1,w2, geom);

    copy(w1,x1,N_blocks);
    copy(w2,z1, N_blocks);

    float *x1f = (float *)x1;
    float *p1f = (float *)y1;
    float *r1f  = (float *)z1;
    float *mmp1f = (float *)t1;

    masterPrintf("Testing rmammpNorm2rxpap\n");
    
    float ar=0.6;
    norm = 0;
    for(int t=0; t < Nt; t++) { 
      for(int z=0; z < Nz; z++) { 
	for(int y=0; y < Ny; y++) { 
	  for(int vec=0; vec < geom.nVecs(); vec++) { 
	    for(int col=0; col < 3; col++) {
	      for( int spin=0; spin < 4; spin++) { 
		for( int reim=0; reim < 2; reim++) { 
		  for(int s=0; s < QPHIX_SOALEN; s++) { 
		    int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
		    x1[block][col][spin][reim][s] += ar*y1[block][col][spin][reim][s];
		    z1[block][col][spin][reim][s] -= ar*t1[block][col][spin][reim][s];
		    double tmpd=(double)z1[block][col][spin][reim][s];
		    norm += tmpd * tmpd;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    CommsUtils::sumDouble(&norm);

    copy(w1,x2, N_blocks);
    copy(w2,z2, N_blocks);

    float *x2f = (float *)x2;
    float *p2f = (float *)p1f; // Unchanged
    float *r2f  = (float *)z2;
    float *mmp2f = (float *)mmp1f; // Unchanged
    norm2 = 0;
    rmammpNorm2rxpap<float,VECLEN>(r2f, ar, mmp2f, norm2, x2f, p2f, len);
    try { 
      for(int t=0; t < Nt; t++) { 
	for(int z=0; z < Nz; z++) { 
	  for(int y=0; y < Ny; y++) { 
	    for(int vec=0; vec < geom.nVecs(); vec++) { 
	      for(int col=0; col < 3; col++) {
		for( int spin=0; spin < 4; spin++) { 
		  for( int reim=0; reim < 2; reim++) { 
		    for(int s=0; s < QPHIX_SOALEN; s++) { 
		      int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
		      assertion( toBool( fabs(x1[block][col][spin][reim][s]-x1[block][col][spin][reim][s]) < 1.0e-6 ) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      masterPrintf("Vector x is OK\n");
    }
    catch( std::exception ) { 
      masterPrintf("Vector x is broken\n");
    }
  
    try { 
      for(int t=0; t < Nt; t++) { 
	for(int z=0; z < Nz; z++) { 
	  for(int y=0; y < Ny; y++) { 
	    for(int vec=0; vec < geom.nVecs(); vec++) { 
	      for(int col=0; col < 3; col++) {
		for( int spin=0; spin < 4; spin++) { 
		  for( int reim=0; reim < 2; reim++) { 
		    for(int s=0; s < QPHIX_SOALEN; s++) { 
		      int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
		      assertion( toBool( fabs(z2[block][col][spin][reim][s]-z1[block][col][spin][reim][s]) < 1.0e-6 ) );
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      masterPrintf("Vector z is OK\n");
    }
    catch( std::exception ) { 
      masterPrintf("Vector z is broken\n");
    }
  
    
    try { 
      assertion( toBool( fabs(norm2 - norm)/ norm < 2.0e-12 ));
    }
    catch(std::exception) { 
      masterPrintf("FAIL:  norm = %4.20e \n norm2= %4.20e\n relative diff =%4.16e\n", norm, norm2, fabs(norm2-norm)/norm); 
    }
  
    // Threaded versions
    for(int bt=1; bt <= N_simt; bt++) { 
      masterPrintf("Testing rmammpNorm2rxpap with %d threads\n", bt);
      copy(w1,x2, N_blocks);
      copy(w2,z2, N_blocks);

      float *x2f = (float *)x2;
      float *p2f = (float *)p1f; // Unchanged
      float *r2f  = (float *)z2;
      float *mmp2f = (float *)mmp1f; // Unchanged
      norm2 = 0;
      rmammpNorm2rxpap<float,VECLEN,QPHIX_SOALEN,true>(z2, ar, t1, norm2, x2, y1, geom, bt);

      try { 
#pragma omp parallel for collapse(6)
	for(int t=0; t < Nt; t++) { 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++) { 
	      for(int vec=0; vec < geom.nVecs(); vec++) { 
		for(int col=0; col < 3; col++) {
		  for( int spin=0; spin < 4; spin++) { 
		    for( int reim=0; reim < 2; reim++) { 
		      for(int s=0; s < QPHIX_SOALEN; s++) { 
			int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
			assertion( toBool( fabs(x2[block][col][spin][reim][s]
						-x1[block][col][spin][reim][s]) < 1.0e-6 ) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	masterPrintf("Vector X is OK\n");

      }
      catch( std::exception ) { 
	masterPrintf("FAIL x is bad\n");
      }

      try { 
#pragma omp parallel for collapse(6)
	for(int t=0; t < Nt; t++) { 
	  for(int z=0; z < Nz; z++) { 
	    for(int y=0; y < Ny; y++) { 
	      for(int vec=0; vec < geom.nVecs(); vec++) { 
		for(int col=0; col < 3; col++) {
		  for( int spin=0; spin < 4; spin++) { 
		    for( int reim=0; reim < 2; reim++) { 
		      for(int s=0; s < QPHIX_SOALEN; s++) { 
			int block = t*geom.getPxyz() + z*geom.getPxy() + y*geom.nVecs()+vec;
			assertion( toBool( fabs(z2[block][col][spin][reim][s]
						-z1[block][col][spin][reim][s]) < 1.0e-6 ) );
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	masterPrintf("Vector R is OK\n");
      }
      catch( std::exception ) { 
	masterPrintf("FAIL: r is bad \n");
      }
      
      try { 
	assertion( toBool( fabs(norm2 - norm)/ norm < 2.0e-12 ));
	masterPrintf("Sum is OK\n");
      }
      catch(std::exception) { 
	masterPrintf("FAIL:  norm = %4.20e \n norm2= %4.20e\n relative diff =%4.16e\n", norm, norm2, fabs(norm2-norm)/norm); 
      }
    }
  
    int titers=1000;
    double cp;
    double tstart = omp_get_wtime();
    
    for(int i=0; i < titers; i++) { 
      rmammpNorm2rxpap<float,VECLEN>(x2f, ar, p2f, norm2, r2f, mmp2f, len);
    }
    double tstop = omp_get_wtime();
    double ttime=tstop-tstart;
    double BW = (6.0*len*sizeof(float)*(double)titers)/ttime;
    
    masterPrintf("rmammpNorm2rxpap: len=%d time=%g iters=%d BW=%g\n", 
		 len, ttime, titers, BW*1.0e-9);

  

    for(int bt =1 ; bt <= N_simt; bt++) { 
      masterPrintf("Timing rmammpNorm2rxpap: blasThreads=%d\n",bt);
            
      tstart = omp_get_wtime();
      
      for(int i=0; i < titers; i++) { 
	rmammpNorm2rxpap<float,VECLEN,QPHIX_SOALEN,true>(z2, ar, t1, norm2, x2, y1, geom, bt);
      }

      double tstop = omp_get_wtime();
      double ttime=tstop-tstart;
      double BW = (6.0*len*sizeof(float)*(double)titers)/ttime;
      masterPrintf("rmammpNorm2rxpap: blasThreads=%d len=%d time=%g iters=%d BW=%g\n", 
	       bt, len, ttime, titers, BW*1.0e-9);
      
    }
  }
#endif

#if 1
  // NORM2SPINOR
  {
    masterPrintf("Testing norm2 spinor\n");
   
    // Hand roll it
    double norm = 0;
    double norm2 = 0;
    float *y1f = (float *)y1;

    for(int i=0; i < len; i++) { 
      norm += ( (double)y1f[i] * (double)y1f[i] );
    }
    CommsUtils::sumDouble(&norm);


    // Unoptimized version
    copy(y1,y2,N_blocks);
    norm2 = norm2Spinor<float,VECLEN>((float *)y2, len);

    try { 
      assertion( toBool( fabs(norm2 - norm)/ norm < 2.0e-12 ));
      masterPrintf("Sum OK\n");
    }
    catch(std::exception) { 
      masterPrintf("FAIL:  norm = %4.20e \n norm2= %4.20e\n relative diff =%4.16e\n", norm, norm2, fabs(norm2-norm)/norm); 
    }


    for(int bt =1 ; bt <=N_simt ; bt++) { 
      masterPrintf("ing norm2 with %d threads per core \n", bt);
      norm2Spinor<float,VECLEN,QPHIX_SOALEN,true>(norm2,y2, geom, bt);
      try { 
	assertion( toBool( fabs(norm2 - norm)/ norm < 2.0e-12 ));
	masterPrintf("Sum OK\n");
      }
      catch(std::exception) { 
	masterPrintf("FAIL:  norm = %4.20e \n norm2= %4.20e\n relative diff =%4.16e\n", norm, norm2, fabs(norm2-norm)/norm); 
      }
    }


    masterPrintf("Timing norm2Spinor\n");
    int titers=4000;
    double n2res;
    double tstart = omp_get_wtime();

    for(int i=0; i < titers; i++) { 
      n2res = norm2Spinor<float,VECLEN>((float *)y2, len);
    }
    double tstop = omp_get_wtime();
    double ttime=tstop-tstart;
    double BW = (1.0*len*sizeof(float)*(double)titers)/ttime;
    
    masterPrintf("norm2Spinor: len=%d time=%g iters=%d BW=%g\n", 
	     len, ttime, titers, BW*1.0e-9);

  
    for(int bt=1; bt <= N_simt; bt++) { 
      masterPrintf("Timing norm2Spinor: blasThreads=%d\n", bt);
      int titers=4000;
      double n2res;
      double tstart = omp_get_wtime();
      
      for(int i=0; i < titers; i++) { 
	norm2Spinor<float,VECLEN,QPHIX_SOALEN,true>(n2res,y2, geom, bt);
      }
      double tstop = omp_get_wtime();
      double ttime=tstop-tstart;
      double BW = (1.0*len*sizeof(float)*(double)titers)/ttime;
      
      masterPrintf("norm2Spinor: blasThreads = %d len=%d time=%g iters=%d BW=%g\n", 
		   bt, len, ttime, titers, BW*1.0e-9);
      
    }
  }
#endif

  masterPrintf("Cleaning up\n");

 
  ALIGNED_FREE(x1);
  ALIGNED_FREE(x2);

  ALIGNED_FREE(y1);
  ALIGNED_FREE(y2);

  ALIGNED_FREE(z1);
  ALIGNED_FREE(z2);

  ALIGNED_FREE(t1);
  ALIGNED_FREE(t2);


  ALIGNED_FREE(w1);
  ALIGNED_FREE(w2);
						   

}
