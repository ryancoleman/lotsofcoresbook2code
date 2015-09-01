#include "unittest.h"
#include "timeMatEqMatHermDouble.h"

static double N_SECS=10;

extern void*  alloc_cache_aligned_3mat(unsigned num_sites, REAL64** x, REAL64 **y, REAL64** z);


// M=M*H kernel (M \in SU3)
void
timeMeqMH_QDP::run(void) 
{

  LatticeColorMatrix x;
  LatticeColorMatrix y;
  LatticeColorMatrix z;
  gaussian(x);
  gaussian(y);

  QDPIO::cout << std::endl << "Timing  QDP++ MH Kernel " <<std::endl;

  StopWatch swatch;
  double n_secs = N_SECS;
  int iters=1;
  double time=0;
  QDPIO::cout << "\t Calibrating for " << n_secs << " seconds " << std::endl;
  do {
    swatch.reset();
    swatch.start();
    
    for(int i=0; i < iters; i++) { 
      z = x*adj(y);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();

    // Average time over nodes
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    if (time < n_secs) {
      iters *=2;
      QDPIO::cout << "." << std::flush;
    }
  }
  while ( time < (double)n_secs );
      
  QDPIO::cout << std::endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    z=x*adj(y);
  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();
  time /= (double)iters;

  double flops=(double)(198*Layout::vol());
  double perf=(flops/time)/(double)(1024*1024);
  QDPIO::cout << "QDP++ MH Kernel: " << perf << " Mflops" << std::endl;

}

// Optimized M=MH kernel (M \in SU3)
void
timeMeqMH::run(void) 
{

  LatticeColorMatrixD3 x;
  LatticeColorMatrixD3 y;

  REAL64* xptr;
  REAL64* yptr;
  REAL64* zptr;

  REAL64* top;

  int n_mat = (all.end() - all.start() + 1);
  top = (REAL64 *)alloc_cache_aligned_3mat(n_mat, &xptr, &yptr, &zptr);

  gaussian(x);
  gaussian(y);

  /* Copy x into x_ptr, y into y_ptr */
  REAL64 *f_x = xptr;
  REAL64 *f_y = yptr;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) {
      for(int col2=0; col2 < 3; col2++) { 
	*f_x = x.elem(site).elem().elem(col1,col2).real();
	*f_y = y.elem(site).elem().elem(col1,col2).real();
	f_x++;
	f_y++;
	*f_x = x.elem(site).elem().elem(col1,col2).imag();
	*f_y = y.elem(site).elem().elem(col1,col2).imag();
	f_x++;
	f_y++;
      }
    }
  }

  QDPIO::cout << std::endl << "Timing SSE D  M=MH  Kernel " <<std::endl;

  StopWatch swatch;
  double n_secs = N_SECS;
  int iters=1;
  double time=0;
  QDPIO::cout << "\t Calibrating for " << n_secs << " seconds " << std::endl;
  do {
    swatch.reset();
    swatch.start();
    
    for(int i=0; i < iters; i++) { 
      ssed_m_eq_mh(zptr, xptr, yptr, n_mat);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();

    // Average time over nodes
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    if (time < n_secs) {
      iters *=2;
      QDPIO::cout << "." << std::flush;
    }
  }
  while ( time < (double)n_secs );
      
  QDPIO::cout << std::endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    ssed_m_eq_mh(zptr, xptr, yptr, n_mat);
  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();
  time /= (double)iters;

  double flops=(double)(198*Layout::vol());
  double perf=(flops/time)/(double)(1024*1024);
  QDPIO::cout << "SSED MH Kernel: " << perf << " Mflops" << std::endl;

  free(top);
}

// M+=M*H kernel (M \in SU3)
void
timeMPeqaMH_QDP::run(void) 
{

  LatticeColorMatrix x;
  LatticeColorMatrix y;
  LatticeColorMatrix z;
  gaussian(x);
  gaussian(y);
  gaussian(z);
  Real a(-1.0);


  QDPIO::cout << std::endl << "Timing  QDP++ M+=MH Kernel " <<std::endl;

  StopWatch swatch;
  double n_secs = N_SECS;
  int iters=1;
  double time=0;
  QDPIO::cout << "\t Calibrating for " << n_secs << " seconds " << std::endl;
  do {
    swatch.reset();
    swatch.start();
    
    for(int i=0; i < iters; i++) { 
      z += x*adj(y);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();

    // Average time over nodes
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    if (time < n_secs) {
      iters *=2;
      QDPIO::cout << "." << std::flush;
    }
  }
  while ( time < (double)n_secs );
      
  QDPIO::cout << std::endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    z+=x*adj(y);
  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();
  time /= (double)iters;

  double flops=(double)(216*Layout::vol());
  double perf=(flops/time)/(double)(1024*1024);
  QDPIO::cout << "QDP++ M+=MH Kernel: " << perf << " Mflops" << std::endl;

}

// Optimized M += aM*H kernel
void
timeMPeqaMH::run(void) 
{

  LatticeColorMatrixD3 x;
  LatticeColorMatrixD3 y;
  Double a(+1.0);

  REAL64* xptr;
  REAL64* yptr;
  REAL64* zptr;  
  REAL64* top;

  int n_mat = (all.end() - all.start() + 1);
  top = (REAL64 *)alloc_cache_aligned_3mat(n_mat, &xptr, &yptr, &zptr);

  gaussian(x);
  gaussian(y);

  /* Copy x into x_ptr, y into y_ptr */
  REAL64 *f_x = xptr;
  REAL64 *f_y = yptr;

  for(int site=all.start(); site <= all.end(); site++) { 
    for(int col1=0; col1 < 3; col1++) {
      for(int col2=0; col2 < 3; col2++) { 
	*f_x = x.elem(site).elem().elem(col1,col2).real();
	*f_y = y.elem(site).elem().elem(col1,col2).real();
	f_x++;
	f_y++;
	*f_x = x.elem(site).elem().elem(col1,col2).imag();
	*f_y = y.elem(site).elem().elem(col1,col2).imag();
	f_x++;
	f_y++;
      }
    }
  }

  REAL64* aptr = &(a.elem().elem().elem().elem());

  QDPIO::cout << std::endl << "Timing SSE D  M+=MH  Kernel " <<std::endl;

  StopWatch swatch;
  double n_secs = N_SECS;
  int iters=1;
  double time=0;
  QDPIO::cout << "\t Calibrating for " << n_secs << " seconds " << std::endl;
  do {
    swatch.reset();
    swatch.start();
    
    for(int i=0; i < iters; i++) { 
      ssed_m_peq_amh(zptr, aptr, xptr, yptr, n_mat);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();

    // Average time over nodes
    QDPInternal::globalSum(time);
    time /= (double)Layout::numNodes();

    if (time < n_secs) {
      iters *=2;
      QDPIO::cout << "." << std::flush;
    }
  }
  while ( time < (double)n_secs );
      
  QDPIO::cout << std::endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << std::endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
    ssed_m_peq_amh(zptr, aptr, xptr, yptr, n_mat);
  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  QDPInternal::globalSum(time);
  time /= (double)Layout::numNodes();
  time /= (double)iters;

  double flops=(double)(216*Layout::vol());
  double perf=(flops/time)/(double)(1024*1024);
  QDPIO::cout << "SSED M+=MH Kernel: " << perf << " Mflops" << std::endl;

  free(top);
}
