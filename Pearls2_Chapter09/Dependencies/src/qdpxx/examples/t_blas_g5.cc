// $Id: t_blas_g5.cc,v 1.3 2005-03-18 11:55:29 bjoo Exp $

#include <iostream>
#include <cstdio>

// Ensure that the template matches don't get used
#define QDP_SCALARSITE_SSE_BLAS_G5_H
#define QDP_SCALARSITE_GENERIC_BLAS_G5_H
#include "qdp.h"

#ifndef QDP_USE_SSE
#warning "Non sse target. Using Generics"
#include "scalarsite_generic/generic_blas_vaxpy3_g5.h"
#include "scalarsite_generic/generic_blas_vaypx3_g5.h"
#include "scalarsite_generic/generic_blas_vadd3_g5.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
#include "scalarsite_generic/generic_blas_vaxpby3_g5.h"
#else

#if BASE_PRECISION == 32

#if defined (__GNUC__)
#warning "SSE 32bit target with GCC. Using Some SSE assembler"
#include "scalarsite_sse/sse_blas_vaxpy3_g5.h"
#include "scalarsite_sse/sse_blas_vaypx3_g5.h"
#include "scalarsite_sse/sse_blas_vadd3_g5.h"
#include "scalarsite_sse/sse_blas_vscal3_g5.h"
#include "scalarsite_sse/sse_blas_vaxpby3_g5.h"
#else
#warning "SSE 32bit target but not GCC. Using Generics"
#include "scalarsite_generic/generic_blas_vaxpy3_g5.h"
#include "scalarsite_generic/generic_blas_vaypx3_g5.h"
#include "scalarsite_generic/generic_blas_vadd3_g5.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
#include "scalarsite_generic/generic_blas_vaxpby3_g5.h"
#endif // GNUC

#else
#warning "SSE 64 BIT target. Using Generics"
#include "scalarsite_generic/generic_blas_vaxpy3_g5.h"
#include "scalarsite_generic/generic_blas_vaypx3_g5.h"
#include "scalarsite_generic/generic_blas_vadd3_g5.h"
#include "scalarsite_generic/generic_blas_vscal_g5.h"
#include "scalarsite_generic/generic_blas_vaxpby3_g5.h"
#endif // BASE_PRECISION

#endif // QDP_USE_SSE
 
using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {8,8,8,8};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  
  Real a=Real(1.5);
  LatticeFermion qx;
  LatticeFermion qy;
  LatticeFermion qz;
  LatticeFermion qz2;
 
  gaussian(qx);
  qy=qx;

  REAL* Out;
  REAL* scalep;
  REAL* InScale;
  REAL* Add;
  int n_4vec;
  Double norm_diff;

  // ax + Py  -- 3Nc Ns flops
  qz  = a*qx + chiralProjectPlus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  axpyz_g5ProjPlus(Out, scalep, InScale, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax + P_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = a*qx + chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = a*qx + chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axpyz_g5ProjPlus(Out, scalep, InScale, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axpyz_g5ProjPlus(Out, scalep, InScale, Add,n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  // ax + P{-}y   3 Nc Ns flops
  qz  = a*qx + chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  axpyz_g5ProjMinus(Out, scalep, InScale, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax + P_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = a*qx + chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = a*qx + chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axpyz_g5ProjMinus(Out, scalep, InScale, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axpyz_g5ProjMinus(Out, scalep, InScale, Add,n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // ax - P{+} y  3Nc Ns flops
  
  qz  = a*qx - chiralProjectPlus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  axmyz_g5ProjPlus(Out, scalep, InScale, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax - P_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = a*qx - chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = a*qx - chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axmyz_g5ProjPlus(Out, scalep, InScale, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axmyz_g5ProjPlus(Out, scalep, InScale, Add,n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // ax - P{-}y   3 Nc Ns flops
  qz  = a*qx - chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  axmyz_g5ProjMinus(Out, scalep, InScale, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax - P_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = a*qx - chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = a*qx - chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axmyz_g5ProjMinus(Out, scalep, InScale, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axmyz_g5ProjMinus(Out, scalep, InScale, Add,n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  // x + a P{+}y  2Nc Ns flops
  qz  = qx + a* chiralProjectPlus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  xpayz_g5ProjPlus(Out, scalep, Add, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x + aP_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx + a*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx + a*chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	xpayz_g5ProjPlus(Out, scalep, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      xpayz_g5ProjPlus(Out, scalep, Add, InScale,n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  // x + a P{-} y,       2Nc Ns flops
  qz  = qx + a*chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  xpayz_g5ProjMinus(Out, scalep, Add, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x + aP_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx + a*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx + a*chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	xpayz_g5ProjMinus(Out, scalep, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      xpayz_g5ProjMinus(Out, scalep, Add, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // x - a P{+} y, 2Nc flops
  qz  = qx - a*chiralProjectPlus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  xmayz_g5ProjPlus(Out, scalep,Add,  InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x -a P_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx - a*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx - a*chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	xmayz_g5ProjPlus(Out, scalep, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      xmayz_g5ProjPlus(Out, scalep, Add, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // x - a P{-}y, 2 Nc Ns flops
  qz  = qx - a*chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  xmayz_g5ProjMinus(Out, scalep, Add, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x - a P_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx - a*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx - a*chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	xmayz_g5ProjMinus(Out, scalep, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      xmayz_g5ProjMinus(Out, scalep, Add, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  // x + P{y}  Nc Ns flops
  qz  = qx +  chiralProjectPlus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  add_g5ProjPlus(Out, Add, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x + P_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx + chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx + chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	add_g5ProjPlus(Out, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      add_g5ProjPlus(Out, Add, InScale,n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  //  x + P{-} y   Nc Ns flops
  qz  = qx + chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  add_g5ProjMinus(Out, Add, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x + P_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx + chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx + chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	add_g5ProjMinus(Out, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      add_g5ProjMinus(Out, Add, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // x - P{+} y    Nc Ns
  
  qz  = qx - chiralProjectPlus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  sub_g5ProjPlus(Out, Add,  InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x - P_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx - chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx - chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	sub_g5ProjPlus(Out, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      sub_g5ProjPlus(Out, Add, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x - P{-} y   Nc Ns flops
  qz  = qx - chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  sub_g5ProjMinus(Out, Add, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "x - P_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = qx - chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = qx - chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	sub_g5ProjMinus(Out, Add, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      sub_g5ProjMinus(Out, Add, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }
    
  // a * P_{+} y  Nc, Ns flops
  qz  = a* chiralProjectPlus(qy);

  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  scal_g5ProjPlus(Out, scalep, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "a * P_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = a* chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = a* chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	scal_g5ProjPlus(Out, scalep, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      scal_g5ProjPlus(Out, scalep, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  // a * P{-} y  Nc Ns flops
  qz  = a* chiralProjectMinus(qy);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  n_4vec = all.end()-all.start()+1;

  scal_g5ProjMinus(Out, scalep, InScale, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "a * P_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = a* chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz  = a* chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	scal_g5ProjMinus(Out, scalep, InScale, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      scal_g5ProjMinus(Out, scalep, InScale, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // AX + B P+ y,   4Nc Ns flops

  REAL* scalep2;
  Real b=Real(-3);
  gaussian(qx);
  gaussian(qy);

  qz  = a*qy + b* chiralProjectPlus(qx);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  scalep2 = (REAL *)&(b.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());

  n_4vec = all.end()-all.start()+1;

  axpbyz_g5ProjPlus(Out, scalep, InScale, scalep2, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax + bP_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = b*qy + a* chiralProjectPlus(qx);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
	qz  = b*qy + a* chiralProjectPlus(qx);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axpbyz_g5ProjPlus(Out, scalep, InScale, scalep2, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axpbyz_g5ProjPlus(Out, scalep, InScale, scalep2, Add, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // AX + B P{-} y,  4Nc Ns flops
  gaussian(qx);
  gaussian(qy);

  qz  = a*qy + b* chiralProjectMinus(qx);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  scalep2 = (REAL *)&(b.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());

  n_4vec = all.end()-all.start()+1;

  axpbyz_g5ProjMinus(Out, scalep, InScale, scalep2, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax + bP_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = b*qy + a* chiralProjectMinus(qx);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
	qz  = b*qy + a* chiralProjectMinus(qx);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axpbyz_g5ProjMinus(Out, scalep, InScale, scalep2, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axpbyz_g5ProjMinus(Out, scalep, InScale, scalep2, Add, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // AX - B P+Y,  4 Nc Ns flops
  gaussian(qx);
  gaussian(qy);

  qz  = a*qy - b* chiralProjectPlus(qx);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  scalep2 = (REAL *)&(b.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());

  n_4vec = all.end()-all.start()+1;

  axmbyz_g5ProjPlus(Out, scalep, InScale, scalep2, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax - bP_{+}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = b*qy - a* chiralProjectPlus(qx);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
	qz  = b*qy - a* chiralProjectPlus(qx);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axmbyz_g5ProjPlus(Out, scalep, InScale, scalep2, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axmbyz_g5ProjPlus(Out, scalep, InScale, scalep2, Add, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  // AX + B P{-} Y, 4Nc Ns flops
  gaussian(qx);
  gaussian(qy);

  qz  = a*qy - b* chiralProjectMinus(qx);
  Out = (REAL *)&(qz2.elem(0).elem(0).elem(0).real());
  scalep = (REAL *)&(a.elem().elem().elem().elem());
  scalep2 = (REAL *)&(b.elem().elem().elem().elem());
  InScale = (REAL *)&(qy.elem(0).elem(0).elem(0).real());
  Add = (REAL *)&(qx.elem(0).elem(0).elem(0).real());

  n_4vec = all.end()-all.start()+1;

  axmbyz_g5ProjMinus(Out, scalep, InScale, scalep2, Add, n_4vec);
 
 
  norm_diff=norm2(qz-qz2);
  {
    QDPIO::cout << "ax - bP_{-}y diff=" << norm_diff << endl;
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	qz  = b*qy - a* chiralProjectMinus(qx);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);

    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
	qz  = b*qy - a* chiralProjectMinus(qx);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " Mflop/s per node" << endl;
  }
  {
    StopWatch swatch;
    double time=0;
    int iter=1;
    while( time < 1.0 ) { 
      iter *=2;
      QDPIO::cout << "Calling " << iter << " times " << endl;
      swatch.reset();
      swatch.start();
      for(int i=0; i < iter; i++) {
	axmbyz_g5ProjMinus(Out, scalep, InScale, scalep2, Add, n_4vec);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      axmbyz_g5ProjMinus(Out, scalep, InScale, scalep2, Add, n_4vec);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

    
  // Time to bolt
  QDP_finalize();

  exit(0);
}
  
