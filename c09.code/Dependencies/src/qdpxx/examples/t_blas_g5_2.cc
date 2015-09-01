// $Id: t_blas_g5_2.cc,v 1.3 2005-06-27 14:13:24 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
 
using namespace QDP;

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {4,4,4,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

  
  Real a=Real(1.5);
  LatticeFermion qx;  qx.moveToFastMemoryHint();
  LatticeFermion qy;  qy.moveToFastMemoryHint();
  LatticeFermion qz;  qz.moveToFastMemoryHint();
  LatticeFermion qz2; qz2.moveToFastMemoryHint();
 
  gaussian(qx);
  qy=qx;


  Double norm_diff;


  //  x += a*P{+}y
  qz = zero;
  qz2 = zero;

  qx = Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz += a*qx;
  qz2 +=  a*chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z += a * P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 +=  a*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 +=  a*chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x += a*P{-}y
  qz = zero;
  qz2 = zero;

  qx = Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz += a*qx;
  qz2 +=  a*chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z += a * P- y: Norm diff = " << norm_diff << endl;
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
	qz2 +=  a*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 +=  a*chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  //  x -= a*P{+}y
  gaussian(qz);
  qz2 = qz;

  qx = Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz -= a*qx;
  qz2 -=  a*chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z -= a * P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 -=  a*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 -=  a*chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x -= a*P{-}y
  gaussian(qz);;
  qz2 = qz;

  qx = Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz -= a*qx;
  qz2 -=  a*chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z -= a * P- y: Norm diff = " << norm_diff << endl;
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
	qz2 -=  a*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 -=  a*chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x += P{+}y
  qz = zero;
  qz2 = zero;

  qx = Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz += qx;
  qz2 +=  chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z +=  P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 += chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 += chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x += P{-}y
  qz = zero;
  qz2 = zero;

  qx = Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz += qx;
  qz2 +=  chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z +=  P- y: Norm diff = " << norm_diff << endl;
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
	qz2 += chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 += chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x -= P{+}y
  gaussian(qz);
  qz2 = qz;

  qx = Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz -= qx;
  qz2 -=  chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z -=  P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 -= chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 -= chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  x -= P{-}y
  gaussian(qz);
  qz2 = qz;

  qx = Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz -= qx;
  qz2 -=  chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z -=  P- y: Norm diff = " << norm_diff << endl;
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
	qz2 -= chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2 -= chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = x + a P{+}y
  gaussian(qz);
  gaussian(qx);
 

  qz = qx + a*Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =qx + a* chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = x + a P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =qx + a* chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  =qx + a* chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = x + a P{-}y
  gaussian(qz);
  gaussian(qx);
 

  qz = qx + a*Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =qx + a* chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = x + a P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =qx + a* chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  =qx + a* chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = x - a P{+}y
  gaussian(qz);
  gaussian(qx);
 

  qz = qx - a*Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =qx - a* chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = x - a P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =qx - a* chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  =qx - a* chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = x - a P{-}y
  gaussian(qz);
  gaussian(qx);
 

  qz = qx - a*Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =qx - a* chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = x - a P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =qx - a* chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  =qx - a* chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(2*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = ax + P{+}y
  gaussian(qz);
  gaussian(qx);
 

  qz = a* qx + Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx + chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x + P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx + chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx + chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = ax + P{-}y
  gaussian(qz);
  gaussian(qx);
 

  qz = a* qx + Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx + chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x + P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx + chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx + chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = ax -  P{+}y
  gaussian(qz);
  gaussian(qx);
 

  qz = a* qx - Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx - chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x - P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx - chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx - chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = ax - P{-}y
  gaussian(qz);
  gaussian(qx);
 

  qz = a* qx - Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx - chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x - P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx - chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx - chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(3*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = a P{+}y
  qz = a* Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a* chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = a P{-}y
  qz = a* Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a* chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }

  //  z = ax + b*P{+}y
  gaussian(qz);
  gaussian(qx);
  Real b=Real(-5.3);


  qz = a* qx + b*Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx + b*chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x + b P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx + b*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx + b*chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }
 
 //  z = ax + b*P{-}y
  gaussian(qz);
  gaussian(qx);
  b=Real(-5.3);


  qz = a* qx + b*Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx + b*chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x + b P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx + b*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx + b*chiralProjectMinus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }


  //  z = ax - b*P{+}y
  gaussian(qz);
  gaussian(qx);
  b=Real(-5.3);


  qz = a* qx - b*Real(0.5)*(qy + GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx - b*chiralProjectPlus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x - b P+ y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx - b*chiralProjectPlus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx - b*chiralProjectPlus(qy);
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    QDPInternal::broadcast(time);
    
    double Nflops = (double)(4*Nc*Ns*Layout::sitesOnNode()*iter);
    QDPIO::cout << "Time taken: " << time << "(us) Perf: " << Nflops/time << " MFlops/node" << endl;
  }
 
 //  z = ax - b*P{-}y
  gaussian(qz);
  gaussian(qx);
  b=Real(-5.3);


  qz = a* qx - b*Real(0.5)*(qy - GammaConst<Ns,Ns*Ns-1>()*qy);
  qz2 =a* qx - b*chiralProjectMinus(qy);
 
  norm_diff=norm2(qz-qz2);
  QDPIO::cout << "z = a x - b P- y: Norm diff = " << norm_diff << endl;
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
	qz2 =a*qx - b*chiralProjectMinus(qy);
      }
      swatch.stop();
      time = swatch.getTimeInSeconds();
      QDPInternal::broadcast(time);
    }
    
    QDPIO::cout << "Timing with " << iter << " iters" << endl;
    swatch.reset();
    swatch.start();
    for(int i=0; i < iter; i++) {
      qz2  = a*qx - b*chiralProjectMinus(qy);
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
  
