// $Id: t_cblas.cc,v 1.4 2007-02-22 15:58:30 bjoo Exp $

#include <iostream>
#include <cstdio>

#include <time.h>

#include "qdp.h"

#include "cblas1.h"
 
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

  LatticeFermion x, y, z1, z2,r;

  Complex alpha=cmplx(Real(1), Real(-2));
  Complex beta =cmplx(Real(-2), Real(-2.1));

  gaussian(x);
  // Do z1 = alpha x with an axpy
  z2 = zero;
  z1 = alpha*x + z2;
  
  // Do z2 = alpha x using *=
  z2 = x;
  z2 *= alpha;

  r = z1 - z2;
  Double diff = norm2(r);
  QDPIO::cout << " z2 *= a diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  // Do z2 = alpha * x using = 
  z2 = alpha*x;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << " z2 = a * x  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  // Do z2 = x * alpha using = 
  z2 = x*alpha;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << " z2 = x * a diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;
  
  // Do z2 = x, z2 = a*z2
  z2 = x;
  z2 = alpha*z2;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z2=x, z2 = a * z2 diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;


  gaussian(x);
  gaussian(y);

  z1 = alpha*x;
  z1 += y;

  z2 = y;
  z2 += alpha*x;

  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z += alpha*x  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = y;
  z2 += x*alpha;

  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z += x*alpha  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;


  z1 = alpha*x;
  z1 = -z1;
  z1 += y;

  z2 = y;
  z2 -= alpha*x;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z -= x*alpha  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = y;
  z2 -= x*alpha;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z -= x*alpha  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  
  z1 = alpha*x;
  z1 += y;

  z2 = alpha*x + y;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = alpha * x + y  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = x * alpha + y;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = x * alpha + y  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = y + alpha*x;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = y + alpha * x  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = y + x*alpha;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = y + x * alpha  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;


  z1 = alpha*x;
  z1 -= y;
  z2 = alpha*x - y;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = alpha * x - y  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = x * alpha - y;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = x * alpha - y  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z1 = y;
  z1 -= alpha*x;

  z2 = y - alpha*x;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = y - alpha*x  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = y - x*alpha;
  r = z1 - z2;
  diff = norm2(r);
  QDPIO::cout << "z = y - x*alpha  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;



  // AX + BY type tests
  gaussian(x);
  gaussian(y);

  // Make z1 = alpha x + beta y
  z1 = alpha*x;
  z2 = beta*y;
  z1 = z1 + z2;

  // Now do z2 with the 4 combinations
  z2 = alpha*x + beta*y;
  r=z2-z1;
  diff = norm2(r);
  QDPIO::cout << "z = ax + by  diff = " << sqrt(diff)/sqrt(norm2(z1)) << endl;

  z2 = x*alpha + beta*y;
  r=z2-z1;  
  diff = norm2(r);
  QDPIO::cout << "z = xa + by  diff = " << sqrt(diff)/sqrt(norm2(z1)) << endl;

  z2 = alpha*x + y*beta;
  r=z2-z1;
  diff = norm2(r);
  QDPIO::cout << "z = ax + yb  diff = " << sqrt(diff)/sqrt(norm2(z1)) << endl;

  z2 = x*alpha + y*beta;
  r=z2-z1;
  diff = norm2(r);
  QDPIO::cout << "z = xa + yb  diff = " << sqrt(diff)/sqrt(norm2(z1)) << endl;

  // Make z1 = alpha x + beta y
  z1 = alpha*x;
  z2 = beta*y;
  z1 = z1 - z2;

  // Now do z2 with the 4 combinations
  z2 = alpha*x - beta*y;
  r=z2-z1;
  diff = norm2(r);
  QDPIO::cout << "z = ax - by  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = x*alpha - beta*y;
  r=z2-z1;  
  diff = norm2(r);
  QDPIO::cout << "z = xa - by  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = alpha*x - y*beta;
  r=z2-z1;
  diff = norm2(r);
  QDPIO::cout << "z = ax - yb  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;

  z2 = x*alpha - y*beta;
  r=z2-z1;
  diff = norm2(r);
  QDPIO::cout << "z = xa - yb  diff = " << sqrt(diff)/sqrt(norm2(z2)) << endl;


  // Timings
   // Test VSCAL
  int icnt;
  double tt;
  gaussian(x);

  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=a*V " << icnt << " times" << endl;
    tt = QDP_CSCALE(y, alpha, x, icnt);
    if (tt > 1)
      break;
  }

  {
    double rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;
    tt *= rescale;
    int Nflops = 2*Ns*Nc;
    QDPIO::cout << "time(V=aV) = " << tt
		<< " micro-secs/site/iteration" 
		<< " , " << Nflops / tt << " Mflops" << endl;
  }


   // Test VAXPY
  gaussian(x);
  gaussian(y);

  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=aV+V " << icnt << " times" << endl;
    tt = QDP_CAXPY(z1, alpha, x, y, icnt);
    if (tt > 1)
      break;
  }
  {
    double rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;
    tt *= rescale;
    int Nflops = 4*Ns*Nc;
    QDPIO::cout << "time(V=aV+V) = " << tt
		<< " micro-secs/site/iteration" 
		<< " , " << Nflops / tt << " Mflops" << endl;
  }


   // Test VAXMY
  gaussian(x);
  gaussian(y);

  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=aV-V " << icnt << " times" << endl;
    tt = QDP_CAXMY(z1, alpha, x, y, icnt);
    if (tt > 1)
      break;
  }

  {
    double rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;
    tt *= rescale;
    int Nflops = 4*Ns*Nc;
    QDPIO::cout << "time(V=aV-V) = " << tt
		<< " micro-secs/site/iteration" 
		<< " , " << Nflops / tt << " Mflops" << endl;
  }

   // Test VAXPBY
  gaussian(x);
  gaussian(y);

  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=aV+bV " << icnt << " times" << endl;
    tt = QDP_CAXPBY(z1, alpha, beta, x, y, icnt);
    if (tt > 1)
      break;
  }
  {
    double rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;
    tt *= rescale;
    int Nflops = 3*2*Ns*Nc;
    QDPIO::cout << "time(V=aV+bV) = " << tt
		<< " micro-secs/site/iteration" 
		<< " , " << Nflops / tt << " Mflops" << endl;
  }


   // Test VAXMBY
  gaussian(x);
  gaussian(y);

  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=aV-bV " << icnt << " times" << endl;
    tt = QDP_CAXMBY(z1, alpha, beta, x, y, icnt);
    if (tt > 1)
      break;
  }

  {
    double rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;
    tt *= rescale;
    int Nflops = 3*2*Ns*Nc;
    QDPIO::cout << "time(V=aV-bV) = " << tt
		<< " micro-secs/site/iteration" 
		<< " , " << Nflops / tt << " Mflops" << endl;
  }

  // Time to bolt
  QDP_finalize();

  exit(0);
}
