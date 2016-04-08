// $Id: t_linalg.cc,v 1.20 2006-09-27 17:26:43 bjoo Exp $

#include <iostream>
#include <cstdio>

#include <time.h>

#include "qdp.h"
#include "linalg.h"

using namespace QDP;
#define TIME_OPS

int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {16,16,16,16};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements
  Layout::setLattSize(nrow);
  Layout::create();

#ifdef QDP_USE_LIBXML2
  XMLFileWriter xml("t_linalg.xml");
  push(xml, "linalgTest");

  push(xml,"lattis");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"nrow", nrow);
  pop(xml);
#endif

  QDPIO::cout << "CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << std::endl;

  LatticeColorMatrix a, b, c;
  gaussian(a);
  gaussian(b);
  gaussian(c);

  int icnt;
  double tt;

#define TIME_OPS 
  // Test M = M
  LatticeColorMatrix m1, m2;
  m1 = zero;
  gaussian(m2);

  for(int i=all.start(); i <= all.end(); i++) { 
    m1.elem(i) -= m2.elem(i);
  }

  LatticeColorMatrix m3=zero;
  m3 -= m2;
  LatticeColorMatrix diff_m;
  diff_m = m3 - m1;
  QDPIO::cout << "Diff M=M = " << norm2(diff_m) << std::endl;
  QDP::StopWatch swatch;
  swatch.reset();
  double time = 0;
  icnt = 1;
  while(time <= 1000000) { 
    swatch.start();
    for(int j=0; j < icnt; j++) {
      for(int i=all.start(); i <= all.end(); i++) { 
	m1.elem(i) -= m2.elem(i);
      }
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    swatch.reset();
    icnt*=2;
  }
  QDPIO::cout << "Call time (old M=M) = " << time / icnt << " us per call" << std::endl;
  
  swatch.reset();
  time = 0;
  icnt = 1;
  while(time <= 1000000) { 
    swatch.start();
    for(int j=0; j < icnt; j++) {
      m1-=m2;
    }
    swatch.stop();
    time = swatch.getTimeInMicroseconds();
    swatch.reset();
    icnt*=2;
  }
  
  QDPIO::cout << "Call time (New M=M)= " << time / icnt << " us per call" << std::endl;



  // Test M=M*M
  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling M=M*M " << icnt << " times" << std::endl;
    tt = QDP_M_eq_M_times_M(c, a, b, icnt);
#if defined(TIME_OPS)
    if (tt > 1)
      break;
#else
    // turn off timings for some testing
    QDPIO::cout << "***WARNING*** : debug mode - timings are bogus" << std::endl;
    break;
#endif
  }

  double rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;

  tt *= rescale;
  int Nflops = Nc*Nc*(4*Nc + (4*Nc-2));
#if defined(TIME_OPS)
  QDPIO::cout << "time(M=M*M) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_eq_M_times_M");
  write(xml,"c", c);
  pop(xml);
#endif

#endif
  
  // Test  M=adj(M)*M
  QDPIO::cout << "calling M=adj(M)*M " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_eq_Ma_times_M(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M=adj(M)*M) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_eq_Ma_times_M");
  write(xml,"a",a);
  write(xml,"b",b);
  write(xml,"c",c);
  pop(xml);
#endif
#endif
  
  // Test  M=M*adj(M)
  QDPIO::cout << "calling M=M*adj(M) " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_eq_M_times_Ma(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M=M*adj(M)) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_eq_M_times_Ma");
  write(xml,"a",a);
  write(xml,"b",b);
  write(xml,"c",c);
  pop(xml);
#endif
#endif

 
  // Test  M=adj(M)*adj(M)
  QDPIO::cout << "calling M=adj(M)*adj(M) " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_eq_Ma_times_Ma(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M=adj(M)*adj(M)) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_eq_Ma_times_Ma");
  write(xml,"a",a);
  write(xml,"b",b);
  write(xml,"c",c);
  pop(xml);
#endif
#endif

  
  // Test  M+= M*M
  QDPIO::cout << "calling M+=M*M " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_peq_M_times_M(c, a, b, icnt);
  Nflops += Nc*Nc * 2;
#if defined(TIME_OPS)
  QDPIO::cout << "time(M+=M*M) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_peq_M_times_M");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


  // Test  M+= adj(M)*M
  QDPIO::cout << "calling M+=adj(M)*M " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_peq_Ma_times_M(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M+=adj(M)*M) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_peq_Ma_times_M");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


  // Test  M+= M*adj(M)
  QDPIO::cout << "calling M+=M*adj(M) " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_peq_M_times_Ma(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M+=M*adj(M)) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_peq_M_times_Ma");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


  // Test  M+= adj(M)*adj(M)
  QDPIO::cout << "calling M+=adj(M)*adj(M) " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_peq_Ma_times_Ma(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M+=adj(M)*adj(M)) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_peq_Ma_times_Ma");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


  // Test  M-= M*M
  QDPIO::cout << "calling M-=M*M " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_meq_M_times_M(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M-=M*M) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_meq_M_times_M");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


  // Test  M-= adj(M)*M
  QDPIO::cout << "calling M-=adj(M)*M " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_meq_Ma_times_M(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M-=adj(M)*M) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_meq_Ma_times_M");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


  // Test  M-= M*adj(M)
  QDPIO::cout << "calling M-=M*adj(M) " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_meq_M_times_Ma(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M-=M*adj(M)) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_meq_M_times_Ma");
  write(xml,"c",c);
  pop(xml)
#endif
#endif


  // Test  M-= adj(M)*adj(M)
  QDPIO::cout << "calling M-=adj(M)*adj(M) " << icnt << " times" << std::endl;
  tt = rescale * QDP_M_meq_Ma_times_Ma(c, a, b, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(M-=adj(M)*adj(M)) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << Nflops / tt << " Mflops" << std::endl;
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_M_meq_Ma_times_Ma");
  write(xml,"c",c);
  pop(xml);
#endif
#endif


//----------------------------------------------------------------------------
  LatticeColorVector lv1,lv2,lv3;
  gaussian(lv1);
  gaussian(lv2);
  gaussian(lv3);

  // Test LatticeColorVector = LatticeColorMatrix * LatticeColorVector
  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=M*V " << icnt << " times" << std::endl;
    tt = QDP_V_eq_M_times_V(lv2, a, lv1, icnt);
#if defined(TIME_OPS)
    if (tt > 1)
      break;
#else
    // turn off timings for some testing
    break;
#endif
  }

  rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;

  tt *= rescale;
#if defined(TIME_OPS)
  QDPIO::cout << "time(V=M*V) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 66 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_V_eq_M_times_V");
  write(xml,"lv2",lv2);
  pop(xml)
#endif
#endif


  // Test LatticeColorVector = LatticeColorMatrix * LatticeColorVector
  QDPIO::cout << "calling V=adj(M)*V " << icnt << " times" << std::endl;
  tt = rescale * QDP_V_eq_Ma_times_V(lv2, a, lv1, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(V=adj(M)*V) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 66 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_V_eq_Ma_times_V");
  write(xml,"lv2",lv2);
  pop(xml);
#endif
#endif


//----------------------------------------------------------------------------
  // Test LatticeColorVector = LatticeColorVector + LatticeColorVector
  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling V=V+V " << icnt << " times" << std::endl;
    tt = QDP_V_eq_V_plus_V(lv3, lv1, lv2, icnt);
#if defined(TIME_OPS)
    if (tt > 1)
      break;
#else
    // turn off timings for some testing
    break;
#endif
  }

  rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;

  tt *= rescale;
#if defined(TIME_OPS)
  QDPIO::cout << "time(V=V+V) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 6 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_V_eq_V_plus_V");
  write(xml,"lv3",lv3);
  pop(xml);
#endif
#endif



//----------------------------------------------------------------------------
  LatticeDiracFermion lf1,lf2,lf3;
  gaussian(lf1);
  gaussian(lf2);

  // Test LatticeDiracFermion = LatticeColorMatrix * LatticeDiracFermion
  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling D=M*D " << icnt << " times" << std::endl;
    tt = QDP_D_eq_M_times_D(lf2, a, lf1, icnt);
#if defined(TIME_OPS)
    if (tt > 1)
      break;
#else
    // turn off timings for some testing
    break;
#endif
  }

  rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;

  tt *= rescale;
#if defined(TIME_OPS)
  QDPIO::cout << "time(D=M*D) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 264 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_D_eq_M_times_D");
  write(xml,"lf2",lf2);
  pop(xml);
#endif
#endif

  // Test LatticeDiracFermion = adj(LatticeColorMatrix) * LatticeDiracFermion
  QDPIO::cout << "calling D=adj(M)*D " << icnt << " times" << std::endl;
  tt = rescale * QDP_D_eq_Ma_times_D(lf2, a, lf1, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(D=adj(M)*D) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 264 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_D_eq_Ma_times_D");
  write(xml,"lf2",lf2);
  pop(xml);
#endif
#endif


//----------------------------------------------------------------------------
  LatticeHalfFermion lh1,lh2,lh3;
  gaussian(lh1);
  gaussian(lh2);

  // Test LatticeHalfFermion = LatticeColorMatrix * LatticeHalfFermion
  for(icnt=1; ; icnt <<= 1)
  {
    QDPIO::cout << "calling H=M*H " << icnt << " times" << std::endl;
    tt = QDP_H_eq_M_times_H(lh2, a, lh1, icnt);
#if defined(TIME_OPS)
    if (tt > 1)
      break;
#else
    // turn off timings for some testing
    break;
#endif
  }

  rescale = 1000*1000 / double(Layout::sitesOnNode()) / icnt;

  tt *= rescale;
#if defined(TIME_OPS)
  QDPIO::cout << "time(H=M*H) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 132 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_H_eq_M_times_H");
  write(xml,"lh2", lh2);
  pop(xml);
#endif
#endif


  // Test LatticeHalfFermion = adj(LatticeColorMatrix) * LatticeHalfFermion
  QDPIO::cout << "calling H=adj(M)*H " << icnt << " times" << std::endl;
  tt = rescale * QDP_H_eq_Ma_times_H(lh2, a, lh1, icnt);
#if defined(TIME_OPS)
  QDPIO::cout << "time(H=adj(M)*H) = " << tt
	      << " micro-secs/site/iteration" 
	      << " , " << 132 / tt << " Mflops" << std::endl;   // check the flop count
#else
#ifdef QDP_USE_LIBXML2
  push(xml,"QDP_H_eq_Ma_times_H");
  write(xml,"lh2", lh2);
  pop(xml);
#endif

#endif


#ifdef QDP_USE_LIBXML2
  pop(xml);
  xml.close();
#endif
  // Time to bolt
  QDP_finalize();

  exit(0);
}
