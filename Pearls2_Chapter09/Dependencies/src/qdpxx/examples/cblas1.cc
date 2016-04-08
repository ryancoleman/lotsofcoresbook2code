// $Id: cblas1.cc,v 1.1 2004-05-09 11:54:43 bjoo Exp $

#include <stdlib.h>
#include <sys/time.h>

#include "qdp.h"
#include "cblas1.h"


using namespace QDP;


double QDP_CSCALE(LatticeFermion& dest, 
		 const Complex& a, 
		 const LatticeFermion& s,
		 int cnt)
{
  StopWatch swatch;

  swatch.start();
  for (; cnt-- > 0; )
    dest =  a * s;
  swatch.stop(); 
  return swatch.getTimeInSeconds();
//  return 2;
}

double QDP_CAXPY(LatticeFermion& dest, 
		const Complex& a,
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt)
{
  StopWatch swatch;
 
  swatch.start();
  for (; cnt-- > 0; )
    dest = a*s1 + s2;
  swatch.stop();
  return swatch.getTimeInSeconds();
//    return 2.0;
}

double QDP_CAXMY(LatticeFermion& dest, 
		const Complex& a,
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt)
{
  StopWatch swatch;
                                                                                
  swatch.start();
  for (; cnt-- > 0; )
    dest = a*s1 - s2;
  swatch.stop();

  return swatch.getTimeInSeconds();
//  return 2.0;
}

double QDP_CAXPBY(LatticeFermion& dest, 
		 const Complex& a,
		 const Complex& b,
		 const LatticeFermion& s1, 
		 const LatticeFermion& s2,
		 int cnt)
{
  StopWatch swatch;
                                                                                
  swatch.start();
  for (; cnt-- > 0; )
    dest = a*s1 + b*s2;
  swatch.stop();

  return swatch.getTimeInSeconds();
//  return 2.0;
}

double QDP_CAXMBY(LatticeFermion& dest, 
		 const Complex& a,
		 const Complex &b,
		 const LatticeFermion& s1, 
		 const LatticeFermion& s2,
		 int cnt)
{
  StopWatch swatch;
                                                                                
  swatch.start();
  for (; cnt-- > 0; )
    dest = a*s1 - b*s2;
  swatch.stop();

  return swatch.getTimeInSeconds();
//  return 2.0;
}


