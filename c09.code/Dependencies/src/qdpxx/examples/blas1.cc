// $Id: blas1.cc,v 1.6 2005-05-25 04:20:32 edwards Exp $

#include <stdlib.h>
#include <sys/time.h>

#include "qdp.h"
#include "blas1.h"


using namespace QDP;


double QDP_SCALE(LatticeFermion& dest, 
		 const Real& a, 
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

double QDP_AXPY(LatticeFermion& dest, 
		const Real& a,
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

double QDP_AXMY(LatticeFermion& dest, 
		const Real& a,
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

double QDP_AXPBY(LatticeFermion& dest, 
		 const Real& a,
		 const Real& b,
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

double QDP_AXMBY(LatticeFermion& dest, 
		 const Real& a,
		 const Real &b,
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

double QDP_VADD(LatticeFermion& dest, 
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt)
{

  StopWatch swatch;
  swatch.start();

  for (; cnt-- > 0; )
    dest = s1 + s2; 

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2.0;
}

double QDP_VSUB(LatticeFermion& dest, 
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest = s1 - s2;

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}

double QDP_NORM2(const LatticeFermion& s1, int cnt)			
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    norm2(s1);

  swatch.stop();

  return swatch.getTimeInSeconds();
  // return 2;
}

double QDP_INNER_PROD(const LatticeFermion& s1, const LatticeFermion& s2, 
		      int cnt)			
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    innerProduct(s1,s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
  // return 2;
}

double QDP_INNER_PROD_REAL(const LatticeFermion& s1, const LatticeFermion& s2,
			   int cnt)			
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    innerProductReal(s1,s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_NORM2(const multi1d<LatticeFermion>& s1, int cnt)	 		
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    norm2(s1);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}

double QDP_INNER_PROD(const multi1d<LatticeFermion>& s1, const multi1d<LatticeFermion>& s2, 
		      int cnt)			
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    innerProduct(s1,s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}

double QDP_INNER_PROD_REAL(const multi1d<LatticeFermion>& s1, const multi1d<LatticeFermion>& s2,
			   int cnt)			
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    innerProductReal(s1,s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}

