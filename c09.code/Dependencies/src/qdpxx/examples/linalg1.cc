// $Id: linalg1.cc,v 1.19 2006-09-26 01:58:36 edwards Exp $

#include <stdlib.h>
#include <sys/time.h>

#include "qdp.h"
#include "linalg.h"


using namespace QDP;


double QDP_M_eq_M_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt)
{
  StopWatch swatch;

  swatch.start();
  for (; cnt-- > 0; )
    dest = s1 * s2;
  swatch.stop(); 
  return swatch.getTimeInSeconds();
//  return 2;
}

double QDP_M_eq_Ma_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt)
{
  StopWatch swatch;
 
  swatch.start();
  for (; cnt-- > 0; )
    dest = adj(s1) * s2;
  swatch.stop();
  return swatch.getTimeInSeconds();
//    return 2.0;
}

double QDP_M_eq_M_times_Ma(LatticeColorMatrix& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeColorMatrix& s2,
			   int cnt)
{
  StopWatch swatch;
                                                                                
  swatch.start();
  for (; cnt-- > 0; )
    dest = s1 * adj(s2);
  swatch.stop();

  return swatch.getTimeInSeconds();
//  return 2.0;
}

double QDP_M_eq_Ma_times_Ma(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt)
{

  StopWatch swatch;
  swatch.start();

  for (; cnt-- > 0; )
    dest = adj(s1) * adj(s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2.0;
}

double QDP_M_peq_M_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest += s1 * s2;

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_peq_M_times_Ma(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest += s1 * adj(s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_peq_Ma_times_M(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest += adj(s1) * s2;

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_peq_Ma_times_Ma(LatticeColorMatrix& dest, 
			     const LatticeColorMatrix& s1, 
			     const LatticeColorMatrix& s2,
			     int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest += adj(s1) * adj(s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_meq_M_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest -= s1 * s2;

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_meq_M_times_Ma(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest -= s1 * adj(s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_meq_Ma_times_M(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest -= adj(s1) * s2;

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_M_meq_Ma_times_Ma(LatticeColorMatrix& dest, 
			     const LatticeColorMatrix& s1, 
			     const LatticeColorMatrix& s2,
			     int cnt)
{
  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest -= adj(s1) * adj(s2);

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2;
}


double QDP_V_eq_M_times_V(LatticeColorVector& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorVector& s2,
			  int cnt)
{

  StopWatch swatch;

  swatch.start();
  for (; cnt-- > 0; )
    dest = s1 * s2;
  swatch.stop();


  return swatch.getTimeInSeconds();
//    return 2.0;
}


double QDP_V_eq_Ma_times_V(LatticeColorVector& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeColorVector& s2,
			   int cnt)
{

  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest = adj(s1) * s2;
  
  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2.0;
}


double QDP_V_eq_V_plus_V(LatticeColorVector& dest, 
			 const LatticeColorVector& s1, 
			 const LatticeColorVector& s2,
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


double QDP_D_eq_M_times_D(LatticeDiracFermion& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeDiracFermion& s2,
			  int cnt)
{
  StopWatch swatch;

  swatch.start();
  for (; cnt-- > 0; )
    dest = s1 * s2;

  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2.0;
}


double QDP_D_eq_Ma_times_D(LatticeDiracFermion& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeDiracFermion& s2,
			   int cnt)
{

  StopWatch swatch;

  swatch.start();

  for (; cnt-- > 0; )
    dest = adj(s1) * s2;

  swatch.stop();


  return swatch.getTimeInSeconds();
//    return 2.0;
}


double QDP_H_eq_M_times_H(LatticeHalfFermion& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeHalfFermion& s2,
			  int cnt)
{

  StopWatch swatch;

  swatch.start();
  for (; cnt-- > 0; )
    dest = s1 * s2;
  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2.0;
}


double QDP_H_eq_Ma_times_H(LatticeHalfFermion& dest, 
		 	   const LatticeColorMatrix& s1, 
			   const LatticeHalfFermion& s2,
			   int cnt)
{
  StopWatch swatch;

  swatch.start();
  for (; cnt-- > 0; )
    dest = adj(s1) * s2;
  swatch.stop();

  return swatch.getTimeInSeconds();
//    return 2.0;
}


