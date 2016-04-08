// -*- C++ -*-
// $Id: linalg.h,v 1.9 2006-09-26 01:58:36 edwards Exp $
//
// Include file for test suite

#include "qdp.h"

using namespace QDP;

double QDP_M_eq_M_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt);

double QDP_M_eq_Ma_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt);

double QDP_M_eq_M_times_Ma(LatticeColorMatrix& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeColorMatrix& s2,
			   int cnt);

double QDP_M_eq_Ma_times_Ma(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt);

double QDP_M_peq_M_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt);

double QDP_M_peq_Ma_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt);

double QDP_M_peq_M_times_Ma(LatticeColorMatrix& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeColorMatrix& s2,
			   int cnt);

double QDP_M_peq_Ma_times_Ma(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt);

double QDP_M_meq_M_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt);

double QDP_M_meq_Ma_times_M(LatticeColorMatrix& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorMatrix& s2,
			  int cnt);

double QDP_M_meq_M_times_Ma(LatticeColorMatrix& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeColorMatrix& s2,
			   int cnt);

double QDP_M_meq_Ma_times_Ma(LatticeColorMatrix& dest, 
			    const LatticeColorMatrix& s1, 
			    const LatticeColorMatrix& s2,
			    int cnt);

double QDP_V_eq_M_times_V(LatticeColorVector& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeColorVector& s2,
			  int cnt);

double QDP_V_eq_Ma_times_V(LatticeColorVector& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeColorVector& s2,
			   int cnt);

double QDP_V_eq_V_plus_V(LatticeColorVector& dest, 
			 const LatticeColorVector& s1, 
			 const LatticeColorVector& s2,
			 int cnt);

double QDP_D_eq_M_times_D(LatticeDiracFermion& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeDiracFermion& s2,
			  int cnt);

double QDP_D_eq_Ma_times_D(LatticeDiracFermion& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeDiracFermion& s2,
			   int cnt);

double QDP_H_eq_M_times_H(LatticeHalfFermion& dest, 
			  const LatticeColorMatrix& s1, 
			  const LatticeHalfFermion& s2,
			  int cnt);

double QDP_H_eq_Ma_times_H(LatticeHalfFermion& dest, 
			   const LatticeColorMatrix& s1, 
			   const LatticeHalfFermion& s2,
			   int cnt);

