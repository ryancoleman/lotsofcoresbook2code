// -*- C++ -*-
// $Id: blas1.h,v 1.6 2005-05-25 04:20:32 edwards Exp $
//
// Include file for test suite

#ifndef BLAS1_H
#define BLAS1_H

#include "qdp.h"

using namespace QDP;


double QDP_SCALE(LatticeFermion& dest, 
		 const Real& a, 
		 const LatticeFermion& s,
		 int cnt);

double QDP_AXPY(LatticeFermion& dest, 
		const Real& a,
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt);

double QDP_AXMY(LatticeFermion& dest, 
		const Real& a,
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt);

double QDP_AXPBY(LatticeFermion& dest, 
		 const Real& a,
		 const Real& b,
		 const LatticeFermion& s1, 
		 const LatticeFermion& s2,
		 int cnt);

double QDP_AXMBY(LatticeFermion& dest, 
		 const Real& a,
		 const Real& b,
		 const LatticeFermion& s1, 
		 const LatticeFermion& s2,
		 int cnt);

double QDP_VADD(LatticeFermion& dest, 
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt);

double QDP_VSUB(LatticeFermion& dest, 
		const LatticeFermion& s1, 
		const LatticeFermion& s2,
		int cnt);

double QDP_NORM2(const LatticeFermion& s1, int cnt);
double QDP_INNER_PROD(const LatticeFermion& s1, const LatticeFermion& s2, 
		      int cnt);
double QDP_INNER_PROD_REAL(const LatticeFermion& s1, const LatticeFermion& s2,
			   int cnt);

double QDP_NORM2(const multi1d<LatticeFermion>& s1, int cnt);
double QDP_INNER_PROD(const multi1d<LatticeFermion>& s1, const multi1d<LatticeFermion>& s2, 
		      int cnt);
double QDP_INNER_PROD_REAL(const multi1d<LatticeFermion>& s1, const multi1d<LatticeFermion>& s2,
			   int cnt);

#endif
