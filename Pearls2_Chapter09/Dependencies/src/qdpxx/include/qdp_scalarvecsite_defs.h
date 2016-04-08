// -*- C++ -*-

/*! \file
 * \brief Type definitions
 */

#ifndef QDP_SCALARVECSITE_DEFS_H
#define QDP_SCALARVECSITE_DEFS_H

namespace QDP {


/*! \addtogroup defs Type definitions
 *
 * User constructed types made from QDP type compositional nesting.
 * The layout is suitable for a scalar-like implementation. Namely,
 * sites are the slowest varying index.
 *
 * @{
 */

#include <qdp_config.h>
#include "qdp_precision.h"

//----------------------------------------------------------------------
//! Gamma matrices are conveniently defined for this Ns
typedef GammaType<Ns> Gamma;

// For now, fix this inner-grid length to 4 or 2 depending on base precision. 
// This causes problems for doubles which I eventually must work out.
// Here, INNER_LOG is the log_2(INNER)
#define INNER_LEN (1 << INNER_LOG)

// Aliases for a scalarvec-like architecture

// Fixed fermion type
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, 4> > LatticeDiracFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, 1> > LatticeStaggeredFermion;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, 4> > LatticeDiracPropagator;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, 1> > LatticeStaggeredPropagator;

// Floating aliases
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, Nc> > > LatticeColorVector;
typedef OLattice< PSpinVector< PScalar< RComplex< ILattice<REAL,INNER_LEN> > >, Ns> > LatticeSpinVector;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc> > > LatticeColorMatrix;
typedef OLattice< PSpinMatrix< PScalar< RComplex< ILattice<REAL,INNER_LEN> > >, Ns> > LatticeSpinMatrix;
typedef OLattice< PSpinMatrix< PScalar< RComplex< ILattice<REAL,INNER_LEN> > >, (Ns>>1) > > LatticeHalfSpinMatrix;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, Ns> > LatticeFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, Ns>>1 > > LatticeHalfFermion;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, Ns> > LatticePropagator;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, (Ns>>1) > > LatticeHalfPropagator;
typedef OLattice< PScalar< PScalar< RComplex< ILattice<REAL,INNER_LEN> > > > > LatticeComplex;

typedef OLattice< PScalar< PSeed < RScalar< ILattice<INTEGER32,INNER_LEN> > > > > LatticeSeed;
typedef OLattice< PScalar< PScalar< RScalar< ILattice<INTEGER32,INNER_LEN> > > > > LatticeInteger;
typedef OLattice< PScalar< PScalar< RScalar< ILattice<REAL,INNER_LEN> > > > > LatticeReal;
typedef OLattice< PScalar< PScalar< RScalar< ILattice<DOUBLE,INNER_LEN> > > > > LatticeDouble;
typedef OLattice< PScalar< PScalar< RScalar< ILattice<LOGICAL,INNER_LEN> > > > > LatticeBoolean;

typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL> >, Nc> > > ColorVector;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL> >, Nc> > > ColorMatrix;
typedef OScalar< PSpinVector< PScalar< RComplex< IScalar<REAL> > >, Ns> > SpinVector;
typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<REAL> > >, Ns> > SpinMatrix;
typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<REAL> > >, (Ns>>1) > > HalfSpinMatrix;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL> >, Nc>, 1> > StaggeredFermion;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL> >, Nc>, Ns> > Fermion;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL> >, Nc>, Ns>>1 > > HalfFermion;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL> >, Nc>, Ns> > Propagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL> >, Nc>, (Ns>>1) > > HalfPropagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL> >, Nc>, 1> > StaggeredPropagator;
typedef OScalar< PScalar< PScalar< RComplex< IScalar<REAL> > > > > Complex;

typedef OScalar< PScalar< PSeed< RScalar< IScalar<INTEGER32> > > > > Seed;
typedef OScalar< PScalar< PScalar< RScalar< IScalar<INTEGER32> > > > > Integer;
typedef OScalar< PScalar< PScalar< RScalar< IScalar<REAL> > > > > Real;
typedef OScalar< PScalar< PScalar< RScalar< IScalar<DOUBLE> > > > > Double;
typedef OScalar< PScalar< PScalar< RScalar< IScalar<LOGICAL> > > > > Boolean;

typedef OScalar< PScalar< PColorVector< RComplex< IScalar<DOUBLE> >, Nc> > > DColorVector;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<DOUBLE> >, Nc> > > DColorMatrix;
typedef OScalar< PSpinVector< PScalar< RComplex< IScalar<DOUBLE> > >, Ns> > DSpinVector;
typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<DOUBLE> > >, Ns> > DSpinMatrix;
typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<DOUBLE> > >, (Ns>>1) > > DHalfSpinMatrix;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<DOUBLE> >, Nc>, 1> > DStaggeredFermion;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<DOUBLE> >, Nc>, Ns> > DFermion;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<DOUBLE> >, Nc>, Ns>>1 > > DHalfFermion;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<DOUBLE> >, Nc>, Ns> > DPropagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<DOUBLE> >, Nc>, (Ns>>1) > > DHalfPropagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<DOUBLE> >, Nc>, 1> > DStaggeredPropagator;
typedef OScalar< PScalar< PScalar< RComplex< IScalar<DOUBLE> > > > > DComplex;

//baryonblocks
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL> >, Nc>, Ns>, Nc>, Ns> > Baryonblock;
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfBaryonblock;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, Ns>, Nc>, Ns> > LatticeBaryonblock;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfBaryonblock;

//fourquark tensors
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, Ns>, Nc>, Ns> > Fourquark;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, Ns>, Nc>, Ns> > LatticeFourquark;
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfFourquark;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfFourquark;

// Other useful names
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL> >, Nc>, Ns> > ColorVectorSpinVector;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL> >, Nc>, Ns> > ColorMatrixSpinMatrix;

// Level below outer for internal convenience
typedef PScalar< PScalar< RScalar< IScalar<REAL> > > > IntReal;
typedef PScalar< PScalar< RScalar< IScalar<REAL32> > > > IntReal32;
typedef PScalar< PScalar< RScalar< IScalar<INTEGER32> > > > IntInteger;
typedef PScalar< PScalar< RScalar< IScalar<REAL64> > > > IntReal64;
typedef PScalar< PScalar< RScalar< IScalar<DOUBLE> > > > IntDouble;
typedef PScalar< PScalar< RScalar< IScalar<LOGICAL> > > > IntBoolean;

// Odd-ball to support random numbers
typedef OScalar< PScalar< PScalar < RScalar< ILattice<REAL,INNER_LEN> > > > > ILatticeReal;
typedef OScalar< PScalar< PSeed < RScalar< ILattice<INTEGER32,INNER_LEN> > > > > ILatticeSeed;

// Floating precision, but specific to a fixed color or spin
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 3>, 4> > LatticeDiracFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 2>, 4> > LatticeDiracFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 1>, 4> > LatticeDiracFermion1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 3>, 1> > LatticeStaggeredFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 2>, 1> > LatticeStaggeredFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 1>, 1> > LatticeStaggeredFermion1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 3>, 4> > LatticeDiracPropagator3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 2>, 4> > LatticeDiracPropagator2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 1>, 4> > LatticeDiracPropagator1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 3>, 1> > LatticeStaggeredPropagator3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 2>, 1> > LatticeStaggeredPropagator2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 1>, 1> > LatticeStaggeredPropagator1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 3>, Ns> > LatticeFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 2>, Ns> > LatticeFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 1>, Ns> > LatticeFermion1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 3>, Ns>>1 > > LatticeHalfFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 2>, Ns>>1 > > LatticeHalfFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 1>, Ns>>1 > > LatticeHalfFermion1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 3>, Ns> > LatticePropagator3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 2>, Ns> > LatticePropagator2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 1>, Ns> > LatticePropagator1;

typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 3> > > LatticeColorMatrix3;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 2> > > LatticeColorMatrix2;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL,INNER_LEN> >, 1> > > LatticeColorMatrix1;

typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 3> > > LatticeColorVector3;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 2> > > LatticeColorVector2;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL,INNER_LEN> >, 1> > > LatticeColorVector1;

typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL> >, 3> > > ColorMatrix3;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL> >, 2> > > ColorMatrix2;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL> >, 1> > > ColorMatrix1;

typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL> >, 3> > > ColorVector3;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL> >, 2> > > ColorVector2;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL> >, 1> > > ColorVector1;

//
// Fixed precision
//
// REAL32 types
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, 4> > LatticeDiracFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, 4> > LatticeDiracFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, 4> > LatticeDiracFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, 4> > LatticeDiracFermionF1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, 1> > LatticeStaggeredFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, 1> > LatticeStaggeredFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, 1> > LatticeStaggeredFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, 1> > LatticeStaggeredFermionF1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, 4> > LatticeDiracPropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, 4> > LatticeDiracPropagatorF3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, 4> > LatticeDiracPropagatorF2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, 4> > LatticeDiracPropagatorF1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, 1> > LatticeStaggeredPropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, 1> > LatticeStaggeredPropagatorF3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, 1> > LatticeStaggeredPropagatorF2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, 1> > LatticeStaggeredPropagatorF1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, Ns> > LatticeFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, Ns> > LatticeFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, Ns> > LatticeFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, Ns> > LatticeFermionF1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, Ns>>1 > > LatticeHalfFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, Ns>>1 > > LatticeHalfFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, Ns>>1 > > LatticeHalfFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, Ns>>1 > > LatticeHalfFermionF1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, Ns> > LatticePropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, (Ns>>1) > > LatticeHalfPropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 3>, Ns> > LatticePropagatorF3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 2>, Ns> > LatticePropagatorF2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 1>, Ns> > LatticePropagatorF1;

typedef OLattice< PSpinMatrix< PScalar< RComplex< ILattice<REAL32,INNER_LEN> > >, Ns> > LatticeSpinMatrixF;
typedef OLattice< PSpinMatrix< PScalar< RComplex< ILattice<REAL32,INNER_LEN> > >, (Ns>>1) > > LatticeHalfSpinMatrixF;
typedef OLattice< PSpinVector< PScalar< RComplex< ILattice<REAL32,INNER_LEN> > >, Ns> > LatticeSpinVectorF;

typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc> > > LatticeColorMatrixF;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 3> > > LatticeColorMatrixF3;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 2> > > LatticeColorMatrixF2;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, 1> > > LatticeColorMatrixF1;

typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, Nc> > > LatticeColorVectorF;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 3> > > LatticeColorVectorF3;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 2> > > LatticeColorVectorF2;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL32,INNER_LEN> >, 1> > > LatticeColorVectorF1;

typedef OLattice< PScalar< PScalar< RComplex< ILattice<REAL32,INNER_LEN> > > > > LatticeComplexF;
typedef OLattice< PScalar< PScalar< RScalar< ILattice<REAL32,INNER_LEN> > > > > LatticeRealF;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, Nc>, 4> > DiracFermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 3>, 4> > DiracFermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 2>, 4> > DiracFermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 1>, 4> > DiracFermionF1;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, Nc>, 1> > StaggeredFermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 3>, 1> > StaggeredFermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 2>, 1> > StaggeredFermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 1>, 1> > StaggeredFermionF1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, Nc>, 4> > DiracPropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 3>, 4> > DiracPropagatorF3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 2>, 4> > DiracPropagatorF2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 1>, 4> > DiracPropagatorF1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, Nc>, 1> > StaggeredPropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 3>, 1> > StaggeredPropagatorF3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 2>, 1> > StaggeredPropagatorF2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 1>, 1> > StaggeredPropagatorF1;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, Nc>, Ns> > FermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 3>, Ns> > FermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 2>, Ns> > FermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 1>, Ns> > FermionF1;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, Nc>, Ns>>1 > > HalfFermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 3>, Ns>>1 > > HalfFermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 2>, Ns>>1 > > HalfFermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL32> >, 1>, Ns>>1 > > HalfFermionF1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, Nc>, Ns> > PropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, Nc>, (Ns>>1) > > HalfPropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 3>, Ns> > PropagatorF3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 2>, Ns> > PropagatorF2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, 1>, Ns> > PropagatorF1;

typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<REAL32> > >, Ns> > SpinMatrixF;
typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<REAL32> > >, (Ns>>1) > > HalfSpinMatrixF;
typedef OScalar< PSpinVector< PScalar< RComplex< IScalar<REAL32> > >, Ns> > SpinVectorF;

typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL32> >, Nc> > > ColorMatrixF;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL32> >, 3> > > ColorMatrixF3;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL32> >, 2> > > ColorMatrixF2;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL32> >, 1> > > ColorMatrixF1;

typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL32> >, Nc> > > ColorVectorF;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL32> >, 3> > > ColorVectorF3;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL32> >, 2> > > ColorVectorF2;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL32> >, 1> > > ColorVectorF1;

typedef OScalar< PScalar< PScalar< RComplex< IScalar<REAL32> > > > > ComplexF;
typedef OScalar< PScalar< PScalar< RScalar< IScalar<REAL32> > > > > RealF;

//baryonblocks 
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, Nc>, Ns>, Nc>, Ns> > BaryonblockF;
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL32> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfBaryonblockF;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, Ns >, Nc>, Ns > > LatticeBaryonblockF;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfBaryonblockF;

//fourquark objects
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, Ns>, Nc>, Ns> > FourquarkF;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, Ns>, Nc>, Ns> > LatticeFourquarkF;
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfFourquarkF;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL32,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfFourquarkF;

// REAL64 types
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, 4> > LatticeDiracFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, 4> > LatticeDiracFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, 4> > LatticeDiracFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, 4> > LatticeDiracFermionD1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, 1> > LatticeStaggeredFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, 1> > LatticeStaggeredFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, 1> > LatticeStaggeredFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, 1> > LatticeStaggeredFermionD1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, 4> > LatticeDiracPropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, 4> > LatticeDiracPropagatorD3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, 4> > LatticeDiracPropagatorD2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, 4> > LatticeDiracPropagatorD1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, 1> > LatticeStaggeredPropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, 1> > LatticeStaggeredPropagatorD3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, 1> > LatticeStaggeredPropagatorD2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, 1> > LatticeStaggeredPropagatorD1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, Ns> > LatticeFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, Ns> > LatticeFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, Ns> > LatticeFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, Ns> > LatticeFermionD1;

typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, Ns>>1 > > LatticeHalfFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, Ns>>1 > > LatticeHalfFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, Ns>>1 > > LatticeHalfFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, Ns>>1 > > LatticeHalfFermionD1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, Ns> > LatticePropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, (Ns>>1) > > LatticePropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 3>, Ns> > LatticePropagatorD3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 2>, Ns> > LatticePropagatorD2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 1>, Ns> > LatticePropagatorD1;

typedef OLattice< PSpinMatrix< PScalar< RComplex< ILattice<REAL64,INNER_LEN> > >, Ns> > LatticeSpinMatrixD;
typedef OLattice< PSpinMatrix< PScalar< RComplex< ILattice<REAL64,INNER_LEN> > >, (Ns>>1) > > LatticeHalfSpinMatrixD;
typedef OLattice< PSpinVector< PScalar< RComplex< ILattice<REAL64,INNER_LEN> > >, Ns> > LatticeSpinVectorD;

typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc> > > LatticeColorMatrixD;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 3> > > LatticeColorMatrixD3;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 2> > > LatticeColorMatrixD2;
typedef OLattice< PScalar< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, 1> > > LatticeColorMatrixD1;

typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, Nc> > > LatticeColorVectorD;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 3> > > LatticeColorVectorD3;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 2> > > LatticeColorVectorD2;
typedef OLattice< PScalar< PColorVector< RComplex< ILattice<REAL64,INNER_LEN> >, 1> > > LatticeColorVectorD1;

typedef OLattice< PScalar< PScalar< RComplex< ILattice<REAL64,INNER_LEN> > > > > LatticeComplexD;
typedef OLattice< PScalar< PScalar< RScalar< ILattice<REAL64,INNER_LEN> > > > > LatticeRealD;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, Nc>, 4> > DiracFermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 3>, 4> > DiracFermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 2>, 4> > DiracFermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 1>, 4> > DiracFermionD1;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, Nc>, 1> > StaggeredFermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 3>, 1> > StaggeredFermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 2>, 1> > StaggeredFermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 1>, 1> > StaggeredFermionD1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, 4> > DiracPropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 3>, 4> > DiracPropagatorD3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 2>, 4> > DiracPropagatorD2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 1>, 4> > DiracPropagatorD1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, 1> > StaggeredPropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 3>, 1> > StaggeredPropagatorD3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 2>, 1> > StaggeredPropagatorD2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 1>, 1> > StaggeredPropagatorD1;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, Nc>, Ns> > FermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 3>, Ns> > FermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 2>, Ns> > FermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 1>, Ns> > FermionD1;

typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, Nc>, Ns>>1 > > HalfFermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 3>, Ns>>1 > > HalfFermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 2>, Ns>>1 > > HalfFermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex< IScalar<REAL64> >, 1>, Ns>>1 > > HalfFermionD1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, Ns> > PropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, (Ns>>1) > > HalfPropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 3>, Ns> > PropagatorD3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 2>, Ns> > PropagatorD2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, 1>, Ns> > PropagatorD1;

typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<REAL64> > >, Ns> > SpinMatrixD;
typedef OScalar< PSpinMatrix< PScalar< RComplex< IScalar<REAL64> > >, (Ns>>1) > > HalfSpinMatrixD;
typedef OScalar< PSpinVector< PScalar< RComplex< IScalar<REAL64> > >, Ns> > SpinVectorD;

typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL64> >, Nc> > > ColorMatrixD;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL64> >, 3> > > ColorMatrixD3;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL64> >, 2> > > ColorMatrixD2;
typedef OScalar< PScalar< PColorMatrix< RComplex< IScalar<REAL64> >, 1> > > ColorMatrixD1;

typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL64> >, Nc> > > ColorVectorD;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL64> >, 3> > > ColorVectorD3;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL64> >, 2> > > ColorVectorD2;
typedef OScalar< PScalar< PColorVector< RComplex< IScalar<REAL64> >, 1> > > ColorVectorD1;

typedef OScalar< PScalar< PScalar< RComplex< IScalar<REAL64> > > > > ComplexD;
typedef OScalar< PScalar< PScalar< RScalar< IScalar<REAL64> > > > > RealD;

//baryonblocks
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, Ns>, Nc>, Ns> > BaryonblockD;
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfBaryonblockD;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, Ns>, Nc>, Ns> > LatticeBaryonblockD;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex< ILattice<REAL64,INNER_LEN> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfBaryonblockD;

//fourquark objects
//fourquark operators
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, Ns>, Nc>, Ns> > FourquarkD;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, Ns>, Nc>, Ns> > LatticeFourquarkD;
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfFourquarkD;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex< IScalar<REAL64> >, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfFourquarkD;

// Equivalent names
typedef Integer  Int;

typedef RealF  Real32;
typedef RealD  Real64;
typedef ComplexF  Complex32;
typedef ComplexD  Complex64;

typedef LatticeInteger  LatticeInt;


/*! @} */   // end of group defs

} // namespace QDP

#endif
