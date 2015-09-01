// -*- C++ -*-

/*! \file
 * \brief Type definitions
 */

#ifndef QDP_SCALARSITE_DEFS_H
#define QDP_SCALARSITE_DEFS_H

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

//! Gamma matrices are conveniently defined for this Ns
typedef GammaTypeDP<Ns> GammaDP;


// Aliases for a scalar architecture

// Fixed fermion type
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, 4> > LatticeDiracFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, 1> > LatticeStaggeredFermion;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, 4> > LatticeDiracPropagator;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, 1> > LatticeStaggeredPropagator;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL>, Nc>, 4> > DiracFermion;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL>, Nc>, 1> > StaggeredFermion;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, 4> > DiracPropagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, 1> > StaggeredPropagator;

// Floating aliases
typedef OLattice< PScalar< PColorVector< RComplex<REAL>, Nc> > > LatticeColorVector;
typedef OLattice< PSpinVector< PScalar< RComplex<REAL> >, Ns> > LatticeSpinVector;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL>, Nc> > > LatticeColorMatrix;
typedef OLattice< PSpinMatrix< PScalar< RComplex<REAL> >, Ns> > LatticeSpinMatrix;
typedef OLattice< PSpinMatrix< PScalar< RComplex<REAL> >, (Ns>>1) > > LatticeHalfSpinMatrix;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns> > LatticeFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, (Ns>>1) > > LatticeHalfFermion;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns> > LatticePropagator;
typedef OLattice< PScalar< PScalar< RComplex<REAL> > > > LatticeComplex;

typedef OLattice< PScalar< PSeed < RScalar<INTEGER32> > > > LatticeSeed;
typedef OLattice< PScalar< PScalar< RScalar<INTEGER32> > > > LatticeInteger;
typedef OLattice< PScalar< PScalar< RScalar<REAL> > > > LatticeReal;
typedef OLattice< PScalar< PScalar< RScalar<DOUBLE> > > > LatticeDouble;
typedef OLattice< PScalar< PScalar< RScalar<LOGICAL> > > > LatticeBoolean;

typedef OScalar< PScalar< PColorVector< RComplex<REAL>, Nc> > > ColorVector;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL>, Nc> > > ColorMatrix;
typedef OScalar< PSpinVector< PScalar< RComplex<REAL> >, Ns> > SpinVector;
typedef OScalar< PSpinMatrix< PScalar< RComplex<REAL> >, Ns> > SpinMatrix;
typedef OScalar< PSpinMatrix< PScalar< RComplex<REAL> >, (Ns>>1) > > HalfSpinMatrix;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns> > Fermion;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL>, Nc>, (Ns>>1) > > HalfFermion;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns> > Propagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, (Ns>>1) > > HalfPropagator;
typedef OScalar< PScalar< PScalar< RComplex<REAL> > > > Complex;

typedef OScalar< PScalar< PSeed< RScalar<INTEGER32> > > > Seed;
typedef OScalar< PScalar< PScalar< RScalar<INTEGER32> > > > Integer;
typedef OScalar< PScalar< PScalar< RScalar<REAL> > > > Real;
typedef OScalar< PScalar< PScalar< RScalar<DOUBLE> > > > Double;
typedef OScalar< PScalar< PScalar< RScalar<LOGICAL> > > > Boolean;

//baryonblocks
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns>, Nc>, Ns> > Baryonblock;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns>, Nc>, Ns> > LatticeBaryonblock;
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfBaryonblock;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfBaryonblock;

//fourquark tensors
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns>, Nc>, Ns> > Fourquark;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns>, Nc>, Ns> > LatticeFourquark;
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfFourquark;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfFourquark;

// Other useful names
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns> > ColorVectorSpinVector;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns> > ColorMatrixSpinMatrix;

// Level below outer for internal convenience
typedef PScalar< PScalar< RScalar<REAL> > > IntReal;
typedef PScalar< PScalar< RScalar<REAL32> > > IntReal32;
typedef PScalar< PScalar< RScalar<INTEGER32> > > IntInteger;
typedef PScalar< PScalar< RScalar<REAL64> > > IntReal64;
typedef PScalar< PScalar< RScalar<DOUBLE> > > IntDouble;
typedef PScalar< PScalar< RScalar<LOGICAL> > > IntBoolean;

// Odd-ball to support random numbers
typedef Real ILatticeReal;
typedef Seed ILatticeSeed;

// Floating aliases but (possibly) in a higher precision. 
// The SINGLE might in fact be the same as DOUBLE
// Fixed fermion type
typedef OLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, 4> > LatticeDDiracFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, 1> > LatticeDStaggeredFermion;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, 4> > LatticeDDiracPropagator;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, 1> > LatticeDStaggeredPropagator;

typedef OScalar< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, 4> > DDiracFermion;
typedef OScalar< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, 1> > DStaggeredFermion;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, 4> > DDiracPropagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, 1> > DStaggeredPropagator;

// Floating aliases
typedef OLattice< PScalar< PColorVector< RComplex<DOUBLE>, Nc> > > LatticeDColorVector;
typedef OLattice< PSpinVector< PScalar< RComplex<DOUBLE> >, Ns> > LatticeDSpinVector;
typedef OLattice< PScalar< PColorMatrix< RComplex<DOUBLE>, Nc> > > LatticeDColorMatrix;
typedef OLattice< PSpinMatrix< PScalar< RComplex<DOUBLE> >, Ns> > LatticeDSpinMatrix;
typedef OLattice< PSpinMatrix< PScalar< RComplex<DOUBLE> >, (Ns>>1) > > LatticeDHalfSpinMatrix;
typedef OLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, Ns> > LatticeDFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, (Ns>>1) > > LatticeDHalfFermion;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, Ns> > LatticeDPropagator;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, (Ns>>1) > > LatticeDHalfPropagator;
typedef OLattice< PScalar< PScalar< RComplex<DOUBLE> > > > LatticeDComplex;

typedef OScalar< PScalar< PColorVector< RComplex<DOUBLE>, Nc> > > DColorVector;
typedef OScalar< PScalar< PColorMatrix< RComplex<DOUBLE>, Nc> > > DColorMatrix;
typedef OScalar< PSpinVector< PScalar< RComplex<DOUBLE> >, Ns> > DSpinVector;
typedef OScalar< PSpinMatrix< PScalar< RComplex<DOUBLE> >, Ns> > DSpinMatrix;
typedef OScalar< PSpinMatrix< PScalar< RComplex<DOUBLE> >, (Ns>>1) > > DHalfSpinMatrix;
typedef OScalar< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, Ns> > DFermion;
typedef OScalar< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, (Ns>>1) > > DHalfFermion;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, Ns> > DPropagator;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, (Ns>>1) > > DHalfPropagator;
typedef OScalar< PScalar< PScalar< RComplex<DOUBLE> > > > DComplex;


// Floating precision, but specific to a fixed color or spin
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, 4> > LatticeDiracFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, 4> > LatticeDiracFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, 4> > LatticeDiracFermion1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, 1> > LatticeStaggeredFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, 1> > LatticeStaggeredFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, 1> > LatticeStaggeredFermion1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 3>, 4> > LatticeDiracPropagator3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 2>, 4> > LatticeDiracPropagator2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 1>, 4> > LatticeDiracPropagator1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 3>, 1> > LatticeStaggerdPropagator3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 2>, 1> > LatticeStaggerdPropagator2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 1>, 1> > LatticeStaggerdPropagator1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, Ns> > LatticeFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, Ns> > LatticeFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, Ns> > LatticeFermion1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, (Ns>>1) > > LatticeHalfFermion3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, (Ns>>1) > > LatticeHalfFermion2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, (Ns>>1) > > LatticeHalfFermion1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 3>, Ns> > LatticePropagator3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 2>, Ns> > LatticePropagator2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 1>, Ns> > LatticePropagator1;

typedef OLattice< PScalar< PColorMatrix< RComplex<REAL>, 3> > > LatticeColorMatrix3;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL>, 2> > > LatticeColorMatrix2;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL>, 1> > > LatticeColorMatrix1;

typedef OLattice< PScalar< PColorVector< RComplex<REAL>, 3> > > LatticeColorVector3;
typedef OLattice< PScalar< PColorVector< RComplex<REAL>, 2> > > LatticeColorVector2;
typedef OLattice< PScalar< PColorVector< RComplex<REAL>, 1> > > LatticeColorVector1;

typedef OScalar< PScalar< PColorMatrix< RComplex<REAL>, 3> > > ColorMatrix3;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL>, 2> > > ColorMatrix2;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL>, 1> > > ColorMatrix1;

typedef OScalar< PScalar< PColorVector< RComplex<REAL>, 3> > > ColorVector3;
typedef OScalar< PScalar< PColorVector< RComplex<REAL>, 2> > > ColorVector2;
typedef OScalar< PScalar< PColorVector< RComplex<REAL>, 1> > > ColorVector1;

//
// Fixed precision
//
// REAL32 types
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, 4> > LatticeDiracFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, 4> > LatticeDiracFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, 4> > LatticeDiracFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, 4> > LatticeDiracFermionF1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, 1> > LatticeStaggeredFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, 1> > LatticeStaggeredFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, 1> > LatticeStaggeredFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, 1> > LatticeStaggeredFermionF1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, 4> > LatticeDiracPropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, 4> > LatticeDiracPropagatorF3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, 4> > LatticeDiracPropagatorF2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, 4> > LatticeDiracPropagatorF1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, 1> > LatticeStaggeredPropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, 1> > LatticeStaggeredPropagatorF3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, 1> > LatticeStaggeredPropagatorF2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, 1> > LatticeStaggeredPropagatorF1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, Ns> > LatticeFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, Ns> > LatticeFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, Ns> > LatticeFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, Ns> > LatticeFermionF1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, (Ns>>1) > > LatticeHalfFermionF;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, (Ns>>1) > > LatticeHalfFermionF3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, (Ns>>1) > > LatticeHalfFermionF2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, (Ns>>1) > > LatticeHalfFermionF1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns> > LatticePropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, (Ns>>1) > > LatticeHalfPropagatorF;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, Ns> > LatticePropagatorF3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, Ns> > LatticePropagatorF2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, Ns> > LatticePropagatorF1;

typedef OLattice< PSpinMatrix< PScalar< RComplex<REAL32> >, Ns> > LatticeSpinMatrixF;
typedef OLattice< PSpinMatrix< PScalar< RComplex<REAL32> >, (Ns>>1) > > LatticeHalfSpinMatrixF;
typedef OLattice< PSpinVector< PScalar< RComplex<REAL32> >, Ns> > LatticeSpinVectorF;

typedef OLattice< PScalar< PColorMatrix< RComplex<REAL32>, Nc> > > LatticeColorMatrixF;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL32>, 3> > > LatticeColorMatrixF3;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL32>, 2> > > LatticeColorMatrixF2;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL32>, 1> > > LatticeColorMatrixF1;

typedef OLattice< PScalar< PColorVector< RComplex<REAL32>, Nc> > > LatticeColorVectorF;
typedef OLattice< PScalar< PColorVector< RComplex<REAL32>, 3> > > LatticeColorVectorF3;
typedef OLattice< PScalar< PColorVector< RComplex<REAL32>, 2> > > LatticeColorVectorF2;
typedef OLattice< PScalar< PColorVector< RComplex<REAL32>, 1> > > LatticeColorVectorF1;

typedef OLattice< PScalar< PScalar< RComplex<REAL32> > > > LatticeComplexF;
typedef OLattice< PScalar< PScalar< RScalar<REAL32> > > > LatticeRealF;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, 4> > DiracFermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 3>, 4> > DiracFermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 2>, 4> > DiracFermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 1>, 4> > DiracFermionF1;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, 1> > StaggeredFermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 3>, 1> > StaggeredFermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 2>, 1> > StaggeredFermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 1>, 1> > StaggeredFermionF1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, 4> > DiracPropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, 4> > DiracPropagatorF3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, 4> > DiracPropagatorF2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, 4> > DiracPropagatorF1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, 1> > StaggeredPropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, 1> > StaggeredPropagatorF3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, 1> > StaggeredPropagatorF2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, 1> > StaggeredPropagatorF1;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, Ns> > FermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 3>, Ns> > FermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 2>, Ns> > FermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 1>, Ns> > FermionF1;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, (Ns>>1) > > HalfFermionF;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 3>, (Ns>>1) > > HalfFermionF3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 2>, (Ns>>1) > > HalfFermionF2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL32>, 1>, (Ns>>1) > > HalfFermionF1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns> > PropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, (Ns>>1) > > HalfPropagatorF;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, Ns> > PropagatorF3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, Ns> > PropagatorF2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, Ns> > PropagatorF1;

typedef OScalar< PSpinMatrix< PScalar< RComplex<REAL32> >, Ns> > SpinMatrixF;
typedef OScalar< PSpinMatrix< PScalar< RComplex<REAL32> >, (Ns>>1) > > HalfSpinMatrixF;
typedef OScalar< PSpinVector< PScalar< RComplex<REAL32> >, Ns> > SpinVectorF;

typedef OScalar< PScalar< PColorMatrix< RComplex<REAL32>, Nc> > > ColorMatrixF;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL32>, 3> > > ColorMatrixF3;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL32>, 2> > > ColorMatrixF2;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL32>, 1> > > ColorMatrixF1;

typedef OScalar< PScalar< PColorVector< RComplex<REAL32>, Nc> > > ColorVectorF;
typedef OScalar< PScalar< PColorVector< RComplex<REAL32>, 3> > > ColorVectorF3;
typedef OScalar< PScalar< PColorVector< RComplex<REAL32>, 2> > > ColorVectorF2;
typedef OScalar< PScalar< PColorVector< RComplex<REAL32>, 1> > > ColorVectorF1;

typedef OScalar< PScalar< PScalar< RComplex<REAL32> > > > ComplexF;
typedef OScalar< PScalar< PScalar< RScalar<REAL32> > > > RealF;

//baryonblocks
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns>, Nc>, Ns> > BaryonblockF;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns>, Nc>, Ns> > LatticeBaryonblockF;
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfBaryonblockF;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfBaryonblockF;

//fourquark objects
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns>, Nc>, Ns> > FourquarkF;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns>, Nc>, Ns> > LatticeFourquarkF;
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfFourquarkF;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfFourquarkF;

// REAL64 types
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, 4> > LatticeDiracFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, 4> > LatticeDiracFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, 4> > LatticeDiracFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, 4> > LatticeDiracFermionD1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, 1> > LatticeStaggeredFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, 1> > LatticeStaggeredFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, 1> > LatticeStaggeredFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, 1> > LatticeStaggeredFermionD1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, 4> > LatticeDiracPropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, 4> > LatticeDiracPropagatorD3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, 4> > LatticeDiracPropagatorD2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, 4> > LatticeDiracPropagatorD1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, 1> > LatticeStaggeredPropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, 1> > LatticeStaggeredPropagatorD3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, 1> > LatticeStaggeredPropagatorD2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, 1> > LatticeStaggeredPropagatorD1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, Ns> > LatticeFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, Ns> > LatticeFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, Ns> > LatticeFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, Ns> > LatticeFermionD1;

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, (Ns>>1) > > LatticeHalfFermionD;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, (Ns>>1) > > LatticeHalfFermionD3;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, (Ns>>1) > > LatticeHalfFermionD2;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, (Ns>>1) > > LatticeHalfFermionD1;

typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns> > LatticePropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, (Ns>>1) > > LatticeHalfPropagatorD;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, Ns> > LatticePropagatorD3;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, Ns> > LatticePropagatorD2;
typedef OLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, Ns> > LatticePropagatorD1;

typedef OLattice< PSpinMatrix< PScalar< RComplex<REAL64> >, Ns> > LatticeSpinMatrixD;
typedef OLattice< PSpinMatrix< PScalar< RComplex<REAL64> >, (Ns>>1) > > LatticeHalfSpinMatrixD;
typedef OLattice< PSpinVector< PScalar< RComplex<REAL64> >, Ns> > LatticeSpinVectorD;

typedef OLattice< PScalar< PColorMatrix< RComplex<REAL64>, Nc> > > LatticeColorMatrixD;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL64>, 3> > > LatticeColorMatrixD3;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL64>, 2> > > LatticeColorMatrixD2;
typedef OLattice< PScalar< PColorMatrix< RComplex<REAL64>, 1> > > LatticeColorMatrixD1;

typedef OLattice< PScalar< PColorVector< RComplex<REAL64>, Nc> > > LatticeColorVectorD;
typedef OLattice< PScalar< PColorVector< RComplex<REAL64>, 3> > > LatticeColorVectorD3;
typedef OLattice< PScalar< PColorVector< RComplex<REAL64>, 2> > > LatticeColorVectorD2;
typedef OLattice< PScalar< PColorVector< RComplex<REAL64>, 1> > > LatticeColorVectorD1;

typedef OLattice< PScalar< PScalar< RComplex<REAL64> > > > LatticeComplexD;
typedef OLattice< PScalar< PScalar< RScalar<REAL64> > > > LatticeRealD;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, 4> > DiracFermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 3>, 4> > DiracFermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 2>, 4> > DiracFermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 1>, 4> > DiracFermionD1;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, 1> > StaggeredFermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 3>, 1> > StaggeredFermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 2>, 1> > StaggeredFermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 1>, 1> > StaggeredFermionD1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, 4> > DiracPropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, 4> > DiracPropagatorD3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, 4> > DiracPropagatorD2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, 4> > DiracPropagatorD1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, 1> > StaggeredPropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, 1> > StaggeredPropagatorD3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, 1> > StaggeredPropagatorD2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, 1> > StaggeredPropagatorD1;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, Ns> > FermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 3>, Ns> > FermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 2>, Ns> > FermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 1>, Ns> > FermionD1;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, (Ns>>1) > > HalfFermionD;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 3>, (Ns>>1) > > HalfFermionD3;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 2>, (Ns>>1) > > HalfFermionD2;
typedef OScalar< PSpinVector< PColorVector< RComplex<REAL64>, 1>, (Ns>>1) > > HalfFermionD1;

typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns> > PropagatorD;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, Ns> > PropagatorD3;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, Ns> > PropagatorD2;
typedef OScalar< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, Ns> > PropagatorD1;

typedef OScalar< PSpinMatrix< PScalar< RComplex<REAL64> >, Ns> > SpinMatrixD;
typedef OScalar< PSpinVector< PScalar< RComplex<REAL64> >, Ns> > SpinVectorD;

typedef OScalar< PScalar< PColorMatrix< RComplex<REAL64>, Nc> > > ColorMatrixD;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL64>, 3> > > ColorMatrixD3;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL64>, 2> > > ColorMatrixD2;
typedef OScalar< PScalar< PColorMatrix< RComplex<REAL64>, 1> > > ColorMatrixD1;

typedef OScalar< PScalar< PColorVector< RComplex<REAL64>, Nc> > > ColorVectorD;
typedef OScalar< PScalar< PColorVector< RComplex<REAL64>, 3> > > ColorVectorD3;
typedef OScalar< PScalar< PColorVector< RComplex<REAL64>, 2> > > ColorVectorD2;
typedef OScalar< PScalar< PColorVector< RComplex<REAL64>, 1> > > ColorVectorD1;

typedef OScalar< PScalar< PScalar< RComplex<REAL64> > > > ComplexD;
typedef OScalar< PScalar< PScalar< RScalar<REAL64> > > > RealD;

//baryonblocks
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns>, Nc>, Ns> > BaryonblockD;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns>, Nc>, Ns> > LatticeBaryonblockD;
typedef OScalar< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfBaryonblockD;
typedef OLattice< PSpinMatrix < PColorVector< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfBaryonblockD;

//fourquark operators
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns>, Nc>, Ns> > FourquarkD;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns>, Nc>, Ns> > LatticeFourquarkD;
typedef OScalar< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > HalfFourquarkD;
typedef OLattice< PSpinMatrix < PColorMatrix< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, (Ns>>1) >, Nc>, (Ns>>1) > > LatticeHalfFourquarkD;

// Equivalent names
typedef Integer  Int;

typedef RealF  Real32;
typedef RealD  Real64;
typedef ComplexF  Complex32;
typedef ComplexD  Complex64;

typedef LatticeInteger  LatticeInt;


//--------------------------------------------------------------------------------
// Aliases for a scalar architecture

// Fixed fermion type
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, 4> > SubLatticeDiracFermion;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, 1> > SubLatticeStaggeredFermion;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, 4> > SubLatticeDiracPropagator;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, 1> > SubLatticeStaggeredPropagator;

// Floating aliases
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL>, Nc> > > SubLatticeColorVector;
typedef OSubLattice< PSpinVector< PScalar< RComplex<REAL> >, Ns> > SubLatticeSpinVector;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL>, Nc> > > SubLatticeColorMatrix;
typedef OSubLattice< PSpinMatrix< PScalar< RComplex<REAL> >, Ns> > SubLatticeSpinMatrix;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns> > SubLatticeFermion;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, (Ns>>1) > > SubLatticeHalfFermion;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, Nc>, Ns> > SubLatticePropagator;
typedef OSubLattice< PScalar< PScalar< RComplex<REAL> > > > SubLatticeComplex;

typedef OSubLattice< PScalar< PSeed < RScalar<INTEGER32> > > > SubLatticeSeed;
typedef OSubLattice< PScalar< PScalar< RScalar<INTEGER32> > > > SubLatticeInteger;
typedef OSubLattice< PScalar< PScalar< RScalar<REAL> > > > SubLatticeReal;
typedef OSubLattice< PScalar< PScalar< RScalar<DOUBLE> > > > SubLatticeDouble;
typedef OSubLattice< PScalar< PScalar< RScalar<LOGICAL> > > > SubLatticeBoolean;

// Floating aliases but (possibly) in a higher precision. 
// The SINGLE might in fact be the same as DOUBLE
// Fixed fermion type
typedef OSubLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, 4> > SubLatticeDDiracFermion;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, 1> > SubLatticeDStaggeredFermion;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, 4> > SubLatticeDDiracPropagator;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, 1> > SubLatticeDStaggeredPropagator;

// Floating aliases
typedef OSubLattice< PScalar< PColorVector< RComplex<DOUBLE>, Nc> > > SubLatticeDColorVector;
typedef OSubLattice< PSpinVector< PScalar< RComplex<DOUBLE> >, Ns> > SubLatticeDSpinVector;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<DOUBLE>, Nc> > > SubLatticeDColorMatrix;
typedef OSubLattice< PSpinMatrix< PScalar< RComplex<DOUBLE> >, Ns> > SubLatticeDSpinMatrix;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, Ns> > SubLatticeDFermion;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<DOUBLE>, Nc>, (Ns>>1) > > SubLatticeDHalfFermion;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<DOUBLE>, Nc>, Ns> > SubLatticeDPropagator;
typedef OSubLattice< PScalar< PScalar< RComplex<DOUBLE> > > > SubLatticeDComplex;

// Floating precision, but specific to a fixed color or spin
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, 4> > SubLatticeDiracFermion3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, 4> > SubLatticeDiracFermion2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, 4> > SubLatticeDiracFermion1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, 1> > SubLatticeStaggeredFermion3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, 1> > SubLatticeStaggeredFermion2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, 1> > SubLatticeStaggeredFermion1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 3>, 4> > SubLatticeDiracPropagator3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 2>, 4> > SubLatticeDiracPropagator2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 1>, 4> > SubLatticeDiracPropagator1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 3>, 1> > SubLatticeStaggerdPropagator3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 2>, 1> > SubLatticeStaggerdPropagator2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 1>, 1> > SubLatticeStaggerdPropagator1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, Ns> > SubLatticeFermion3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, Ns> > SubLatticeFermion2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, Ns> > SubLatticeFermion1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 3>, (Ns>>1) > > SubLatticeHalfFermion3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 2>, (Ns>>1) > > SubLatticeHalfFermion2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL>, 1>, (Ns>>1) > > SubLatticeHalfFermion1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 3>, Ns> > SubLatticePropagator3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 2>, Ns> > SubLatticePropagator2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL>, 1>, Ns> > SubLatticePropagator1;

typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL>, 3> > > SubLatticeColorMatrix3;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL>, 2> > > SubLatticeColorMatrix2;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL>, 1> > > SubLatticeColorMatrix1;

typedef OSubLattice< PScalar< PColorVector< RComplex<REAL>, 3> > > SubLatticeColorVector3;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL>, 2> > > SubLatticeColorVector2;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL>, 1> > > SubLatticeColorVector1;

//
// Fixed precision
//
// REAL32 types
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, 4> > SubLatticeDiracFermionF;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, 4> > SubLatticeDiracFermionF3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, 4> > SubLatticeDiracFermionF2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, 4> > SubLatticeDiracFermionF1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, 1> > SubLatticeStaggeredFermionF;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, 1> > SubLatticeStaggeredFermionF3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, 1> > SubLatticeStaggeredFermionF2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, 1> > SubLatticeStaggeredFermionF1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, 4> > SubLatticeDiracPropagatorF;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, 4> > SubLatticeDiracPropagatorF3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, 4> > SubLatticeDiracPropagatorF2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, 4> > SubLatticeDiracPropagatorF1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, 1> > SubLatticeStaggeredPropagatorF;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, 1> > SubLatticeStaggeredPropagatorF3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, 1> > SubLatticeStaggeredPropagatorF2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, 1> > SubLatticeStaggeredPropagatorF1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, Ns> > SubLatticeFermionF;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, Ns> > SubLatticeFermionF3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, Ns> > SubLatticeFermionF2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, Ns> > SubLatticeFermionF1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, Nc>, (Ns>>1) > > SubLatticeHalfFermionF;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 3>, (Ns>>1) > > SubLatticeHalfFermionF3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 2>, (Ns>>1) > > SubLatticeHalfFermionF2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL32>, 1>, (Ns>>1) > > SubLatticeHalfFermionF1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, Nc>, Ns> > SubLatticePropagatorF;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 3>, Ns> > SubLatticePropagatorF3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 2>, Ns> > SubLatticePropagatorF2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL32>, 1>, Ns> > SubLatticePropagatorF1;

typedef OSubLattice< PSpinMatrix< PScalar< RComplex<REAL32> >, Ns> > SubLatticeSpinMatrixF;
typedef OSubLattice< PSpinVector< PScalar< RComplex<REAL32> >, Ns> > SubLatticeSpinVectorF;

typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL32>, Nc> > > SubLatticeColorMatrixF;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL32>, 3> > > SubLatticeColorMatrixF3;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL32>, 2> > > SubLatticeColorMatrixF2;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL32>, 1> > > SubLatticeColorMatrixF1;

typedef OSubLattice< PScalar< PColorVector< RComplex<REAL32>, Nc> > > SubLatticeColorVectorF;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL32>, 3> > > SubLatticeColorVectorF3;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL32>, 2> > > SubLatticeColorVectorF2;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL32>, 1> > > SubLatticeColorVectorF1;

typedef OSubLattice< PScalar< PScalar< RComplex<REAL32> > > > SubLatticeComplexF;
typedef OSubLattice< PScalar< PScalar< RScalar<REAL32> > > > SubLatticeRealF;

// REAL64 types
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, 4> > SubLatticeDiracFermionD;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, 4> > SubLatticeDiracFermionD3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, 4> > SubLatticeDiracFermionD2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, 4> > SubLatticeDiracFermionD1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, 1> > SubLatticeStaggeredFermionD;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, 1> > SubLatticeStaggeredFermionD3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, 1> > SubLatticeStaggeredFermionD2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, 1> > SubLatticeStaggeredFermionD1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, 4> > SubLatticeDiracPropagatorD;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, 4> > SubLatticeDiracPropagatorD3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, 4> > SubLatticeDiracPropagatorD2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, 4> > SubLatticeDiracPropagatorD1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, 1> > SubLatticeStaggeredPropagatorD;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, 1> > SubLatticeStaggeredPropagatorD3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, 1> > SubLatticeStaggeredPropagatorD2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, 1> > SubLatticeStaggeredPropagatorD1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, Ns> > SubLatticeFermionD;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, Ns> > SubLatticeFermionD3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, Ns> > SubLatticeFermionD2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, Ns> > SubLatticeFermionD1;

typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, Nc>, (Ns>>1) > > SubLatticeHalfFermionD;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 3>, (Ns>>1) > > SubLatticeHalfFermionD3;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 2>, (Ns>>1) > > SubLatticeHalfFermionD2;
typedef OSubLattice< PSpinVector< PColorVector< RComplex<REAL64>, 1>, (Ns>>1) > > SubLatticeHalfFermionD1;

typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, Nc>, Ns> > SubLatticePropagatorD;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 3>, Ns> > SubLatticePropagatorD3;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 2>, Ns> > SubLatticePropagatorD2;
typedef OSubLattice< PSpinMatrix< PColorMatrix< RComplex<REAL64>, 1>, Ns> > SubLatticePropagatorD1;

typedef OSubLattice< PSpinMatrix< PScalar< RComplex<REAL64> >, Ns> > SubLatticeSpinMatrixD;
typedef OSubLattice< PSpinVector< PScalar< RComplex<REAL64> >, Ns> > SubLatticeSpinVectorD;

typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL64>, Nc> > > SubLatticeColorMatrixD;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL64>, 3> > > SubLatticeColorMatrixD3;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL64>, 2> > > SubLatticeColorMatrixD2;
typedef OSubLattice< PScalar< PColorMatrix< RComplex<REAL64>, 1> > > SubLatticeColorMatrixD1;

typedef OSubLattice< PScalar< PColorVector< RComplex<REAL64>, Nc> > > SubLatticeColorVectorD;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL64>, 3> > > SubLatticeColorVectorD3;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL64>, 2> > > SubLatticeColorVectorD2;
typedef OSubLattice< PScalar< PColorVector< RComplex<REAL64>, 1> > > SubLatticeColorVectorD1;

typedef OSubLattice< PScalar< PScalar< RComplex<REAL64> > > > SubLatticeComplexD;
typedef OSubLattice< PScalar< PScalar< RScalar<REAL64> > > > SubLatticeRealD;




/*! @} */   // end of group defs

} // namespace QDP

#endif
