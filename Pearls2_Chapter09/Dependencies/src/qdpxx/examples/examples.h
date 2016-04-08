// -*- C++ -*-
// $Id: examples.h,v 1.8 2007-02-24 01:00:29 bjoo Exp $
//
// Include file for test suite

#include "qdp.h"

using namespace QDP;

#if defined(QDP_DEBUG_MEMORY)
#define START_CODE() QDP::Allocator::theQDPAllocator::Instance().pushFunc(__func__, __LINE__)
#define END_CODE()   QDP::Allocator::theQDPAllocator::Instance().popFunc()

#else
#define START_CODE() QDP_PUSH_PROFILE(QDP::getProfileLevel())
#define END_CODE()   QDP_POP_PROFILE()

#endif


enum Reunitarize {REUNITARIZE, REUNITARIZE_ERROR, REUNITARIZE_LABEL};


void reunit(LatticeColorMatrix& xa);
void reunit(LatticeColorMatrix& xa, LatticeBoolean& bad, int& numbad, enum Reunitarize ruflag);


#ifdef QDP_USE_LIBXML2
void junk(XMLWriter&, LatticeColorMatrix& b3, const LatticeColorMatrix& b1, const LatticeColorMatrix& b2, const Subset& s);
#endif

void MesPlq(const multi1d<LatticeColorMatrix>& u, Double& w_plaq, Double& s_plaq, 
	    Double& t_plaq, Double& link);
void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi1d< multi1d<Real> >& meson_propagator, 
	    const multi1d<int>& t_source, int j_decay);
void baryon(const LatticePropagator& quark_propagator, 
	    multi1d< multi1d<Complex> >& barprop, 
	    const multi1d<int>& t_source, int j_decay, int bc_spec);
void dslash_2d_plus(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi,
	    int cb);
void dslash(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi,
	    int isign, int cb);

void dslash2(LatticeFermion& chi, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& psi,
	    int isign, int cb);

#ifdef QDP_USE_LIBXML2
void FormFac(const multi1d<LatticeColorMatrix>& u, const LatticePropagator& quark_propagator,
	     const LatticePropagator& seq_quark_prop, const multi1d<int>& t_source, 
	     int t_sink, int j_decay, XMLWriter& xml);
#endif

void expm12(LatticeColorMatrix& a);

void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g);

void taproj(LatticeColorMatrix& a);

void polylp(const multi1d<LatticeColorMatrix>& u, DComplex& poly_loop, int mu);
