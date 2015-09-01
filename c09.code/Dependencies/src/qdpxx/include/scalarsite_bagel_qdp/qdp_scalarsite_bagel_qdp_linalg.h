// -*- C++ -*-
// $Id: qdp_scalarsite_bagel_qdp_linalg.h,v 1.5 2008-05-09 17:45:20 bjoo Exp $

/*! @file
 * @brief Qcdoc optimizations
 *
 * Qcdoc version of optimized basic operations
 */

#ifndef QDP_SCALARSITE_BAGEL_QDP_LINALG_H
#define QDP_SCALARSITE_BAGEL_QDP_LINALG_H

namespace QDP {

/*! @defgroup optimizations  Optimizations
 *
 * Optimizations for basic QDP operations
 *
 * @{
 */

// Use this def just to safe some typing later on in the file



#include "bagel_qdp.h"

#if 1
typedef RComplex<BAGELQDPFloat>  RComplexFloat;


template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M*M) subset = s " << endl;
#endif
  
   if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_mm(resptr, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_mm(resptr, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                              UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(A*M) subset = s " << endl;
#endif
  
   if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_am(resptr, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_am(resptr, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                               Reference<
	                                 QDPType<
                                               PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                                       OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                                         > 
                                       >,
 
	                                UnaryNode< OpIdentity,
                                        Reference<
                                           QDPType<
                                               PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                                       OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                                           > 
                                         > 
                                        > 
                                      
	                   >,

	                   OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                     >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M*A) subset = s " << endl;
#endif
  
   if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_ma(resptr, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_ma(resptr, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}


template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(A*A) subset = s " << endl;
#endif
  BAGELQDPFloat one_minus_i[2] QDP_ALIGN16;
  one_minus_i[0] = (BAGELQDPFloat)1;
  one_minus_i[1] = (BAGELQDPFloat)(-1);

   if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_aa(resptr, lptr, rptr, num_sites, (unsigned long)one_minus_i);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_aa(resptr, lptr, rptr, one_site, (unsigned long)one_minus_i);

     }
   }

}


  // += 
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M*M) subset = s " << endl;
#endif

  BAGELQDPFloat plus_one[2] QDP_ALIGN16;
  plus_one[0] = (BAGELQDPFloat)1;
  plus_one[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_su3_mm_peq(resptr, plus_one, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_mm_peq(resptr, plus_one, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(A*M) subset = s " << endl;
#endif

  BAGELQDPFloat plus_one[2] QDP_ALIGN16;
  plus_one[0] = (BAGELQDPFloat)1;
  plus_one[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_am_peq(resptr, plus_one, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_am_peq(resptr, plus_one, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	      Reference< QDPType<
                  PScalar<PColorMatrix<RComplexFloat, 3> >, 
	          OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
	        > 
	      >,
 
	      UnaryNode< OpIdentity,
	      Reference< QDPType<
	          PScalar<PColorMatrix<RComplexFloat, 3> >, 
	          OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
	        > 
	      > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M*A) subset = s " << endl;
#endif
  BAGELQDPFloat plus_one[2] QDP_ALIGN16;
  plus_one[0] = (BAGELQDPFloat)1;
  plus_one[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_ma_peq(resptr, plus_one, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_ma_peq(resptr, plus_one, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}


template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(A*A) subset = s " << endl;
#endif
  BAGELQDPFloat one_minus_i[2] QDP_ALIGN16;
  one_minus_i[0] = (BAGELQDPFloat)1;
  one_minus_i[1] = (BAGELQDPFloat)(-1);

  BAGELQDPFloat plus_one[2] QDP_ALIGN16;
  plus_one[0] = (BAGELQDPFloat)1;
  plus_one[1] = (BAGELQDPFloat)0;

  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_aa_peq(resptr,  plus_one, lptr, rptr, num_sites, (unsigned long)one_minus_i);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_aa_peq(resptr, plus_one, lptr, rptr, one_site, (unsigned long)one_minus_i);

     }
   }

}



  //  -= 
  
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M*M) subset = s " << endl;
#endif

  BAGELQDPFloat minus_one[2] QDP_ALIGN16;
  minus_one[0] = (BAGELQDPFloat)-1;
  minus_one[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_su3_mm_peq(resptr, minus_one, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_mm_peq(resptr, minus_one, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                              UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(A*M) subset = s " << endl;
#endif

  BAGELQDPFloat minus_one[2] QDP_ALIGN16;
  minus_one[0] = (BAGELQDPFloat)-1;
  minus_one[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_am_peq(resptr, minus_one, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_am_peq(resptr, minus_one, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                               Reference<
	                                 QDPType<
                                               PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                                       OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                                         > 
                                       >,
 
	                                UnaryNode< OpIdentity,
                                        Reference<
                                           QDPType<
                                               PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                                       OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                                           > 
                                         > 
                                        > 
                                      
	                   >,

	                   OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                     >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M*A) subset = s " << endl;
#endif
  BAGELQDPFloat minus_one[2] QDP_ALIGN16;
  minus_one[0] = (BAGELQDPFloat)-1;
  minus_one[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_ma_peq(resptr, minus_one, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_ma_peq(resptr, minus_one, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}


template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >, 
	      UnaryNode<OpIdentity, Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

#ifdef DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(A*A) subset = s " << endl;
#endif
  BAGELQDPFloat one_minus_i[2] QDP_ALIGN16;
  one_minus_i[0] = (BAGELQDPFloat)1;
  one_minus_i[1] = (BAGELQDPFloat)(-1);

  BAGELQDPFloat minus_one[2] QDP_ALIGN16;
  minus_one[0] = (BAGELQDPFloat)-1;
  minus_one[1] = (BAGELQDPFloat)0;

  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));


     qdp_su3_aa_peq(resptr, minus_one, lptr, rptr, num_sites, (unsigned long)one_minus_i);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_aa_peq(resptr, minus_one, lptr, rptr, one_site, (unsigned long)one_minus_i);

     }
   }

}




  // += *
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 

       const QDPExpr<
                    BinaryNode< OpMultiply, 
	              BinaryNode<OpMultiply, 
                        Reference<
                          QDPType<
                            PScalar<PScalar<RScalar<BAGELQDPFloat> > >,
	                    OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > >
                          > 
                        >,
                        Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
	             >, 
	             Reference<
	               QDPType<
                         PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                 OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                       > 
	             >
                    >,
	            OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
              >& rhs,

	      const Subset& s) {

  typedef OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > > F;

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;

  typedef BinaryNode<OpMultiply, 
                        Reference<
                          QDPType<
                            PScalar<PScalar<RScalar<BAGELQDPFloat> > >,
	                    OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > >
                          > 
                        >,
                        Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
    > BN;


  const BN& node = static_cast<const BN&>(rhs.expression().left());
  const F& scal = static_cast<const F&>(node.left());
  const C& l = static_cast<const C&>(node.right());
  const C& r = static_cast<const C&>(rhs.expression().right());

#if DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M += alpha*M*M 2 ) subset = s " << endl;
#endif

  BAGELQDPFloat scalar[2] QDP_ALIGN16;
  scalar[0] = (BAGELQDPFloat)(scal.elem().elem().elem().elem());
  scalar[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_su3_mm_peq(resptr, scalar, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_mm_peq(resptr, scalar, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}

  // +-= *
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 

       const QDPExpr<
                    BinaryNode< OpMultiply, 
	              BinaryNode<OpMultiply, 
                        Reference<
                          QDPType<
                            PScalar<PScalar<RScalar<BAGELQDPFloat> > >,
	                    OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > >
                          > 
                        >,
                        Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
	             >, 
	             Reference<
	               QDPType<
                         PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                 OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                       > 
	             >
                    >,
	            OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
              >& rhs,

	      const Subset& s) {

  typedef OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > > F;

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;

  typedef BinaryNode<OpMultiply, 
                        Reference<
                          QDPType<
                            PScalar<PScalar<RScalar<BAGELQDPFloat> > >,
	                    OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > >
                          > 
                        >,
                        Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
    > BN;


  const BN& node = static_cast<const BN&>(rhs.expression().left());
  const F& scal = static_cast<const F&>(node.left());
  const C& l = static_cast<const C&>(node.right());
  const C& r = static_cast<const C&>(rhs.expression().right());

#if DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M -= alpha*M*M 2 ) subset = s " << endl;
#endif

  BAGELQDPFloat scalar[2] QDP_ALIGN16;
  scalar[0] = -(BAGELQDPFloat)(scal.elem().elem().elem().elem());
  scalar[1] = (BAGELQDPFloat)0;
  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;

     BAGELQDPFloat *resptr = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(start).elem().elem(0,0).real()));
     BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_su3_mm_peq(resptr, scalar, lptr, rptr, num_sites, (unsigned long)0);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       unsigned long one_site = 1;
       BAGELQDPFloat *resptr = &(d.elem(i).elem().elem(0,0).real());

       BAGELQDPFloat *lptr   = const_cast<BAGELQDPFloat*>(&(l.elem(i).elem().elem(0,0).real()));
       BAGELQDPFloat *rptr   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_su3_mm_peq(resptr, scalar, lptr, rptr, one_site, (unsigned long)0);

     }
   }

}


  // +=a * M
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 

	      const QDPExpr<
	              BinaryNode< OpMultiply, 
	                Reference<
                          QDPType<
                            PScalar<PScalar<RScalar<BAGELQDPFloat> > >,
	                    OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > >
                          > 
                        >,
                        Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
	              >,
	            OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
              >& rhs,

	      const Subset& s) {

  typedef OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > > F;

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;



  const F& scal = static_cast<const F&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

#if DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M += alpha*M ) subset = s " << endl;
#endif

  BAGELQDPFloat *scalar   = const_cast<BAGELQDPFloat*>(&(scal.elem().elem().elem().elem()));  
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;
     
     BAGELQDPFloat *y = &(d.elem(start).elem().elem(0,0).real());

     BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_vaxpy3(y, scalar, x, y, 3*num_sites);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       BAGELQDPFloat *y = &(d.elem(i).elem().elem(0,0).real());
       BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_vaxpy3(y, scalar, x,y, 3);

     }
   }

}

  // -=a * M
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 

	      const QDPExpr<
	              BinaryNode< OpMultiply, 
	                Reference<
                          QDPType<
                            PScalar<PScalar<RScalar<BAGELQDPFloat> > >,
	                    OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > >
                          > 
                        >,
                        Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
	              >,
	            OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
              >& rhs,

	      const Subset& s) {

  typedef OScalar<PScalar<PScalar<RScalar<BAGELQDPFloat> > > > F;

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;



  const F& scal = static_cast<const F&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

#if DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M -= alpha*M ) subset = s " << endl;
#endif

  BAGELQDPFloat scalar = -scal.elem().elem().elem().elem();
  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;
     
     BAGELQDPFloat *y = &(d.elem(start).elem().elem(0,0).real());

     BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_vaxpy3(y, &scalar, x, y, 3*num_sites);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       BAGELQDPFloat *y = &(d.elem(i).elem().elem(0,0).real());
       BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_vaxpy3(y, &scalar, x,y, 3);

     }
   }

}

  // += M
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAddAssign& op, 

	      const QDPExpr< 
	              UnaryNode<OpIdentity,
	                Reference<
                          QDPType<
	      PScalar<PColorMatrix<RComplexFloat, 3> >,
                            OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >
                          > 
                        >
	              >,
	              OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                    >& rhs,

	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  const C& r = static_cast<const C&>(rhs.expression().child());

#if DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M += M ) subset = s " << endl;
#endif

  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;
     
     BAGELQDPFloat *y = &(d.elem(start).elem().elem(0,0).real());
     BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_vadd3(y, x, y, 3*num_sites);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       BAGELQDPFloat *y = &(d.elem(i).elem().elem(0,0).real());
       BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_vadd3(y,x,y, 3);

     }
   }

}

  // -=a * M
template<>
inline
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpSubtractAssign& op, 

	      const QDPExpr< 
	              UnaryNode< OpIdentity, 
	                Reference<
                          QDPType<
                            PScalar<PColorMatrix<RComplexFloat, 3> >, 
	                    OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                          > 
                        >
                      >,
 	              OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > 
                   >& rhs,

	      const Subset& s) {

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;

  const C& r = static_cast<const C&>(rhs.expression().child());

#if DEBUG_BAGELQDP_LINALG
  QDPIO::cout << "evaluate(M -= M ) subset = s " << endl;
#endif

  if( s.hasOrderedRep() ) { 
     // Do whole subset
     unsigned int start = s.start();
     unsigned int end   = s.end();

     unsigned long num_sites = end - start + 1;
     
     BAGELQDPFloat *y = &(d.elem(start).elem().elem(0,0).real());

     BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(start).elem().elem(0,0).real()));

     qdp_vsub3(y, y,x, 3*num_sites);

   }
   else { 
     // Do site by site
     const int* tab = s.siteTable().slice();
     unsigned int num_sites = s.numSiteTable();
     for(unsigned int j=0; j < num_sites; j++) {
       int i = tab[j];
       BAGELQDPFloat *y = &(d.elem(i).elem().elem(0,0).real());
       BAGELQDPFloat *x   = const_cast<BAGELQDPFloat*>(&(r.elem(i).elem().elem(0,0).real()));

       qdp_vsub3(y,y,x, 3);

     }
   }

}


#endif


#if defined(DEBUG_BAGELQDP_LINALG)
#undef DEBUG_BAGELQDP_LINALG
#endif

} // namespace QDP;

#endif
