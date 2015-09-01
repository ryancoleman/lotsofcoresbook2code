// $Id: qdp_scalarsite_sse.cc,v 1.3 2009-09-15 20:49:12 bjoo Exp $

/*! @file
 * @brief Intel SSE optimizations
 * 
 * SSE optimizations of basic operations
 */


#include "qdp.h"


// These SSE asm instructions are only supported under GCC/G++
#if defined(__GNUC__)
#include "qdp_sse_intrin.h"
namespace QDP {



#if 1
//-------------------------------------------------------------------
// Specialization to optimize the case   
//    LatticeColorMatrix[ Subset] = LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_M_times_M" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 

      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *dm = (su3_matrixf *)&(d.elem(i).elem().elem(0,0).real());

      intrin_sse_mult_su3_nn(lm, rm, dm);
			     
    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *dm = (su3_matrixf *)&(d.elem(i).elem().elem(0,0).real());

      intrin_sse_mult_su3_nn(lm, rm, dm);

    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_aM_times_M" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *dm = (su3_matrixf *)&(d.elem(i).elem().elem(0,0).real());

      intrin_sse_mult_su3_an(lm, rm, dm);

    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {

      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *dm = (su3_matrixf *)&(d.elem(i).elem().elem(0,0).real());

      intrin_sse_mult_su3_an(lm, rm, dm);

    }
  }

}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_M_times_aM" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *dm = (su3_matrixf *)&(d.elem(i).elem().elem(0,0).real());

      intrin_sse_mult_su3_na(lm, rm, dm);

    }
  }
  else { 

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *dm = (su3_matrixf *)&(d.elem(i).elem().elem(0,0).real());

      intrin_sse_mult_su3_na(lm, rm, dm);

    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_Ma_times_Ma" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  PColorMatrix<RComplexFloat,3> tmp;

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(rm, lm, tmpm);
      
      // Take the adj(r*l) = adj(l)*adj(r)
      d.elem(i).elem().elem(0,0).real() =  tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() = -tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() =  tmp.elem(1,0).real();
      d.elem(i).elem().elem(0,1).imag() = -tmp.elem(1,0).imag();
      d.elem(i).elem().elem(0,2).real() =  tmp.elem(2,0).real();
      d.elem(i).elem().elem(0,2).imag() = -tmp.elem(2,0).imag();
      
      d.elem(i).elem().elem(1,0).real() =  tmp.elem(0,1).real();
      d.elem(i).elem().elem(1,0).imag() = -tmp.elem(0,1).imag();
      d.elem(i).elem().elem(1,1).real() =  tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() = -tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() =  tmp.elem(2,1).real();
      d.elem(i).elem().elem(1,2).imag() = -tmp.elem(2,1).imag();
      
      d.elem(i).elem().elem(2,0).real() =  tmp.elem(0,2).real();
      d.elem(i).elem().elem(2,0).imag() = -tmp.elem(0,2).imag();
      d.elem(i).elem().elem(2,1).real() =  tmp.elem(1,2).real();
      d.elem(i).elem().elem(2,1).imag() = -tmp.elem(1,2).imag();
      d.elem(i).elem().elem(2,2).real() =  tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() = -tmp.elem(2,2).imag();
    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(rm, lm, tmpm);
      
      // Take the adj(r*l) = adj(l)*adj(r)
      d.elem(i).elem().elem(0,0).real() =  tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() = -tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() =  tmp.elem(1,0).real();
      d.elem(i).elem().elem(0,1).imag() = -tmp.elem(1,0).imag();
      d.elem(i).elem().elem(0,2).real() =  tmp.elem(2,0).real();
      d.elem(i).elem().elem(0,2).imag() = -tmp.elem(2,0).imag();
      
      d.elem(i).elem().elem(1,0).real() =  tmp.elem(0,1).real();
      d.elem(i).elem().elem(1,0).imag() = -tmp.elem(0,1).imag();
      d.elem(i).elem().elem(1,1).real() =  tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() = -tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() =  tmp.elem(2,1).real();
      d.elem(i).elem().elem(1,2).imag() = -tmp.elem(2,1).imag();
      
      d.elem(i).elem().elem(2,0).real() =  tmp.elem(0,2).real();
      d.elem(i).elem().elem(2,0).imag() = -tmp.elem(0,2).imag();
      d.elem(i).elem().elem(2,1).real() =  tmp.elem(1,2).real();
      d.elem(i).elem().elem(2,1).imag() = -tmp.elem(1,2).imag();
      d.elem(i).elem().elem(2,2).real() =  tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() = -tmp.elem(2,2).imag();
    }
  }
}

//-------------------------------------------------------------------

// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_M_times_M" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  PColorMatrix<RComplexFloat,3> tmp;

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(lm, rm, tmpm);

      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_aM_times_M" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

  PColorMatrix<RComplexFloat,3> tmp;

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_an(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_an(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_M_times_aM" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  PColorMatrix<RComplexFloat,3> tmp;

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_na(lm, rm, tmpm);

      
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
  else { 

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_na(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] += adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_peq_Ma_times_Ma" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  PColorMatrix<RComplexFloat,3> tmp;

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(rm, lm, tmpm);
      
      // Take the adj(r*l) = adj(l)*adj(r)
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(2,0).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(2,1).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(0,2).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(1,2).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();

    }
  }
  else { 

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); ++j) {
      int i = tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(rm, lm, tmpm);
      
      // Take the adj(r*l) = adj(l)*adj(r)
      d.elem(i).elem().elem(0,0).real() += tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() += tmp.elem(1,0).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(0,2).real() += tmp.elem(2,0).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(2,0).imag();
      
      d.elem(i).elem().elem(1,0).real() += tmp.elem(0,1).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(1,1).real() += tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() += tmp.elem(2,1).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(2,1).imag();
      
      d.elem(i).elem().elem(2,0).real() += tmp.elem(0,2).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(0,2).imag();
      d.elem(i).elem().elem(2,1).real() += tmp.elem(1,2).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(1,2).imag();
      d.elem(i).elem().elem(2,2).real() += tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
}

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix = LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< TCol, OLattice< TCol > > > >,
                 OLattice< TCol > >& rhs, 
	      const Subset& s) 
{
  typedef OLattice<TCol> C;
  const C& l = static_cast<const C&>(rhs.expression().child());


  if( s.hasOrderedRep() ) {

    const int start = s.start();
    const int end = s.end();

    REAL32* d_ptr =&(d.elem(start).elem().elem(0,0).real());
    const REAL32* r_ptr =&(l.elem(start).elem().elem(0,0).real());
  
    const unsigned int total_reals = (end-start+1)*3*3*2;
    const unsigned int total_v4sf = total_reals/4;
    const unsigned int remainder = total_reals%4;
  
    float* d_ptr_v4sf = (float *)d_ptr;
    float* r_ptr_v4sf = (float *)r_ptr;

  
    for(unsigned int i = 0 ; i < total_v4sf; i++, d_ptr_v4sf +=4, r_ptr_v4sf+=4 ) { 
      _mm_store_ps( d_ptr_v4sf, _mm_load_ps(r_ptr_v4sf));
    }
  
    
    r_ptr = (REAL32 *)r_ptr_v4sf;
    d_ptr = (REAL32 *)d_ptr_v4sf;
    for(unsigned int i=0; i < remainder; i++, r_ptr++, d_ptr++) { 
      *d_ptr = *r_ptr;
    }
  }
  else {
    // Unordered case 
    const int* tab = s.siteTable().slice();
    
    // Loop through the sites
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      // Do the copy in the dumb way -- this could become quite complex
      // Depending on whether the individual matrices are aligned or not.
      d.elem(i).elem() = l.elem(i).elem();

    }

  }
}

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix += LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpAddAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< TCol, OLattice< TCol > > > >,
                 OLattice< TCol > >& rhs, 
	      const Subset& s) 
{
  typedef OLattice<TCol> C;
  const C& l = static_cast<const C&>(rhs.expression().child());

  if( s.hasOrderedRep() ) { 
    const int start = s.start();
    const int end = s.end();

    REAL32* d_ptr =&(d.elem(start).elem().elem(0,0).real());
    const REAL32* r_ptr =&(l.elem(start).elem().elem(0,0).real());
    
    const unsigned int total_reals = (end-start+1)*3*3*2;
    const unsigned int total_v4sf = total_reals/4;
    const unsigned int remainder = total_reals%4;
    
    float* d_ptr_v4sf = (float *)d_ptr;
    float* r_ptr_v4sf = (float *)r_ptr;
    
    
    for(unsigned int i = 0 ; i < total_v4sf; i++, d_ptr_v4sf +=4, r_ptr_v4sf+=4 ) { 
      _mm_store_ps( d_ptr_v4sf, _mm_add_ps( _mm_load_ps(d_ptr_v4sf), 
					    _mm_load_ps(r_ptr_v4sf) ) );
    }
    
    
    r_ptr = (REAL32 *)r_ptr_v4sf;
    d_ptr = (REAL32 *)d_ptr_v4sf;
    for(unsigned int i=0; i < remainder; i++, r_ptr++, d_ptr++) { 
      *d_ptr += *r_ptr;
    }
  }
  else {
    // Unordered case 
    const int* tab = s.siteTable().slice();
    
    // Loop through the sites
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      // Do the copy in the dumb way -- this could become quite complex
      // Depending on whether the individual matrices are aligned or not.
      d.elem(i).elem() += l.elem(i).elem();

    }

  }
}

//-------------------------------------------------------------------
// Specialization to optimize the case
//   LatticeColorMatrix -= LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<
	         UnaryNode<OpIdentity, Reference< QDPType< TCol, OLattice< TCol > > > >,
                 OLattice< TCol > >& rhs, 
	      const Subset& s) 
{
  typedef OLattice<TCol> C;
  const C& l = static_cast<const C&>(rhs.expression().child());
  if (s.hasOrderedRep()) { 
    const int start = s.start();
    const int end = s.end();
    
    REAL32* d_ptr =&(d.elem(start).elem().elem(0,0).real());
    const REAL32* r_ptr =&(l.elem(start).elem().elem(0,0).real());
    
    const unsigned int total_reals = (end-start+1)*3*3*2;
    const unsigned int total_v4sf = total_reals/4;
    const unsigned int remainder = total_reals%4;
    
    float* d_ptr_v4sf = (float *)d_ptr;
    float* r_ptr_v4sf = (float *)r_ptr;
    
    
    for(unsigned int i = 0 ; i < total_v4sf; i++, d_ptr_v4sf +=4, r_ptr_v4sf+=4 ) { 
      _mm_store_ps( d_ptr_v4sf, _mm_sub_ps( _mm_load_ps( d_ptr_v4sf), 
					    _mm_load_ps( r_ptr_v4sf) ) );
    }
    
    
    r_ptr = (REAL32 *)r_ptr_v4sf;
    d_ptr = (REAL32 *)d_ptr_v4sf;
    for(unsigned int i=0; i < remainder; i++, r_ptr++, d_ptr++) { 
      *d_ptr -= *r_ptr;
    }
  }
  else {   
    // Unordered case 
    const int* tab = s.siteTable().slice();
    
    // Loop through the sites
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i = tab[j];

      // Do the copy in the dumb way -- this could become quite complex
      // Depending on whether the individual matrices are aligned or not.
      d.elem(i).elem() -= l.elem(i).elem();

    }

  }
}

// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_M_times_M" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  PColorMatrix<RComplexFloat,3> tmp;


  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
  else { 

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= adj(LatticeColorMatrix) * LatticeColorMatrix
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiply, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >, 
	                    Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_aM_times_M" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right());

  PColorMatrix<RComplexFloat,3> tmp;
  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_an(lm, rm, tmpm);

      
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
  else { 

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_an(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= LatticeColorMatrix * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiplyAdj, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_M_times_aM" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  PColorMatrix<RComplexFloat,3> tmp;

  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_na(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_na(lm, rm, tmpm);
      
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() -= tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(0,1).imag() -= tmp.elem(0,1).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(0,2).imag() -= tmp.elem(0,2).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(1,0).imag() -= tmp.elem(1,0).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() -= tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(1,2).imag() -= tmp.elem(1,2).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(2,0).imag() -= tmp.elem(2,0).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(2,1).imag() -= tmp.elem(2,1).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() -= tmp.elem(2,2).imag();
    }
  }
}


// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] -= adj(LatticeColorMatrix) * adj(LatticeColorMatrix)
template<>
void evaluate(OLattice< TCol >& d, 
	      const OpSubtractAssign& op, 
	      const QDPExpr<BinaryNode<OpAdjMultiplyAdj, 
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > >,
	                    UnaryNode<OpIdentity, Reference<QDPType< TCol, OLattice< TCol > > > > >,
	                    OLattice< TCol > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_meq_Ma_times_Ma" << endl;

  typedef OLattice< TCol >    C;

  const C& l = static_cast<const C&>(rhs.expression().left().child());
  const C& r = static_cast<const C&>(rhs.expression().right().child());

  PColorMatrix<RComplexFloat,3> tmp;
  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(rm, lm, tmpm);

      // Take the adj(r*l) = adj(l)*adj(r)
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(2,0).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(2,1).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(0,2).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(1,2).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();

    }
  }
  else { 
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      su3_matrixf *lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      su3_matrixf *rm = (su3_matrixf *)&(r.elem(i).elem().elem(0,0).real());
      su3_matrixf *tmpm = (su3_matrixf *)&(tmp.elem(0,0).real());

      intrin_sse_mult_su3_nn(rm, lm, tmpm);
      
      
      // Take the adj(r*l) = adj(l)*adj(r)
      d.elem(i).elem().elem(0,0).real() -= tmp.elem(0,0).real();
      d.elem(i).elem().elem(0,0).imag() += tmp.elem(0,0).imag();
      d.elem(i).elem().elem(0,1).real() -= tmp.elem(1,0).real();
      d.elem(i).elem().elem(0,1).imag() += tmp.elem(1,0).imag();
      d.elem(i).elem().elem(0,2).real() -= tmp.elem(2,0).real();
      d.elem(i).elem().elem(0,2).imag() += tmp.elem(2,0).imag();
      
      d.elem(i).elem().elem(1,0).real() -= tmp.elem(0,1).real();
      d.elem(i).elem().elem(1,0).imag() += tmp.elem(0,1).imag();
      d.elem(i).elem().elem(1,1).real() -= tmp.elem(1,1).real();
      d.elem(i).elem().elem(1,1).imag() += tmp.elem(1,1).imag();
      d.elem(i).elem().elem(1,2).real() -= tmp.elem(2,1).real();
      d.elem(i).elem().elem(1,2).imag() += tmp.elem(2,1).imag();
      
      d.elem(i).elem().elem(2,0).real() -= tmp.elem(0,2).real();
      d.elem(i).elem().elem(2,0).imag() += tmp.elem(0,2).imag();
      d.elem(i).elem().elem(2,1).real() -= tmp.elem(1,2).real();
      d.elem(i).elem().elem(2,1).imag() += tmp.elem(1,2).imag();
      d.elem(i).elem().elem(2,2).real() -= tmp.elem(2,2).real();
      d.elem(i).elem().elem(2,2).imag() += tmp.elem(2,2).imag();
    }
  }
}


//-------------------------------------------------------------------

// Specialization to optimize the case   
//    LatticeHalfFermion = LatticeColorMatrix * LatticeHalfFermion
// NOTE: let this be a subroutine to save space
template<>
void evaluate(OLattice< TVec2 >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	                    Reference<QDPType< TCol, OLattice< TCol > > >, 
	                    Reference<QDPType< TVec2, OLattice< TVec2 > > > >,
	                    OLattice< TVec2 > >& rhs,
	      const Subset& s)
{
#if defined(QDP_SCALARSITE_DEBUG)
  cout << "specialized QDP_H_M_times_H" << endl;
#endif

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >       C;
  typedef OLattice<PSpinVector<PColorVector<RComplexFloat, 3>, 2> > H;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const H& r = static_cast<const H&>(rhs.expression().right());



  if( s.hasOrderedRep() ) { 
    for(int i=s.start(); i <= s.end(); i++) { 

      su3_matrixf* lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      half_wilson_vectorf* rh = (half_wilson_vectorf*)&(r.elem(i).elem(0).elem(0).real());
      half_wilson_vectorf* dh = (half_wilson_vectorf*)&(d.elem(i).elem(0).elem(0).real());

      intrin_sse_mult_su3_mat_hwvec(lm,rh,dh);

    }
  }
  else { 

    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      int i=tab[j];
      su3_matrixf* lm = (su3_matrixf *)&(l.elem(i).elem().elem(0,0).real());
      half_wilson_vectorf* rh = (half_wilson_vectorf*)&(r.elem(i).elem(0).elem(0).real());
      half_wilson_vectorf* dh = (half_wilson_vectorf*)&(d.elem(i).elem(0).elem(0).real());

      intrin_sse_mult_su3_mat_hwvec(lm,rh,dh);
    }
  }
}
#endif


//-------------------------------------------------------------------
// GNUC vector type

//#define DEBUG_BLAS

// AXPY and AXMY routines
void vaxpy3(REAL32 *Out,REAL32 *scalep,REAL32 *InScale, REAL32 *Add,int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vaxpy3" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

  v4sf vscalep = _mm_load_ss(scalep);
  vscalep = _mm_shuffle_ps(vscalep, vscalep, 0);


  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_add_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+ 0)), _mm_load_ps(Add+ 0)));
    _mm_store_ps(Out+ 4, _mm_add_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+ 4)), _mm_load_ps(Add+ 4)));
    _mm_store_ps(Out+ 8, _mm_add_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+ 8)), _mm_load_ps(Add+ 8)));
    _mm_store_ps(Out+12, _mm_add_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+12)), _mm_load_ps(Add+12)));
    _mm_store_ps(Out+16, _mm_add_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+16)), _mm_load_ps(Add+16)));
    _mm_store_ps(Out+20, _mm_add_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+20)), _mm_load_ps(Add+20)));

    Out += 24; InScale += 24; Add += 24;
  }
}


void vaxmy3(REAL32 *Out,REAL32 *scalep,REAL32 *InScale, REAL32 *Sub,int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vaxmy3" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

//  v4sf va = load_v4sf((float *)&a);
  v4sf vscalep = _mm_load_ss(scalep);
  vscalep = _mm_shuffle_ps( vscalep, vscalep, 0);

  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_sub_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+ 0)), _mm_load_ps(Sub+ 0)));
    _mm_store_ps(Out+ 4, _mm_sub_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+ 4)), _mm_load_ps(Sub+ 4)));
    _mm_store_ps(Out+ 8, _mm_sub_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+ 8)), _mm_load_ps(Sub+ 8)));
    _mm_store_ps(Out+12, _mm_sub_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+12)), _mm_load_ps(Sub+12)));
    _mm_store_ps(Out+16, _mm_sub_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+16)), _mm_load_ps(Sub+16)));
    _mm_store_ps(Out+20, _mm_sub_ps(_mm_mul_ps(vscalep, _mm_load_ps(InScale+20)), _mm_load_ps(Sub+20)));

    Out += 24; InScale += 24; Sub += 24;
  }
}


void vadd(REAL32 *Out, REAL32 *In1, REAL32 *In2, int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vadd" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_add_ps(_mm_load_ps(In1+ 0), _mm_load_ps(In2+ 0)));
    _mm_store_ps(Out+ 4, _mm_add_ps(_mm_load_ps(In1+ 4), _mm_load_ps(In2+ 4)));
    _mm_store_ps(Out+ 8, _mm_add_ps(_mm_load_ps(In1+ 8), _mm_load_ps(In2+ 8)));
    _mm_store_ps(Out+12, _mm_add_ps(_mm_load_ps(In1+12), _mm_load_ps(In2+12)));
    _mm_store_ps(Out+16, _mm_add_ps(_mm_load_ps(In1+16), _mm_load_ps(In2+16)));
    _mm_store_ps(Out+20, _mm_add_ps(_mm_load_ps(In1+20), _mm_load_ps(In2+20)));

    Out += 24; In1 += 24; In2 += 24;
  }
}


void vsub(REAL32 *Out, REAL32 *In1, REAL32 *In2, int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vsub" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_sub_ps(_mm_load_ps(In1+ 0), _mm_load_ps(In2+ 0)));
    _mm_store_ps(Out+ 4, _mm_sub_ps(_mm_load_ps(In1+ 4), _mm_load_ps(In2+ 4)));
    _mm_store_ps(Out+ 8, _mm_sub_ps(_mm_load_ps(In1+ 8), _mm_load_ps(In2+ 8)));
    _mm_store_ps(Out+12, _mm_sub_ps(_mm_load_ps(In1+12), _mm_load_ps(In2+12)));
    _mm_store_ps(Out+16, _mm_sub_ps(_mm_load_ps(In1+16), _mm_load_ps(In2+16)));
    _mm_store_ps(Out+20, _mm_sub_ps(_mm_load_ps(In1+20), _mm_load_ps(In2+20)));

    Out += 24; In1 += 24; In2 += 24;
  }
}

void vscal(REAL32 *Out, REAL32 *scalep, REAL32 *In, int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vadd" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

//  v4sf va = load_v4sf((float *)&a);
  v4sf vscalep = _mm_load_ss(scalep);
  vscalep = _mm_shuffle_ps(vscalep, vscalep, 0);

  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_mul_ps(vscalep, _mm_load_ps(In+ 0)));
    _mm_store_ps(Out+ 4, _mm_mul_ps(vscalep, _mm_load_ps(In+ 4)));
    _mm_store_ps(Out+ 8, _mm_mul_ps(vscalep, _mm_load_ps(In+ 8)));
    _mm_store_ps(Out+12, _mm_mul_ps(vscalep, _mm_load_ps(In+12)));
    _mm_store_ps(Out+16, _mm_mul_ps(vscalep, _mm_load_ps(In+16)));
    _mm_store_ps(Out+20, _mm_mul_ps(vscalep, _mm_load_ps(In+20)));

    Out += 24; In += 24;
  }
}  


void vaxpby3(REAL32 *Out, REAL32 *a, REAL32 *x, REAL32 *b, REAL32 *y, int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vaxpby3: a*x+b*y" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

//  v4sf va = load_v4sf((float *)&a);
  v4sf va = _mm_load_ss(a);
  v4sf vb = _mm_load_ss(b);
  va = _mm_shuffle_ps(va, va, 0);
  vb = _mm_shuffle_ps(vb, vb, 0);


  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_add_ps(_mm_mul_ps(va, _mm_load_ps(x+ 0)), _mm_mul_ps(vb, _mm_load_ps(y+ 0))));
    _mm_store_ps(Out+ 4, _mm_add_ps(_mm_mul_ps(va, _mm_load_ps(x+ 4)), _mm_mul_ps(vb, _mm_load_ps(y+ 4))));
    _mm_store_ps(Out+ 8, _mm_add_ps(_mm_mul_ps(va, _mm_load_ps(x+ 8)), _mm_mul_ps(vb, _mm_load_ps(y+ 8))));
    _mm_store_ps(Out+12, _mm_add_ps(_mm_mul_ps(va, _mm_load_ps(x+12)), _mm_mul_ps(vb, _mm_load_ps(y+12))));
    _mm_store_ps(Out+16, _mm_add_ps(_mm_mul_ps(va, _mm_load_ps(x+16)), _mm_mul_ps(vb, _mm_load_ps(y+16))));
    _mm_store_ps(Out+20, _mm_add_ps(_mm_mul_ps(va, _mm_load_ps(x+20)), _mm_mul_ps(vb, _mm_load_ps(y+20))));

    Out += 24; x += 24; y += 24;
  }
}


void vaxmby3(REAL32 *Out, REAL32 *a, REAL32 *x, REAL32 *b, REAL32 *y, int n_4vec)
{
#ifdef DEBUG_BLAS
  QDPIO::cout << "SSE_TEST: vaxmby3: a*x-b*y" << endl;
#endif

//  int n_loops = n_4vec >> 2;   // only works on multiple of length 4 vectors
  int n_loops = n_4vec;   // only works on multiple of length 24 vectors

//  v4sf va = load_v4sf((float *)&a);
  v4sf va = _mm_load_ss(a);
  v4sf vb = _mm_load_ss(b);
  va = _mm_shuffle_ps( va, va, 0);
  vb = _mm_shuffle_ps( vb, vb, 0);

  for (; n_loops-- > 0; )
  {
    _mm_store_ps(Out+ 0, _mm_sub_ps(_mm_mul_ps(va, _mm_load_ps(x+ 0)), _mm_mul_ps(vb, _mm_load_ps(y+ 0))));
    _mm_store_ps(Out+ 4, _mm_sub_ps(_mm_mul_ps(va, _mm_load_ps(x+ 4)), _mm_mul_ps(vb, _mm_load_ps(y+ 4))));
    _mm_store_ps(Out+ 8, _mm_sub_ps(_mm_mul_ps(va, _mm_load_ps(x+ 8)), _mm_mul_ps(vb, _mm_load_ps(y+ 8))));
    _mm_store_ps(Out+12, _mm_sub_ps(_mm_mul_ps(va, _mm_load_ps(x+12)), _mm_mul_ps(vb, _mm_load_ps(y+12))));
    _mm_store_ps(Out+16, _mm_sub_ps(_mm_mul_ps(va, _mm_load_ps(x+16)), _mm_mul_ps(vb, _mm_load_ps(y+16))));
    _mm_store_ps(Out+20, _mm_sub_ps(_mm_mul_ps(va, _mm_load_ps(x+20)), _mm_mul_ps(vb, _mm_load_ps(y+20))));

    Out += 24; x += 24; y += 24;
  }
}



/*! Sum squared routine. USE THIS ONLY FOR arrays that
 * are a multiple of 24 or 48 floats */

void local_sumsq_24_48(REAL64 *Out, REAL32 *In, int n_4vec) 
{

  __m128d vsum = _mm_setzero_pd();
  __m128d vsum2 = _mm_setzero_pd();
  __m128d lower_2, upper_2;
  __m128d dat_sq, dat_sq2;
  __m128  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  
  REAL32* num=In;
  int loop_end;



  loop_end = n_4vec-1;
  tmp1 = _mm_load_ps(num); num+=4;

  for(int site=0; site < loop_end;  site++)  {    

    tmp3 = _mm_load_ps(num); num+=4;                   // Load 4

    // First 4 numbers
    tmp2 = _mm_shuffle_ps(tmp1, tmp1, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp1);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
    tmp4 = _mm_load_ps(num); num+=4;                   // Load 4
      
    // 2nd 4 numbers
    tmp2 = _mm_shuffle_ps(tmp3, tmp3, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp3);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
    tmp5 = _mm_load_ps(num); num+=4;                   // Load 4


    // 3rd 4 numbers      
    tmp2 = _mm_shuffle_ps(tmp4, tmp4, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp4);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
    tmp6 = _mm_load_ps(num); num+=4;                   // Load 4
  
    // 4th 4 numbers
    tmp2 = _mm_shuffle_ps(tmp5, tmp5, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp5);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
    tmp7 = _mm_load_ps(num); num+=4;                   // Load 4
  
    // 5th 4 numbers
    tmp2 = _mm_shuffle_ps(tmp6, tmp6, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp6);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
    tmp1 = _mm_load_ps(num); num+=4;                   // Load 4
 
    // 6th 4 numbers
    tmp2 = _mm_shuffle_ps(tmp7, tmp7, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp7);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
      
        
  }
  /* Last one */
  {     
    tmp3 = _mm_load_ps(num); num+=4;                   // Load 4

    // First 4 numbers
    tmp2 = _mm_shuffle_ps(tmp1, tmp1, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp1);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
    
    tmp4 = _mm_load_ps(num); num+=4;                   // Load 4
      
    // 2nd 4 numbers
    tmp2 = _mm_shuffle_ps(tmp3, tmp3, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp3);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
      
    tmp5 = _mm_load_ps(num); num+=4;                   // Load 4
      
    // 3rd 4 numbers
    tmp2 = _mm_shuffle_ps(tmp4, tmp4, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp4);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
    
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);

    tmp6 = _mm_load_ps(num); num+=4;                   // Load 4
    
    // 4th 4 numbers
    tmp2 = _mm_shuffle_ps(tmp5, tmp5, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp5);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
    
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);

    tmp7 = _mm_load_ps(num); num+=4;                   // Load 4
    
    // 5th 4 numbers
    tmp2 = _mm_shuffle_ps(tmp6, tmp6, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp6);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
    
    // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
      
    // 6th 4 numbers
    tmp2 = _mm_shuffle_ps(tmp7, tmp7, 0x0e);    // flip numbers into tmp2
    lower_2 = _mm_cvtps_pd(tmp7);               // convert to double
    upper_2 = _mm_cvtps_pd(tmp2);
      
      // f*f
    dat_sq = _mm_mul_pd(lower_2, lower_2);        
    vsum = _mm_add_pd(vsum, dat_sq);
    dat_sq2 = _mm_mul_pd(upper_2, upper_2);
    vsum2 = _mm_add_pd(vsum2, dat_sq2);
  }


  vsum = _mm_add_pd(vsum, vsum2);
  /* Now sum horizontally on a and b */
  /* Move high word of vsum to low of lower_2 */
  /* Cross operation: SRC = vsum = | sum_1 | sum_0 |
                      DEST= vsum = | sum 1 | sum_0 |

     Bit field 0 != 0 => DEST[63-0] = DEST[127-64] = sum_1
     Bit field 1 == 0 => DEST[127:64] ← SRC[63:0] = sum_0 ;
     So output register: sum_0 | sum_1
  */
  lower_2 = _mm_shuffle_pd(vsum, vsum, 0x01);

  /* Now sum up:    lower_2 = sum_0 | sum_1
                +   vsum    = sum_1 | sum_0
                  ========================
                              Sum   | Sum 
  */

  upper_2 = _mm_add_pd(lower_2, vsum);
  
  /* Sum a should now have 2 copies of the sum */
  /* I should store the result */
  _mm_storeh_pd((double *)Out, upper_2);

  //  *Out=result;
}


void local_vcdot(REAL64 *Out_re, REAL64 *Out_im, REAL32 *V1, REAL32 *V2, int n_3vec)
{
  double result_re;
  double result_im;

  
  double v1_0r;
  double v1_0i;
  double v1_1r;
  double v1_1i;
  double v1_2r;
  double v1_2i;

  double v2_0r;
  double v2_0i;
  double v2_1r;
  double v2_1i;
  double v2_2r;
  double v2_2i;

  int counter=0;
  unsigned long vecptr1=0;
  unsigned long vecptr2=0;
  result_re= 0;
  result_im= 0;

  
  if( n_3vec > 0 ) { 

    v1_0r = (REAL64)V1[vecptr1++];
    v2_0r = (REAL64)V2[vecptr2++];
    
    v1_0i = (REAL64)V1[vecptr1++];
    v2_0i = (REAL64)V2[vecptr2++];
    
    v1_1r = (REAL64)V1[vecptr1++];
    v2_1r = (REAL64)V2[vecptr2++];

    for(counter=0; counter < n_3vec-1; counter++) {


      
      result_re = result_re + v1_0r*v2_0r;
      v1_1i =(REAL64)V1[vecptr1++];
      result_im = result_im - v1_0i*v2_0r;
      v2_1i = (REAL64)V2[vecptr2++];
      result_im = result_im + v1_0r*v2_0i;
      v1_2r = (REAL64)V1[vecptr1++];
      result_re = result_re + v1_0i*v2_0i;
      v2_2r = (REAL64)V2[vecptr2++];
      
      result_re = result_re + v1_1r*v2_1r;
      v1_2i = (REAL64)V1[vecptr1++];
      result_im = result_im - v1_1i*v2_1r;
      v2_2i = (REAL64)V2[vecptr2++];
      result_im = result_im + v1_1r*v2_1i;
      v1_0r = (REAL64)V1[vecptr1++];
      result_re = result_re + v1_1i*v2_1i;
      v2_0r = (REAL64)V2[vecptr2++];

      result_re = result_re + v1_2r*v2_2r;
      v1_0i = (REAL64)V1[vecptr1++];
      result_im = result_im - v1_2i*v2_2r;
      v2_0i = (REAL64)V2[vecptr2++];
      result_im = result_im + v1_2r*v2_2i;
      v1_1r = (REAL64)V1[vecptr1++];
      result_re = result_re + v1_2i*v2_2i;
      v2_1r = (REAL64)V2[vecptr2++];

    }

    // Last one plus drain...
    result_re = result_re + v1_0r*v2_0r;
    v1_1i =(REAL64)V1[vecptr1++];
    result_im = result_im - v1_0i*v2_0r;
    v2_1i = (REAL64)V2[vecptr2++];
    result_im = result_im + v1_0r*v2_0i;
    v1_2r = (REAL64)V1[vecptr1++];
    result_re = result_re + v1_0i*v2_0i;
    v2_2r = (REAL64)V2[vecptr2++];
      
    result_re = result_re + v1_1r*v2_1r;
    v1_2i = (REAL64)V1[vecptr1++];
    result_im = result_im - v1_1i*v2_1r;
    v2_2i = (REAL64)V2[vecptr2++];
 
    result_im = result_im + v1_1r*v2_1i;
    result_re = result_re + v1_1i*v2_1i;
    

    result_re = result_re + v1_2r*v2_2r;
    result_im = result_im - v1_2i*v2_2r;
    result_im = result_im + v1_2r*v2_2i;
    result_re = result_re + v1_2i*v2_2i;    

  }
  
  *Out_re=(REAL64)result_re;
  *Out_im=(REAL64)result_im;
}


void local_vcdot_real(REAL64 *Out, REAL32 *V1, REAL32 *V2, int n_3vec)
{
  REAL64 result;
  
  REAL64 v1_0r;
  REAL64 v1_0i;
  REAL64 v1_1r;
  REAL64 v1_1i;
  REAL64 v1_2r;
  REAL64 v1_2i;

  REAL64 v2_0r;
  REAL64 v2_0i;
  REAL64 v2_1r;
  REAL64 v2_1i;
  REAL64 v2_2r;
  REAL64 v2_2i;

  int counter=0;
  unsigned long vecptr1=0;
  unsigned long vecptr2=0;
  result= 0;

  
  if( n_3vec > 0 ) { 

    // Prefetch 
    v1_0r = (REAL64)V1[vecptr1++];
    v2_0r = (REAL64)V2[vecptr2++];

    v1_0i = (REAL64)V1[vecptr1++];
    v2_0i = (REAL64)V2[vecptr2++];

    v1_1r = (REAL64)V1[vecptr1++];
    v2_1r = (REAL64)V2[vecptr2++];

    v1_1i =(REAL64)V1[vecptr1++];
    v2_1i = (REAL64)V2[vecptr2++];
    
    v1_2r = (REAL64)V1[vecptr1++];
    v2_2r = (REAL64)V2[vecptr2++];
    
    v1_2i = (REAL64)V1[vecptr1++];
    v2_2i = (REAL64)V2[vecptr2++];

    for(counter=0; counter < n_3vec-1; counter++) {
      result = result + v1_0r*v2_0r;
      v1_0r = (REAL64)V1[vecptr1++];
      v2_0r = (REAL64)V2[vecptr2++];    

      result = result + v1_0i*v2_0i;
      v1_0i = (REAL64)V1[vecptr1++];
      v2_0i = (REAL64)V2[vecptr2++];
      
      result = result + v1_1r*v2_1r;
      v1_1r = (REAL64)V1[vecptr1++];
      v2_1r = (REAL64)V2[vecptr2++];

      result = result + v1_1i*v2_1i;
      v1_1i =(REAL64)V1[vecptr1++];
      v2_1i = (REAL64)V2[vecptr2++];
      
      result = result + v1_2r*v2_2r;
      v1_2r = (REAL64)V1[vecptr1++];
      v2_2r = (REAL64)V2[vecptr2++];

      result = result + v1_2i*v2_2i;
      v1_2i = (REAL64)V1[vecptr1++];
      v2_2i = (REAL64)V2[vecptr2++];


    }

    // Last one plus drain...
    result = result + v1_0r*v2_0r;
    result = result + v1_0i*v2_0i;
    result = result + v1_1r*v2_1r;
    result = result + v1_1i*v2_1i;
    result = result + v1_2r*v2_2r;
    result = result + v1_2i*v2_2i;    

  }
  
  *Out=(REAL64)result;
}




} // namespace QDP;

#endif  // defined(__GNUC__)
