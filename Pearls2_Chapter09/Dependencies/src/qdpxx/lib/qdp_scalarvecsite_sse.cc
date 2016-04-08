
/*! @file
 * @brief Intel SSE optimizations
 * 
 * SSE optimizations of basic operations
 */


#include "qdp.h"

// These SSE asm instructions are only supported under GCC/G++
#if defined(__GNUC__) && __GNUC_MINOR__ >= 2  &&  QDP_USE_SSE == 1

namespace QDP {

// Specialization to optimize the case   
//    LatticeColorMatrix[Subset] = LatticeColorMatrix * LatticeColorMatrix
template<>
void evaluate(OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >& d, 
	      const OpAssign& op, 
	      const QDPExpr<BinaryNode<OpMultiply, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > >, 
	      Reference<QDPType<PScalar<PColorMatrix<RComplexFloat, 3> >, 
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > > > >,
	      OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > > >& rhs,
	      const Subset& s)
{
//  cout << "call single site QDP_M_eq_M_times_M" << std::endl;

  typedef OLattice<PScalar<PColorMatrix<RComplexFloat, 3> > >  C;

  const C& l = static_cast<const C&>(rhs.expression().left());
  const C& r = static_cast<const C&>(rhs.expression().right());

  if( s.hasOrderedRep()) { 
    const int istart = s.start() >> INNER_LOG;
    const int iend   = s.end()   >> INNER_LOG;

    for(int i=istart; i <= iend; ++i) 
      {
	float *dd = (float*)&(d.elem(i).elem());
	float *ll = (float*)&(l.elem(i).elem());
	float *rr = (float*)&(r.elem(i).elem());
	
	_inline_ssevec_mult_su3_nn(dd,ll,rr,0);
	_inline_ssevec_mult_su3_nn(dd,ll,rr,1);
	_inline_ssevec_mult_su3_nn(dd,ll,rr,2);
      }
  }
  else { 
    // Do unoptimised - hope this still workl
    const int *tab = s.siteTable().slice();
    for(int j=0; j < s.numSiteTable(); j++) { 
      i = tab[j]; 
      d.elem(i).elem() = l.elem(i).elem() * r.elem(i).elem();
    }
  }

}

} // namespace QDP;

#endif  // defined(__GNUC__)
