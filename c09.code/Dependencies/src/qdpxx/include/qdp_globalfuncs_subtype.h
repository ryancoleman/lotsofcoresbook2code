// -*- C++ -*-

#ifndef QDP_GLOBALFUNCS_SUBTYPE_H
#define QDP_GLOBALFUNCS_SUBTYPE_H

namespace QDP
{


  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPSubType<T1,C1>& s1, const QDPType<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }
  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPType<T1,C1>& s1, const QDPSubType<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }
#if 0
  template<class T1, class C1, class T2, class C2>
  inline typename BinaryReturn<C1, C2, FnInnerProduct>::Type_t
  innerProduct(const QDPSubType<T1,C1>& s1, const QDPSubType<T2,C2>& s2)
  {
    return sum(localInnerProduct(s1,s2));
  }
#endif

  template<class T1,class C1,class T2,class C2>
  typename QDPSubTypeTrait< typename BinaryReturn<C1,C2,FnLocalInnerProduct>::Type_t >::Type_t
  localInnerProduct(const QDPSubType<T1,C1> & l,const QDPType<T2,C2> & r)
  {
    if (!l.getOwnsMemory())
      QDP_error_exit("localInnerProduct with subtype view called");

    typename QDPSubTypeTrait< typename BinaryReturn<C1,C2,FnLocalInnerProduct>::Type_t >::Type_t ret;
    ret.setSubset( l.subset() );

    //QDP_info("localInnerProduct %d sites",l.subset().numSiteTable());

    const int *tab = l.subset().siteTable().slice();
    for(int j=0; j < l.subset().numSiteTable(); ++j)
      {
	int i = tab[j];
	FnLocalInnerProduct op;
	ret.getF()[j] = op( l.getF()[j] , r.elem(i) );
      }

    return ret;
  }

  template<class T1,class C1,class T2,class C2>
  typename QDPSubTypeTrait< typename BinaryReturn<C1,C2,FnLocalInnerProduct>::Type_t >::Type_t
  localInnerProduct(const QDPType<T1,C1> & l,const QDPSubType<T2,C2> & r)
  {
    if (!r.getOwnsMemory())
      QDP_error_exit("localInnerProduct with subtype view called");

    typename QDPSubTypeTrait< typename BinaryReturn<C1,C2,FnLocalInnerProduct>::Type_t >::Type_t ret;
    ret.setSubset( r.subset() );

    const int *tab = r.subset().siteTable().slice();
    for(int j=0; j < r.subset().numSiteTable(); ++j)
      {
	int i = tab[j];
	FnLocalInnerProduct op;
	ret.getF()[j] = op( l.elem(i) , r.getF()[j] );
      }

    return ret;
  }


#if 0
  template<class T1,class C1,class T2,class C2>
  typename QDPSubTypeTrait< typename BinaryReturn<C1,C2,FnLocalInnerProduct>::Type_t >::Type_t
  localInnerProduct(const QDPSubType<T1,C1> & l,const QDPSubType<T2,C2> & r)
  {
    if (!l.getOwnsMemory())
      QDP_error_exit("localInnerProduct with subtype view called");
    if (!r.getOwnsMemory())
      QDP_error_exit("localInnerProduct with subtype view called");
    if (r.subset().numSiteTable() != l.subset().numSiteTable())
      QDP_error_exit("localInnerProduct with incompatible subset sizes");

    typename QDPSubTypeTrait< typename BinaryReturn<C1,C2,FnLocalInnerProduct>::Type_t >::Type_t ret;
    ret.setSubset( r.subset() );

    for(int j=0; j < r.subset().numSiteTable(); ++j)
      {
	FnLocalInnerProduct op;
	ret.getF()[j] = op( l.getF()[j] , r.getF()[j] );
      }

    return ret;
  }
#endif



template<class T>
typename UnaryReturn<OLattice<T>, FnSum>::Type_t
sum( const OSubLattice<T>& s1 )
{
  typename UnaryReturn<OLattice<T>, FnSum>::Type_t  d;

  // Must initialize to zero since we do not know if the loop will be entered
  zero_rep(d.elem());

  for(int j=0; j < s1.subset().numSiteTable(); ++j) 
    {
      d.elem() += s1.getF()[j];
    }

  // Do a global sum on the result
  QDPInternal::globalSum(d);

  return d;
 }


  template<class T> 
  inline
  void zero_rep_F( T* dest, const Subset& s)
  {
    for(int j=0; j < s.numSiteTable(); ++j) 
      {
	zero_rep( dest[j] );
      }
  }


  //! dest  = 0 
  template<class T>
  void zero_rep(OSubLattice<T> dd) 
  {
    if (dd.getOwnsMemory()) {
      zero_rep_F(dd.getF(),dd.subset());
    } else {
      OLattice<T> tmp(dd.getF(),1.0);
      zero_rep(tmp,dd.subset());
    }
  }


} // namespace

#endif
