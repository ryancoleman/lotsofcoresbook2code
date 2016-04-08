// -*- C++ -*-
/*! @file
 * @brief QDPType after a subset
 *
 * Subclass of QDPType used for subset operations
 */

#ifndef QDP_QDPSUBTYPE_H
#define QDP_QDPSUBTYPE_H

namespace QDP {


//! QDPSubType - type representing a field living on a subset
/*! 
 * This class is meant to be an auxilliary class used only for
 * things like lvalues - left hand side of expressions, arguments
 * to calls that modify the source (like RNG), etc.
 */
template<class T, class C> 
class QDPSubType
{
  //! This is a type name like OSubLattice<T> or OSubScalar<T>
  typedef typename QDPSubTypeTrait<C>::Type_t CC;

public:
  //! Default constructor 
  QDPSubType() {}

  //! Copy constructor
  QDPSubType(const QDPSubType&) {}

  //! Destructor
  ~QDPSubType() {}


  //---------------------------------------------------------
  // Operators

  inline
  void assign(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  inline
  void assign(const Zero&)
  {
    if (getOwnsMemory()) {
      zero_rep_F(getF(),subset());
    } else {
      C tmp(getF(),1.0);
      zero_rep(tmp,subset());
    }
  }

  template<class T1,class C1>
  inline
  void assign(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpAssign(),PETE_identity(rhs),subset());
    }
  }


  template<class T1,class C1>
  inline
  void assign(const QDPSubType<T1,C1>& rhs)
  {
    if (getOwnsMemory() && rhs.getOwnsMemory()) {
      //std::cout << "own = own\n";
      if (subset().numSiteTable() != rhs.subset().numSiteTable())
	QDP_error_exit("assignment with incompatible subset sizes");
      for(int j=0; j < subset().numSiteTable(); ++j)
	getF()[j] = rhs.getF()[j];
    } else

#if 0
    if (!getOwnsMemory() && rhs.getOwnsMemory()) {
      //std::cout << "view = own\n";
      if (subset().numSiteTable() != rhs.subset().numSiteTable())
	QDP_error_exit("assignment with incompatible subset sizes");
      const int *tab = subset().siteTable().slice();
      for(int j=0; j < subset().numSiteTable(); ++j) {
	int i = tab[j];
	getF()[i] = rhs.getF()[j];
      }
    }
    if (getOwnsMemory() && !rhs.getOwnsMemory()) {
      //std::cout << "own = view\n";
      if (subset().numSiteTable() != rhs.subset().numSiteTable())
	QDP_error_exit("assignment with incompatible subset sizes");
      const int *tab = rhs.subset().siteTable().slice();
      for(int j=0; j < rhs.subset().numSiteTable(); ++j) {
	int i = tab[j];
	getF()[j] = rhs.getF()[i];
      }
    }
    if (!getOwnsMemory() && !rhs.getOwnsMemory())
#endif
      QDP_error_exit("assignment of two view subtypes is not supported");
  }


  template<class T1,class C1>
  inline
  void assign(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpAssign(),rhs,subset());
    }
  }

  inline
  void operator+=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpAddAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpAddAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator+=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpAddAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpAddAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator+=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpAddAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpAddAssign(),rhs,subset());
    }
  }

  inline
  void operator-=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpSubtractAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpSubtractAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator-=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpSubtractAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpSubtractAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator-=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpSubtractAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpSubtractAssign(),rhs,subset());
    }
  }

  inline
  void operator*=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpMultiplyAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpMultiplyAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator*=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpMultiplyAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpMultiplyAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator*=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpMultiplyAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpMultiplyAssign(),rhs,subset());
    }
  }

  inline
  void operator/=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpDivideAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpDivideAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator/=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpDivideAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpDivideAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator/=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpDivideAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpDivideAssign(),rhs,subset());
    }
  }

  inline
  void operator%=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpModAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpModAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator%=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpModAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpModAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator%=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpModAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpModAssign(),rhs,subset());
    }
  }

  inline
  void operator|=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseOrAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseOrAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator|=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseOrAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseOrAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator|=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseOrAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseOrAssign(),PETE_identity(rhs),subset());
    }
  }

  inline
  void operator&=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseAndAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseAndAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator&=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseAndAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseAndAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator&=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseAndAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseAndAssign(),rhs,subset());
    }
  }

  inline
  void operator^=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseXorAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseXorAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator^=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseXorAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseXorAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator^=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpBitwiseXorAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpBitwiseXorAssign(),rhs,subset());
    }
  }

  inline
  void operator<<=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpLeftShiftAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpLeftShiftAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator<<=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpLeftShiftAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpLeftShiftAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator<<=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpLeftShiftAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpLeftShiftAssign(),rhs,subset());
    }
  }


  inline
  void operator>>=(const typename WordType<C>::Type_t& rhs)
  {
    typedef typename SimpleScalar<typename WordType<C>::Type_t>::Type_t  Scalar_t;
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpRightShiftAssign(),PETE_identity(Scalar_t(rhs)),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpRightShiftAssign(),PETE_identity(Scalar_t(rhs)),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator>>=(const QDPType<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpRightShiftAssign(),PETE_identity(rhs),subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpRightShiftAssign(),PETE_identity(rhs),subset());
    }
  }

  template<class T1,class C1>
  inline
  void operator>>=(const QDPExpr<T1,C1>& rhs)
  {
    if (getOwnsMemory()) {
      evaluate_F(getF(),OpRightShiftAssign(),rhs,subset());
    } else {
      C tmp(getF(),1.0);
      evaluate(tmp,OpRightShiftAssign(),rhs,subset());
    }
  }

private:
  //! Hide default operator=
  inline
  C& operator=(const QDPSubType& rhs) {}

public:
  //C& field() {return static_cast<CC*>(this)->field();}
  bool getOwnsMemory() { return static_cast<CC*>(this)->getOwnsMemory(); }
  bool getOwnsMemory() const { return static_cast<const CC*>(this)->getOwnsMemory(); }
  T* getF() {return static_cast<CC*>(this)->getF();}
  T* getF() const {return static_cast<CC const *>(this)->getF();}
  const Subset& subset() const {return static_cast<const CC*>(this)->subset();}

};


} // namespace QDP

#endif



