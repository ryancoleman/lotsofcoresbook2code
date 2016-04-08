// -*- C++ -*-
/*! \file
 * \brief Outer grid classes after a subset
 */


#ifndef QDP_OUTERSUBTYPE_H
#define QDP_OUTERSUBTYPE_H

namespace QDP {

//! OScalar class narrowed to a subset
/*! 
 * Only used for lvalues
 */

template<class T>
class OSubScalar: public QDPSubType<T, OScalar<T> >
{
  typedef OScalar<T> C;

public:
  //OSubScalar(OScalar<T>& a, const Subset& ss): F(a.getF()), s(ss) {}
  OSubScalar(OScalar<T>& a, const Subset& ss): F(a.getF()), s(&(const_cast<Subset&>(ss))) {}
  OSubScalar(const OSubScalar& a): F(a.F), s(a.s) {}
  ~OSubScalar() {}

  //---------------------------------------------------------
  // Operators
  // NOTE: all assignment-like operators except operator= are
  // inherited from QDPSubType

  inline
  void operator=(const typename WordType<T>::Type_t& rhs)
    {
      this->assign(rhs);
    }

  inline
  void operator=(const Zero& rhs)
    {
      this->assign(rhs);
    }

  template<class T1,class C1>
  inline
  void operator=(const QDPType<T1,C1>& rhs)
    {
      this->assign(rhs);
    }

  template<class T1,class C1>
  inline
  void operator=(const QDPExpr<T1,C1>& rhs)
    {
      this->assign(rhs);
    }


  inline
  void operator=(const OSubScalar& rhs)
    {
      this->assign(rhs);
    }


private:
  // Hide default constructor
  OSubScalar() {}

public:
  T* getF() {return F;}
  T* getF() const {return F;}
  const Subset& subset() const {return *s;}
  bool getOwnsMemory() const { return ownsMemory; }
  bool getOwnsMemory() { return ownsMemory; }

private:
  T* F;
  bool ownsMemory;
  Subset* s;

  //const Subset& s;
};



//-------------------------------------------------------------------------------------
//! OLattice class narrowed to a subset
/*! 
 * Only used for lvalues
 */
template<class T> 
class OSubLattice: public QDPSubType<T, OLattice<T> >
{
  typedef OLattice<T> C;

public:
  //OSubLattice(OLattice<T>& a, const Subset& ss): F(a.getF()), s(ss), ownsMemory(false) {}
  OSubLattice(OLattice<T>& a, const Subset& ss): F(a.getF()), s(&(const_cast<Subset&>(ss))), ownsMemory(false) {}
  OSubLattice(const OSubLattice& a) {
    ownsMemory = a.ownsMemory; 
    s = a.s;
    if (a.ownsMemory) {
      //std::cout << "OSubLattice copy ctor, must copy\n";
      alloc_mem();
      for( int i = 0 ; i < s->numSiteTable() ; ++i )
	F[i] = a.F[i];
    } else {
      F = a.F;
    }
  }


#if 1
  OSubLattice(): s(NULL), F(NULL), ownsMemory(false) {}


  void setSubset( const Subset& ss ) {
    if (F && !ownsMemory)
      QDP_error_exit("You try to set the subset on an OSubLattice that is a view of an OLattice!");
    if (F) {
      free_mem();
    }
    s = &(const_cast<Subset&>(ss));
    alloc_mem();
    ownsMemory = true;
  }


  OSubLattice(const Subset& ss , OLattice<T>& a): s(&(const_cast<Subset&>(ss))), ownsMemory(true) {
    alloc_mem();
    const int *tab = s->siteTable().slice();
    for( int j = 0 ; j < s->numSiteTable() ; ++j ) {
      int i = tab[j];
      F[j] = a.elem(i);
    }
  }


  void alloc_mem() {
    //QDP_info("OSubLattice alloc for %d sites",s->numSiteTable());
    try 
      {
	F = (T*)QDP::Allocator::theQDPAllocator::Instance().allocate(sizeof(T)*s->numSiteTable(),QDP::Allocator::DEFAULT);
      }
      catch(std::bad_alloc) 
      {
	QDPIO::cerr << "Allocation failed in OSubLattice" << std::endl;
	QDP::Allocator::theQDPAllocator::Instance().dump();
	QDP_abort(1);
      }
  }

  ~OSubLattice() {
    if (ownsMemory) {
      free_mem();
    }
  }

  void free_mem() {
    QDP::Allocator::theQDPAllocator::Instance().free(F);
  }
#endif




  //---------------------------------------------------------
  // Operators
  // NOTE: all assignment-like operators except operator= are
  // inherited from QDPType

  inline
  void operator=(const typename WordType<T>::Type_t& rhs)
    {
      this->assign(rhs);
    }

  inline
  void operator=(const Zero& rhs)
    {
      this->assign(rhs);
    }

  template<class T1,class C1>
  inline
  void operator=(const QDPType<T1,C1>& rhs)
    {
      this->assign(rhs);
    }

  template<class T1,class C1>
  inline
  void operator=(const QDPExpr<T1,C1>& rhs)
    {
      this->assign(rhs);
    }


  inline
  void operator=(const OSubLattice& rhs)
    {
      this->assign(rhs);
    }


private:
  // Hide default constructor
  //OSubLattice() {}

public:
  bool getOwnsMemory() const { return ownsMemory; }
  bool getOwnsMemory() { return ownsMemory; }
  T* getF() {return F;}
  T* getF() const {return F;}
  //const Subset& subset() const {return s;}
  const Subset& subset() const {return *s;}

private:
  T* F;
  bool ownsMemory;
  Subset* s;

  //C&      F;
  //const Subset& s;
};


//-----------------------------------------------------------------------------
// Traits class for returning the subset-ted class name of a outer grid class
//-----------------------------------------------------------------------------

template<class T>
struct QDPSubTypeTrait<OScalar<T> > 
{
  typedef OSubScalar<T>  Type_t;
};


template<class T>
struct QDPSubTypeTrait<OLattice<T> > 
{
  typedef OSubLattice<T>  Type_t;
};


//-----------------------------------------------------------------------------
// Traits classes to support operations of simple scalars (floating constants, 
// etc.) on QDPTypes
//-----------------------------------------------------------------------------

template<class T>
struct WordType<OSubScalar<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};


template<class T>
struct WordType<OSubLattice<T> > 
{
  typedef typename WordType<T>::Type_t  Type_t;
};


// ------------------------------------------------------------
// Get Single Precision Types of OuterSubType templates
// ------------------------------------------------------------
template<class T>
struct SinglePrecType<OSubScalar<T> > 
{
  typedef OSubScalar<typename SinglePrecType<T>::Type_t>  Type_t;
};


template<class T>
struct SinglePrecType<OSubLattice<T> > 
{
  typedef OSubLattice<typename SinglePrecType<T>::Type_t>  Type_t;
};


// ------------------------------------------------------------
// Get Single Precision Types of OuterSubType templates
// ------------------------------------------------------------
template<class T>
struct DoublePrecType<OSubScalar<T> > 
{
  typedef OSubScalar<typename DoublePrecType<T>::Type_t>  Type_t;
};

template<class T>
struct DoublePrecType<OSubLattice<T> > 
{
  typedef OSubLattice<typename DoublePrecType<T>::Type_t>  Type_t;
};


//-----------------------------------------------------------------------------
// Scalar Operations
//-----------------------------------------------------------------------------

//! dest = 0
template<class T> 
void zero_rep(OScalar<T>& dest, const Subset& s) 
{
  zero_rep(dest.field().elem());
}

//! dest = 0
template<class T>
void zero_rep(OSubScalar<T> dest) 
{
  zero_rep(dest.field().elem());
}

//! dest = (mask) ? s1 : dest
template<class T1, class T2> 
void copymask(OSubScalar<T2> dest, const OScalar<T1>& mask, 
	      const OScalar<T2>& s1) 
{
  copymask(dest.field().elem(), mask.elem(), s1.elem());
}


//-----------------------------------------------------------------------------
// Random numbers
//! dest  = random  
/*! Implementation is in the specific files */
template<class T>
void random(OSubScalar<T> d);

//! dest  = gaussian
template<class T>
void gaussian(OSubScalar<T> dd)
{
  OLattice<T>& d = dd.field();
  const Subset& s = dd.subset();

  OScalar<T>  r1, r2;

  random(r1(s));
  random(r2(s));

  fill_gaussian(d.elem(), r1.elem(), r2.elem());
}

} // namespace QDP

#endif
