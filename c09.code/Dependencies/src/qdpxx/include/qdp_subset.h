// -*- C++ -*-
/*! @file
 * @brief Sets and subsets
 */

#ifndef QDP_SUBSET_H
#define QDP_SUBSET_H

namespace QDP {

/*! @defgroup subsets Sets and Subsets
 *
 * Sets are the objects that facilitate colorings of the lattice.
 * Subsets are groups of sites that are all of one color.
 * Subsets (and in a few cases Sets) can be used to restrict operations
 * to only a particular coloring of the lattice.
 *
 * @{
 */

//-----------------------------------------------------------------------
//! SetMap
class SetFunc
{
public:
  // Virtual destructor to stop compiler warnings - no cleanup needed
  virtual ~SetFunc() {}
  virtual int operator() (const multi1d<int>& coordinate) const = 0;
  virtual int numSubsets() const = 0;
};

//-----------------------------------------------------------------------
// Forward declaration
class Set;

//-----------------------------------------------------------------------
//! Subsets - controls how lattices are looped
class Subset 
{
public:
  //! There can be an empty constructor
  Subset() {}

  //! Copy constructor
  Subset(const Subset& s):
    ordRep(s.ordRep), startSite(s.startSite), endSite(s.endSite), 
    sub_index(s.sub_index), sitetable(s.sitetable), set(s.set)
    {}

  // Simple constructor
  void make(const Subset& s);

  //! Destructor for a subset
  virtual ~Subset() {}

  //! The = operator
  Subset& operator=(const Subset& s);

  //! Access the coloring for this subset
  int color() const {return sub_index;}

protected:
  // Simple constructor
  void make(bool rep, int start, int end, multi1d<int>* ind, int cb, Set* set);

private:
  bool ordRep;
  int startSite;
  int endSite;
  int sub_index;

  //! Site lookup table
  multi1d<int>* sitetable;

  //! Original set
  Set *set;

public:
  inline bool hasOrderedRep() const {return ordRep;}
  inline int start() const {return startSite;}
  inline int end() const {return endSite;}

  const multi1d<int>& siteTable() const {return *sitetable;}
  inline int numSiteTable() const {return sitetable->size();}

  //! The super-set of this subset
  const Set& getSet() const { return *set; }

  friend class Set;
};


//-----------------------------------------------------------------------
//! Set - collection of subsets controlling which sites are involved in an operation
class Set 
{
public:
  //! There can be an empty constructor
  Set() {}

  //! Constructor from a function object
  Set(const SetFunc& fn) {make(fn);}

  //! Constructor from a function object
  void make(const SetFunc& fn);

  //! Index operator selects a subset from a set
  const Subset& operator[](int subset_index) const {return sub[subset_index];}

  //! Return number of subsets
  int numSubsets() const {return sub.size();}

  //! Destructor for a set
  virtual ~Set() {}

  //! The = operator
  Set& operator=(const Set& s);

protected:
  //! A set is composed of an array of subsets
  multi1d<Subset> sub;

  //! Index or color array of lattice
  multi1d<int> lat_color;

  //! Array of sitetable arrays
  multi1d<multi1d<int> > sitetables;

public:
  //! The coloring of the lattice sites
  const multi1d<int>& latticeColoring() const {return lat_color;}
};



//-----------------------------------------------------------------------
//! Default all subset
extern Subset all;


//! Experimental 3d checkerboarding for temp_precond
extern Set rb3;

//! Default 2-checkerboard (red/black) subset
extern Set rb;

//! Default 2^{Nd+1}-checkerboard subset. Useful for pure gauge updating.
extern Set mcb;

//! Default even subset
extern Subset even;

//! Default odd subset
extern Subset odd;


/** @} */ // end of group subsetss

} // namespace QDP

#endif
