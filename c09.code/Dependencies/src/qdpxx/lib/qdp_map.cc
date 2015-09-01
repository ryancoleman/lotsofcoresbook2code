
/*! @file
 * @brief Support routines for Maps
 * 
 * Common routines for implementation of Maps
 */

#include "qdp.h"
#include "qdp_util.h"

namespace QDP {

//! Definition of shift function object
ArrayBiDirectionalMap  shift;

//! Function object used for constructing the default nearest neighbor map
struct NearestNeighborMapFunc : public ArrayMapFunc
{
  NearestNeighborMapFunc() {}

  // Virtual destructor - no cleanup needed
  virtual ~NearestNeighborMapFunc() {} 

  virtual multi1d<int> operator() (const multi1d<int>& coord, int sign, int dir) const
    {
      multi1d<int> lc = coord;

      const multi1d<int>& nrow = Layout::lattSize();
      lc[dir] = (coord[dir] + sgnum(sign) + 4*nrow[dir]) % nrow[dir];

      return lc;
    }

  virtual int numArray() const {return Nd;}

private:
  int sgnum(int x) const {return (x > 0) ? 1 : -1;}
}; 


//! Initializer for maps
void initDefaultMaps()
{
  // Initialize the nearest neighbor map
  NearestNeighborMapFunc bbb;

  shift.make(bbb);
}


//----------------------------------------------------------------------------
// ArrayMap

// This class is is used for binding the direction index of an ArrayMapFunc
// so as to construct a MapFunc
struct PackageArrayMapFunc : public MapFunc
{
  PackageArrayMapFunc(const ArrayMapFunc& mm, int dd): pmap(mm), dir(dd) {}

  // Virtual Destructor - no cleanup needed
  virtual ~PackageArrayMapFunc() {}
  virtual multi1d<int> operator() (const multi1d<int>& coord, int isign) const
    {
      return pmap(coord, isign, dir);
    }

private:
  const ArrayMapFunc& pmap;
  const int dir;
}; 


//! Initializer for an array of maps
void ArrayMap::make(const ArrayMapFunc& func)
{
  // We are allowed to declare a mapsa, but not allocate one.
  // There is an empty constructor for Map. Hence, the resize will
  // actually allocate the space.
  mapsa.resize(func.numArray());

  // Loop over each direction making the Map
  for(int dir=0; dir < func.numArray(); ++dir)
  {
    PackageArrayMapFunc  my_local_map(func,dir);

    mapsa[dir].make(my_local_map);
  }
}


//----------------------------------------------------------------------------
// BiDirectionalMap

// This class is is used for binding the direction index of an BiDirectionalMapFunc
// so as to construct a MapFunc
struct PackageBiDirectionalMapFunc : public MapFunc
{
  PackageBiDirectionalMapFunc(const MapFunc& mm, int mmult): pmap(mm), mult(mmult) {}

  // Virtual destructor -- no real cleanup needed
  virtual ~PackageBiDirectionalMapFunc() {}
  virtual multi1d<int> operator() (const multi1d<int>& coord, int isign) const
    {
      return pmap(coord, mult*isign);
    }

private:
  const MapFunc& pmap;
  const int mult;
}; 


//! Initializer for a generic bi-directional maps
void BiDirectionalMap::make(const MapFunc& func)
{
  // We are allowed to declare a bimaps, but not allocate one.
  // There is an empty constructor for Map. Hence, the resize will
  // actually allocate the space.
  bimaps.resize(2);

  // Construct maps for each sign
  PackageBiDirectionalMapFunc  my_neg_map(func,-1);
  bimaps[0].make(my_neg_map);

  PackageBiDirectionalMapFunc  my_pos_map(func,+1);
  bimaps[1].make(my_pos_map);
}


//----------------------------------------------------------------------------
// ArrayBiDirectionalMap

// This class is is used for binding the direction index of an ArrayBiDirectionalMapFunc
// so as to construct a MapFunc
struct PackageArrayBiDirectionalMapFunc : public MapFunc
{
  PackageArrayBiDirectionalMapFunc(const ArrayMapFunc& mm, int mmult, int dd) : 
    pmap(mm), mult(mmult), dir(dd) {}

  virtual multi1d<int> operator() (const multi1d<int>& coord, int isign) const
    {
      return pmap(coord, mult*isign, dir);
    }

private:
  const ArrayMapFunc& pmap;
  const int mult;
  const int dir;
}; 


//! Initializer for an array of bi-directional maps
void ArrayBiDirectionalMap::make(const ArrayMapFunc& func)
{
  // We are allowed to declare a mapsa, but not allocate one.
  // There is an empty constructor for Map. Hence, the resize will
  // actually allocate the space.
  bimapsa.resize(2,func.numArray());

  // Loop over each direction making the Map
  for(int dir=0; dir < func.numArray(); ++dir)
  {
    // Construct maps for each sign
    PackageArrayBiDirectionalMapFunc  my_neg_map(func,-1,dir);
    bimapsa(0,dir).make(my_neg_map);

    PackageArrayBiDirectionalMapFunc  my_pos_map(func,+1,dir);
    bimapsa(1,dir).make(my_pos_map);
  }
}


//-----------------------------------------------------------------------------


} // namespace QDP;
