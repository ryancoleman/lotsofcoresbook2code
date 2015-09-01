
/*! @file
 * @brief Scalar-like architecture specific routines
 * 
 * Routines common to all scalar-like architectures including scalar and parscalar
 */


#include "qdp.h"
#include "qdp_util.h"

namespace QDP {

//-----------------------------------------------------------------------------
namespace Layout
{
  //! coord[mu]  <- mu  : fill with lattice coord in mu direction
  /* Assumes no inner grid */
  LatticeInteger latticeCoordinate(int mu)
  {
    const int nodeSites = Layout::sitesOnNode();
    const int nodeNumber = Layout::nodeNumber();
    LatticeInteger d;

    if (mu < 0 || mu >= Nd)
      QDP_error_exit("dimension out of bounds");

    /* This pragma is from Jacques... Should be OK, since each i is independent
     * no danger of concurrent writes */
#pragma omp parallel for
    for(int i=0; i < nodeSites; ++i) 
    {
      Integer cc = Layout::siteCoords(nodeNumber,i)[mu];
      d.elem(i) = cc.elem();
    }

    return d;
  }
}


//-----------------------------------------------------------------------------
// IO routine solely for debugging. Only defined here
template<class T>
std::ostream& operator<<(std::ostream& s, const multi1d<T>& s1)
{
  for(int i=0; i < s1.size(); ++i)
    s << " " << s1[i];

  return s;
}


//-----------------------------------------------------------------------------
//! Constructor from a function object
void Set::make(const SetFunc& fun)
{
  int nsubset_indices = fun.numSubsets();
  const int nodeSites = Layout::sitesOnNode();
  const int nodeNumber = Layout::nodeNumber();

#if QDP_DEBUG >= 2
  QDP_info("Set a subset: nsubset = %d",nsubset_indices);
#endif

  // This actually allocates the subsets
  sub.resize(nsubset_indices);

  // Create the space of the colorings of the lattice
  lat_color.resize(nodeSites);

  // Create the array holding the array of sitetable info
  sitetables.resize(nsubset_indices);

  // Loop over linear sites determining their color
  /* This OMP pragma added by Jacques. Should be OK since in the end 
     each value of linear is independent */
#pragma omp parallel for
  for(int linear=0; linear < nodeSites; ++linear)
  {
    multi1d<int> coord = Layout::siteCoords(nodeNumber, linear);

    int node   = Layout::nodeNumber(coord);
    int lin    = Layout::linearSiteIndex(coord);
    int icolor = fun(coord);

#if QDP_DEBUG >= 3
    std::cerr<<"linear="<<linear<<" coord="<<coord<<" node="<<node<<" col="<<icolor << std::endl;
#endif

    // Sanity checks
    if (node != nodeNumber)
      QDP_error_exit("Set: found site with node outside current node!");

    if (lin != linear)
      QDP_error_exit("Set: inconsistent linear sites");

    if (icolor < 0 || icolor >= nsubset_indices)
      QDP_error_exit("Set: coloring is outside legal range: color[%d]=%d",linear,icolor);

    // The coloring of this linear site
    lat_color[linear] = icolor;
  }
  
  
  /*
   * Loop over the lexicographic sites.
   * This implementation of the Set will always use a
   * sitetable.
   */

  /* NB: This OMP parallel added by Jacques.
   * Should be OK because each subset is independent, 
   * tho the number of subsets may be small... so scope for improvement
   * from threading may be limited */

#pragma omp parallel for
  for(int cb=0; cb < nsubset_indices; ++cb)
  {
    // Always construct the sitetables. 

    // First loop and see how many sites are needed
    int num_sitetable = 0;

    /* FIXME: This is a 'histogram' -- not yet threaded */
    for(int linear=0; linear < nodeSites; ++linear) {
      if (lat_color[linear] == cb) {
	++num_sitetable;
      }
    }
    // Now take the inverse of the lattice coloring to produce
    // the site list
    multi1d<int>& sitetable = sitetables[cb];
    sitetable.resize(num_sitetable);


    // Site ordering stuff for later
    bool ordRep;
    int start, end;

    // Handle the case that there are no sites
    if (num_sitetable > 0)
    {
      // For later sanity, initialize this to something 
      for(int i=0; i < num_sitetable; ++i)
	sitetable[i] = -1;

      for(int linear=0, j=0; linear < nodeSites; ++linear)
	if (lat_color[linear] == cb)
	  sitetable[j++] = linear;


      // Check *if* this coloring is contiguous and find the start
      // and ending sites
      ordRep = true;
      start = sitetable[0];   // this is the beginning
      end = sitetable[sitetable.size()-1];  // the absolute last site
      
      // Now look for a hole
      for(int prev=sitetable[0], i=0; i < sitetable.size(); ++i)
	if (sitetable[i] != prev++)
	{
#if QDP_DEBUG >= 2
	  QDP_info("Set(%d): sitetable[%d]=%d",cb,i,sitetable[i]);
#endif
	
	  // Found a hold. The rep is not ordered.
	  ordRep = false;
	  start = end = -1;
	  break;
	}
    }
    else  // num_sitetable == 0
    {
      ordRep = false;
      start = end = -1;
    }

    sub[cb].make(ordRep, start, end, &(sitetables[cb]), cb, this);

#if QDP_DEBUG >= 2
    QDP_info("Subset(%d)",cb);
#endif
  }
}
	  

} // namespace QDP;
