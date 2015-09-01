
/*! @file
 * @brief Scalarvec layout routines
 * 
 * Layout routines for scalarvec implementation
 * QDP data parallel interface
 *
 * Layout
 *
 * This routine provides various layouts, including
 *    lexicographic
 *    2-checkerboard  (even/odd-checkerboarding of sites)
 *    32-style checkerboard (even/odd-checkerboarding of hypercubes)
 */

#include "qdp.h"
#include "qdp_util.h"

namespace QDP 
{

//-----------------------------------------------------------------------------
// Layout stuff specific to a scalarvec architecture
  namespace Layout
  {
    //-----------------------------------------------------
    //! Local data specific to a scalarvec architecture
    /*! 
     * NOTE: the disadvantage to using a struct to keep things together is
     * that subsequent groupings of namespaces can not just add onto the
     * current namespace. 
     */
    struct LocalLayout_t
    {
      //! Total lattice volume
      int vol;

      //! Lattice size
      multi1d<int> nrow;

      //! Subgrid lattice volume
      int subgrid_vol;

      //! Outergrid lattice volume
      int olattice_vol;

      //! Logical node coordinates
      multi1d<int> logical_coord;

      //! Logical system size
      multi1d<int> logical_size;

    } _layout;


    //-----------------------------------------------------
    // Functions

    //! Main destruction routine
    void destroy() {}

    //! Set virtual grid (problem grid) lattice size
    void setLattSize(const multi1d<int>& nrows) {_layout.nrow = nrows;}

    //! Set SMP flag -- true if using smp/multiprocessor mode on a node
    /*! For now, this is ignored */
    void setSMPFlag(bool flag) {}

    //! Set number of processors in a multi-threaded implementation
    /*! For now, this is ignored */
    void setNumProc(int N) {}

    //! Virtual grid (problem grid) lattice size
    const multi1d<int>& lattSize() {return _layout.nrow;}

    //! Total lattice volume
    int vol() {return _layout.vol;}

    //! Subgrid lattice volume
    int sitesOnNode() {return _layout.subgrid_vol;}

    //! Subgrid lattice volume
    int outerSitesOnNode() {return _layout.olattice_vol;}

    //! Returns whether this is the primary node
    /*! Always true on a scalarvec platform */
    bool primaryNode() {return true;}

    //! Subgrid (grid on each node) lattice size
    const multi1d<int>& subgridLattSize() {return _layout.nrow;}

    //! Returns the node number of this node
    int nodeNumber() {return 0;}

    //! Returns the logical node number for the corresponding lattice coordinate
    int nodeNumber(const multi1d<int>& coord) {return 0;}

    //! Returns the number of nodes
    int numNodes() {return 1;}

    //! Returns the logical node coordinates for this node
    const multi1d<int>& nodeCoord() {return _layout.logical_coord;}

    //! Returns the logical size of this machine
    const multi1d<int>& logicalSize() {return _layout.logical_size;}

    //! Returns the node number given some logical node coordinate
    /*! This is not meant to be speedy */
    int getNodeNumberFrom(const multi1d<int>& node_coord) {return 0;}

    //! Returns the logical node coordinates given some node number
    /*! This is not meant to be speedy */
    multi1d<int> getLogicalCoordFrom(int node) 
    {
      multi1d<int> node_coord(Nd);
      node_coord = 0;
      return node_coord;
    }

    //! Initializer for layout
    void init() {}

    //! The linearized site index for the corresponding lexicographic site
    int linearSiteIndex(int lexicosite)
    {
      return linearSiteIndex(crtesn(lexicosite, lattSize()));
    }

    //! Initializer for all the layout defaults
    void initDefaults()
    {
      // Default set and subsets
      initDefaultSets();

      // Default maps
      initDefaultMaps();

      // Initialize RNG
      RNG::initDefaultRNG();

      // Set default profile level
      setProfileLevel(getProgramProfileLevel());
    }

    //! Initializer for layout
    void create()
    {
      if ( ! QDP_isInitialized() )
	QDP_error_exit("QDP is not initialized");

      if (_layout.nrow.size() != Nd)
	QDP_error_exit("dimension of lattice size not the same as the default");

      _layout.vol=1;
      for(int i=0; i < Nd; ++i) 
	_layout.vol *= _layout.nrow[i];
      _layout.subgrid_vol  = _layout.vol;
      _layout.olattice_vol = _layout.vol >> INNER_LOG;
  
      _layout.logical_coord.resize(Nd);
      _layout.logical_size.resize(Nd);

      _layout.logical_coord = 0;
      _layout.logical_size = 1;

#if QDP_DEBUG >= 2
      fprintf(stderr,"vol=%d\n",_layout.vol);
#endif

      // This implementation requires there be a multiple of INNER_LEN sites 
      if (_layout.vol % INNER_LEN != 0)
	QDP_error_exit("Layout::create() - this scalarvec implementation requires there be a multiple of %d sites", INNER_LEN);


      // Sanity check - check the layout functions make sense
      for(int i=0; i < vol(); ++i) 
      {
	int ii = Layout::linearSiteIndex(Layout::siteCoords(Layout::nodeNumber(),i));
	if (i != ii)
	  QDP_error_exit("Layout::create - Layout problems, the layout functions do not work correctly with this lattice size");
      }

      // Initialize various defaults
      initDefaults();

      QDP_info("Finished lattice layout");
    }
  }


//-----------------------------------------------------------------------------
#if QDP_USE_LEXICO_LAYOUT == 1

#warning "Using a lexicographic layout"

  namespace Layout
  {
    //! Reconstruct the lattice coordinate from the node and site number
    /*! 
     * This is the inverse of the nodeNumber and linearSiteIndex functions.
     * The API requires this function to be here.
     */
    multi1d<int> siteCoords(int node, int linearsite) // ignore node
    {
      return crtesn(linearsite, lattSize());
    }

    //! The linearized site index for the corresponding coordinate
    /*! This layout is a simple lexicographic lattice ordering */
    int linearSiteIndex(const multi1d<int>& coord)
    {
      return local_site(coord, lattSize());
    }
  }


//-----------------------------------------------------------------------------

#elif QDP_USE_CB2_LAYOUT == 1

#warning "Using a 2 checkerboard (red/black) layout"

  namespace Layout
  {
    //! Reconstruct the lattice coordinate from the node and site number
    /*! 
     * This is the inverse of the nodeNumber and linearSiteIndex functions.
     * The API requires this function to be here.
     */
    multi1d<int> siteCoords(int node, int linearsite) // ignore node
    {
      int vol_cb = vol() >> 1;
      multi1d<int> cb_nrow = lattSize();
      cb_nrow[0] >>= 1;

      int cb = linearsite / vol_cb;
      multi1d<int> coord = crtesn(linearsite % vol_cb, cb_nrow);

      int cbb = cb;
      for(int m=1; m<coord.size(); ++m)
	cbb += coord[m];
      cbb = cbb & 1;

      coord[0] = 2*coord[0] + cbb;

      return coord;
    }

    //! The linearized site index for the corresponding coordinate
    /*! This layout is appropriate for a 2 checkerboard (red/black) lattice */
    int linearSiteIndex(const multi1d<int>& coord)
    {
      int vol_cb = vol() >> 1;
      multi1d<int> cb_nrow = lattSize();
      cb_nrow[0] >>= 1;

      multi1d<int> cb_coord = coord;

      cb_coord[0] >>= 1;    // Number of checkerboards
    
      int cb = 0;
      for(int m=0; m<coord.size(); ++m)
	cb += coord[m];
      cb = cb & 1;

      return local_site(cb_coord, cb_nrow) + cb*vol_cb;
    }
  }

//-----------------------------------------------------------------------------

#elif QDP_USE_CB32_LAYOUT == 1

#warning "Using a 32 checkerboard layout"

  namespace Layout
  {
    //! Reconstruct the lattice coordinate from the node and site number
    /*! 
     * This is the inverse of the nodeNumber and linearSiteIndex functions.
     * The API requires this function to be here.
     */
    multi1d<int> siteCoords(int node, int linearsite) // ignore node
    {
      int vol_cb = vol() >> (Nd+1);
      multi1d<int> cb_nrow(Nd);
      cb_nrow[0] = lattSize()[0] >> 2;
      for(int i=1; i < Nd; ++i) 
	cb_nrow[i] = lattSize()[i] >> 1;

      int subl = linearsite / vol_cb;
      multi1d<int> coord = crtesn(linearsite % vol_cb, cb_nrow);

      int cb = 0;
      for(int m=1; m<Nd; ++m)
	cb += coord[m];
      cb &= 1;

      coord[0] <<= 2;
      for(int m=1; m<Nd; ++m)
	coord[m] <<= 1;

      subl ^= (cb << Nd);
      for(int m=0; m<Nd; ++m)
	coord[m] ^= (subl & (1 << m)) >> m;
      coord[0] ^= (subl & (1 << Nd)) >> (Nd-1);   // this gets the hypercube cb

      return coord;
    }

    //! The linearized site index for the corresponding coordinate
    /*! This layout is appropriate for a 32-style checkerboard lattice */
    int linearSiteIndex(const multi1d<int>& coord)
    {
      int vol_cb = vol() >> (Nd+1);
      multi1d<int> cb_nrow(Nd);
      cb_nrow[0] = lattSize()[0] >> 2;
      for(int i=1; i < Nd; ++i) 
	cb_nrow[i] = lattSize()[i] >> 1;

      int subl = coord[Nd-1] & 1;
      for(int m=Nd-2; m >= 0; --m)
	subl = (subl << 1) + (coord[m] & 1);

      int cb = 0;
      for(int m=0; m < Nd; ++m)
	cb += coord[m] >> 1;

      subl += (cb & 1) << Nd;   // Final color or checkerboard

      // Construct the checkerboard lattice coord
      multi1d<int> cb_coord(Nd);

      cb_coord[0] = coord[0] >> 2;
      for(int m=1; m < Nd; ++m)
	cb_coord[m] = coord[m] >> 1;

      return local_site(cb_coord, cb_nrow) + subl*vol_cb;
    }
  }

#else

#error "no appropriate layout defined"

#endif

//-----------------------------------------------------------------------------


} // namespace QDP;
