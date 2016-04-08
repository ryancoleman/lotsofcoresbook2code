
/*! @file
 * @brief Parscalarvec layout routines
 * 
 * Layout routines for parscalarvec implementation
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

#include "qmp.h"

namespace QDP 
{

//-----------------------------------------------------------------------------
// IO routine solely for debugging. Only defined here
  template<class T>
  ostream& operator<<(ostream& s, const multi1d<T>& s1)
  {
    for(int i=0; i < s1.size(); ++i)
      s << " " << s1[i];

    return s;
  }


//-----------------------------------------------------------------------------
  namespace Layout
  {
    //-----------------------------------------------------
    //! Local data specific to all architectures
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

      //! Outer-subgrid lattice volume
      int olattice_vol;

      //! Subgrid lattice size
      multi1d<int> subgrid_nrow;

      //! Logical node coordinates
      multi1d<int> logical_coord;

      //! Logical system size
      multi1d<int> logical_size;

      //! Node rank
      int node_rank;

      //! Total number of nodes
      int num_nodes;
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
    bool primaryNode() {return (_layout.node_rank == 0) ? true : false;}

    //! Subgrid (grid on each node) lattice size
    const multi1d<int>& subgridLattSize() {return _layout.subgrid_nrow;}

    //! Returns the node number of this node
    int nodeNumber() {return _layout.node_rank;}

    //! Returns the number of nodes
    int numNodes() {return _layout.num_nodes;}

    //! Returns the logical node coordinates for this node
    const multi1d<int>& nodeCoord() {return _layout.logical_coord;}

    //! Returns the logical size of this machine
    const multi1d<int>& logicalSize() {return _layout.logical_size;}

    //! Returns the node number given some logical node coordinate
    /*! This is not meant to be speedy */
    int getNodeNumberFrom(const multi1d<int>& node_coord) 
    {
      return QMP_get_node_number_from(node_coord.slice());
    }

    //! Returns the logical node coordinates given some node number
    /*! This is not meant to be speedy */
    multi1d<int> getLogicalCoordFrom(int node) 
    {
      multi1d<int> node_coord(Nd);
      const int* node_crd = QMP_get_logical_coordinates_from(node);  // QMP mallocs here

      for(int i=0; i < Nd; ++i)
	node_coord[i] = node_crd[i];

      // Hackery as free cannot take a const void*, so grab
      // node_crd with a non const pointer
      int* non_const_node_crd = const_cast<int*>(node_crd);
      free(non_const_node_crd);   // free up QMP's memory
      return node_coord;
    }

    //! Return the smallest lattice size per node allowed
    multi1d<int> minimalLayoutMapping();

    //! Initializer for layout
    void init()
    {
      _layout.num_nodes = QMP_get_number_of_nodes();
      _layout.node_rank = QMP_get_node_number();
    }



    //! The linearized site index for the corresponding lexicographic site
    int linearSiteIndex(int site)
    { 
      multi1d<int> coord = crtesn(site, lattSize());
    
      return linearSiteIndex(coord);
    }


    //! Initializer for all the layout defaults
    void initDefaults()
    {
#if QDP_DEBUG >= 2
      QDP_info("Create default subsets");
#endif
      // Default set and subsets
      initDefaultSets();

      // Default maps
      initDefaultMaps();

      // Initialize RNG
      RNG::initDefaultRNG();

      // Set default profile level
      setProfileLevel(getProgramProfileLevel());
    }


    //! Main lattice creation routine
    void create()
    {
      if ( ! QDP_isInitialized() )
	QDP_error_exit("QDP is not initialized");

      if (_layout.nrow.size() != Nd)
	QDP_error_exit("dimension of lattice size not the same as the default");

      _layout.vol=1;
      for(int i=0; i < Nd; ++i) 
	_layout.vol *= _layout.nrow[i];
  
#if QDP_DEBUG >= 2
      QDP_info("vol=%d",_layout.vol);
#endif

#if QDP_DEBUG >= 2
      QDP_info("Initialize layout");
#endif

      // This implementation requires there be a multiple of INNER_LEN sites 
      if (_layout.vol % INNER_LEN != 0)
	QDP_error_exit("Layout::create() - this scalarvec implementation requires there be a multiple of %d sites", INNER_LEN);


      // Simple check - we insist here that the total volume is divisible
      // by the number of processors. Will also insist that the problem
      // size is regular on each node
      if (_layout.vol % numNodes() != 0)
	QDP_error_exit("Layout::create - problem size not divisible by number of processors");

    
      // Return the smallest lattice size per node allowed
      multi1d<int> min_dim = minimalLayoutMapping();

      // Another sanity check - if the specified lattice size is not a multiple
      // of the minimal size, then barf
      for(int i=0; i < Nd; ++i)
	if (_layout.nrow[i] % min_dim[i] != 0)
	  QDP_error_exit("Layout::create - Error: problem size not multiple of smallest size allowed for this type of layout");


      // Now, layout the machine. Note, if the logical_machine size was set previously
      // it will be used
      multi1d<int> nrow(Nd);
      for(int i=0; i < Nd; ++i)
	nrow[i] = _layout.nrow[i] / min_dim[i];

      int* nrow_slice=const_cast<int*>(nrow.slice());
      QMP_layout_grid(nrow_slice, Nd);


      // Pull out useful stuff
      const int* phys_size = QMP_get_logical_dimensions();
      const int* phys_coord = QMP_get_logical_coordinates();

      _layout.subgrid_nrow.resize(Nd);
      _layout.logical_coord.resize(Nd);
      _layout.logical_size.resize(Nd);

      _layout.subgrid_vol = 1;

      for(int i=0; i < Nd; ++i)
      {
	_layout.logical_coord[i] = phys_coord[i];
	_layout.logical_size[i] = phys_size[i];
	_layout.subgrid_nrow[i] = _layout.nrow[i] / _layout.logical_size[i];

	_layout.subgrid_vol *= _layout.subgrid_nrow[i];
      }

      _layout.olattice_vol = _layout.subgrid_vol >> INNER_LOG; // outer-subgrid volume


      // Diagnostics
      QDPIO::cout << "Lattice initialized:\n";
      QDPIO::cout << "  problem size =";
      for(int i=0; i < Nd; ++i)
	QDPIO::cout << " " << _layout.nrow[i];
      QDPIO::cout << std::endl;

      QDPIO::cout << "  logical machine size =";
      for(int i=0; i < Nd; ++i)
	QDPIO::cout << " " << _layout.logical_size[i];
      QDPIO::cout << std::endl;

      QDPIO::cout << "  logical node coord =";
      for(int i=0; i < Nd; ++i)
	QDPIO::cout << " " << _layout.logical_coord[i];
      QDPIO::cout << std::endl;

      QDPIO::cout << "  subgrid size =";
      for(int i=0; i < Nd; ++i)
	QDPIO::cout << " " << _layout.subgrid_nrow[i];
      QDPIO::cout << std::endl;

      QDPIO::cout << "  total volume = " << _layout.vol << std::endl;
      QDPIO::cout << "  subgrid volume = " << _layout.subgrid_vol << std::endl;

      // This implementation requires there be a multiple of INNER_LEN sites 
      // on a node
      if (_layout.subgrid_vol % INNER_LEN != 0)
	QDP_error_exit("Layout::create() - this parscalarvec implementation requires there be a multiple of %d sites", INNER_LEN);


      // Sanity check - check the QMP node number functions
      for(int node=0; node < Layout::numNodes(); ++node)
      { 
	multi1d<int> coord = Layout::getLogicalCoordFrom(node);
	int node2 = Layout::getNodeNumberFrom(coord);

	if (node != node2)
	  QDP_error_exit("Layout::create - Layout problems, the QMP logical to physical node map functions do not work correctly with this lattice size");
      }

      // Sanity check - check the layout functions make sense
      for(int site=0; site < vol(); ++site) 
      {
	multi1d<int> coord1 = crtesn(site, lattSize());

	int linear = linearSiteIndex(coord1);
	int node   = nodeNumber(coord1);

	multi1d<int> coord2 = siteCoords(node,linear);
	int j = local_site(coord2, lattSize());
      
#if QDP_DEBUG >= 2
	QDP_info("site= %d   coord= %d %d %d %d   linear= %d node=%d   crd=%d %d %d %d   j= %d",
		 site,coord1[0],coord1[1],coord1[2],coord1[3],
		 linear,node,
		 coord2[0],coord2[1],coord2[2],coord2[3],
		 j);
#endif
	if (site != j)
	  QDP_error_exit("Layout::create - Layout problems, the layout functions do not work correctly with this lattice size");
      }

      // Initialize various defaults
      initDefaults();

      QDPIO::cout << "Finished lattice layout" << std::endl;
    }
  }


//-----------------------------------------------------------------------------
#if QDP_USE_LEXICO_LAYOUT == 1

#warning "Using a lexicographic layout"

  namespace Layout
  {
    //! The linearized site index for the corresponding coordinate
    /*! This layout is a simple lexicographic lattice ordering */
    int linearSiteIndex(const multi1d<int>& coord)
    {
      multi1d<int> tmp_coord(Nd);

      for(int i=0; i < coord.size(); ++i)
	tmp_coord[i] = coord[i] % Layout::subgridLattSize()[i];
    
      return local_site(tmp_coord, Layout::subgridLattSize());
    }


    //! The node number for the corresponding lattice coordinate
    /*! This layout is a simple lexicographic lattice ordering */
    int nodeNumber(const multi1d<int>& coord)
    {
      multi1d<int> tmp_coord(Nd);

      for(int i=0; i < coord.size(); ++i)
	tmp_coord[i] = coord[i] / Layout::subgridLattSize()[i];
    
      return Layout::getNodeNumberFrom(tmp_coord);
    }


    //! Returns the lattice site for some input node and linear index
    /*! This layout is a simple lexicographic lattice ordering */
    multi1d<int> siteCoords(int node, int linear)
    {
      multi1d<int> coord = getLogicalCoordFrom(node);

      // Get the base (origins) of the absolute lattice coord
      coord *= Layout::subgridLattSize();
    
      // Find the coordinate within a node and accumulate
      // This is a lexicographic ordering
      coord += crtesn(linear, Layout::subgridLattSize());

      return coord;
    }


    //! Return the smallest lattice size per node allowed
    /*! This layout is a simple lexicographic lattice ordering */
    multi1d<int> minimalLayoutMapping()
    {
      multi1d<int> dim(Nd);
      dim = 1;

      return dim;
    }
  }

//-----------------------------------------------------------------------------

#elif QDP_USE_CB2_LAYOUT == 1

#warning "Using a 2 checkerboard (red/black) layout"

  namespace Layout
  {
    //! The linearized site index for the corresponding coordinate
    /*! This layout is appropriate for a 2 checkerboard (red/black) lattice */
    int linearSiteIndex(const multi1d<int>& coord)
    {
      int subgrid_vol_cb = Layout::sitesOnNode() >> 1;
      multi1d<int> subgrid_cb_nrow = Layout::subgridLattSize();
      subgrid_cb_nrow[0] >>= 1;

      int cb = 0;
      for(int m=0; m < Nd; ++m)
	cb += coord[m];
      cb &= 1;

      multi1d<int> subgrid_cb_coord(Nd);
      subgrid_cb_coord[0] = (coord[0] >> 1) % subgrid_cb_nrow[0];
      for(int i=1; i < Nd; ++i)
	subgrid_cb_coord[i] = coord[i] % subgrid_cb_nrow[i];
    
      return local_site(subgrid_cb_coord, subgrid_cb_nrow) + cb*subgrid_vol_cb;
    }


    //! The node number for the corresponding lattice coordinate
    /*! 
     * This layout is appropriate for a 2 checkerboard (red/black) lattice,
     * but to find the nodeNumber this function resembles a simple lexicographic 
     * layout
     */
    int nodeNumber(const multi1d<int>& coord)
    {
      multi1d<int> tmp_coord(Nd);

      for(int i=0; i < coord.size(); ++i)
	tmp_coord[i] = coord[i] / Layout::subgridLattSize()[i];
    
      return Layout::getNodeNumberFrom(tmp_coord);
    }


    //! Reconstruct the lattice coordinate from the node and site number
    /*! 
     * This is the inverse of the nodeNumber and linearSiteIndex functions.
     * The API requires this function to be here.
     */
    multi1d<int> siteCoords(int node, int linearsite) // ignore node
    {
      int subgrid_vol_cb = Layout::sitesOnNode() >> 1;
      multi1d<int> subgrid_cb_nrow = Layout::subgridLattSize();
      subgrid_cb_nrow[0] >>= 1;

      // Get the base (origins) of the absolute lattice coord
      multi1d<int> coord = getLogicalCoordFrom(node);
      coord *= Layout::subgridLattSize();
    
      int cb = linearsite / subgrid_vol_cb;
      multi1d<int> tmp_coord = crtesn(linearsite % subgrid_vol_cb, subgrid_cb_nrow);

      // Add on position within the node
      // NOTE: the cb for the x-coord is not yet determined
      coord[0] += 2*tmp_coord[0];
      for(int m=1; m < Nd; ++m)
	coord[m] += tmp_coord[m];

      // Determine cb including global node cb
      int cbb = cb;
      for(int m=1; m < Nd; ++m)
	cbb += coord[m];
      coord[0] += (cbb & 1);

      return coord;
    }


    //! Return the smallest lattice size per node allowed
    /*! This layout is appropriate for a 2 checkerboard (red/black) lattice */
    multi1d<int> minimalLayoutMapping()
    {
      multi1d<int> dim(Nd);
      dim = 1;
      dim[0] = 2;       // must have multiple length 2 for cb

      return dim;
    }
  }


//-----------------------------------------------------------------------------

#elif QDP_USE_CB32_LAYOUT == 1

#warning "Using a 32 checkerboard layout"

#error "THIS BIT STILL UNDER CONSTRUCTION"

  namespace Layout
  {
    //! The linearized site index for the corresponding coordinate
    /*! This layout is appropriate for a 32-style checkerboard lattice */
    int linearSiteIndex(const multi1d<int>& coord)
    {
      int subgrid_vol_cb = Layout::sitesOnNode() >> (Nd+1);
      multi1d<int> subgrid_cb_nrow = Layout::subgridLattSize();
      subgrid_cb_nrow[0] >>= 2;
      for(int i=1; i < Nd; ++i) 
	subgrid_cb_nrow[i] >>= 1;

      int subl = coord[Nd-1] & 1;
      for(int m=Nd-2; m >= 0; --m)
	subl = (subl << 1) + (coord[m] & 1);

      int cb = 0;
      for(int m=0; m < Nd; ++m)
	cb += coord[m] >> 1;

      subl += (cb & 1) << Nd;   // Final color or checkerboard

      // Construct the checkerboard lattice coord
      multi1d<int> subgrid_cb_coord(Nd);

      subgrid_cb_coord[0] = (coord[0] >> 2) % subgrid_cb_nrow[0];
      for(int i=1; i < Nd; ++i)
	subgrid_cb_coord[i] = (coord[i] >> 1) % subgrid_cb_nrow[i];
    
      return local_site(subgrid_cb_coord, subgrid_cb_nrow) + subl*subgrid_vol_cb;
    }


    //! The node number for the corresponding lattice coordinate
    /*! 
     * This layout is appropriate for a 32 checkerboard (red/black) lattice,
     * but to find the nodeNumber this function resembles a simple lexicographic 
     * layout
     */
    int nodeNumber(const multi1d<int>& coord)
    {
      multi1d<int> tmp_coord(Nd);

      for(int i=0; i < coord.size(); ++i)
	tmp_coord[i] = coord[i] / Layout::subgridLattSize()[i];
    
      return Layout::getNodeNumberFrom(tmp_coord);
    }


    //! Reconstruct the lattice coordinate from the node and site number
    /*! 
     * This is the inverse of the nodeNumber and linearSiteIndex functions.
     * The API requires this function to be here.
     */
    multi1d<int> siteCoords(int node, int linearsite) // ignore node
    {
      int subgrid_vol_cb = Layout::sitesOnNode() >> (Nd+1);
      multi1d<int> subgrid_cb_nrow = Layout::subgridLattSize();
      subgrid_cb_nrow[0] >>= 2;
      for(int i=1; i < Nd; ++i) 
	subgrid_cb_nrow[i] >>= 1;

      // Get the base (origins) of the absolute lattice coord
      multi1d<int> coord = Layout::getLogicalCoordFrom(node);
      coord *= Layout::subgridLattSize();
    
      int subl = linearsite / subgrid_vol_cb;
      multi1d<int> tmp_coord = crtesn(linearsite % subgrid_vol_cb, subgrid_cb_nrow);

      // Add on position within the node
      // NOTE: the cb for the x-coord is not yet determined
      coord[0] += tmp_coord[0] << 2;
      for(int m=1; m < Nd; ++m)
	coord[m] += tmp_coord[m] << 1;

      int cb = 0;
      for(int m=1; m<Nd; ++m)
	cb += coord[m];
      cb &= 1;

      subl ^= (cb << Nd);
      for(int m=0; m<Nd; ++m)
	coord[m] ^= (subl & (1 << m)) >> m;
      coord[0] ^= (subl & (1 << Nd)) >> (Nd-1);   // this gets the hypercube cb

      return coord;
    }


    //! Return the smallest lattice size per node allowed
    /*! This layout is appropriate for a 32-style checkerboard lattice */
    multi1d<int> minimalLayoutMapping()
    {
      multi1d<int> dim(Nd);

      dim[0] = 4;       // must have multiple length 2 for cb and hypercube
      for(int m=1; m < Nd; ++m)
	dim[m] = 2;     //  must have multiple length 2 for hypercube

      return dim;
    }
  }

#else

#error "no appropriate layout defined"

#endif

//-----------------------------------------------------------------------------


} // namespace QDP;
