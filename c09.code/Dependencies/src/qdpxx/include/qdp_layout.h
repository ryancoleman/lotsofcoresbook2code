// -*- C++ -*-

/*! @file
* @brief Lattice layout
*
* Lattice layout namespace and operations
*/

#ifndef QDP_LAYOUT_H
#define QDP_LAYOUT_H

namespace QDP {

	/*! @defgroup layout  Layout 
	*
	* Namespace holding info on problem size and machine info
	*
	* @{
	*/

	//! Layout namespace holding info on problem size and machine info
	/*! 
	* This is a namespace instead of a class since it is completely static -
	* no such object should be created 
	*
	* The functions here should be common to all architectures
	*/
	namespace Layout
	{
		//! Initialize some fundamental pieces of the layout
		/*! This routine is used to boostrap the   create   function below */
		void init();

		//! Main lattice creation routine
		void create();

		//! Main destruction routine
		void destroy();

		//! Set lattice size -- problem size
		void setLattSize(const multi1d<int>& nrows);

		//! Set SMP flag -- true if using smp/multiprocessor mode on a node
		void setSMPFlag(bool);

		//! Set number of processors in a multi-threaded implementation
		void setNumProc(int N);

		//! Returns the logical node number for the corresponding lattice coordinate
		/*! The API requires this function to be here */
		int nodeNumber(const multi1d<int>& coord) QDP_CONST;

		//returns local lexicographical site coordinate from linear index:
		multi1d<int> localLexiCoordFromLinear(const int& linearr) QDP_CONST;

		//! The linearized site index within a node for the corresponding lattice coordinate
		/*! The API requires this function to be here */
		int linearSiteIndex(const multi1d<int>& coord) QDP_CONST;

		//! Reconstruct the lattice coordinate from the node and site number
		/*! 
		* This is the inverse of the nodeNumber and linearSiteIndex functions.
		* The API requires this function to be here.
		*/
		multi1d<int> siteCoords(int node, int index) QDP_CONST;
  
		extern "C" { 
			/* Export this to "C" */
			void QDPXX_getSiteCoords(int coord[], int node, int linear) QDP_CONST;
			int QDPXX_getLinearSiteIndex(const int coord[]);
			int QDPXX_nodeNumber(const int coord[]);

		};

		//! Returns the node number of this node
		int nodeNumber() QDP_CONST;

		//! Returns the number of nodes
		int numNodes() QDP_CONST;
	

		//! Virtual grid (problem grid) lattice size
		const multi1d<int>& lattSize() QDP_CONST;

		//! Total lattice volume
		int vol() QDP_CONST;

		//! Number of sites on node
		int sitesOnNode() QDP_CONST;

		//! Returns whether this is the primary node
		bool primaryNode() QDP_CONST;

		//! The linearized site index for the corresponding lexicographic site
		int linearSiteIndex(int lexicosite) QDP_CONST;

		//! Returns the logical node coordinates for this node
		const multi1d<int>& nodeCoord() QDP_CONST;

		//! Returns the logical node coordinates for the corresponding lattice coordinate
		multi1d<int> nodeCoord(const multi1d<int>& coord);

		//! Subgrid (grid on each node) lattice size
		const multi1d<int>& subgridLattSize() QDP_CONST;

		//! Returns the logical size of this machine
		const multi1d<int>& logicalSize() QDP_CONST;

		//! Returns the node number given some logical node coordinate
		/*! This is not meant to be speedy */
		int getNodeNumberFrom(const multi1d<int>& node_coord);

		//! Returns the logical node coordinates given some node number
		/*! This is not meant to be speedy */
		multi1d<int> getLogicalCoordFrom(int node);


		//! Check if I/O grid is defined
		bool isIOGridDefined(void) QDP_CONST;

		//! set the IO node grid
		void setIONodeGrid(const multi1d<int>& io_grid);

		//! number of I/O nodes
		int numIONodeGrid(void) QDP_CONST;
	
		//! Get the I/O Node grid
		const multi1d<int>& getIONodeGrid() QDP_CONST;


	}

	//! Declaration of shift function object
	extern ArrayBiDirectionalMap  shift;


	/*! @} */   // end of group layout

} // namespace QDP

#endif
