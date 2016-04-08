/*! @file
 * @brief Parscalar specific routines
 * 
 * Routines for parscalar implementation
 */


#include "qdp.h"
#include "qdp_util.h"
#include "qmp.h"


namespace QDP {

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
//! Initializer for generic map constructor
  void Map::make(const MapFunc& func)
  {
#if QDP_DEBUG >= 3
    QDP_info("Map::make");
#endif
    const int nodeSites = Layout::sitesOnNode();

    //--------------------------------------
    // Setup the communication index arrays
    goffsets.resize(nodeSites);
    srcnode.resize(nodeSites);
    dstnode.resize(nodeSites);

    const int my_node = Layout::nodeNumber();

    // Loop over the sites on this node
    /* This for loop added by Jacques. Should be OK.
       Loop works out forward and backward neighbours and writes
       one value per value of linear. So no write conflicts */

#pragma omp parallel for
    for(int linear=0; linear < nodeSites; ++linear)
    {
      // Get the true lattice coord of this linear site index
      multi1d<int> coord = Layout::siteCoords(my_node, linear);

      // Source neighbor for this destination site
      multi1d<int> fcoord = func(coord,+1);

      // Destination neighbor receiving data from this site
      // This functions as the inverse map
      multi1d<int> bcoord = func(coord,-1);

      int fnode = Layout::nodeNumber(fcoord);
      int bnode = Layout::nodeNumber(bcoord);

      // Source linear site and node
      goffsets[linear] = Layout::linearSiteIndex(fcoord);
      srcnode[linear]  = fnode;

      // Destination node
      dstnode[linear]  = bnode;

#if QDP_DEBUG >= 3
      QDP_info("linear=%d  coord=%d %d %d %d   fcoord=%d %d %d %d   bcoord=%d %d %d %d   goffsets=%d", 
	       linear, 
	       coord[0], coord[1], coord[2], coord[3],
	       fcoord[0], fcoord[1], fcoord[2], fcoord[3],
	       bcoord[0], bcoord[1], bcoord[2], bcoord[3],
	       goffsets[linear]);
#endif
    }	// end for-loop

    // Return a list of the unique nodes in the list
    // NOTE: my_node may be included as a unique node, so one extra
    multi1d<int> srcenodes_tmp = uniquify_list(srcnode);
    multi1d<int> destnodes_tmp = uniquify_list(dstnode);

    // To simplify life, remove a possible occurance of my_node
    // This may mean that the list could be empty afterwards, so that
    // means mark offnodeP as false
    int cnt_srcenodes = 0;
    for(int i=0; i < srcenodes_tmp.size(); ++i)
      if (srcenodes_tmp[i] != my_node)
				++cnt_srcenodes;

    int cnt_destnodes = 0;
    for(int i=0; i < destnodes_tmp.size(); ++i)
      if (destnodes_tmp[i] != my_node)
				++cnt_destnodes;

#if QDP_DEBUG >= 3
    // Debugging
    for(int i=0; i < srcenodes_tmp.size(); ++i)
      QDP_info("srcenodes_tmp(%d) = %d",i,srcenodes_tmp[i]);

    for(int i=0; i < destnodes_tmp.size(); ++i)
      QDP_info("destnodes_tmp(%d) = %d",i,destnodes_tmp[i]);

    QDP_info("cnt_srcenodes = %d", cnt_srcenodes);
    QDP_info("cnt_destnodes = %d", cnt_destnodes);
#endif

    // A sanity check - both counts must be either 0 or non-zero
    if (cnt_srcenodes > 0 && cnt_destnodes == 0)
      QDP_error_exit("Map: some bizarre error - no dest nodes but have srce nodes");

    if (cnt_srcenodes == 0 && cnt_destnodes > 0)
      QDP_error_exit("Map: some bizarre error - no srce nodes but have dest nodes");

    // If no srce/dest nodes, then we know no off-node communications
    offnodeP = (cnt_srcenodes > 0) ? true : false;
  
    //
    // The rest of the routine is devoted to supporting off-node communications
    // If there is not any communications, then return
    //
    if (! offnodeP)
    {
#if QDP_DEBUG >= 3
      QDP_info("no off-node communications: exiting Map::make");
#endif
      return;
    }

    // Finally setup the srce and dest nodes now without my_node
    srcenodes.resize(cnt_srcenodes);
    destnodes.resize(cnt_destnodes);
  
    for(int i=0, j=0; i < srcenodes_tmp.size(); ++i)
      if (srcenodes_tmp[i] != my_node)
				srcenodes[j++] = srcenodes_tmp[i];

    for(int i=0, j=0; i < destnodes_tmp.size(); ++i)
      if (destnodes_tmp[i] != my_node)
				destnodes[j++] = destnodes_tmp[i];

#if QDP_DEBUG >= 3
    // Debugging
    for(int i=0; i < srcenodes.size(); ++i)
      QDP_info("srcenodes(%d) = %d",i,srcenodes(i));

    for(int i=0; i < destnodes.size(); ++i)
      QDP_info("destnodes(%d) = %d",i,destnodes(i));
#endif

    // Run through the lists and find the number of each unique node
    srcenodes_num.resize(srcenodes.size());
    destnodes_num.resize(destnodes.size());

    srcenodes_num = 0;
    destnodes_num = 0;

		// fast so no threading necessary
    for(int linear=0; linear < nodeSites; ++linear)
    {
			int this_node = srcnode[linear];
      if (this_node != my_node)
				for(int i=0; i < srcenodes_num.size(); ++i)
				{
					if (srcenodes[i] == this_node)
					{
						srcenodes_num[i]++;
						break;
					}
				}

      int that_node = dstnode[linear];
      if (that_node != my_node)
			for(int i=0; i < destnodes_num.size(); ++i)
			{
				if (destnodes[i] == that_node)
				{
					destnodes_num[i]++;
					break;
				}
			}
    }	// end for linear
  

#if QDP_DEBUG >= 3
    for(int i=0; i < destnodes.size(); ++i)
    {
      QDP_info("srcenodes(%d) = %d",i,srcenodes(i));
      QDP_info("destnodes(%d) = %d",i,destnodes(i));
    }

    for(int i=0; i < destnodes_num.size(); ++i)
    {
      QDP_info("srcenodes_num(%d) = %d",i,srcenodes_num(i));
      QDP_info("destnodes_num(%d) = %d",i,destnodes_num(i));
    }
#endif

    // Implementation limitation in the Map::operator(). Only support
    // a node sending data all to one node or no sending at all (offNodeP == false).
    if (srcenodes.size() != 1)
      QDP_error_exit("Map: for now only allow 1 destination node");
      
    if (destnodes.size() != 1)
      QDP_error_exit("Map: for now only allow receives from 1 node");


    // Now make a small scatter array for the dest_buf so that when data
    // is sent, it is put in an order the gather can pick it up
    // If we allow multiple dest nodes, then soffsets here needs to be
    // an array of arrays
    soffsets.resize(destnodes_num[0]);
  
    // Loop through sites on my *destination* node - here I assume all nodes have
    // the same number of sites. Mimic the gather pattern needed on that node and
    // set my scatter array to scatter into the correct site order
    for(int i=0, si=0; i < nodeSites; ++i) 
    {
      // Get the true lattice coord of this linear site index
      multi1d<int> coord = Layout::siteCoords(destnodes[0], i);
      multi1d<int> fcoord = func(coord,+1);
      int fnode = Layout::nodeNumber(fcoord);
      int fline = Layout::linearSiteIndex(fcoord);

      if (fnode == my_node)
				soffsets[si++] = fline;
    }

#if QDP_DEBUG >= 3
    // Debugging
    for(int i=0; i < soffsets.size(); ++i)
      QDP_info("soffsets(%d) = %d",i,soffsets(i));
#endif


#if QDP_DEBUG >= 3
    QDP_info("exiting Map::make");
#endif
  }


//------------------------------------------------------------------------
// Message passing convenience routines
//------------------------------------------------------------------------

  namespace QDPInternal
  {
    //! Broadcast a string from primary node to all other nodes
    void broadcast_str(std::string& result)
    {
      char *dd_tmp;
      int lleng;

      // Only primary node can grab string
      if (Layout::primaryNode()) 
      {
	lleng = result.length();
      }

      // First must broadcast size of string
      QDPInternal::broadcast(lleng);

      // Now every node can alloc space for string
      dd_tmp = new(std::nothrow) char[lleng];
      if( dd_tmp == 0x0 ) { 
	QDP_error_exit("Unable to allocate dd_tmp\n");
      }

      if (Layout::primaryNode())
	memcpy(dd_tmp, result.c_str(), lleng);
  
      // Now broadcast char array out to all nodes
      QDPInternal::broadcast((void *)dd_tmp, lleng);

      // All nodes can now grab char array and make a string, but only
      // need this on non-primary nodes
      if (! Layout::primaryNode())
      {
	result.assign(dd_tmp, lleng);
      }

      // Clean-up and boogie
      delete[] dd_tmp;
    }

    //! Is this a grid architecture
    bool gridArch()
    { 
      return (QMP_get_msg_passing_type() == QMP_GRID) ? true : false;
    }

    //! Send a clear-to-send
    void clearToSend(void *buffer, int count, int node)
    { 
      // On non-grid machines, use a clear-to-send like protocol
      if (! QDPInternal::gridArch())
      {
	if(Layout::nodeNumber() == 0 && node != 0)
	  QDPInternal::sendToWait(buffer, node, count);
	if(Layout::nodeNumber() == node && node != 0)
	  QDPInternal::recvFromWait(buffer, 0, count);
      }
    }

    //! Route to another node (blocking)
    void route(void *buffer, int srce_node, int dest_node, int count)
    { 
#if QDP_DEBUG >= 2
      QDP_info("starting a route, count=%d, srcenode=%d destnode=%d", 
	       count,srce_node,dest_node);
#endif

//    QMP_route(buffer, count, srce_node, dest_node);
      DML_route_bytes((char*)buffer, count, srce_node, dest_node);

#if QDP_DEBUG >= 2
      QDP_info("finished a route");
#endif
    }

    //! Send to another node (wait)
    void 
    sendToWait(void *send_buf, int dest_node, int count)
    {
#if QDP_DEBUG >= 2
      QDP_info("starting a sendToWait, count=%d, destnode=%d", count,dest_node);
#endif

      QMP_msgmem_t request_msg = QMP_declare_msgmem(send_buf, count);
      QMP_msghandle_t request_mh = QMP_declare_send_to(request_msg, dest_node, 0);

      if (QMP_start(request_mh) != QMP_SUCCESS)
	QDP_error_exit("sendToWait failed\n");

      QMP_wait(request_mh);

      QMP_free_msghandle(request_mh);
      QMP_free_msgmem(request_msg);

#if QDP_DEBUG >= 2
      QDP_info("finished a sendToWait");
#endif
    }

    //! Receive from another node (wait)
    void 
    recvFromWait(void *recv_buf, int srce_node, int count)
    {
#if QDP_DEBUG >= 2
      QDP_info("starting a recvFromWait, count=%d, srcenode=%d", count, srce_node);
#endif

      QMP_msgmem_t request_msg = QMP_declare_msgmem(recv_buf, count);
      QMP_msghandle_t request_mh = QMP_declare_receive_from(request_msg, srce_node, 0);

      if (QMP_start(request_mh) != QMP_SUCCESS)
	QDP_error_exit("recvFromWait failed\n");

      QMP_wait(request_mh);

      QMP_free_msghandle(request_mh);
      QMP_free_msgmem(request_msg);

#if QDP_DEBUG >= 2
      QDP_info("finished a recvFromWait");
#endif
    }

  };


//-----------------------------------------------------------------------------
// Write a lattice quantity
  void writeOLattice(BinaryWriter& bin, 
		     const char* output, size_t size, size_t nmemb)
  {
    const int xinc = Layout::subgridLattSize()[0];

    size_t sizemem = size*nmemb;
    size_t tot_size = sizemem*xinc;
    char *recv_buf = new(std::nothrow) char[tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf\n");
    }

    // Find the location of each site and send to primary node
    int old_node = 0;

    for(int site=0; site < Layout::vol(); site += xinc)
    {
      // first site in each segment uniquely identifies the node
      int node = Layout::nodeNumber(crtesn(site, Layout::lattSize()));

      // Send nodes must wait for a ready signal from the master node
      // to prevent message pileups on the master node
      if (node != old_node)
      {
	// On non-grid machines, use a clear-to-send like protocol
	QDPInternal::clearToSend(recv_buf,sizeof(int),node);
	old_node = node;
      }
    
      // Copy to buffer: be really careful since max(linear) could vary among nodes
      if (Layout::nodeNumber() == node)
      {
	for(int i=0; i < xinc; ++i)
	{
	  int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));
	  memcpy(recv_buf+i*sizemem, output+linear*sizemem, sizemem);
	}
      }

      // Send result to primary node. Avoid sending prim-node sending to itself
      if (node != 0)
      {
#if 1
	// All nodes participate
	QDPInternal::route((void *)recv_buf, node, 0, tot_size);
#else
	if (Layout::primaryNode())
	  QDPInternal::recvFromWait((void *)recv_buf, node, tot_size);

	if (Layout::nodeNumber() == node)
	  QDPInternal::sendToWait((void *)recv_buf, 0, tot_size);
#endif
      }

      bin.writeArrayPrimaryNode(recv_buf, size, nmemb*xinc);
    }

    delete[] recv_buf;
  }


// Write a single site of a lattice quantity
  void writeOLattice(BinaryWriter& bin, 
		     const char* output, size_t size, size_t nmemb,
		     const multi1d<int>& coord)
  {
    size_t tot_size = size*nmemb;
    char *recv_buf = new(std::nothrow) char[tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recvbuf\n");
    }


    // Send site to primary node
    int node   = Layout::nodeNumber(coord);
    int linear = Layout::linearSiteIndex(coord);

    // Send nodes must wait for a ready signal from the master node
    // to prevent message pileups on the master node
    QDPInternal::clearToSend(recv_buf,sizeof(int),node);
  
    // Copy to buffer: be really careful since max(linear) could vary among nodes
    if (Layout::nodeNumber() == node)
      memcpy(recv_buf, output+linear*tot_size, tot_size);
  
    // Send result to primary node. Avoid sending prim-node sending to itself
    if (node != 0)
    {
#if 1
      // All nodes participate
      QDPInternal::route((void *)recv_buf, node, 0, tot_size);
#else
      if (Layout::primaryNode())
	QDPInternal::recvFromWait((void *)recv_buf, node, tot_size);

      if (Layout::nodeNumber() == node)
	QDPInternal::sendToWait((void *)recv_buf, 0, tot_size);
#endif
    }

    bin.writeArray(recv_buf, size, nmemb);

    delete[] recv_buf;
  }


//-----------------------------------------------------------------------------
// Write a lattice quantity
  void writeOLattice(BinaryWriter& bin, 
		     const char* output, size_t size, size_t nmemb,
		     const Subset& sub)
  {
    // Single node code
    const Set& set    = sub.getSet();
    const multi1d<int>& lat_color = set.latticeColoring();
    const int color = sub.color();

    const int xinc = Layout::subgridLattSize()[0];

    size_t sizemem = size*nmemb;
    size_t max_tot_size = sizemem*xinc;
    char *recv_buf = new(std::nothrow) char[max_tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf\n");
    }

    char *recv_buf_size = new(std::nothrow) char[sizeof(int)];
    if( recv_buf_size == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf_size\n");
    }

    // Find the location of each site and send to primary node
    int old_node = 0;

    for(int site=0; site < Layout::vol(); site += xinc)
    {
      // This algorithm is cumbersome. We do not keep the coordinate function for the subset.
      // So, we have to ask the sending node how many sites are to be transferred,
      // and rely on each node to send whatever number of sites live in the desired
      // subgridLattSize strip.

      // first site in each segment uniquely identifies the node
      int node = Layout::nodeNumber(crtesn(site, Layout::lattSize()));

      // Send nodes must wait for a ready signal from the master node
      // to prevent message pileups on the master node
      if (node != old_node)
      {
	// On non-grid machines, use a clear-to-send like protocol
	QDPInternal::clearToSend(recv_buf,sizeof(int),node);
	old_node = node;
      }
    
      // Copy to buffer: be really careful since max(linear) could vary among nodes
      int site_cnt = 0;
      if (Layout::nodeNumber() == node)
      {
	for(int i=0; i < xinc; ++i)
	{
	  int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));
	  if (lat_color[linear] == color)
	  {
	    memcpy(recv_buf+site_cnt*sizemem, output+linear*sizemem, sizemem);
	    site_cnt++;
	  }
	}
	memcpy(recv_buf_size, (void *)&site_cnt, sizeof(int));
      }

      // Send result to primary node. Avoid sending prim-node sending to itself
      if (node != 0)
      {
#if 0
	// All nodes participate
	// First send the byte size for this 
	QDP_error_exit("Do not support route in writeOLattice(sub)");

#else
	// We are using the point-to-point version
	if (Layout::primaryNode())
	{
	  QDPInternal::recvFromWait((void *)recv_buf_size, node, sizeof(int));
	  memcpy((void *)&site_cnt, recv_buf_size, sizeof(int));
	  QDPInternal::recvFromWait((void *)recv_buf, node, site_cnt*sizemem);
	}

	if (Layout::nodeNumber() == node)
	{
	  QDPInternal::sendToWait((void *)recv_buf_size, 0, sizeof(int));
	  QDPInternal::sendToWait((void *)recv_buf, 0, site_cnt*sizemem);
	}
#endif
      }

      bin.writeArrayPrimaryNode(recv_buf, size, nmemb*site_cnt);
    }

    delete[] recv_buf_size;
    delete[] recv_buf;
  }



//! Read a lattice quantity
/*! This code assumes no inner grid */
  void readOLattice(BinaryReader& bin, 
		    char* input, size_t size, size_t nmemb)
  {
    const int xinc = Layout::subgridLattSize()[0];

    size_t sizemem = size*nmemb;
    size_t tot_size = sizemem*xinc;
    char *recv_buf = new(std::nothrow) char[tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recvbuf\n");
    }

    // Find the location of each site and send to primary node
    for(int site=0; site < Layout::vol(); site += xinc)
    {
      // first site in each segment uniquely identifies the node
      int node = Layout::nodeNumber(crtesn(site, Layout::lattSize()));

      // Only on primary node read the data
      bin.readArrayPrimaryNode(recv_buf, size, nmemb*xinc);

      // Send result to destination node. Avoid sending prim-node sending to itself
      if (node != 0)
      {
#if 1
	// All nodes participate
	QDPInternal::route((void *)recv_buf, 0, node, tot_size);
#else
	if (Layout::primaryNode())
	  QDPInternal::sendToWait((void *)recv_buf, node, tot_size);

	if (Layout::nodeNumber() == node)
	  QDPInternal::recvFromWait((void *)recv_buf, 0, tot_size);
#endif
      }

      if (Layout::nodeNumber() == node)
      {
	for(int i=0; i < xinc; ++i)
	{
	  int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));

	  memcpy(input+linear*sizemem, recv_buf+i*sizemem, sizemem);
	}
      }
    }

    delete[] recv_buf;
  }

//! Read a single site worth of a lattice quantity
/*! This code assumes no inner grid */
  void readOLattice(BinaryReader& bin, 
		    char* input, size_t size, size_t nmemb,
		    const multi1d<int>& coord)
  {
    size_t tot_size = size*nmemb;
    char *recv_buf = new(std::nothrow) char[tot_size];
    if( recv_buf == 0x0 ) {
      QDP_error_exit("Unable to allocate recv_buf\n");
    }


    // Find the location of each site and send to primary node
    int node   = Layout::nodeNumber(coord);
    int linear = Layout::linearSiteIndex(coord);

    // Only on primary node read the data
    bin.readArrayPrimaryNode(recv_buf, size, nmemb);

    // Send result to destination node. Avoid sending prim-node sending to itself
    if (node != 0)
    {
#if 1
      // All nodes participate
      QDPInternal::route((void *)recv_buf, 0, node, tot_size);
#else
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)recv_buf, node, tot_size);

      if (Layout::nodeNumber() == node)
	QDPInternal::recvFromWait((void *)recv_buf, 0, tot_size);
#endif
    }

    if (Layout::nodeNumber() == node)
      memcpy(input+linear*tot_size, recv_buf, tot_size);

    delete[] recv_buf;
  }

//! Read a lattice quantity
/*! This code assumes no inner grid */
  void readOLattice(BinaryReader& bin, 
		    char* input, size_t size, size_t nmemb,
		    const Subset& sub)
  {
    // Single node code
    const Set& set    = sub.getSet();
    const multi1d<int>& lat_color = set.latticeColoring();
    const int color = sub.color();

    const int xinc = Layout::subgridLattSize()[0];

    size_t sizemem = size*nmemb;
    size_t max_tot_size = sizemem*xinc;
    char *recv_buf = new(std::nothrow) char[max_tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf\n");
    }

    char *recv_buf_size = new(std::nothrow) char[sizeof(int)];
    if( recv_buf_size == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf_size\n");
    }

    // Find the location of each site and send to primary node
    for(int site=0; site < Layout::vol(); site += xinc)
    {
      // This algorithm is cumbersome. We do not keep the coordinate function for the subset.
      // So, we have to ask the sending node how many sites are to be transferred,
      // and rely on each node to send whatever number of sites live in the desired
      // subgridLattSize strip.

      // first site in each segment uniquely identifies the node
      int node = Layout::nodeNumber(crtesn(site, Layout::lattSize()));

      // Find the amount of data to read. Unfortunately, have to ask the remote node
      // Place the result in a send buffer
      int site_cnt = 0;
      if (Layout::nodeNumber() == node)
      {
	for(int i=0; i < xinc; ++i)
	{
	  int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));
	  if (lat_color[linear] == color)
	  {
	    site_cnt++;
	  }
	}
	memcpy(recv_buf_size, (void *)&site_cnt, sizeof(int));
      }

      if (node != 0)
      {
	// Send the data size to the primary node.
	// We are using the point-to-point version
	if (Layout::primaryNode())
	{
	  QDPInternal::recvFromWait((void *)recv_buf_size, node, sizeof(int));
	  memcpy((void *)&site_cnt, recv_buf_size, sizeof(int));
	}

	if (Layout::nodeNumber() == node)
	{
	  QDPInternal::sendToWait((void *)recv_buf_size, 0, sizeof(int));
	}
      }

      // Only on primary node read the data
      bin.readArrayPrimaryNode(recv_buf, size, nmemb*site_cnt);

      // Send result to destination node. Avoid sending prim-node sending to itself
      if (node != 0)
      {
#if 0
	// All nodes participate
	// First send the byte size for this 
	QDP_error_exit("Do not support route in readOLattice(sub)");

#else
	// We are using the point-to-point version
	if (Layout::primaryNode())
	  QDPInternal::sendToWait((void *)recv_buf, node, site_cnt*sizemem);

	if (Layout::nodeNumber() == node)
	  QDPInternal::recvFromWait((void *)recv_buf, 0, site_cnt*sizemem);
#endif
      }

      if (Layout::nodeNumber() == node)
      {
	for(int i=0,j=0; i < xinc; ++i)
	{
	  int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));
	  if (lat_color[linear] == color)
	  {
	    memcpy(input+linear*sizemem, recv_buf+j*sizemem, sizemem);
	    j++;
	  }
	}
      }
    }

    delete[] recv_buf_size;
    delete[] recv_buf;
  }


  // **************************************************************
  namespace LatticeTimeSliceIO 
  {
    void readOLatticeSlice(BinaryReader& bin, char* input, 
			   size_t size, size_t nmemb,
			   int start_lexico, int stop_lexico)
    {
      const int xinc = Layout::subgridLattSize()[0];

      if ((stop_lexico % xinc) != 0)
      {
	QDPIO::cerr << __func__ << ": erorr: stop_lexico= " << stop_lexico << "  xinc= " << xinc << std::endl;
	QDP_abort(1);
      }

      size_t sizemem = size*nmemb;
      size_t tot_size = sizemem*xinc;
      char *recv_buf = new(std::nothrow) char[tot_size];
      if( recv_buf == 0x0 ) { 
	QDP_error_exit("Unable to allocate recvbuf\n");}

      // Find the location of each site and send to primary node
      for (int site=start_lexico; site < stop_lexico; site += xinc)
      {
	// first site in each segment uniquely identifies the node
	int node = Layout::nodeNumber(crtesn(site, Layout::lattSize()));

	// Only on primary node read the data
	bin.readArrayPrimaryNode(recv_buf, size, nmemb*xinc);

	// Send result to destination node. Avoid sending prim-node sending to itself
	if (node != 0)
	{
#if 1
	  // All nodes participate
	  QDPInternal::route((void *)recv_buf, 0, node, tot_size);
#else
	  if (Layout::primaryNode())
	    QDPInternal::sendToWait((void *)recv_buf, node, tot_size);
	  if (Layout::nodeNumber() == node)
	    QDPInternal::recvFromWait((void *)recv_buf, 0, tot_size);
#endif
	}

	if (Layout::nodeNumber() == node)
	{
	  for(int i=0; i < xinc; ++i)
	  {
	    int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));
	    memcpy(input+linear*sizemem, recv_buf+i*sizemem, sizemem);
	  }
	}
      }

      delete[] recv_buf;
    }

 
    // Write a time slice of a lattice quantity (time must be most slowly varying)
    void writeOLatticeSlice(BinaryWriter& bin, const char* output, 
			    size_t size, size_t nmemb,
			    int start_lexico, int stop_lexico)
    {
      const int xinc = Layout::subgridLattSize()[0];

      if ((stop_lexico % xinc) != 0)
      {
	QDPIO::cerr << __func__ << ": erorr: stop_lexico= " << stop_lexico << "  xinc= " << xinc << std::endl;
	QDP_abort(1);
      }

      size_t sizemem = size*nmemb;
      size_t tot_size = sizemem*xinc;
      char *recv_buf = new(std::nothrow) char[tot_size];
      if( recv_buf == 0x0 ) { 
	QDP_error_exit("Unable to allocate recv_buf\n");}

      // Find the location of each site and send to primary node
      int old_node = 0;

      for (int site=start_lexico; site < stop_lexico; site += xinc)
      {
	// first site in each segment uniquely identifies the node
	int node = Layout::nodeNumber(crtesn(site, Layout::lattSize()));

	// Send nodes must wait for a ready signal from the master node
	// to prevent message pileups on the master node
	if (node != old_node){
	  // On non-grid machines, use a clear-to-send like protocol
	  QDPInternal::clearToSend(recv_buf,sizeof(int),node);
	  old_node = node;}
    
	// Copy to buffer: be really careful since max(linear) could vary among nodes
	if (Layout::nodeNumber() == node){
	  for(int i=0; i < xinc; ++i){
	    int linear = Layout::linearSiteIndex(crtesn(site+i, Layout::lattSize()));
	    memcpy(recv_buf+i*sizemem, output+linear*sizemem, sizemem);
	  }
	}

	// Send result to primary node. Avoid sending prim-node sending to itself
	if (node != 0)
	{
#if 1
	  // All nodes participate
	  QDPInternal::route((void *)recv_buf, node, 0, tot_size);
#else
	  if (Layout::primaryNode())
	    QDPInternal::recvFromWait((void *)recv_buf, node, tot_size);
	  if (Layout::nodeNumber() == node)
	    QDPInternal::sendToWait((void *)recv_buf, 0, tot_size);
#endif
	}

	bin.writeArrayPrimaryNode(recv_buf, size, nmemb*xinc);
      }
      delete[] recv_buf;
    }
  }


//-----------------------------------------------------------------------
// Compute simple NERSC-like checksum of a gauge field
/*
 * \ingroup io
 *
 * \param u          gauge configuration ( Read )
 *
 * \return checksum
 */    

  n_uint32_t computeChecksum(const multi1d<LatticeColorMatrix>& u,
			     int mat_size)
  {
    size_t size = sizeof(REAL32);
    size_t su3_size = size*mat_size;
    n_uint32_t checksum = 0;   // checksum

    multi1d<multi1d<ColorMatrix> > sa(Nd);   // extract gauge fields
    const int nodeSites = Layout::sitesOnNode();

    for(int dd=0; dd<Nd; dd++)        /* dir */
    {
      sa[dd].resize(nodeSites);
      QDP_extract(sa[dd], u[dd], all);
    }

    char  *chk_buf = new(std::nothrow) char[su3_size];
    if( chk_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate chk_buf\n");
    }

    for(int linear=0; linear < nodeSites; ++linear)
    {
      for(int dd=0; dd<Nd; dd++)        /* dir */
      {
	switch (mat_size)
	{
	case 12:
	{
	  REAL32 su3[2][3][2];

	  for(int kk=0; kk<Nc; kk++)      /* color */
	    for(int ii=0; ii<2; ii++)    /* color */
	    {
	      Complex sitecomp = peekColor(sa[dd][linear],ii,kk);
	      su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
	      su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	    }

	  memcpy(chk_buf, &(su3[0][0][0]), su3_size);
	}
	break;

	case 18:
	{
	  REAL32 su3[3][3][2];

	  for(int kk=0; kk<Nc; kk++)      /* color */
	    for(int ii=0; ii<Nc; ii++)    /* color */
	    {
	      Complex sitecomp = peekColor(sa[dd][linear],ii,kk);
	      su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
	      su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	    }

	  memcpy(chk_buf, &(su3[0][0][0]), su3_size);
	}
	break;

	default:
	  QDPIO::cerr << __func__ << ": unexpected size" << std::endl;
	  QDP_abort(1);
	}

	// Compute checksum
	n_uint32_t* chk_ptr = (n_uint32_t*)chk_buf;
	for(unsigned int i=0; i < mat_size*size/sizeof(n_uint32_t); ++i)
	  checksum += chk_ptr[i];
      }
    }

    delete[] chk_buf;

    // Get all nodes to contribute
    QDPInternal::globalSumArray((unsigned int*)&checksum, 1);   // g++ requires me to narrow the type to unsigned int

    return checksum;
  }


//-----------------------------------------------------------------------
// Read a QCD archive file
// Read a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param cfg_in     binary writer object ( Modify )
 * \param u          gauge configuration ( Modify )
 */    

  void readArchiv(BinaryReader& cfg_in, multi1d<LatticeColorMatrix>& u, 
		  n_uint32_t& checksum, int mat_size, int float_size)
  {
    size_t size = float_size;
    size_t su3_size = size*mat_size;
    size_t tot_size = su3_size*Nd;
    const int nodeSites = Layout::sitesOnNode();

    char  *input = new(std::nothrow) char[tot_size*nodeSites];  // keep another copy in input buffers
    if( input == 0x0 ) { 
      QDP_error_exit("Unable to allocate input\n");
    }

    char  *recv_buf = new(std::nothrow) char[tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf\n");
    }

    checksum = 0;

    // Find the location of each site and send to primary node
    for(int site=0; site < Layout::vol(); ++site)
    {
      multi1d<int> coord = crtesn(site, Layout::lattSize());

      int node   = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);

      // Only on primary node read the data
      cfg_in.readArrayPrimaryNode(recv_buf, size, mat_size*Nd);

      if (Layout::primaryNode()) 
      {
	// Compute checksum
	n_uint32_t* chk_ptr = (n_uint32_t*)recv_buf;
	for(unsigned int i=0; i < mat_size*Nd*size/sizeof(n_uint32_t); ++i)
	  checksum += chk_ptr[i];
      }

      // Send result to destination node. Avoid sending prim-node sending to itself
      if (node != 0)
      {
#if 1
	// All nodes participate
	QDPInternal::route((void *)recv_buf, 0, node, tot_size);
#else
	if (Layout::primaryNode())
	  QDPInternal::sendToWait((void *)recv_buf, node, tot_size);

	if (Layout::nodeNumber() == node)
	  QDPInternal::recvFromWait((void *)recv_buf, 0, tot_size);
#endif
      }

      if (Layout::nodeNumber() == node)
	memcpy(input+linear*tot_size, recv_buf, tot_size);
    }

    delete[] recv_buf;

    QDPInternal::broadcast(checksum);

    // Reconstruct the gauge field
    ColorMatrix  sitefield;
    REAL su3[3][3][2];

    for(int linear=0; linear < nodeSites; ++linear)
    {
      for(int dd=0; dd<Nd; dd++)        /* dir */
      {
	// Transfer the data from input into SU3
	if (float_size == 4) 
	{
	  REAL* su3_p = (REAL *)su3;
	  REAL32* input_p = (REAL32 *)( input+su3_size*(dd+Nd*linear) );
	  for(int cp_index=0; cp_index < mat_size; cp_index++) {
	    su3_p[cp_index] = (REAL)(input_p[cp_index]);
	  }
	}
	else if (float_size == 8) 
	{
	  // IEEE64BIT case
	  REAL *su3_p = (REAL *)su3;
	  REAL64 *input_p = (REAL64 *)( input+su3_size*(dd+Nd*linear) );
	  for(int cp_index=0; cp_index < mat_size; cp_index++) { 
	    su3_p[cp_index] = (REAL)input_p[cp_index];
	  }
	}
	else { 
	  QDPIO::cerr << __func__ << ": Unknown mat size" << std::endl;
	  QDP_abort(1);
	}

	/* Reconstruct the third column  if necessary */
	if (mat_size == 12) 
	{
	  su3[2][0][0] = su3[0][1][0]*su3[1][2][0] - su3[0][1][1]*su3[1][2][1]
	    - su3[0][2][0]*su3[1][1][0] + su3[0][2][1]*su3[1][1][1];
	  su3[2][0][1] = su3[0][2][0]*su3[1][1][1] + su3[0][2][1]*su3[1][1][0]
	    - su3[0][1][0]*su3[1][2][1] - su3[0][1][1]*su3[1][2][0];

	  su3[2][1][0] = su3[0][2][0]*su3[1][0][0] - su3[0][2][1]*su3[1][0][1]
	    - su3[0][0][0]*su3[1][2][0] + su3[0][0][1]*su3[1][2][1];
	  su3[2][1][1] = su3[0][0][0]*su3[1][2][1] + su3[0][0][1]*su3[1][2][0]
	    - su3[0][2][0]*su3[1][0][1] - su3[0][2][1]*su3[1][0][0];
          
	  su3[2][2][0] = su3[0][0][0]*su3[1][1][0] - su3[0][0][1]*su3[1][1][1]
	    - su3[0][1][0]*su3[1][0][0] + su3[0][1][1]*su3[1][0][1];
	  su3[2][2][1] = su3[0][1][0]*su3[1][0][1] + su3[0][1][1]*su3[1][0][0]
	    - su3[0][0][0]*su3[1][1][1] - su3[0][0][1]*su3[1][1][0];
	}

	/* Copy into the big array */
	for(int kk=0; kk<Nc; kk++)      /* color */
	{
	  for(int ii=0; ii<Nc; ii++)    /* color */
	  {
	    Complex sitecomp = cmplx(Real(su3[ii][kk][0]), Real(su3[ii][kk][1]));
	    pokeColor(sitefield,sitecomp,ii,kk);
	  }
	}
      
	u[dd].elem(linear) = sitefield.elem();
      }
    }
  
    delete[] input;
  }


//-----------------------------------------------------------------------
// Write a QCD archive file
// Write a QCD (NERSC) Archive format gauge field
/*
 * \ingroup io
 *
 * \param cfg_out    binary writer object ( Modify )
 * \param u          gauge configuration ( Read )
 */    
  void writeArchiv(BinaryWriter& cfg_out, const multi1d<LatticeColorMatrix>& u,
		   int mat_size)
  {
    size_t size = sizeof(REAL32);
    size_t su3_size = size*mat_size;
    size_t tot_size = su3_size*Nd;
    char *recv_buf = new(std::nothrow) char[tot_size];
    if( recv_buf == 0x0 ) { 
      QDP_error_exit("Unable to allocate recv_buf\n");
    }

    const int nodeSites = Layout::sitesOnNode();

    multi1d<multi1d<ColorMatrix> > sa(Nd);   // extract gauge fields

    for(int dd=0; dd<Nd; dd++)        /* dir */
    {
      sa[dd].resize(nodeSites);
      QDP_extract(sa[dd], u[dd], all);
    }

    // Find the location of each site and send to primary node
    for(int site=0; site < Layout::vol(); ++site)
    {
      multi1d<int> coord = crtesn(site, Layout::lattSize());

      int node   = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);

      // Copy to buffer: be really careful since max(linear) could vary among nodes
      if (Layout::nodeNumber() == node)
      {
	char *recv_buf_tmp = recv_buf;

	for(int dd=0; dd<Nd; dd++)        /* dir */
	{
	  if ( mat_size == 12 ) 
	  {
	    REAL32 su3[2][3][2];

	    for(int kk=0; kk<Nc; kk++)      /* color */
	      for(int ii=0; ii<2; ii++)    /* color */
	      {
		Complex sitecomp = peekColor(sa[dd][linear],ii,kk);
		su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
		su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	      }

	    memcpy(recv_buf_tmp, &(su3[0][0][0]), su3_size);
	  }
	  else
	  {
	    REAL32 su3[3][3][2];

	    for(int kk=0; kk<Nc; kk++)      /* color */
	      for(int ii=0; ii<Nc; ii++)    /* color */
	      {
		Complex sitecomp = peekColor(sa[dd][linear],ii,kk);
		su3[ii][kk][0] = toFloat(Real(real(sitecomp)));
		su3[ii][kk][1] = toFloat(Real(imag(sitecomp)));
	      }

	    memcpy(recv_buf_tmp, &(su3[0][0][0]), su3_size);
	  }

	  recv_buf_tmp += su3_size;
	}
      }

      // Send result to primary node. Avoid sending prim-node sending to itself
      if (node != 0)
      {
#if 1
	// All nodes participate
	QDPInternal::route((void *)recv_buf, node, 0, tot_size);
#else
	if (Layout::primaryNode())
	  QDPInternal::recvFromWait((void *)recv_buf, node, tot_size);

	if (Layout::nodeNumber() == node)
	  QDPInternal::sendToWait((void *)recv_buf, 0, tot_size);
#endif
      }

      cfg_out.writeArrayPrimaryNode(recv_buf, size, mat_size*Nd);
    }

    delete[] recv_buf;

    if (cfg_out.fail())
    {
      QDPIO::cerr << __func__ << ": error writing configuration" << std::endl;
      QDP_abort(1);
    }
  }

} // namespace QDP;
