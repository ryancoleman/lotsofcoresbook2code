//
/*! @file
 * @brief IO support via QIO
 */

#include "qdp.h"

namespace QDP 
{

  //-----------------------------------------
  static int get_node_number(const int coord[])
  {
    multi1d<int> crd(Nd);
    crd = coord;   // an array copy
    int node = Layout::nodeNumber(crd);
    return node;
  }

  static int get_node_index(const int coord[])
  {
    multi1d<int> crd(Nd);
    crd = coord;   // an array copy
    int linear = Layout::linearSiteIndex(crd);
    return linear;
  }

  static void get_coords(int coord[], int node, int linear)
  {
    multi1d<int> crd = Layout::siteCoords(node, linear);
    for(int i=0; i < Nd; ++i)
      coord[i] = crd[i];
  }

  static int get_sites_on_node(int node) 
  {
    return Layout::sitesOnNode();
  }

  // Setting up the QIO Filesystem
  //! Master IO node is always IO node 0
  int master_io_node(void)
  {
	return 0;
  }
	
	
   //! io_node(node) returns the I/O node for node 'node'
   /*! The plan is:
	     a) If there is only 1 I/O node (e.g. scalar) that is the io_node for all
	     b) otherwise if I/O geom is not defined, each node is their own I/O node.
		 c) otherwise block layout by I/O geom to compute the I/O node
	*/
   int io_node(int node) 
   {
	    // If no io grid is defined, each node is its own I/O node
	  
	   if( ! Layout::isIOGridDefined() ) {
			return node;
		}
		
		// A grid is defined. A shortcut for when there is only 1 
	   if ( Layout::numIONodeGrid() == 1 ) {
			return master_io_node();
		}
		
	    // Compute my I/O node. Block lat size into I/O grid
	    const multi1d<int>& proc_size = Layout::logicalSize();
	    const multi1d<int>& io_geom = Layout::getIONodeGrid();

	    multi1d<int> block_sizes(Nd);
	    for(int mu=0; mu < Nd; mu++) { 
		block_sizes[mu] = proc_size[mu]/io_geom[mu];
						
	        // Pick up slack if CPU dimension not divisible
		// Last block will stretch beyond processor grid, but this is OK
		// since no node outside the proc grid will be asked for
		if( proc_size[mu] % io_geom[mu] != 0 ) block_sizes[mu]++;
	    }
		
	    // My node coords -- always in the processor grid
	    multi1d<int> node_coords=Layout::getLogicalCoordFrom(node);
		
	    // Coords of I/O node, basically the origin of the block the 
	    // current node is in
		multi1d<int> io_node_coords(Nd);
		for(int mu=0; mu < Nd; mu++) { 
			// Integer division: will truncate to origin of block.
			io_node_coords[mu] = node_coords[mu] / block_sizes[mu];
		        io_node_coords[mu] *= block_sizes[mu];
		}
			
	    // Now just convert the io_node_coords to a node number and we're done
	   return Layout::getNodeNumberFrom(io_node_coords);
	}

	//! io_node for multifile...
	/*! code gets confused if multifile is set, but filesystem is not 'multfile' */
	int io_node_multifile(int node) 
	{
		return node;
	}
	
  //-----------------------------------------------------------------------------
  // QDP QIO support
  QDPFileReader::QDPFileReader() {iop=false;}

  QDPFileReader::QDPFileReader(XMLReader& xml, 
			       const std::string& path,
			       QDP_serialparallel_t serpar)
  {open(xml,path,serpar);}

  void QDPFileReader::open(XMLReader& file_xml, 
			   const std::string& path, 
			   QDP_serialparallel_t serpar)
  {
    QIO_Layout layout;
		  
    int latsize[Nd];

    for(int m=0; m < Nd; ++m)
      latsize[m] = Layout::lattSize()[m];

    layout.node_number = &get_node_number;
    layout.node_index  = &get_node_index;
    layout.get_coords  = &get_coords;
    layout.num_sites = &get_sites_on_node;
    layout.latsize = latsize;
    layout.latdim = Nd; 
    layout.volume = Layout::vol(); 
    layout.sites_on_node = Layout::sitesOnNode(); 
    layout.this_node = Layout::nodeNumber(); 
    layout.number_of_nodes = Layout::numNodes(); 

    QIO_Filesystem fs;
    fs.my_io_node = &io_node;
    fs.master_io_node = &master_io_node;
	  
    // Initialize string objects 
    QIO_String *xml_c  = QIO_string_create();

    QIO_Iflag iflag;
    if( serpar == QDPIO_PARALLEL ) {   
      iflag.serpar = QIO_PARALLEL;
    }
    else { 
      iflag.serpar = QIO_SERIAL;
    }

    iflag.volfmt = QIO_UNKNOWN;

    // Call QIO read
    // At this moment, serpar (which is an enum in QDP++) is ignored here.
    if ((qio_in = QIO_open_read(xml_c, path.c_str(), &layout, &fs, &iflag)) == NULL)
    {
      iostate = QDPIO_badbit;  // not helpful

      QDPIO::cerr << "QDPFileReader: failed to open file " << path << std::endl;
      QDP_abort(1);  // just bail, otherwise xml stuff below fails.
    }
    else
    {
      iostate = QDPIO_goodbit;
    }

    // Use string to initialize XMLReader
    std::istringstream ss;
    if (Layout::primaryNode())
    {
      std::string foo = QIO_string_ptr(xml_c);
      ss.str(foo);
    }
    file_xml.open(ss);

    QIO_string_destroy(xml_c);

    iop=true;
  }


  void QDPFileReader::close()
  {
    if (is_open()) 
    {
      //int status = QIO_close_read(qio_in);
      QIO_close_read(qio_in);
    }

    iop = false;
    iostate = QDPIO_badbit;
  }

  bool QDPFileReader::is_open() {return iop;}

  bool QDPFileReader::eof() const {return false;}

  bool QDPFileReader::bad() const {return iostate;}

  void QDPFileReader::clear(QDP_iostate_t state)
  {
    iostate = state;
  }

  QDPFileReader::~QDPFileReader() {close();}

  //! Close a QDPFileReader
  void close(QDPFileReader& qsw)
  {
    qsw.close();
  }

  //! Is a QDPFileReader open
  bool is_open(QDPFileReader& qsw)
  {
    return qsw.is_open();
  }

  // Reads a record header only (state of reader left intact
  // so subsequent read to get header and data still works)

  void QDPFileReader::read(XMLReader& rec_xml)
  {
    QIO_RecordInfo rec_info;
    QIO_String* xml_c = QIO_string_create();
    int status;
  
    status = QIO_read_record_info(qio_in, &rec_info, xml_c);
    if( status != QIO_SUCCESS) { 
      QDPIO::cerr << "Failed to read the Record Info" << std::endl;
      QDP_abort(1);
    }
  
    std::istringstream ss;
    if (Layout::primaryNode()) {
      std::string foo = QIO_string_ptr(xml_c);
      ss.str(foo);
    }
    rec_xml.open(ss);
  
    QIO_string_destroy(xml_c);
  }


  // Reads a BinaryBufferReader object
  /*!
    \param rec_xml The (user) record metadata.
    \param sl The data
  */
  void QDPFileReader::read(XMLReader& rec_xml, BinaryBufferReader& s1)
  {
    QIO_RecordInfo rec_info;
    QIO_String* xml_c = QIO_string_create();
    int status;
  
    status = QIO_read_record_info(qio_in, &rec_info, xml_c);
    if( status != QIO_SUCCESS) { 
      QDPIO::cerr << "Failed to read the Record Info" << std::endl;
      QDP_abort(1);
    }
  
    QDPIO::cout << "BinaryBufferRead" << std::endl;
    std::string from_disk;
    from_disk.resize(QIO_get_datacount(&rec_info));
    status = QIO_read_record_data(qio_in,
				  &(QDPOScalarFactoryPut<char> ),
				  from_disk.size()*sizeof(char),
				  sizeof(char),
				  (void *)&(from_disk[0]));
    if (status != QIO_SUCCESS) { 
      QDPIO::cerr << "Failed to read data" << std::endl;
      clear(QDPIO_badbit);
      QDP_abort(1);
    }
    QDPIO::cout << "QIO_read_finished" << std::endl;
      
    // Cast appropriately
//    for(int i=0; i < from_disk.size(); i++) { 
    s1.open(from_disk);
//    }
  
    std::istringstream ss;
    if (Layout::primaryNode()) {
      std::string foo = QIO_string_ptr(xml_c);
      ss.str(foo);
    }
    rec_xml.open(ss);
  
    QIO_string_destroy(xml_c);
  }


  //-----------------------------------------------------------------------------
  // QDP QIO support (writers)
  QDPFileWriter::QDPFileWriter() {iop=false;}

  QDPFileWriter::QDPFileWriter(XMLBufferWriter& xml, 
			       const std::string& path,
			       QDP_volfmt_t qdp_volfmt,
			       QDP_serialparallel_t qdp_serpar,
			       QDP_filemode_t qdp_mode) 
  {
    open(xml,path,qdp_volfmt,qdp_serpar,qdp_mode, std::string());
  }

  QDPFileWriter::QDPFileWriter(XMLBufferWriter& xml, 
			       const std::string& path,
			       QDP_volfmt_t qdp_volfmt,
			       QDP_serialparallel_t qdp_serpar,
			       QDP_filemode_t qdp_mode,
			       const std::string& data_LFN) 
  {
    open(xml,path,qdp_volfmt,qdp_serpar,qdp_mode, data_LFN);
  }

  // filemode not specified
  void QDPFileWriter::open(XMLBufferWriter& file_xml, 
			   const std::string& path,
			   QDP_volfmt_t qdp_volfmt,
			   QDP_serialparallel_t qdp_serpar)
  {
    open(file_xml,path,qdp_volfmt,qdp_serpar,QDPIO_OPEN, std::string());
  }

  void QDPFileWriter::open(XMLBufferWriter& file_xml, 
			   const std::string& path,
			   QDP_volfmt_t qdp_volfmt,
			   QDP_serialparallel_t qdp_serpar,
			   const std::string& data_LFN) 
  {
    open(file_xml,path,qdp_volfmt,qdp_serpar,QDPIO_OPEN, data_LFN);
  }


  // filemode not specified
  QDPFileWriter::QDPFileWriter(XMLBufferWriter& xml, 
			       const std::string& path,
			       QDP_volfmt_t qdp_volfmt,
			       QDP_serialparallel_t qdp_serpar) 
  {
    open(xml,path,qdp_volfmt,qdp_serpar,QDPIO_OPEN, std::string());
  }

  // filemode not specified
  QDPFileWriter::QDPFileWriter(XMLBufferWriter& xml, 
			       const std::string& path,
			       QDP_volfmt_t qdp_volfmt,
			       QDP_serialparallel_t qdp_serpar,
			       const std::string& data_LFN) 
  {
    open(xml,path,qdp_volfmt,qdp_serpar,QDPIO_OPEN, data_LFN);
  }

  void QDPFileWriter::open(XMLBufferWriter& file_xml, 
			   const std::string& path,
			   QDP_volfmt_t qdp_volfmt,
			   QDP_serialparallel_t qdp_serpar,
			   QDP_filemode_t qdp_mode, 
			   const std::string& data_LFN) 
  {

    QIO_Layout layout;
    int latsize[Nd];

    for(int m=0; m < Nd; ++m)
      latsize[m] = Layout::lattSize()[m];

    layout.node_number = &get_node_number;
    layout.node_index  = &get_node_index;
    layout.get_coords  = &get_coords;
    layout.num_sites = &get_sites_on_node;
    layout.latsize = latsize;
    layout.latdim = Nd; 
    layout.volume = Layout::vol(); 
    layout.sites_on_node = Layout::sitesOnNode(); 
    layout.this_node = Layout::nodeNumber(); 
    layout.number_of_nodes = Layout::numNodes(); 

    // Copy metadata string into simple qio string container
    QIO_String* xml_c = QIO_string_create();
    QIO_string_set(xml_c, file_xml.str().c_str());

    if (xml_c == NULL)
    {
      QDPIO::cerr << "QDPFileWriter - error in creating QIO string" << std::endl;
      iostate = QDPIO_badbit;
    }
    else
    {
      iostate = QDPIO_goodbit;
    }

	QIO_Filesystem fs;
	fs.my_io_node = &io_node;
	fs.master_io_node = &master_io_node;
	  
    // Wrappers over simple ints
    int volfmt;
    switch(qdp_volfmt)
    {
    case QDPIO_SINGLEFILE:
      volfmt = QIO_SINGLEFILE;
      break;
    
    case QDPIO_MULTIFILE:
      volfmt = QIO_MULTIFILE;
	  fs.my_io_node = &io_node_multifile;
      break;

    case QDPIO_PARTFILE:
      volfmt = QIO_PARTFILE;
      
      break;

    default: 
      QDPIO::cerr << "Unknown value for qdp_volfmt " << qdp_volfmt << std::endl;
      QDP_abort(1);
      return;
    }
  
    // Wrappers over simple ints
    int mode;
    switch(qdp_mode)
    {
    case QDPIO_CREATE:
      mode = QIO_CREAT;
      break;
    
    case QDPIO_OPEN:
      mode = QIO_TRUNC;
      break;

    case QDPIO_APPEND:
      mode = QIO_APPEND;
      break;

    default: 
      QDPIO::cerr << "Unknown value for qdp_mode " << qdp_mode << std::endl;
      QDP_abort(1);
      return;
    }
  
    // QIO write
    // For now, serpar (which is an enum in QDP) is ignored here
    QIO_Oflag oflag;
    if( qdp_serpar == QDPIO_SERIAL ) { 
      oflag.serpar = QIO_SERIAL;
    }
    else { 
      oflag.serpar = QIO_PARALLEL;
    }
    oflag.mode   = mode;
    oflag.ildgstyle = QIO_ILDGLAT;
    if( data_LFN.length() == 0 ) { 
      oflag.ildgLFN = NULL;
    }
    else {
      oflag.ildgLFN = QIO_string_create();
      QIO_string_set(oflag.ildgLFN, data_LFN.c_str());
    }

	  
    // This is the QIO Way - older way 
    if ((qio_out = QIO_open_write(xml_c, path.c_str(), 
				  volfmt, 
				  &layout, 
				  &fs, &oflag)) == NULL )
    {
      iostate = QDPIO_badbit;  // not helpful

      QDPIO::cerr << "QDPFileWriter: failed to open file " << path << std::endl;
      QDP_abort(1);  // just bail. Not sure I want this. This is not stream semantics
    }
    else
    {
      iostate = QDPIO_goodbit;
    }


    // Free memory -- this is OK< as it should'of been copied
    QIO_string_destroy(oflag.ildgLFN);
    // Cleanup
    QIO_string_destroy(xml_c);

    iop=true;
  }

  void QDPFileWriter::close()
  {
    if (is_open()) 
    {
      // int status = QIO_close_write(qio_out);
      QIO_close_write(qio_out);
    }

    iop = false;
    iostate = QDPIO_badbit;
  }

  bool QDPFileWriter::is_open() {return iop;}

  bool QDPFileWriter::bad() const {return iostate;}

  void QDPFileWriter::clear(QDP_iostate_t state)
  {
    iostate = state;
  }

  QDPFileWriter::~QDPFileWriter() {close();}


  void close(QDPFileWriter& qsw)
  {
    qsw.close();
  }


  bool is_open(QDPFileWriter& qsw)
  {
    return qsw.is_open();
  }


  //! Writes a BinaryBufferWriter object
  /*!
    \param rec_xml The (user) record metadata.
    \param sl The data
  */
  void QDPFileWriter::write(XMLBufferWriter& rec_xml, BinaryBufferWriter& s1)
  {
    std::string ss = s1.str();
    char *typestr=(char *)"char";
    char *signtype=(char *)"U";

    QIO_RecordInfo* info = QIO_create_record_info(QIO_GLOBAL, NULL, NULL, 0,
						  typestr,
						  signtype,
						  0, 0, 
						  sizeof(char), ss.size());


    // Copy metadata string into simple qio string container
    QIO_String* xml_c = QIO_string_create();
    if (xml_c == NULL)
    {
      QDPIO::cerr << "QDPFileWriter::write - error in creating XML string" << std::endl;
      QDP_abort(1);
    }

    if (Layout::primaryNode())
      QIO_string_set(xml_c, rec_xml.str().c_str());

    // Big call to qio
    if (QIO_write(get(), info, xml_c,
		  &(QDPOScalarFactoryGet<char>),
		  ss.size()*sizeof(char), 
		  sizeof(char), 
		  (void *)ss.c_str()) != QIO_SUCCESS)
    {
      QDPIO::cerr << "QDPFileWriter: error in write" << std::endl;
      clear(QDPIO_badbit);
      QDP_abort(1);
    }

    // Cleanup
    QIO_string_destroy(xml_c);
    QIO_destroy_record_info(info);
  }


} // namespace QDP;
