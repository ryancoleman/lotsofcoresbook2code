// -*- C++ -*-
/*! \file
 *  \brief A Map Object that works lazily from Disk
 */


#ifndef __qdp_map_obj_disk_h__
#define __qdp_map_obj_disk_h__

#include "qdp_map_obj.h"
#include <unordered_map>

namespace QDP
{

  namespace MapObjDiskEnv { 
    typedef unsigned int file_version_t;
   
    //! Get the file magic
    std::string getFileMagic();

    //! Get the meta-data from a file
    std::string getMetaData(const std::string& filename);

    //! Check if this will be a new file
    bool checkForNewFile(const std::string& filename, std::ios_base::openmode mode);
  };




  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectDisk : public MapObject<K,V>
  {
  public:
    //! Empty constructor
    MapObjectDisk() : file_version(1), state(INIT), level(0) {}

    //! Finalizes object
    ~MapObjectDisk();

    //! Set debugging level
    void setDebug(int level);

    //! Get debugging level
    int getDebug() const {return level;}

    //! Open a file
    void open(const std::string& file, std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out);

    //! Check if a DB file exists before opening.
    bool fileExists(const std::string& file) const {
      return (! MapObjDiskEnv::checkForNewFile(file, std::ios_base::in));
    }

    //! Close the file
    void close();

    /**
     * Insert a pair of data and key into the database
     * @param key a key
     * @param val a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    int insert(const K& key, const V& val);

    /**
     * Get data for a given key
     * @param key user supplied key
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int get(const K& key, V& val) const;


    /**
     * Flush database in memory to disk
     */
    void flush();


    /**
     * Does this key exist in the store
     * @param key a key object
     * @return true if the answer is yes
     */
    bool exist(const K& key) const;

    /** 
     * The number of elements
     */
    unsigned int size() const {return static_cast<unsigned long>(src_map.size());}

    /**
     * Return all available keys to user
     * @param keys user suppled an empty vector which is populated
     * by keys after this call.
     */
    void keys(std::vector<K>& keys_) const;
    
    /**
     * Insert user data into the  metadata database
     *
     * @param user_data user supplied data
     * @return returns 0 if success, else failure
     */
    int insertUserdata(const std::string& user_data);
    
    /**
     * Get user user data from the metadata database
     *
     * @param user_data user supplied buffer to store user data
     * @return returns 0 if success. Otherwise failure.
     */
    int getUserdata(std::string& user_data) const;

  public:
    typedef std::iostream::pos_type pos_type;  // position in buffer
    typedef std::iostream::off_type off_type;  // offset in buffer

  private:
    //! Sometime ago, chose 16 bytes for stream positions. Yup, really big.
    union priv_pos_type_t
    {
      unsigned char     c[16];
      uint64_t          p;
    };

    //! Type for the map
    typedef std::unordered_map<std::string, priv_pos_type_t> MapType_t;

    //! State 
    enum State {INIT, UNCHANGED, MODIFIED};

    //! State of db
    State state;

    //! Debugging
    int level;

    //! File related stuff. Unsigned int is as close to uint32 as I can get
    MapObjDiskEnv::file_version_t file_version;
    
    //! Map of objects
    mutable MapType_t src_map;
    
    //! The parameters
    std::string filename;
    
    //! Metadata
    std::string user_data;

    //! Reader and writer interfaces
    mutable BinaryFileReaderWriter streamer;
    
    //! Convert to known size
    priv_pos_type_t convertToPrivate(const pos_type& input) const;

    //! Convert from known size
    pos_type convertFromPrivate(const priv_pos_type_t& input) const;

    //! Open a new DB, and will write map
    void openWrite(const std::string& file, std::ios_base::openmode mode);

    //! Open an existing DB, and read map
    void openRead(const std::string& file, std::ios_base::openmode mode);

    // Internal Utility: Create/Skip past header
    void writeSkipHeader(void);
    
    //! Internal Utility: Read/Check header 
    priv_pos_type_t readCheckHeader(void);
    
    //! Internal Utility: Dump the map to disk
    void writeMapBinary(void);  
    
    //! Internal Utility: Read the map from disk
    void readMapBinary(const priv_pos_type_t& md_start);
    
    //! Internal Utility: Close File after write mode
    void closeWrite(void);
    
    //! Sink State for errors:
    void errorState(const std::string err) const {
      throw err;
    }
  };


  
  /* ****************** IMPLEMENTATIONS ********************** */
  
  //! Set debugging level
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::setDebug(int level_)
  {
    level = level_;
  }


  //! Convert to known size
  template<typename K, typename V>
  typename MapObjectDisk<K,V>::priv_pos_type_t
  MapObjectDisk<K,V>::convertToPrivate(const pos_type& input) const
  {
    priv_pos_type_t f;
    f.p = static_cast<uint64_t>(input);
    return f;
  }

  //! Convert from known size
  template<typename K, typename V>
  typename MapObjectDisk<K,V>::pos_type 
  MapObjectDisk<K,V>::convertFromPrivate(const priv_pos_type_t& input) const
  {
    return static_cast<pos_type>(input.p);
  }


  //! Open a file
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::open(const std::string& file, std::ios_base::openmode mode)
  {
    if ( MapObjDiskEnv::checkForNewFile(file, mode) )
    {
      openWrite(file, mode);
    }
    else
    {
      openRead(file, mode);
    }
  }


  //! Open a new DB, and will write map
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::openWrite(const std::string& file, std::ios_base::openmode mode)
  {
    switch(state) { 
    case INIT: 
    {
      filename = file;

      QDPIO::cout << "MapObjectDisk: opening file " << filename
		  << " for writing" << std::endl;
      
      streamer.open(filename, mode);
          
      if (level >= 2) {
        QDPIO::cout << "sizeof(unsigned char) = " << sizeof(unsigned char) << std::endl;
        QDPIO::cout << "sizeof(int) = " << sizeof(int) << std::endl;
        QDPIO::cout << "sizeof(pos_type) = " << sizeof(pos_type) << std::endl;
        QDPIO::cout << "sizeof(priv_pos_type_t) = " << sizeof(priv_pos_type_t) << std::endl;
        QDPIO::cout << "sizeof(file_version)t) = " << sizeof(MapObjDiskEnv::file_version_t) << std::endl;
      }

      if (level >= 2) {
	QDPIO::cout << "Writing file magic: len= " << MapObjDiskEnv::getFileMagic().length() << std::endl;
      }

      // Write string
      streamer.writeDesc(MapObjDiskEnv::getFileMagic());
      
      if (level >= 2) {
	QDPIO::cout << "Wrote magic. Current Position: " << streamer.currentPosition() << std::endl;
      }
      
      write(streamer, (MapObjDiskEnv::file_version_t)file_version);
    
      if (level >= 2) {
	QDPIO::cout << "Wrote Version. Current Position is: " << streamer.currentPosition() << std::endl;
      }
      
      if (level >= 2) {
	QDPIO::cout << "Writing User Data string=" << user_data << std::endl;
      }
      writeDesc(streamer, user_data);
      
      if (level >= 2) {
	QDPIO::cout << "Wrote User Data string. Current Position is: " << streamer.currentPosition() << std::endl;
      }
      
      priv_pos_type_t dummypos = convertToPrivate(streamer.currentPosition());
    
      if (level >= 2) {
	int user_len = user_data.length();
	QDPInternal::broadcast(user_len);

	QDPIO::cout << "Sanity Check 1" << std::endl; ;
	uint64_t cur_pos = convertToPrivate(streamer.currentPosition()).p;
	uint64_t exp_pos = 
	  MapObjDiskEnv::getFileMagic().length()+sizeof(int)
	  +user_len+sizeof(int)
	  +sizeof(MapObjDiskEnv::file_version_t);

	QDPIO::cout << "cur pos=" << (size_t)(cur_pos) << " expected " << (size_t)(exp_pos) << std::endl;

	if ( cur_pos != exp_pos ) {
	  QDPIO::cout << "ERROR: Sanity Check 1 failed." << std::endl;
	  QDPIO::cout << "cur pos=" << (size_t)(cur_pos) << " expected " << (size_t)(exp_pos) << std::endl;
	  QDP_abort(1);
	}
      }
    
      /* Write a dummy link - make room for it */
      streamer.writeArray((char *)&dummypos, sizeof(priv_pos_type_t), 1);
      
      if (level >= 2) {
	QDPIO::cout << "Wrote dummy link: Current Position " << streamer.currentPosition() << std::endl;
	int user_len = user_data.length();
	QDPInternal::broadcast(user_len);

	QDPIO::cout << "Sanity Check 2" << std::endl;
	uint64_t cur_pos = convertToPrivate(streamer.currentPosition()).p;
	uint64_t exp_pos = 
	  MapObjDiskEnv::getFileMagic().length()+sizeof(int)
	  +user_len+sizeof(int)
	  +sizeof(MapObjDiskEnv::file_version_t)
	  +sizeof(priv_pos_type_t);

	if ( cur_pos != exp_pos ) {
	  QDPIO::cout << "Cur pos = " << (size_t)(cur_pos) << std::endl;
	  QDPIO::cout << "Expected: " << (size_t)(exp_pos) << std::endl;
	  QDPIO::cout << "ERROR: Sanity Check 2 failed." << std::endl;
	  QDP_abort(1);
	}
	QDPIO::cout << "Finished sanity Check 2" << std::endl;
      }
      
      // Advance state machine state
      state = MODIFIED;
      break;      
    }

    default:
      errorState("MapOjectDisk: openWrite called from invalid state");
      break;
    }
  
    return;
  }

 

  //! Open an existing DB, and read map
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::openRead(const std::string& file, std::ios_base::openmode mode)
  {  
    switch (state) { 
    case INIT:
    {
      filename = file;

      QDPIO::cout << "MapObjectDisk: opening file " << filename
		  << " for reading" << std::endl;
      
      // Open the reader
      streamer.open(filename, mode);
	
      QDPIO::cout << "MapObjectDisk: reading and checking header" << std::endl;

      priv_pos_type_t md_start = readCheckHeader();
	
      // Seek to metadata
      QDPIO::cout << "MapObjectDisk: reading key/fileposition data" << std::endl;
	
      /* Read the map in (metadata) */
      readMapBinary(md_start);
	
      /* And we are done */
      state = UNCHANGED;
    }
    break;
    default:
      errorState("MapObjectDisk: openRead() called from invalid state");
      break;
    }
    
    return;
  }


  
  //! Close
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::close() 
  {
    switch(state) { 
    case UNCHANGED:
      if( streamer.is_open() ) { 
	streamer.close();
      }
      break;
    case MODIFIED:
      closeWrite(); // This finalizes files for us
      if( streamer.is_open() ) { 
	streamer.close();
      }
      break;
    case INIT:
      break;
    default:
      errorState("close: destructor called from invalid state");
      break;
    }

    state = INIT;
  }
  

  //! Destructor
  template<typename K, typename V>
  MapObjectDisk<K,V>::~MapObjectDisk() 
  {
    close();
  }
  

  //! Destructor
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::flush() 
  {
    switch(state) { 
    case MODIFIED:
      closeWrite();  // not optimal
      state = UNCHANGED;
      break;
    case UNCHANGED:
      break;
    case INIT:
      break;
    default:
      break;
    }
  }
  

  //! Dump keys
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::keys(std::vector<K>& keys_) const 
  {
    if( streamer.is_open() ) 
    {
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter) 
      { 
	BinaryBufferReader bin(iter->first);
	K key;
	read(bin, key);
	keys_.push_back(key);
      }
    }
  }
    

  /**
   * Insert user data into the  metadata database
   */
  template<typename K, typename V>
  int 
  MapObjectDisk<K,V>::insertUserdata(const std::string& _user_data)
  {
    int ret = 0;
    switch(state) { 
    case INIT:
      user_data = _user_data;
      break;

    case UNCHANGED:
    case MODIFIED:
      ret = 1;
      break;

    default:
      errorState("MapObjectDisk::insertUserdata() called from invalid state");
      break;
    }

    return ret;
  }

    
  /**
   * Get user user data from the metadata database
   *
   * @param user_data user supplied buffer to store user data
   * @return returns 0 if success. Otherwise failure.
   */
  template<typename K, typename V>
  int
  MapObjectDisk<K,V>::getUserdata(std::string& _user_data) const
  {
    int ret = 0;

    switch(state) { 
    case INIT:
    {
      ret = 1;
      break;
    }
    case UNCHANGED:
    case MODIFIED:
    {
      _user_data = user_data;
      break;
    }
    default:
      errorState("MapObjectDisk::getUserdata called from unknown state");
      break;
    }

    return ret;
  }


  /*! 
   * Insert a value into the Map.
   */
  template<typename K, typename V>
  int 
  MapObjectDisk<K,V>::insert(const K& key, const V& val) 
  {
    int ret = 0;

    switch (state)  { 
    case MODIFIED :
    case UNCHANGED : {
      //  Find key
      BinaryBufferWriter bin;
      write(bin, key);
      typename MapType_t::const_iterator key_ptr = src_map.find(bin.str());

      if (key_ptr != src_map.end()) { 
	// Key does exist
	pos_type wpos = convertFromPrivate(key_ptr->second);
	if (level >= 2) {
	  QDPIO::cout << "Found key to update. Position is " << wpos << std::endl;
	}
	
	streamer.seek(wpos);
	
	if (level >= 2) {
	  QDPIO::cout << "Sought write position. Current Position: " << streamer.currentPosition() << std::endl;
	}
	streamer.resetChecksum(); // Reset checksum. It gets calculated on write.
	write(streamer, val);
	
	if (level >= 2) {	
	  QDPIO::cout << "Wrote value to disk. Current Position: " << streamer.currentPosition() << std::endl;
	}
	write(streamer, streamer.getChecksum()); // Write Checksum
	streamer.flush(); // Sync the file
	
	if (level >= 2) {
	  QDPIO::cout << "Wrote checksum " << streamer.getChecksum() << " to disk. Current Position: " << streamer.currentPosition() << std::endl;
	}

	// Done
	state = MODIFIED;
      }
      else {
	// Key does not exist

	// Make note of current writer position
	priv_pos_type_t pos = convertToPrivate(streamer.currentPosition());
      
	// Insert pos into map
	src_map.insert(std::make_pair(bin.str(),pos));
     
	streamer.resetChecksum();

	// Add position to the map
	StopWatch swatch;
	swatch.reset();
	swatch.start();

	write(streamer, val); // DO write
	swatch.stop();
     
	// Get diagnostics.
	if (level >= 1) {
	  priv_pos_type_t end_pos = convertToPrivate(streamer.currentPosition());
	  double MiBWritten = (double)(end_pos.p - pos.p)/(double)(1024*1024);
	  double time = swatch.getTimeInSeconds();

	  QDPIO::cout << " wrote: " << MiBWritten << " MiB. Time: " << time << " sec. Write Bandwidth: " << MiBWritten/time<<std::endl;
        }

	if (level >= 2) {
	  QDPIO::cout << "Wrote value to disk. Current Position: " << streamer.currentPosition() << std::endl;
	}

	write(streamer, streamer.getChecksum()); // Write Checksum
	streamer.flush();
	
	if (level >= 2) {
	  QDPIO::cout << "Wrote checksum " << streamer.getChecksum() << " to disk. Current Position: " << streamer.currentPosition() << std::endl;
	}

	// Done
	state = MODIFIED;
      }
      break;
    }
    default:
      ret = 1;
      break;
    }

    return ret;
  }



  /*! 
   * Lookup an item in the map.
   */
  template<typename K, typename V>
  int 
  MapObjectDisk<K,V>::get(const K& key, V& val) const
  { 
    int ret = 0;

    switch(state) { 
    case UNCHANGED: // Deliberate fallthrough
    case MODIFIED: {
      BinaryBufferWriter bin;
      write(bin, key);
      typename MapType_t::const_iterator key_ptr = src_map.find(bin.str());

      if (key_ptr != src_map.end())
      {
	// If key exists find file offset
	priv_pos_type_t pos = key_ptr->second;

	// Do the seek and time it 
	StopWatch swatch;

	swatch.reset();
	swatch.start();
	streamer.seek(convertFromPrivate(pos));
	swatch.stop();
	double seek_time = swatch.getTimeInSeconds();

	// Reset the checkums
	streamer.resetChecksum();

	// Grab start pos: We've just seeked it
	priv_pos_type_t start_pos = pos;

	// Time the read
	swatch.reset();
	swatch.start();
	read(streamer, val);
	swatch.stop();

	double read_time = swatch.getTimeInSeconds();
	priv_pos_type_t end_pos = convertToPrivate(streamer.currentPosition());

	// Print data
	if (level >= 1) { 
	  double MiBRead = (double)(end_pos.p - start_pos.p)/(double)(1024*1024);
	  QDPIO::cout << " seek time: " << seek_time 
	  	      << " sec. read time: " << read_time 
		      << "  " << MiBRead <<" MiB, " << MiBRead/read_time << " MiB/sec" << std::endl;
        }


	if (level >= 2) { 
	  QDPIO::cout << "Read record. Current position: " << streamer.currentPosition() << std::endl;
	}

	QDPUtil::n_uint32_t calc_checksum=streamer.getChecksum();
	QDPUtil::n_uint32_t read_checksum;
	read(streamer, read_checksum);

	if (level >= 2) {
	  QDPIO::cout << " Record checksum: " << read_checksum << "  Current Position: " << streamer.currentPosition() << std::endl;
	}

	if( read_checksum != calc_checksum ) { 
	  QDPIO::cout << "Mismatched Checksums: Expected: " << calc_checksum << " but read " << read_checksum << std::endl;
	  QDP_abort(1);
	}

	if (level >= 2) {
	  QDPIO::cout << "  Checksum OK!" << std::endl;
	}
      }
      else {
	ret = 1;
      }
      break;
    }
    default:
      ret = 1;
      break;
    }

    return ret;
  }
  
  
  /**
   * Does this key exist in the store
   * @param key a key object
   * @return true if the answer is yes
   */
  template<typename K, typename V>
  bool 
  MapObjectDisk<K,V>::exist(const K& key) const 
  {
    BinaryBufferWriter bin;
    write(bin, key);
    return (src_map.find(bin.str()) == src_map.end()) ? false : true;
  }
  
  

  /***************** UTILITY ******************/


  //! Skip past header
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::writeSkipHeader(void) 
  { 
    switch(state) { 
    case MODIFIED: {
      if ( streamer.is_open() ) 
      { 
	int user_len = user_data.length();
	QDPInternal::broadcast(user_len);
	
	streamer.seek( MapObjDiskEnv::getFileMagic().length() + sizeof(int)
		       + user_len + sizeof(int)
		       + sizeof(MapObjDiskEnv::file_version_t) );
      }
      else { 
	QDPIO::cerr << "Attempting writeSkipHeader, not in write mode" <<std::endl;
	QDP_abort(1);
      }
    }
    break;
    default:
      errorState("MapObjectDisk: writeSkipHeader() called not in MODIFIED state");
      break;
    }
  }
  
  //! Check the header 
  template<typename K, typename V>
  typename MapObjectDisk<K,V>::priv_pos_type_t
  MapObjectDisk<K,V>::readCheckHeader(void) 
  {
    priv_pos_type_t md_position;
    bzero(&md_position.c, sizeof(priv_pos_type_t));

    if( streamer.is_open() ) 
    {
      if (level >= 2) {
	QDPIO::cout << "Rewinding File" << std::endl;
      }
      
      streamer.rewind();
      
      std::string read_magic;
      streamer.readDesc(read_magic);
      
      // Check magic
      if (read_magic != MapObjDiskEnv::getFileMagic()) { 
	QDPIO::cerr << "Magic String Wrong: Expected: " << MapObjDiskEnv::getFileMagic() 
		    << " but read: " << read_magic << std::endl;
	QDP_abort(1);
      }

      if (level >= 2) {
	QDPIO::cout << "Read File Magic. Current Position: " << streamer.currentPosition() << std::endl;
      }
      
      MapObjDiskEnv::file_version_t read_version;
      read(streamer, read_version);
      
      if (level >= 2) {
	QDPIO::cout << "Read File Verion. Current Position: " << streamer.currentPosition() << std::endl;
      }
      
      // Check version
      QDPIO::cout << "MapObjectDisk: file has version: " << read_version << std::endl;
      
      QDP::readDesc(streamer, user_data);
      if (level >= 2) {
	QDPIO::cout << "User data. String=" << user_data << ". Current Position: " << streamer.currentPosition() << std::endl;
      }
      
      // Read MD location
      streamer.readArray((char *)&md_position, sizeof(priv_pos_type_t), 1);

      if (level >= 2) {
	QDPIO::cout << "Read MD Location. Current position: " << streamer.currentPosition() << std::endl;
      }
      
      if (level >= 2) {
	QDPIO::cout << "Metadata starts at position: " << convertFromPrivate(md_position) << std::endl;
      }
	
    }
    else { 
      QDPIO::cerr << "readCheckHeader needs reader mode to be opened. It is not" << std::endl;
      QDP_abort(1);
    }

    return md_position;
  }

  //! Dump the map
  // Private utility function -- no one else should use.
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::writeMapBinary(void)  
  {
    unsigned int map_size = src_map.size();

    streamer.resetChecksum();
    write(streamer, map_size);
    if (level >= 2) {
      QDPIO::cout << "Wrote map size: " << map_size << " entries.  Current position : " << streamer.currentPosition() << std::endl;
    }
    
    typename MapType_t::const_iterator iter;
    for(iter  = src_map.begin();
	iter != src_map.end();
	++iter) 
    { 
      priv_pos_type_t pos=iter->second;
      
      writeDesc(streamer, iter->first); 
      streamer.writeArray((char *)&pos,sizeof(priv_pos_type_t),1);
      
      if (level >= 2) {
	QDPIO::cout << "Wrote Key/Position pair:  Current Position: " << streamer.currentPosition() << std::endl;
      }
    }
    write(streamer, streamer.getChecksum());
    QDPIO::cout << "Wrote Checksum On Map: " << streamer.getChecksum() << std::endl;
    streamer.flush();
  }
  
  //! read the map 
  // assume positioned at start of map data
  // Private utility function -- no one else should use.
  template<typename K, typename V>
  void 
  MapObjectDisk<K,V>::readMapBinary(const priv_pos_type_t& md_start)
  {
    streamer.seek(convertFromPrivate(md_start));
    streamer.resetChecksum();

    if (level >= 2) {
      QDPIO::cout << "Sought start of metadata. Current position: " << streamer.currentPosition() << std::endl;
    }
    
    unsigned int num_records;
    read(streamer, num_records);

    if (level >= 2) {
      QDPIO::cout << "Read num of entries: " << num_records << " records. Current Position: " << streamer.currentPosition() << std::endl;
    }
    
    for(unsigned int i=0; i < num_records; i++) 
    { 
      priv_pos_type_t rpos;
      std::string key_str;
      readDesc(streamer, key_str);
      
      streamer.readArray((char *)&rpos, sizeof(priv_pos_type_t),1);
      
      if (level >= 2) {
	QDPIO::cout << "Read Key/Position pair. Current position: " << streamer.currentPosition() << std::endl;
      }
      // Add position to the map
      src_map.insert(std::make_pair(key_str,rpos));
    }
    QDPUtil::n_uint32_t calc_checksum = streamer.getChecksum();
    QDPUtil::n_uint32_t read_checksum;
    read(streamer, read_checksum);

    if (level >= 2) {
      QDPIO::cout << "Read Map checksum: " << read_checksum << "  Current Position: " << streamer.currentPosition();
    }
    if( read_checksum != calc_checksum ) { 
      QDPIO::cout << "Mismatched Checksums: Expected: " << calc_checksum << " but read " << read_checksum << std::endl;
      QDP_abort(1);
    }

    if (level >= 2) {
      QDPIO::cout << " Map Checksum OK!" << std::endl;
    }
  }



  /*!
   * This is a utility function to sync the in memory offset map
   * with the one on the disk, and then close the file for writing.
   * Should not be called by user.
   */
  template<typename K, typename V>
  void
  MapObjectDisk<K,V>::closeWrite(void) 
  {
    switch(state) { 
    case MODIFIED:
    {
      if (level >= 2) {
	QDPIO::cout << "Beginning closeWrite: current position: " << streamer.currentPosition() << std::endl;
      }

      // Go to end of file
      streamer.seekEnd(0);

      // Take note of current position
      priv_pos_type_t metadata_start = convertToPrivate(streamer.currentPosition());
	
      if (level >= 2) {
	QDPIO::cout << "CloseWrite: Metadata starts at position: " << convertFromPrivate(metadata_start) << std::endl;
      }

      // Dump metadata
      writeMapBinary();
	
      // Rewind and Skip header 
      streamer.rewind();
      if (level >= 2) {
	QDPIO::cout << "Rewound file. Current Position: " << streamer.currentPosition() << std::endl;
      }
      writeSkipHeader();
      if (level >= 2) {
	QDPIO::cout << "Skipped Header. Current Position: " << streamer.currentPosition() << std::endl;
      }
      // write start position of metadata
      streamer.writeArray((const char *)&metadata_start,sizeof(priv_pos_type_t),1);
	
      if (level >= 2) {
	QDPIO::cout << "Wrote link to metadata. Current Position: " << streamer.currentPosition() << std::endl;
      }
	
      // skip to end and close
      streamer.seekEnd(0);
      streamer.flush();
	
      QDPIO::cout << "MapObjectDisk: Closed file " << filename<< " for write access" <<  std::endl;
    }
    break;
    default:
      errorState("MapObjectDisk: closeWrite() called in an invalid state");
      break;
    }
  }


} // namespace Chroma

#endif
