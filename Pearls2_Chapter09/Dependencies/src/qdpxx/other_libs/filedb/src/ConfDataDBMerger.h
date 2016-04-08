// -*- C++ -*-
/*----------------------------------------------------------------------------
 * Copyright (c) 2007      Jefferson Science Associates, LLC               
 *                         Under U.S. DOE Contract No. DE-AC05-06OR23177  
 *                                                                      
 *                         Thomas Jefferson National Accelerator Facility
 *
 *                         Jefferson Lab 
 *                         Scientific Computing Group,
 *                         12000 Jefferson Ave.,      
 *                         Newport News, VA 23606 
 *----------------------------------------------------------------------------
 *  
 * Description:
 *     Process all Configurations stored in databases and insert all information
 *     into a secondary database storage
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *   $Log: ConfDataDBMerger.h,v $
 *   Revision 1.5  2009-09-22 14:56:44  edwards
 *   Commented out some debugging.
 *
 *   Revision 1.4  2009/03/04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.3  2009/03/02 23:27:26  chen
 *   Test DBMerge Code
 *
 *   Revision 1.2  2009/02/28 21:03:17  edwards
 *   Rerranged some code in merge. Declare objects when they are used.
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#ifndef _FILEDB_CONF_DATA_DB_MERGER_H
#define _FILEDB_CONF_DATA_DB_MERGER_H

#include <vector>
#include "ConfigInfo.h"
#include "ConfDataStoreDB.h"
#include "AllConfStoreDB.h"
#include "HashFunc.h"

#include <unordered_map>

namespace FILEDB
{
  template <typename K, typename D>
  class ConfDataDBMerger
  {
  private:
    /**
     * In memory hash table for key and data
     * default == operator and default hash function are used
     */
    typedef std::unordered_map< std::string, std::vector< std::string > > MemMap_t;


    /**
     * Database holding all Keys and Data
     */
    AllConfStoreDB< K , D >* dbh_;


    /**
     * Internal all configuration information
     */
    std::vector<ConfigInfo> configs_;


    /**
     * maximum memory we are going to use
     */
    unsigned int max_memory_;


    /**
     * Populate initial keys and data from the first configuration
     *
     * @param info configuration information for this set of keys and data
     * @param keys all keys
     * @param data all data
     * @param chunksize each iteration chunk size through the vectors
     * @param chunkloop current chunk loop number
     * @param numconfigs total number of configs
     * @param bag hash map holding chunksize of keys and data for all configs
     */
    static void initBag (const ConfigInfo& info,
			 std::vector< std::string >& keys,
			 std::vector< std::string >& data,
			 int chunksize,
			 int chunkloop,
			 int numconfigs,
			 MemMap_t& bag)
    {
      unsigned int i, ielem, felem;

      // get initial element
      ielem = chunkloop * chunksize;

      // get final element index
      felem = (chunkloop + 1) * chunksize;

      if (felem > keys.size())
	felem = keys.size();
      
      for (i = ielem; i < felem; i++) {
#if 0
	std::cerr << "Element " << i << std::endl;
#endif
	std::vector< std::string > thr;
	
	thr.resize (numconfigs);

	if (bag.find (keys[i])  != bag.end()) {
	  std::cerr << "Cannot be right, this key should not be here" 
		    << std::endl;
	}
	bag.insert (std::make_pair (keys[i], thr));
	
	// insert data into the vector
	bag.find (keys[i])->second[info.index()] = data[i];
#if 0
	std::cerr << "Init index " << info.index() << std::endl;
#endif
      }    
    }


    /**
     * Populate keys and data from a configuration
     *
     * @param c index points to which configuration we are working on
     * @param keys all keys
     * @param chunksize each iteration chunk size through the vectors
     * @param chunkloop current chunk loop number
     * @param numconfigs total number of configs
     * @param bag hash map holding chunksize of keys and data for all configs
     */
    static void populateBag (const ConfigInfo& info,
			     std::vector< std::string >& keys,
			     int chunksize,
			     int chunkloop,
			     int numconfigs,
			     unsigned int dbasecachesize,
			     MemMap_t& bag)
    {
      unsigned int i, ielem, felem;

      // check whether this config info has valid urlname
      if (info.urlname().length() > 0) {
	ConfDataStoreDB< K, D >* sdb = 0;
	sdb = new ConfDataStoreDB< K, D > ();
	sdb->setCacheSize (dbasecachesize);
	if (sdb->open  (info.urlname(), O_RDONLY, 0400) != 0) {
	  std::cerr << __FILE__ << " cannot open individual database " << info.urlname() << std::endl;
	  ::exit (123);
	}

	// get initial position of interested keys 
	ielem = chunkloop * chunksize;
	
	// get final index of intertested keys
	felem = (chunkloop + 1) * chunksize;
	if (felem > keys.size())
	  felem = keys.size();

#if 0
	std::cerr << "initial item = " << ielem << " final item = " << felem << std::endl;
#endif
	// retrieve data for each key in this chunk
	for (i = ielem; i < felem; i++) {
	  std::string  tdata;
	  
	  if (sdb->getBinary(keys[i], tdata) != 0) {
	    std::cerr << "This configuration " << info.configNumber() 
		      << " does not have this key " << i << std::endl;
	    ::exit (143);
	  }
	
	  if (bag.find (keys[i])  == bag.end() ) {
	    std::cerr << "Cannot be right, this key should be in the bag" << std::endl;
	    ::exit (143);
	  }
	  // insert data into the ensemble
#if 0
	  std::cerr << "Update index " << info.index() << std::endl;
#endif
	  bag.find (keys[i])->second[info.index()] = tdata;
	}

	// free memory
	delete sdb;
      }
      else
	std::cerr << "Skip this configuration " << info << std::endl;
    }



    /**
     * Insert a memory hash map loaded with keys and data into 
     * into final big database
     */
    void insertBag (MemMap_t& bag)
    {
      int num = 0;
      // iterator is a template also
      MemMap_t::iterator ite;
      
      for (ite = bag.begin(); ite != bag.end (); ite++) {
	// make compiler happy
	std::string tkey((*ite).first);

	dbh_->insertBinaryData (tkey, (*ite).second);
	num++;
      }
    }


    /**
     * Calculate number of elements each iteration to use.
     * This number of elements will be used to accumulate keys and data
     * stored in a hash map in memory (limited by the following memory size)
     * for all configurations.
     * 
     * This routine also returns cache size for each database when the
     * database is opened
     *
     * @param key a repentative key value usesd to calculate size
     * @param data a repentative data value usesd to calculate size
     * @param numelems tptal of number of (key, data) pairs
     * @param numconfs number of configurations
     * @param cachesize database cache size when a small database is opened
     *
     */
    static int iterateSize (std::string& key,
			    std::string& data,
			    unsigned int numelems,
			    int numconfs,
			    unsigned int max_memory,
			    unsigned int& cachesize)
    {
      int chunk;

      // calculate size for one key and its data for all configuration
      unsigned int asize = (key.length() + data.length()) * numconfs;

      // number of elements each iteration to use
      chunk = max_memory/asize > numelems ? numelems : max_memory/asize;

      // calculate cachesize for each single database
      cachesize = chunk * (key.length() + data.length());
      
      return chunk;
    }


    /**
     * Update configuration data information to configuration information
     * manager database
     */
    void updateConfigInfo (void)
    {
      std::vector<ConfigInfo>::iterator ite;
      std::string data;
      ConfigInfo cinfo;

      for (ite = configs_.begin(); ite != configs_.end(); ite++) {
	if ((*ite).urlname().length() != 0) {
	  cinfo = *ite;
	  cinfo.insert(1);
	  cinfo.modifiedTime ((int)current_time());
	  dbh_->setConfigInfo (cinfo);
	}
      }
    }

  public:
    /**
     * Constructor
     * @param configs a vector of configuration information
     * @param maximum_memory how much memory (in bytes) we are going to use 
     * to merge data in memory before writing to the combined database
     *
     */
    ConfDataDBMerger (const std::vector<ConfigInfo>& configs,
		      unsigned int max_memory)
      :dbh_ (0), configs_ (configs), max_memory_ (max_memory)
    {
      // empty
    }

    /**
     * Constructor
     * Take a vector of database names, configuration numbers
     */
    ConfDataDBMerger (const std::vector<int>& config_numbers,
		      const std::vector<std::string>& dbase_names,
		      unsigned int max_memory)
      :dbh_(0), configs_ (), max_memory_(max_memory)
    {
      ConfigInfo info;
      int i;
    
      // all these vectors should have the same size
      if (dbase_names.size () == config_numbers.size ()) {
	std::vector<std::string>::const_iterator site;
	std::vector<int>::const_iterator cite;

	cite = config_numbers.begin();
	site = dbase_names.begin();
	i = 0;
	while (site != dbase_names.end() && cite != config_numbers.end()) {
	  info.configNumber (*cite);
	  info.index (i);
	  info.type (0);
	  if ((*site) != "N/A")
	    info.urlname (*site);
	  configs_.push_back (info);
	  
	  site++;
	  cite++;
	  i++;
	}
      }
      else {
	std::cerr << "Incorrect parameter to construct ConfDataManager: vector size mismatch" << std::endl;
	::exit (143);
      }
    }

    /**
     * Destructor
     */
    ~ConfDataDBMerger (void)
    {
      // dbh_ and informan_ are managed by others
    }


    /**
     * Set database handle
     */
    void setDatabaseHandle (AllConfStoreDB< K, D >* dbh)
    {
      dbh_ = dbh;
    }


    
    /**
     * Merge all data from these databases into one single database
     * which should be empty in the very beginning
     */
    void merge (void)
    {
      // find everything in the database that have valid database name in it
      unsigned int i;
      for (i = 0; i < configs_.size(); i++) {
	if (configs_[i].urlname().length () != 0) 
	  break;
      }
      if (i >= configs_.size()) {
	std::cerr << "Fatal: no database name provided for any configuration" << std::endl;
	::exit (1);
      }

      // open the first valid configuration
      ConfDataStoreDB< K, D >* fsconfdb = 0;
      fsconfdb = new ConfDataStoreDB< K, D > ();
      if (fsconfdb->open (configs_[i].urlname(), O_RDONLY, 0400) != 0) {
	std::cerr << __FILE__ << " cannot open individual database " << configs_[i].urlname() << std::endl;
	::exit (1);
      }
    
      std::vector< std::string > allkeys;
      std::vector< std::string > data0;

      // get time of the starting of all merge
      double it = current_time ();

      // get keys and data in binary form
      fsconfdb->binaryKeysAndData (allkeys, data0);

      if (allkeys.size() == 0 || data0.size() == 0)
      {
	std::cerr << __func__ << ": empty key or data size in db " << configs_[i].urlname() << std::endl;
	::exit (1);
      }

#if 0
      std::cerr << "Number of keys and data = " << allkeys.size() << std::endl;
      std::cerr << "Data size = " << data0.size() << std::endl;
#endif

      // calculate the number of elements we use in the above two vectors each
      // time we accumulate all configurations into memory before we
      // write out to the final database.
      unsigned int dcachesize;
      
      int chunk = iterateSize (allkeys[0], data0[0], 
			       allkeys.size(), configs_.size(),
			       max_memory_,
			       dcachesize);

#if 1
      std::cerr << "Keep " << chunk << " in memory for " << configs_.size() 
		<< " of configurations" << std::endl;
#endif
      // calculate how many loops are we going to run.
      // Take care the case when size is multiple of chunk
      unsigned int numloops = 
	(allkeys.size () == allkeys.size()/chunk * chunk) ? allkeys.size()/chunk : allkeys.size()/chunk + 1;

      MemMap_t bag;
      for (int loop = 0; loop < numloops; loop++) {
	// clear the bag
	bag.clear ();

	// populate bag using the data from the first database
	initBag (configs_[i], allkeys, data0, chunk, loop, configs_.size(), bag);

	// retrieve all information from each configuration
	for (int c = 0; c < configs_.size(); c++) {
	  // populate bag for each configurations
	  if (c != i)
	    populateBag (configs_[c], allkeys, chunk, loop, 
			 configs_.size(), 4*dcachesize, bag);
	}
	// we write this bag into the database
	this->insertBag (bag);

      }
      double ft = current_time ();
      
      std::cerr << "Writing " << configs_.size() << " configurations with " 
		<< allkeys.size() << " data pairs takes " << ft - it 
		<< " seconds " << std::endl;
      
      // Update Cinfiguration Information
      updateConfigInfo ();

      // free memory
      delete fsconfdb;
    }

  };
}
#endif
