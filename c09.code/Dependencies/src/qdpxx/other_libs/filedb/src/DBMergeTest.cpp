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
 *     Testing Merge databases each of which is a configuration 
 *     into a large database
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *    
 *
 * Revision History:
 *   $Log: DBMergeTest.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/03/02 23:27:26  chen
 *   Test DBMerge Code
 *
 *
 *
 */
#include <cstdio>
#include <cctype>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "DBString.h"
#include "VFloatData.h"
#include "AllConfStoreDB.h"
#include "ConfDataDBMerger.h"

using namespace std;
using namespace FILEDB;

static void
usage (int argc, char** argv)
{
  cerr << "Usage: " << argv[0] << " -f dbasename -l list -m memorysize" << endl;
}

static void
help (int argc, char** argv)
{
  usage (argc, argv);
  cerr << "The list is a file containing two columns: confignum dbasename" << endl;
  cerr << "where confignum are configuration numbers, and " << endl;
  cerr << "dbasename are database name for the configuration number." << endl;
  cerr << "If the dbasename is empty, this configuration has not been finished" << endl;
  cerr << "memorysize is the maximum memory (bytes) we are going to use (<= 4GB)" << endl;
}

static int
parse_options (int argc, char** argv, 
	       string& dbase, string& listfile, unsigned int& msize)
{
  int cnt = 0;

  int i = 0;
  while (i < argc) {
    if (::strcmp (argv[i], "-f") == 0) {
      if (i == argc - 1) {
	// this is the last one
	cerr << "-f has to be followed by a database name" << endl;
	return -1;
      }
      if (argv[i + 1][0] == '-') {
	// another option
	cerr << "-f cannot be followed by another option" << endl;
	return -1;
      }
      dbase = argv[++i];
      cnt++;
    }
    else if (::strcmp (argv[i], "-l") == 0) {
      if (i == argc - 1) {
	// this is the last one
	cerr << "-l has to be followed by a list file name" << endl;
	return -1;
      }
      if (argv[i + 1][0] == '-') {
	// another option
	cerr << "-l cannot be followed by another option" << endl;
	return -1;
      }
      listfile = argv[++i];
      cnt++;
    }
    else if (::strcmp (argv[i], "-m") == 0) {
      if (i == argc - 1) {
	// this is the last one
	cerr << "-m has to be followed memory size" << endl;
	return -1;
      }
      if (argv[i + 1][0] == '-') {
	// another option
	cerr << "-m cannot be followed by another option" << endl;
	return -1;
      }   
      // get string stream from this line
      istringstream istr (argv[++i]);
      if (istr >> msize)
	cnt++;
      else {
	cerr << "-m has to be followed by a unsigned integer" << endl;
	return -1;
      }
    }
    else if (::strcmp (argv[i], "-h") == 0) {
      help (argc, argv);
      exit (1);
    }
    else
      i++;
  }
  if (cnt != 3) {
    return -1;
  }
  return 0;
}

/**
 * Parse list file and populate configuration information
 */
int
config_infos (const string& lfile,
	      vector<int>& cfgs, vector<string>& sdbs)
{
  ifstream ifs (lfile.c_str());

  if (!ifs) {
    cerr << "Cannot open list file " << lfile << endl;
    return -1;
  }

  int cfg;
  string dname;
  char line[256];
  while (ifs.good()) {
    // get whole line
    ifs.getline (line, sizeof(line));
    
    // get string stream from this line
    istringstream istr (line);
    if (istr >> cfg) {
      if (istr >> dname) {
	sdbs.push_back (dname);
	cfgs.push_back (cfg);
      }
      else {
	cfgs.push_back (cfg);
	sdbs.push_back ("N/A");
      }
    }
  }
  ifs.close ();

  return 0;
}

int
main (int argc, char** argv)
{
  string dbase, lfile;
  unsigned int mem_size;

  if (parse_options (argc, argv, dbase, lfile, mem_size) != 0) {
    usage (argc, argv);
    return -1;
  }
  
  // construct configuration list and database list
  std::vector<int> cfgs;
  std::vector<string> sdbs;

  if (config_infos (lfile, cfgs, sdbs) != 0) 
    return -1;

  AllConfStoreDB<StringKey, VFloatData> database (cfgs, sdbs);
  database.setPageSize (8192);
  database.setNumberBuckets(64);
  database.enablePageMove ();
  
  if (database.open (dbase, O_RDWR|O_CREAT|O_TRUNC, 0666) != 0) {
    cerr << "Cannot open combined database " << dbase << endl;
    return -1;
  }

  ConfDataDBMerger<StringKey, VFloatData> merger (cfgs, sdbs, mem_size);

  // set database handle
  merger.setDatabaseHandle (&database);
    
  // let us start merging
  merger.merge ();

  // update configuration information is done inside merge code

  // Display all data and keys
  std::vector<StringKey> keys;
  std::vector< std::vector<VFloatData> > data;

  database.keysAndData (keys, data);

  database.close ();

  // Display Keys and data
  for (int i = 0; i < keys.size(); i++) {
    cerr << "Key " << i << endl;
    cerr << (string)keys[i] << endl;
    cerr << "value " << i << endl;

    for (int k = 0; k < data[i].size(); k++) {
      cerr << "config " << k << endl;
      for (int m = 0; m < ((vector<float>)(data[i][k])).size(); m++) 
	cerr << ((vector<float>)(data[i][k]))[m];
      cerr << endl;
    }
  }	     	     

  return 0;
}

