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
 *     Simple Test Program to Test ConfDataStoreDB Class
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: ConfReadTest.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#include <fstream>
#include <iostream>
#include "DBString.h"
#include "ConfDataStoreDB.h"

using namespace std;
using namespace FILEDB;

int
main (int argc, char** argv)
{
  if (argc < 5) {
    cerr << "Usage: " << argv[0] << " cachesize rearrange(0|1) dbasename stringfile" << endl;
    return -1;
  }
  
  unsigned int cachesize = atoi (argv[1]);
  int rearrange = atoi (argv[2]);
  std::string dbase (argv[3]);
  std::string strfile(argv[4]);

  // open string file stream
  ifstream sf(argv[4]);
  if (sf.bad()) {
    cerr << "Cannot open file " << strfile << endl;
    return -1;
  }
  
  // Open database
  ConfDataStoreDB<StringKey, UserData> dbtest;

  dbtest.setCacheSize (cachesize);
  if  (rearrange)
    dbtest.enablePageMove ();

  if (dbtest.open (dbase, O_RDONLY, 0664) != 0) {
    cerr << "cannot open database " << dbase << endl;
    sf.close ();
    return -1;
  }

  // get user data
  std::string userdata;
  if (dbtest.getUserdata (userdata) != 0) {
    cerr << " Cannot get user data" << endl;
    sf.close ();
    dbtest.close ();
    return -1;
  }
  cerr << userdata << endl;

  while (sf.good()) {
    string t1, t2;

    sf >> t1;
    sf >> t2;

    StringKey key(t1);
    UserData  data;

    if (dbtest.get(key, data) != 0) {
      cerr << "Cannot find data  " << endl;
      sf.close ();
      dbtest.close ();
      return -1;
    }
    else {
      if ((string)data != t2) {
	cerr << "Retreive data error" << endl;
	sf.close ();
	dbtest.close ();
	return -1;
      }
    }
  }
  sf.close ();
  dbtest.close ();

  return 0;
}
