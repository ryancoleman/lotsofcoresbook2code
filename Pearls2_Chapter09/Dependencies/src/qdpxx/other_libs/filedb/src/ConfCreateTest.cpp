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
 *   $Log: ConfCreateTest.cpp,v $
 *   Revision 1.4  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.3  2009/02/25 18:16:07  edwards
 *   Changed userdata to go through a string.
 *
 *   Revision 1.2  2009/02/24 21:43:15  chen
 *   Add O_CREAT
 *
 *   Revision 1.1  2009/02/20 20:44:47  chen
 *   initial import
 *
 *
 *
 */
#include <fstream>
#include <iostream>
#include "DBString.h"
#include "ConfDataStoreDB.h"

#define USER_STRING "Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCDHello threre I am the new user string for QCD  Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD"

using namespace std;
using namespace FILEDB;

int
main (int argc, char** argv)
{
  if (argc < 6) {
    cerr << "Usage: " << argv[0] << " pagesize numpages rearrange(0|1) dbasename stringfile" << endl;
    return -1;
  }
  
  unsigned int pagesize = atoi (argv[1]);
  unsigned int nbuckets = atoi (argv[2]);
  int rearrange = atoi (argv[3]);
  std::string dbase (argv[4]);
  std::string strfile(argv[5]);

  // open string file stream
  ifstream sf(argv[5]);
  if (sf.bad()) {
    cerr << "Cannot open file " << strfile << endl;
    return -1;
  }
  
  // Open database
  ConfDataStoreDB<StringKey, UserData> dbtest;

  // 8 GB cache
  dbtest.setCacheSizeMB (8*1024);
  dbtest.setPageSize (pagesize);
  dbtest.setNumberBuckets (nbuckets);
  if  (rearrange)
    dbtest.enablePageMove ();
  std::string userdata(USER_STRING);
  dbtest.setMaxUserInfoLen (userdata.size());
  dbtest.setMaxNumberConfigs (2);

  if (dbtest.open (dbase, O_RDWR | O_TRUNC | O_CREAT, 0664) != 0) {
    cerr << "cannot open database " << dbase << endl;
    sf.close ();
    return -1;
  }

  if (dbtest.insertUserdata (userdata) != 0) {
    cerr << "Cannot insert user string" << endl;
    sf.close ();
    dbtest.close ();
    return -1;
  }

  while (sf.good()) {
    string t1, t2;

    sf >> t1;
    sf >> t2;

    StringKey key(t1);
    UserData  data(t2);

    if (dbtest.insert(key, data) != 0) {
      cerr << "Insert data error " << endl;
      sf.close ();
      dbtest.close ();
      return -1;
    }
  }
  sf.close ();
  dbtest.close ();

  return 0;
}
