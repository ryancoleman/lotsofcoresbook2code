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
 *   $Log: SConfCreateTest.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/03/02 23:27:26  chen
 *   Test DBMerge Code
 *
 *
 *
 */
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include "DBString.h"
#include "VFloatData.h"
#include "ConfDataStoreDB.h"

#define USER_STRING "Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCDHello threre I am the new user string for QCD  Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD"

using namespace std;
using namespace FILEDB;

int
main (int argc, char** argv)
{
  if (argc < 8) {
    cerr << "Usage: " << argv[0] << " pagesize numpages rearrange(0|1) dbasename keyfile vectorlength confignum" << endl;
    return -1;
  }
  
  unsigned int pagesize = atoi (argv[1]);
  unsigned int nbuckets = atoi (argv[2]);
  int rearrange = atoi (argv[3]);
  std::string dbase (argv[4]);
  std::string keyfile(argv[5]);
  unsigned int vectorlen = atoi(argv[6]);
  unsigned int confnum = atoi(argv[7]);

  // open string file stream
  ifstream sf(argv[5]);
  if (sf.bad()) {
    cerr << "Cannot open file " << keyfile << endl;
    return -1;
  }
  
  // Open database
  ConfDataStoreDB<StringKey, VFloatData> dbtest;

  dbtest.setCacheSize (100000000);
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


  VFloatData vfdata(vectorlen, (float)confnum);
  int i = 0;
  while (sf.good()) {
    string t1;

    sf >> t1;
    StringKey key(t1);

    if (t1.length() == 0)
      break;

    if (dbtest.insert(key, vfdata) != 0) {
      cerr << "Insert data error " << endl;
      sf.close ();
      dbtest.close ();
      return -1;
    }
    i++;
  }
  sf.close ();
  dbtest.close ();

  return 0;
}
