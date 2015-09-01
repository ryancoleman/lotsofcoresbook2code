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
 *   $Log: AllConfCreateTest.cpp,v $
 *   Revision 1.3  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.2  2009/03/01 01:03:56  chen
 *   Add missing O_CREAT flag
 *
 *   Revision 1.1  2009/02/20 20:44:47  chen
 *   initial import
 *
 *
 *
 */
#include <cstdio>
#include <fstream>
#include <iostream>
#include "DBString.h"
#include "AllConfStoreDB.h"

#define USER_STRING "Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCDHello threre I am the new user string for QCD  Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello threre I am the new user string for QCD Hello There I am Seasoned User"


using namespace std;
using namespace FILEDB;

#define NUM_KEYS 1000

int
main (int argc, char** argv)
{
  char keystr[80], datastr[128];
  if (argc < 6) {
    cerr << "Usage: " << argv[0] << " pagesize numpages numconfigs rearrange(0|1) dbasename " << endl;
    return -1;
  }
  
  unsigned int pagesize = atoi (argv[1]);
  unsigned int nbuckets = atoi (argv[2]);
  unsigned int numconfigs = atoi (argv[3]);
  int rearrange = atoi (argv[4]);
  std::string dbase (argv[5]);

  vector<int> configs;

  for (int i = 0; i < numconfigs; i++)
    configs.push_back(i);
  
  // Open database
  AllConfStoreDB<StringKey, UserData> dbtest(configs);

  dbtest.setCacheSize (100000000);
  dbtest.setPageSize (pagesize);
  dbtest.setNumberBuckets (nbuckets);
  if  (rearrange)
    dbtest.enablePageMove ();
  dbtest.setMaxUserInfoLen (5000);
  dbtest.setMaxNumberConfigs (numconfigs);  

  if (dbtest.open (dbase, O_RDWR | O_TRUNC | O_CREAT, 0664) != 0) {
    cerr << "cannot open database " << dbase << endl;
    return -1;
  }

  std::string userdata(USER_STRING);
  if (dbtest.insertUserdata (userdata) != 0) {
    cerr << "Cannot insert user string" << endl;
    dbtest.close ();
    return -1;
  }

  for (int i = 0; i < NUM_KEYS; i++) {
    ::sprintf (keystr, "Key Test Loop %d", i);
    StringKey key(keystr);
    char** datav = new char*[numconfigs];
    for (int k = 0; k < numconfigs; k++) {
      for (int m = 0; m < 128; m++)
	datastr[m] = (i + k + m) % 127;

      datav[k] = new char[128];
      memcpy (datav[k], datastr, 128);
    }
    
    vector<UserData> tdata;
    for (int k = 0; k < numconfigs; k++) {
      UserData tmp(datav[k], 128);
      tdata.push_back (tmp);
    }

    if (dbtest.insert (key, tdata) != 0) {
      cerr << "Cannot insert vector data at " << i << endl;
      dbtest.close ();
      return -1;
    }
    delete []datav;
  }

  dbtest.close ();

  return 0;
}
