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
 *   $Log: AllConfReadTest.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
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

using namespace std;
using namespace FILEDB;

#define NUM_KEYS 1000

int
main (int argc, char** argv)
{
  char keystr[80], datastr[128];
  if (argc < 5) {
    cerr << "Usage: " << argv[0] << " cachesize numconfigs rearrange(0|1) dbasename " << endl;
    return -1;
  }
  
  unsigned int cachesize = atoi (argv[1]);
  unsigned int numconfigs = atoi (argv[2]);
  int rearrange = atoi (argv[3]);
  std::string dbase (argv[4]);

  // Open database
  AllConfStoreDB<StringKey, UserData> dbtest;

  dbtest.setCacheSize (cachesize);
  if  (rearrange)
    dbtest.enablePageMove ();

  if (dbtest.open (dbase, O_RDONLY, 0400) != 0) {
    cerr << "cannot open database " << dbase << endl;
    return -1;
  }

  std::string userdata;
  if (dbtest.getUserdata (userdata) != 0) {
    cerr << "Cannot get user string" << endl;
    dbtest.close ();
    return -1;
  }

  cerr << userdata << endl;

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

    vector<UserData> tdatav;
    for (int k = 0; k < numconfigs; k++) {
      UserData tmp(datav[k], 128);
      tdatav.push_back (tmp);
    }

    vector<UserData> rdatav;
    if (dbtest.get (key, rdatav) != 0) {
      cerr << "Cannot get vector data at " << i << endl;
      dbtest.close ();
      return -1;
    }
    else {
      for (int k = 0; k < numconfigs; k++) {
	if (tdatav[k] != rdatav[k]) {
	  cerr << "Retrieved data is wrong at " << k << " element loop " << i << endl;
	  dbtest.close ();
	  return -1;	
	}
      }
    }
  }

  dbtest.close ();

  return 0;
}
