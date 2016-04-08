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
 *     
 *     Simple Class to instantiate template class to assist on checking 
 *     code errors
 *     
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: AllConfStoreDB.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:47  chen
 *   initial import
 *
 *
 *
 */
#ifdef _FILEDB_COMPILE_TEST

#include "DBString.h"
#include "AllConfStoreDB.h"

using namespace std;

namespace FILEDB
{
  static void init_allconfstoredb_test (void)
  {
    int i;
    std::vector<int> allconfigs;

    for (i = 0; i < 1000; i++)
      allconfigs[i] = i;

    AllConfStoreDB<StringKey, UserData> dbtest (allconfigs);

    dbtest.setCacheSize (10000000);
    dbtest.setPageSize  (8192);
    dbtest.setNumberBuckets (1000);
    dbtest.enablePageMove ();
    dbtest.setMaxUserInfoLen (5000);
    dbtest.setMaxNumberConfigs (1000);
    
    dbtest.open ("hashtest", O_RDWR, 0644);

    std::vector<UserData> data0;
    StringKey key;

    dbtest.insert (key, data0);

    std::string key1("mytest");
    std::vector<std::string> data1;
    dbtest.insertBinaryData (key1, data1);

    std::vector<UserData> data2;    
    dbtest.get (key, data2);

    std::string key2("mykey2"); 
    std::vector<std::string> data3;   
    
    dbtest.getBinaryData (key2, data3);

    UserData data4("hello");
    dbtest.update (key, data4, 100, 90);

    dbtest.exist (key);
    
    dbtest.flush ();
    dbtest.close ();


    
  }
}
#endif
