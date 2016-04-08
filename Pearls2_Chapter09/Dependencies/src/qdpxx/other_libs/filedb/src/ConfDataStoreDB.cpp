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
 *   $Log: ConfDataStoreDB.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#ifdef _FILEDB_COMPILE_TEST

#include "DBString.h"
#include "ConfDataStoreDB.h"

using namespace std;

namespace FILEDB
{
  static void init_confdatastoredb_test (void)
  {
    ConfDataStoreDB<StringKey, UserData> dbtest;

    dbtest.setCacheSize (10000000);
    dbtest.setPageSize  (8192);
    dbtest.setNumberBuckets (100);
    dbtest.enablePageMove ();
    dbtest.setMaxUserInfoLen (5000);
    dbtest.setMaxNumberConfigs (1000);
    
    dbtest.open ("hashtest", O_RDWR, 0644);

    std::string key("mykey");
    std::string data("mydata");
    dbtest.insertBinary (key, data);

    std::string data1;
    dbtest.getBinary (key, data1);

    StringKey key2("mykey2");
    UserData  data2 ("mydata2");
    
    dbtest.insert (key2, data2);

    UserData data3;
    dbtest.get (key2, data3);
    
    dbtest.flush ();
    dbtest.close ();


    
  }
}
#endif

