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
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *
 * Revision History:
 *   $Log: ConfDataDBMerger.cpp,v $
 *   Revision 1.3  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.2  2009/03/02 23:27:26  chen
 *   Test DBMerge Code
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *
 *
 */
#ifdef _FILEDB_COMPILE_TEST

#include <DBString.h>
#include "AllConfStoreDB.h"
#include "ConfDataDBMerger.h"

using namespace std;

namespace FILEDB
{
  static void test_merge_code (void)
  {
    vector<ConfigInfo> configs;
    vector<string> dnames;
    vector<int> cfigs;
    
    ConfDataDBMerger<StringKey, UserData> m1 (configs, 100000000);
    ConfDataDBMerger<StringKey, UserData> m2 (cfigs, dnames, 100000000);

    AllConfStoreDB<StringKey, UserData> dbh;

    m2.setDatabaseHandle (&dbh);
    m1.setDatabaseHandle (&dbh);

    m2.merge ();
  }
}
#endif

