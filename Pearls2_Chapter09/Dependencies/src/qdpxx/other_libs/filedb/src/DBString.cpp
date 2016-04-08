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
 *     Implementation for String and Char* as Key and Data
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: DBString.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *   Revision 1.4  2008/08/12 17:37:36  chen
 *   Using new interface
 *
 *   Revision 1.3  2008/08/02 20:05:21  edwards
 *   Changed readObject and writeObject semantics to instead take input
 *   and output buffers from the caller instead of a private internal
 *   buffer within the object itself.
 *
 *   Revision 1.2  2008/06/02 20:10:33  chen
 *   First try to put into production
 *
 *   Revision 1.1.1.1  2008/05/29 19:31:39  chen
 *   Initial CVS import of FFDB
 *
 *
 */

#include <cstring>
#include <sstream>
#include "DBString.h"

using namespace std;

namespace FILEDB
{
  /**************************************************************************
   *      Implementation  StringKey Class                                   *
   **************************************************************************/
  StringKey::StringKey (void)
    :DBKey(), str_()
  {
    // empty
  }

  StringKey::StringKey (const StringKey& key)
    :DBKey (), str_ (key.str_)
  {
    // empty
  }

  StringKey::StringKey (const std::string& str)
    :DBKey (), str_ (str)
  {
    // empty
  }

  StringKey::StringKey (const char* buf)
    :DBKey(), str_ (buf)
  {
    // empty
  }

  StringKey&
  StringKey::operator = (const StringKey& key)
  {
    if (this != &key) {
      str_ = key.str_;
    }
    return *this;
  }
    

  StringKey&
  StringKey::operator = (const std::string& str)
  {
    str_ = str;

    return *this;
  }

  int
  StringKey::operator == (const StringKey& key2)
  {
    return (str_.compare (key2.str_) == 0);
  }

  int
  StringKey::operator != (const StringKey& key2)
  {
    return (str_.compare (key2.str_) != 0);
  }

  StringKey::~StringKey (void)
  {
    // empty
  }

  StringKey::operator std::string (void)
  {
    return str_;
  }

  // A string object is preceded by a magic number and an ID number
  void
  StringKey::writeObject (std::string& output) const throw (SerializeException)
  {
    ostringstream ostrm;
    unsigned short magicnum = hton_short(StringKey::MAGIC);
    unsigned short sid = hton_short(serialID());

    // binary write magic id and serial id
    try {
      ostrm.write ((const char *)&magicnum, sizeof(unsigned short));
      if (ostrm.bad()) {
	string msg("writing magic number error.");
	throw SerializeException ("StringKey", msg);
      }

      ostrm.write ((const char *)&sid, sizeof(unsigned short));
      if (ostrm.bad()) {
	string msg("writing serial id error.");
	throw SerializeException ("StringKey", msg);
      }

      // finally write out our string
      ostrm.write ((const char *)&(str_[0]), str_.length());
      if (ostrm.bad()) {
	string msg("writing internal string error.");
	throw SerializeException ("StringKey", msg);
      }
    }
    catch (ios_base::failure& e) {
      throw SerializeException ("StringKey", e.what());
    }
    
    output = ostrm.str();
  }

  // again we need to retrive magic and serial id for security reason
  void
  StringKey::readObject (const std::string& input) throw (SerializeException)
  {
    istringstream istrm (input);
    unsigned short magicnum, sid;

    // read magic number first
    try {
      istrm.read ((char *)&magicnum, sizeof(unsigned short));
      if (!istrm.good()) {
	string msg = string("Read magic number error.");
	throw SerializeException ("StringKey", msg);
      }

      // check magic number
      magicnum = ntoh_short (magicnum);
      if (magicnum != StringKey::MAGIC) {
	string msg = string("Magic number mismatch ") + itostr (magicnum) + string(" != expected ") + itostr (StringKey::MAGIC);
	throw SerializeException ("StringKey", msg);
      }
      
      // read serial id and check id
      istrm.read ((char *)&sid, sizeof(unsigned short)); 
      if (!istrm.good()) {
	string msg = string ("Read serial id error.");
	throw SerializeException ("StringKey", msg);     
      }
      
      // check serial id
      sid = ntoh_short (sid);
      if (sid != serialID()) {
	string msg = string ("Serial ID number mismatch ") + itostr (sid) + string (" != expected ") + itostr (serialID());
	throw SerializeException ("StringKey", msg);
      }

      // finally read out the real string
      str_ = input.substr (istrm.tellg());
    }
    catch (ios_base::failure& e) {
      string msg = string ("ReadObject Error: ") + e.what();
      throw SerializeException ("StringKey", msg);
    }
  }


  /**************************************************************************
   *      Implementation  userData  Class                                   *
   **************************************************************************/
  UserData::UserData (void)
    :DBData(), str_()
  {
    // empty
  }

  UserData::UserData (const UserData& key)
    :DBData(), str_ (key.str_)
  {
    // empty
  }

  UserData::UserData (const std::string& str)
    :DBData(), str_ (str)
  {
    // empty
  }

  UserData::UserData (const char* buf)
    :DBData(), str_ (buf)
  {
    // empty
  }

  UserData::UserData (const char* buf, unsigned int len)
    :DBData(), str_ (buf, len)
  {
    // empty
  }

  UserData&
  UserData::operator = (const UserData& key)
  {
    if (this != &key) {
      str_ = key.str_;
    }
    return *this;
  }
    

  UserData&
  UserData::operator = (const std::string& str)
  {
    str_ = str;

    return *this;
  }

  int
  UserData::operator == (const UserData& data2)
  {
    return (str_.compare (data2.str_) == 0);
  }

  int
  UserData::operator != (const UserData& data2)
  {
    return (str_.compare (data2.str_) != 0);
  }

  UserData::~UserData (void)
  {
    // empty
  }

  UserData::operator std::string (void)
  {
    return str_;
  }

  void
  UserData::writeObject (std::string& output) const throw (SerializeException)
  {
    ostringstream ostrm;
    unsigned short magicnum = hton_short(UserData::MAGIC);
    unsigned short sid = hton_short(serialID());

    // binary write magic id and serial id
    try {
      ostrm.write ((const char *)&magicnum, sizeof(unsigned short));
      if (ostrm.bad()) {
	string msg("writing magic number error.");
	throw SerializeException ("UserData", msg);
      }

      ostrm.write ((const char *)&sid, sizeof(unsigned short));
      if (ostrm.bad()) {
	string msg("writing serial id error.");
	throw SerializeException ("UserData", msg);
      }

      // finally write out our string
      ostrm.write ((const char *)&(str_[0]), str_.length());
      if (ostrm.bad()) {
	string msg("writing internal string error.");
	throw SerializeException ("UserData", msg);
      }
    }
    catch (ios_base::failure& e) {
      throw SerializeException ("UserData", e.what());
    }
    
    output = ostrm.str();
  }

  void
  UserData::readObject (const std::string& input) throw (SerializeException)
  {
    istringstream istrm (input);
    unsigned short magicnum, sid;

    // read magic number first
    try {
      istrm.read ((char *)&magicnum, sizeof(unsigned short));
      if (!istrm.good()) {
	string msg = string("Read magic number error.");
	throw SerializeException ("UserData", msg);
      }
      // check magic number
      magicnum = ntoh_short (magicnum);
      if (magicnum != UserData::MAGIC) {
	string msg = string("Magic number mismatch ") + itostr (magicnum) + string(" != expected ") + itostr (UserData::MAGIC);
	throw SerializeException ("UserData", msg);
      }
      
      // read serial id and check id
      istrm.read ((char *)&sid, sizeof(unsigned short)); 
      if (!istrm.good()) {
	string msg = string ("Read serial id error.");
	throw SerializeException ("UserData", msg);     
      }
      
      // check serial id
      sid = ntoh_short (sid);
      if (sid != serialID()) {
	string msg = string ("Serial ID number mismatch ") + itostr (sid) + string (" != expected ") + itostr (serialID());
	throw SerializeException ("UserData", msg);
      }

      // finally read out the real string
      str_ = input.substr (istrm.tellg());
    }
    catch (ios_base::failure& e) {
      string msg = string ("ReadObject Error: ") + e.what();
      throw SerializeException ("UserData", msg);
    }
  }


}
