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
 *     Implementation for Vector Float as database data Class for test
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: VFloatData.cpp,v $
 *   Revision 1.3  2009-04-21 18:52:24  chen
 *   Minor change
 *
 *   Revision 1.2  2009/03/04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/03/02 23:27:26  chen
 *   Test DBMerge Code
 *
 *
 */
#include <stdlib.h>
#include <sstream>
#include "VFloatData.h"

using namespace std;

namespace FILEDB
{
  /**
   * Constructor
   */
  VFloatData::VFloatData (void)
    :DBData(), data_()
  {

  }

  VFloatData::VFloatData (unsigned int size, float defval)
    :DBData(), data_ (size)
  {
    for (unsigned int i = 0; i < size; i++)
      data_[i] = defval;
  }
   
  VFloatData::VFloatData (const VFloatData& d)
    :DBData(), data_(d.data_)
  {
    // empty
  }

  /**
   * Assignment Operator
   */
  VFloatData &
  VFloatData::operator = (const VFloatData& d)
  {
    if (this != &d) {
      data_ = d.data_;
    }
    return *this;
  }

  /**
   * Comparison operator
   */
  int
  VFloatData::operator == (const VFloatData& d) 
  {
    if (data_ == d.data_)
      return 1;

    return 0;
  }

  int
  VFloatData::operator != (const VFloatData& d) 
  {
    if (data_ == d.data_)
      return 0;

    return 1;
  }



  /**
   * Convertion operator
   */
  VFloatData::operator const std::vector<float>& (void)
  {
    return data_;
  }

  /**
   * Return number of elements in the vector
   */
  const int
  VFloatData::numberOfElements (void) const
  {
    return data_.size();
  }

  /**
   * Convert this object to a buffer
   *
   * @param buffer application allocated memory buffer to hold the flat object
   * @param len the refrence to the initial size of the buffer. When the
   * routine is returned, the len contains the object actual length
   */
  void 
  VFloatData::writeObject (std::string& output) const throw (SerializeException)
  {
    float* vcfbuf;
    unsigned short id;
    unsigned short pad;
    unsigned int   num;

    num = hton_int(data_.size ());
    id = hton_short(serialID());
    vcfbuf = new float[data_.size()];

    // assign complex data values to the array
    for (int i = 0; i < data_.size(); i++) 
      vcfbuf[i] = hton_int((int)data_[i]);

    // output stream
    ostringstream ostrm;
    try {
      ostrm.write ((const char *)&id, sizeof(unsigned short));
      ostrm.write ((const char *)&pad, sizeof(unsigned short));
      ostrm.write ((const char *)&num, sizeof(unsigned int));
      ostrm.write ((const char *)vcfbuf, data_.size() * sizeof(float));
      if (ostrm.bad()) 
	throw SerializeException ("VFloatData",
				  "writing vector complex array error");
    }
    catch (ios_base::failure& e) {
      string msg = string ("writing vector complex array error: ") + e.what();
      throw SerializeException ("VFloatData", msg);
    }

    output = ostrm.str();

    delete []vcfbuf;
  }

  /**
   * Convert a buffer into an object
   *
   * @param buffer an input buffer containing flattened object
   * @param len the input buffer length
   * @param object a new populated object
   */
  void 
  VFloatData::readObject (const std::string& input)  throw (SerializeException)
  {
    float* vcfbuf;
    unsigned short id, pad;
    unsigned int   num;
    unsigned int hlen = 2 * sizeof(unsigned short) + sizeof(unsigned int);
    unsigned int idlen = sizeof(unsigned short);
    
    const char* buf = input.data();
    ::memcpy (&id, buf, idlen);
    ::memcpy (&num, &buf[2*idlen], sizeof(unsigned int));

    id = ntoh_short (id);
    num = ntoh_int (num);
    
    if (id != serialID()) 
      throw SerializeException ("VFloatData", "Serial ID error");

    if (num * sizeof(float) != input.length() - hlen) 
      throw SerializeException ("VFloatData", "Serial Binary Data Size error");

    // copy all float values
    vcfbuf = new float[num];
    ::memcpy (vcfbuf, &buf[hlen], input.length() - hlen);

    // resize this vector
    data_.resize(num);

    // assign vector elements
    for (int i = 0; i < num; i++) 
      data_[i] = ntoh_int((int)vcfbuf[i]);

    delete []vcfbuf;
  }
}
