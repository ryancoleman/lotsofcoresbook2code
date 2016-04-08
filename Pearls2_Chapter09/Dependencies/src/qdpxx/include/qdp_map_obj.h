// -*- C++ -*-
/*! \file
 * \brief Wrapper over maps
 */

#ifndef __qdp_map_obj_h__
#define __qdp_map_obj_h__

#include <map>
#include <vector>

#include "qdp.h"

namespace QDP
{

  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObject
  {
  public:
    //! Virtual Destructor
    virtual	
    ~MapObject() {}

    /**
     * Does this key exist in the store
     * @param key a key object
     * @return true if the answer is yes
     */
    virtual
    bool exist(const K& key) const = 0;
			
    //! Insert
    virtual
    int insert(const K& key, const V& val) = 0;
 
    //! Other accessor
    virtual
    int get(const K& key, V& val) const = 0;

    //! Flush out state of object
    virtual
    void flush() = 0;

    //! Size of Map
    virtual
    unsigned int size() const = 0;

    /**
     * Return all available keys to user
     * @param keys user suppled an empty vector which is populated
     * by keys after this call.
     */
    virtual
    void keys(std::vector<K>& keys_) const = 0;

    //! Insert user data into the metadata database
    virtual
    int insertUserdata(const std::string& user_data) = 0;
    
    //! Get user user data from the metadata database
    virtual
    int getUserdata(std::string& user_data) const = 0;
  };

} // namespace QDP

#endif
