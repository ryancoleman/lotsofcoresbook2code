// -*- C++ -*-
/*! \file
 *  \brief A memory based map object 
 */

#ifndef __qdp_map_obj_memory_h__
#define __qdp_map_obj_memory_h__

#include "qdp_map_obj.h"
#include <vector>
#include <unordered_map>

namespace QDP
{


  //----------------------------------------------------------------------------
  //! Iterator for a MapObject
  template<typename K, typename V>
  class MapObjectMemoryIterator : public std::iterator<std::input_iterator_tag, 
						 std::unordered_map<std::string, std::pair<K,V> > >
  {
  public:
    //! Map type convenience
    typedef std::unordered_map<std::string, std::pair<K,V> > MapType_t;

    //! Take an input iterator coming from MapObject
    MapObjectMemoryIterator(const typename MapType_t::const_iterator& iter_) : src_iter(iter_) {}

    //! Copy constructor
    MapObjectMemoryIterator(const MapObjectMemoryIterator& p) : src_iter(p.src_iter) {}

    //! Prefix increment
    MapObjectMemoryIterator& operator++() {++src_iter; return *this;}
 
    //! Postfix increment
    MapObjectMemoryIterator operator++(int) 
    {
      MapObjectMemoryIterator<K,V> tmp(*this); 
      operator++(); 
      return tmp;
    }

    //! Equality
    bool operator==(const MapObjectMemoryIterator& rhs) const {return src_iter == rhs.src_iter;}

    //! Rather obvious
    bool operator!=(const MapObjectMemoryIterator& rhs) const {return src_iter != rhs.src_iter;}

    // Coercion
    const std::pair<K,V>& operator*() const {return src_iter->second;}
    
    // Coercion
    const std::pair<K,V>* operator->() const {return &(src_iter->second);}
    
  private:
    //! Hold the iterator
    typename MapType_t::const_iterator src_iter;
    
  };


  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectMemory : public MapObject<K,V>
  {
  public:
    //! Re-export these
    typedef K key_type;
    typedef V mapped_type;
    typedef std::pair<K,V> value_type;
    typedef MapObjectMemoryIterator<K,V> const_iterator;

    //! Map type convenience
    typedef std::unordered_map<std::string, value_type> MapType_t;
    
    //! Default constructor
    MapObjectMemory() {}

    //! Destructor
    ~MapObjectMemory() {}

    //! Insert user data into the metadata database
    int insertUserdata(const std::string& user_data_) 
    {
      user_data = user_data_;
      return 0;
    }
    
    //! Get user user data from the metadata database
    int getUserdata(std::string& user_data_) const 
    {
      user_data_ = user_data;
      return 0;
    }

    //! Insert
    int insert(const K& key, const V& val) 
    {
      return this->insert(std::make_pair(key,val));
    }
			
    //! Insert from iterators
    template<typename InputIterator>
    int insert(InputIterator first, InputIterator last)
    {
      int ret = 0;
      for(InputIterator iter = first; iter != last; ++iter)
	this->insert(*iter);
      return ret;
    }
			
    //! Insert
    int insert(const value_type& vv) 
    { 
      int ret = 0;
      BinaryBufferWriter bin;
      write(bin, vv.first);

      const std::string bin_key(bin.str());
      typename MapType_t::iterator iter = src_map.find(bin_key);
      if (iter != src_map.end())
      {
	iter->second = vv;
      }
      else
      {
	src_map.insert(std::make_pair(bin_key,vv));
      }
      return ret;
    }
			
    //! Getter
    int get(const K& key, V& val) const
    {
      int ret = 0;
      BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::const_iterator iter = src_map.find(bin.str());
      if (iter == src_map.end()) {
	ret = 1;
      }
      else {
	val = iter->second.second;
      }
      return ret;
    }
			
    //! Erase a key-value
    void erase(const K& key) 
    {
      BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::const_iterator iter = src_map.find(bin.str());
      if (iter != src_map.end())
      {
	src_map.erase(iter);
      }
    }
			
    //! Clear the object
    void clear() {src_map.clear();}

    //! Flush out state of object
    void flush() {}

    //! Exists?
    bool exist(const K& key) const
    {
      BinaryBufferWriter bin;
      write(bin, key);
      return (src_map.find(bin.str()) == src_map.end()) ? false : true;
    }
			
    //! The number of elements
//    size_t size() const {return src_map.size();}
    unsigned int size() const {return static_cast<unsigned long>(src_map.size());}

    //! Dump keys
    void keys(std::vector<K>& _keys) const
    {
      _keys.resize(0);
      for(typename MapType_t::const_iterator iter  = src_map.begin(); iter != src_map.end(); ++iter)
      {
	_keys.push_back(iter->second.first);
      }
    }

    //! Dump keys and values
    virtual void keysAndValues(std::vector<K>& _keys, std::vector<V>& _vals) const 
    {
      _keys.resize(0);
      _vals.resize(0);
      for(typename MapType_t::const_iterator iter  = src_map.begin(); iter != src_map.end(); ++iter)
      {
	_keys.push_back(iter->second.first);
	_vals.push_back(iter->second.second);
      }
    }

    // Extensions to the basic MapObject
    //! Getter
    const V& operator[](const K& key) const 
    {
      BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::const_iterator iter = src_map.find(bin.str());
      if (iter == src_map.end())
      {
	std::cerr << "MapObject: key not found" << std::endl;
	exit(1);
      }

      return iter->second.second;
    }
			
    //! Setter
    V& operator[](const K& key) 
    {
      BinaryBufferWriter bin;
      write(bin, key);

      typename MapType_t::iterator iter = src_map.find(bin.str());
      if (iter == src_map.end())
      {
	std::cerr << "MapObject: key not found" << std::endl;
	exit(1);
      }

      return iter->second.second;
    }

    //! Begin loop over keys
    virtual const_iterator begin() const {return MapObjectMemoryIterator<K,V>(src_map.begin());}

    //! End loop over keys
    virtual const_iterator end() const {return MapObjectMemoryIterator<K,V>(src_map.end());}


  protected:  
    //! Map of objects
    mutable MapType_t src_map;

    //! Metadata
    std::string user_data;
  };



  //----------------------------------------------------------------------------
  //! Read a MapObject via xml
  template<typename K, typename V>
  inline
  void read(XMLReader& xml, const std::string& s, MapObjectMemory<K,V>& input)
  {
    XMLReader arraytop(xml, s);

    std::ostringstream error_message;
    std::string elemName = "elem";
  
    int array_size;
    try {
      array_size = arraytop.count(elemName);
    }
    catch( const std::string& e) { 
      error_message << "Exception occurred while counting " << elemName
		    << " during array read " << s << std::endl;
      arraytop.close();
      throw error_message.str();
    }
      
    // Get the elements one by one
    for(int i=0; i < array_size; i++) 
    {
      std::ostringstream element_xpath;

      // Create the query for the element 
      element_xpath << elemName << "[" << (i+1) << "]";

      // recursively try and read the element.
      try {
	XMLReader xml_elem(arraytop, element_xpath.str());

	K key;
	V val;

	read(xml_elem, std::string("Key"), key);
	read(xml_elem, std::string("Val"), val);

	input.insert(key, val);
      } 
      catch (const std::string& e) 
      {
	error_message << "Failed to match element " << i
		      << " of array  " << s << "  with query " << element_xpath.str()
		      << std::endl
		      << "Query returned error: " << e;
	arraytop.close();
	throw error_message.str();
      }
    }

    // Arraytop should self destruct but just to be sure.
    arraytop.close();
  }


  //----------------------------------------------------------------------------
  //! Write a MapObject in xml
  template<typename K, typename V>
  inline
  void write(XMLWriter& xml, const std::string& path, const MapObjectMemory<K,V>& param)
  {
    push(xml, path);

    for(typename MapObjectMemory<K,V>::const_iterator iter = param.begin(); iter != param.end(); ++iter)
    {
      push(xml, "elem");
      write(xml, "Key", iter->first);
      write(xml, "Val", iter->second);
      pop(xml);
    }

    pop(xml);
  }

} // namespace QDP

#endif
