#ifndef QPHIX_MEMMAP_H
#define QPHIX_MEMMAP_H

#include "qphix/print_utils.h"
#include <utility>
#include <map>


namespace QPhiX
{
  typedef std::pair<unsigned long, unsigned long> MemMapValue_t;
  
  extern std::map<unsigned long, unsigned long>  buffer_alloc_map;
  extern std::map<unsigned long, unsigned long>  aligned_alloc_map;
  
  
  void init_maps()
  {
    buffer_alloc_map.clear();
    aligned_alloc_map.clear();
  }
  
  void print_maps()
  {
    masterPrintf("Buffer Map:\n");
    if( !buffer_alloc_map.empty() ) { 
      std::map<unsigned long, unsigned long>::iterator map_iterator = buffer_alloc_map.begin();
      while( map_iterator != buffer_alloc_map.end() ) { 
	masterPrintf( " key = %lu  value = %lu\n", map_iterator->first, map_iterator->second );
	map_iterator++;
      }
    }
    else { 
      masterPrintf(" MAP IS EMPTY\n");
    }

    masterPrintf("Aligned Map:\n");
    if( !aligned_alloc_map.empty() ) { 
      std::map<unsigned long, unsigned long>::iterator map_iterator = aligned_alloc_map.begin();
      while( map_iterator != aligned_alloc_map.end() ) { 
	masterPrintf( " key = %lu  value = %lu\n", map_iterator->first, map_iterator->second );
	map_iterator++;
      }
    }
    else { 
      masterPrintf(" MAP IS EMPTY\n");
    }
  }
}; // namespace  
#endif
