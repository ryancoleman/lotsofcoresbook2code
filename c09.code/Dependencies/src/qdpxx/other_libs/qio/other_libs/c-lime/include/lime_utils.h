#ifndef LIME_UTILS_H
#define LIME_UTILS_H

#include "lime_fixed_types.h"
 
  /** \brief Utility function for padding headers and data to 4byte boundaries
   *
   *  \param nbytes  is the number of bytes in the data
   *  \returns the number of bytes of padding required to make the data 4byte
   *           padded
   */
  int lime_padding(size_t nbytes);

  /** \brief Utility function to convert a native long long to big endian or
   *  to convert a big_endian long long to native. (Its the same operation.)
   *
   *  \param the long long to be converted
   *
   *  \returns the input long long (converted if necessary)
   *
   *  On a big endian machine, this routine does nothing.
   *  On a little endian machine, the routine swaps byte order
   *  Here endianness refers solely to byte order */

  n_uint64_t big_endian_long_long(n_uint64_t ull);

  /** \brief Utility function to convert a native long to big endian or
   *  to convert a big_endian long to native. (Its the same operation.)
   *
   *  \param the long to be converted
   *
   *  \returns the input long (converted if necessary)
   *
   *  On a big endian machine, this routine does nothing.
   *  On a little endian machine, the routine swaps byte order
   *  Here endianness refers solely to byte order */

  n_uint32_t big_endian_long(n_uint32_t ul);

  /** \brief Utility function to convert a native short to big endian or
   *  to convert a big_endian short to native. (Its the same operation.)
   *
   *  \param the long to be converted
   *
   *  \returns the input short (converted if necessary)
   *
   *  On a big endian machine, this routine does nothing.
   *  On a little endian machine, the routine swaps byte order
   *  Here endianness refers solely to byte order */

  n_uint16_t big_endian_short(n_uint16_t us);

#endif
