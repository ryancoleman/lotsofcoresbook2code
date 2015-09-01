//
//
//  Determine the byte order of a platform
//
//   returns
//     1  = big-endian    (alpha, intel linux, etc.)
//     2  = little-endian (sun, ibm, hp, etc)
//

#include <iostream>
#include <cstdlib>
#include "qdp_byteorder.h"

namespace QDPUtil
{
  //! Is the native byte order big endian?
  bool big_endian()
  {
    union {
      int  l;
      char c[sizeof(int)];
    } u;
    u.l = 1;
    return (u.c[sizeof(int) - 1] == 1);
  }


  //! Byte-swap an array of data each of size nmemb
  void byte_swap(void *ptr, size_t size, size_t nmemb)
  {
    unsigned int j;

    char char_in[16];		/* characters used in byte swapping */

    char *in_ptr;
    double *double_ptr;		/* Pointer used in the double routines */

    switch (size)
    {
    case 4:  /* n_uint32_t */
    {
      n_uint32_t *w = (n_uint32_t *)ptr;
      n_uint32_t old, recent;

      for(j=0; j<nmemb; j++)
      {
	old = w[j];
	recent = old >> 24 & 0x000000ff;
	recent |= old >> 8 & 0x0000ff00;
	recent |= old << 8 & 0x00ff0000;
	recent |= old << 24 & 0xff000000;
	w[j] = recent;
      }
    }
    break;

    case 1:  /* n_uint8_t: byte - do nothing */
      break;

    case 8:  /* n_uint64_t */
    {
      for(j = 0, double_ptr = (double *) ptr;
	  j < nmemb;
	  j++, double_ptr++)
      {
	in_ptr = (char *) double_ptr; /* Set the character pointer to
					 point to the start of the double */

	/*
	 *  Assign all the byte variables to a character
	 */
	char_in[0] = in_ptr[0];
	char_in[1] = in_ptr[1];
	char_in[2] = in_ptr[2];
	char_in[3] = in_ptr[3];
	char_in[4] = in_ptr[4];
	char_in[5] = in_ptr[5];
	char_in[6] = in_ptr[6];
	char_in[7] = in_ptr[7];

	/*
	 *  Now just swap the order
	 */
	in_ptr[0] = char_in[7];
	in_ptr[1] = char_in[6];
	in_ptr[2] = char_in[5];
	in_ptr[3] = char_in[4];
	in_ptr[4] = char_in[3];
	in_ptr[5] = char_in[2];
	in_ptr[6] = char_in[1];
	in_ptr[7] = char_in[0];
      }
    }
    break;

    case 16:  /* Long Long */
    {
      for(j = 0, double_ptr = (double *) ptr;
	  j < nmemb;
	  j++, double_ptr+=2)
      {

	in_ptr = (char *) double_ptr; /* Set the character pointer to
					 point to the start of the double */

	/*
	 *  Assign all the byte variables to a character
	 */
	char_in[0] = in_ptr[0];
	char_in[1] = in_ptr[1];
	char_in[2] = in_ptr[2];
	char_in[3] = in_ptr[3];
	char_in[4] = in_ptr[4];
	char_in[5] = in_ptr[5];
	char_in[6] = in_ptr[6];
	char_in[7] = in_ptr[7];
	char_in[8] = in_ptr[8];
	char_in[9] = in_ptr[9];
	char_in[10] = in_ptr[10];
	char_in[11] = in_ptr[11];
	char_in[12] = in_ptr[12];
	char_in[13] = in_ptr[13];
	char_in[14] = in_ptr[14];
	char_in[15] = in_ptr[15];

	/*
	 *  Now just swap the order
	 */
	in_ptr[0] = char_in[15];
	in_ptr[1] = char_in[14];
	in_ptr[2] = char_in[13];
	in_ptr[3] = char_in[12];
	in_ptr[4] = char_in[11];
	in_ptr[5] = char_in[10];
	in_ptr[6] = char_in[9];
	in_ptr[7] = char_in[8];

	in_ptr[8] = char_in[7];
	in_ptr[9] = char_in[6];
	in_ptr[10] = char_in[5];
	in_ptr[11] = char_in[4];
	in_ptr[12] = char_in[3];
	in_ptr[13] = char_in[2];
	in_ptr[14] = char_in[1];
	in_ptr[15] = char_in[0];

      }
    }
    break;

    case 2:  /* n_uint16_t */
    {
      n_uint16_t *w = (n_uint16_t *)ptr;
      n_uint16_t old, recent;

      for(j=0; j<nmemb; j++)
      {
	old = w[j];
	recent = old >> 8 & 0x00ff;
	recent |= old << 8 & 0xff00;
	w[j] = recent;
      }
    }
    break;

    default:
      std::cerr << __func__ << ": unsupported word size = " << size << "\n";
      exit(1);
    }
  }


  //! fread on a binary file written in big-endian order
  size_t bfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t n;

    n = fread(ptr, size, nmemb, stream);

    if (! big_endian())
    {
      /* little-endian */
      /* Swap */
      byte_swap(ptr, size, n);
    }

    return n;
  }


  //! fwrite to a binary file in big-endian order
  size_t bfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
  {
    size_t n;

    if (big_endian())
    {
      /* big-endian */
      /* Write */
      n = fwrite(ptr, size, nmemb, stream);
    }
    else
    {
      /* little-endian */
      /* Swap and write and swap */
      byte_swap(ptr, size, nmemb);
      n = fwrite(ptr, size, nmemb, stream);
      byte_swap(ptr, size, nmemb);
    }

    return n;
  }

} // namespace QDPUtil
