#include <lime_config.h>
#include <stdlib.h>
#include <stdio.h>
#include <lime_defs.h>
#include <lime_utils.h>
#include <lime_fixed_types.h>

int lime_padding(size_t nbytes)
{
  int ret_val;
    
  if (nbytes > 0) {

    /* nbytes > 0 -- we may need padding */
    if( nbytes % 8 != 0 ) { 

      /* nbytes not divisible by 8 */
      /* padding is 8 - remainder of nbytes on division by 4 */
      ret_val = 8 - (nbytes % 8);
    }
    else { 

      /* nbytes is divisible by 8 -- no padding needed */
      ret_val = 0;
    }
  }
  else { 
    /* nbytes is 0, no padding needed */
    ret_val = 0;
  }

  return ret_val;
}

int lime_big_endian(void)
{
  union {
    int  l;
    char c[sizeof(int)];
  } u;
  u.l = 1;

  return (u.c[sizeof(int)-1] == 1 ? 1 : 0); 
}

void lime_byte_swap(void *ptr, size_t size, size_t nmemb)
{
  unsigned int j;
  char myname[] = "lime_byte_swap";
                                                                                
  char char_in[8];              /* characters used in byte swapping */
                                                                                
  char *in_ptr;
  double *double_ptr;           /* Pointer used in the double routines */
                                                                                
  switch (size)
  {
  case 4:  /* n_uint32_t */
  {
    n_uint32_t *w = (n_uint32_t *)ptr;
    register n_uint32_t old, recent;
                                                                                
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
        j++, double_ptr++){
                                                                                
                                                                                
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
                                                                                
  case 2:  /* n_uint16_t */
  {
    n_uint16_t *w = (n_uint16_t *)ptr;
    register n_uint16_t old, recent;
                                                                                
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
    fprintf(stderr,"%s: unsupported word size = %lu\n",
	    myname,(unsigned long)size);
    exit(1);
  }
}


n_uint64_t big_endian_long_long(n_uint64_t ull)
{
   n_uint64_t ret_val = ull;

   if ( lime_big_endian() == 0 ) {
	/* We are little endian. We need to convert */
        lime_byte_swap((void *)&ret_val, sizeof(n_uint64_t), 1);
 
   }
   return ret_val;
}
 
n_uint32_t big_endian_long(n_uint32_t ul)
{
   n_uint32_t ret_val = ul;

   if ( lime_big_endian() == 0 ) {
	/* We are little endian. We need to convert */
        lime_byte_swap((void *)&ret_val, sizeof(n_uint32_t), 1);
 
   }
   return ret_val;
}
 
n_uint16_t big_endian_short(n_uint16_t us)
{
   n_uint16_t ret_val = us;
   if( lime_big_endian() == 0 ) { 
       /* We are little endian. Need to convert */
       lime_byte_swap((void *)&ret_val, sizeof(n_uint16_t), 1);
   }
   return ret_val;
}
