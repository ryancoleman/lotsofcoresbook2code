#include <lime_config.h>
#include <lime.h>
#include <lime_fixed_types.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Load LIME record header structure from parameters */

LimeRecordHeader *limeCreateHeader(int MB_flag,
				   int ME_flag,
				   char *type,
				   n_uint64_t reclen)
{

  LimeRecordHeader *ret_val;
  int type_length;

  ret_val = (LimeRecordHeader *)malloc(sizeof(LimeRecordHeader));
  if( ret_val == (LimeRecordHeader *)NULL ) { 
    return NULL;
  }

  /* Get the length of the type */
  type_length=strlen(type);

  ret_val->type = (char *)malloc(sizeof(char)*(type_length + 1));
  if( ret_val->type == (char *)NULL ) { 
    /* Couldn't get room for type --- destroy ret val and return null*/
    free(ret_val);
    return(NULL);
  }

  /* Copy the type string */
  strcpy(ret_val->type, type);

  /* Fill out the rest of the fields */
  ret_val->ME_flag = ME_flag;
  ret_val->MB_flag = MB_flag;
  ret_val->data_length = reclen;
  ret_val->lime_version = LIME_VERSION;


  return ret_val;
}

void limeDestroyHeader(LimeRecordHeader *h)
{
  if ( h != (LimeRecordHeader *) NULL) { 
    if (h->type != (char *)NULL ) { 
      free(h->type);
      h->type = (char *) NULL;
    }
    free(h);
  }
}
  
  
  

  
			       
