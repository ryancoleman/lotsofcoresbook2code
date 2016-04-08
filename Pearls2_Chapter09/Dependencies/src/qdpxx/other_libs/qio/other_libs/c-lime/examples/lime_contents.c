/* Display contents of a LIME formated file */
/* Balint Joo 2003 */
/* C. DeTar 10/26/04 reformatted output */

#include <lime_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <lime.h>
#include <lime_fixed_types.h>
#define MAX_BYTES 64000

/* Scan for non-ASCII characters */
/* Return true if all characters are ASCII */
int all_ascii(char *buf, size_t length){
  size_t i;

  for(i = 0; i < length; i++)
    if(0x80 & buf[i])return 0;
  return 1;
}

int main(int argc, char *argv[]) 
{
  char* data_buf;
  
  LimeReader *reader;
  FILE *fp;
  int status;
  n_uint64_t nbytes, read_bytes;
  int msg,rec,first;
  char *lime_type;
  size_t bytes_pad;
  int MB_flag, ME_flag;
  
  if( argc != 2 ) { 
    fprintf(stderr, "Usage: %s <lime_file>\n", argv[0]);
    return EXIT_FAILURE;
  }
  

  fp = DCAPL(fopen)(argv[1], "r");
  if(fp == (FILE *)NULL) { 
    fprintf(stderr,"Unable to open file %s for reading\n", argv[1]);
    return EXIT_FAILURE;
  }

  reader = limeCreateReader(fp);
  if( reader == (LimeReader *)NULL ) { 
    fprintf(stderr, "Unable to open LimeReader\n");
    return EXIT_FAILURE;
  }

  msg = 0; first = 1; rec = 0;
  while( (status = limeReaderNextRecord(reader)) != LIME_EOF ){
    
    if( status != LIME_SUCCESS ) { 
      fprintf(stderr, "limeReaderNextRecord returned status = %d\n", 
	      status);
      return EXIT_FAILURE;
    }

    nbytes    = limeReaderBytes(reader);
    lime_type = limeReaderType(reader);
    bytes_pad = limeReaderPadBytes(reader);
    MB_flag   = limeReaderMBFlag(reader);
    ME_flag   = limeReaderMEFlag(reader);
    
    if (MB_flag == 1 || first)
      {
	first = 0;
	rec = 0;
	msg++;
      }

    rec++;

    printf("\n\n");
    printf("Message:        %d\n", msg);
    printf("Record:         %d\n", rec);
    printf("Type:           %s\n", lime_type);
    printf("Data Length:    %llu\n", (unsigned long long)nbytes);
    printf("Padding Length: %lu\n", (unsigned long)bytes_pad);
    printf("MB flag:        %d\n", MB_flag);
    printf("ME flag:        %d\n", ME_flag);
    
    /* TO DO: Buffer the input */
    if(nbytes < MAX_BYTES){
      data_buf = (char *)malloc((size_t)nbytes+1);
      if( data_buf == (char *)NULL) { 
	fprintf(stderr, "Couldn't malloc data buf\n");
	return EXIT_FAILURE;
      }
      
      read_bytes = nbytes;
      status = limeReaderReadData((void *)data_buf, &read_bytes, reader);
      
      if( status < 0 ) { 
	if( status != LIME_EOR ) { 
	  fprintf(stderr, "LIME read error occurred: status= %d  %llu bytes wanted, %llu read\n", 
		  status, (unsigned long long)nbytes, 
		  (unsigned long long)read_bytes);
	  return EXIT_FAILURE;
	}
      }
      
      data_buf[nbytes]='\0';
      if(!all_ascii(data_buf, nbytes))
	printf("Data:           [Binary data]\n");
      else
	printf("Data:           \"%s\" \n", data_buf);
      
      free(data_buf);
    }
    else{
	printf("Data:           [Long record skipped]\n");
    }

  }

  limeDestroyReader(reader);
  DCAP(fclose)(fp);

  return EXIT_SUCCESS;

}   
    
  
