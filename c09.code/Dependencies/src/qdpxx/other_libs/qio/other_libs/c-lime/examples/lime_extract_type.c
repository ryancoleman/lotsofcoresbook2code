/* Extract the first LIME record of a given type */
/* C. DeTar 2005 */
/* Usage ...


      lime_extract_type <lime_file> <lime_type>
*/

#include <lime_config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lime.h>
#include <lime_fixed_types.h>
/*#define MAXBUF 1048576*/
#define MAXBUF 8

n_uint64_t mino(n_uint64_t i, n_uint64_t j){
  return i < j ? i : j;
}

int main(int argc, char *argv[]) 
{
  char buf[MAXBUF];
  LimeReader *reader;
  FILE *fp;
  int status;
  n_uint64_t nbytes, bytes_left, bytes_to_copy, read_bytes;
  int rec, msg;
  char *lime_type;
  char *lime_type_target;
  //size_t bytes_pad;
  int MB_flag/*, ME_flag*/;
  
  if( argc < 3 ) { 
    fprintf(stderr, "Usage: %s <lime_file> <lime_type>\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* Open file */

  fp = DCAPL(fopen)(argv[1], "r");
  if(fp == (FILE *)NULL) { 
    fprintf(stderr,"Unable to open file %s for reading\n", argv[1]);
    return EXIT_FAILURE;
  }

  /* Tarte LIME type */

  lime_type_target = argv[2];

  /* Open LIME reader */
  reader = limeCreateReader(fp);
  if( reader == (LimeReader *)NULL ) { 
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(EXIT_FAILURE);
  }

  /* Loop over records */
  rec = 0;
  msg = 0;
  while( (status = limeReaderNextRecord(reader)) != LIME_EOF ){
    
    if( status != LIME_SUCCESS ) { 
      fprintf(stderr, "limeReaderNextRecord returned status = %d\n", 
	      status);
      exit(EXIT_FAILURE);
    }

    nbytes    = limeReaderBytes(reader);
    lime_type = limeReaderType(reader);
    //bytes_pad = limeReaderPadBytes(reader);
    MB_flag   = limeReaderMBFlag(reader);
    //ME_flag   = limeReaderMEFlag(reader);

    /* Update message and record numbers */
    if(MB_flag == 1){
      msg++;
      rec = 0;
    }

    rec++;

    /* Skip to next record until target record is reached */
    if (strcmp(lime_type,lime_type_target) != 0) continue;
    
    /* Buffered copy */
    bytes_left = nbytes;
    while(bytes_left > (n_uint64_t)0){
      bytes_to_copy = mino((n_uint64_t)MAXBUF,bytes_left);
      read_bytes = bytes_to_copy;
      status = limeReaderReadData((void *)buf, &read_bytes, reader);
    
      if( status < 0 && status != LIME_EOR ) { 
	fprintf(stderr, "LIME read error occurred: status= %d", status);
	return EXIT_FAILURE;
      }
      if (read_bytes != bytes_to_copy) {
	fprintf(stderr, "Read error %lld bytes wanted,%lld read\n", 
		(unsigned long long)nbytes, (unsigned long long)read_bytes);
	return EXIT_FAILURE;
      }
    
      /* Print to stdout */
      
      fwrite(buf, bytes_to_copy, 1, stdout);
      bytes_left -= bytes_to_copy;
    }

    /* Quit at this record */
    break;
  }

  limeDestroyReader(reader);
  DCAP(fclose)(fp);

  return EXIT_SUCCESS;
}   
    
  
