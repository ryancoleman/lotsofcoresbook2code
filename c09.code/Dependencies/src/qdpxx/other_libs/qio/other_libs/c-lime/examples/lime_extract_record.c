/* Extract a single LIME record */
/* Balint Joo 2003 */
/* Usage ...

      lime_extract_record <lime_file> <msgno> <recno> <output_file>

   where msgno is the message number and recno is the record number
   (1-based enumeration)

*/

#include <lime_config.h>
#include <stdio.h>
#include <stdlib.h>
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
  FILE *fp,*fpout;
  int status;
  n_uint64_t nbytes, bytes_left, bytes_to_copy, read_bytes;
  int rec_seek,msg_seek;
  int rec, msg;
  //char *lime_type;
  //size_t bytes_pad;
  int MB_flag/*, ME_flag*/;
  
  if( argc < 5 ) { 
    fprintf(stderr, "Usage: %s <lime_file> <msgno> <recno> <output_file>\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* Open file */

  fp = DCAPL(fopen)(argv[1], "r");
  if(fp == (FILE *)NULL) { 
    fprintf(stderr,"Unable to open file %s for reading\n", argv[1]);
    return EXIT_FAILURE;
  }

  fpout = DCAPL(fopen)(argv[4], "w");
  if(fpout == (FILE *)NULL) { 
    fprintf(stderr,"Unable to open file %s for writing\n", argv[4]);
    return EXIT_FAILURE;
  }

  /* Decode message and record number from command line */

  msg_seek = atoi(argv[2]);
  if (msg_seek <= 0) {
    fprintf(stderr,"Invalid message number = %d\n", msg_seek);
    return EXIT_FAILURE;
  }

  rec_seek = atoi(argv[3]);
  if (rec_seek <= 0) {
    fprintf(stderr,"Invalid record number = %d\n", rec_seek);
    return EXIT_FAILURE;
  }

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
    //lime_type = limeReaderType(reader);
    //bytes_pad = limeReaderPadBytes(reader);
    MB_flag   = limeReaderMBFlag(reader);
    //ME_flag   = limeReaderMEFlag(reader);

    /* Update message and record numbers */
    if(MB_flag == 1){
      msg++;
      rec = 0;
    }

    rec++;

#if 0
    printf("\n\n");
    //printf("Type:           %s\n",   lime_type);
    printf("Data Length:    %ld\n",  nbytes);
    //printf("Padding Length: %d\n",   bytes_pad);
    printf("MB flag:        %d\n",   MB_flag);
    //printf("ME flag:        %d\n",   ME_flag);
#endif


    /* Skip to next record until target record is reached */
    if (msg != msg_seek || rec != rec_seek) continue;
    
    /* Buffered copy */
    bytes_left = nbytes;
    while(bytes_left > (n_uint64_t)0){
      bytes_to_copy = mino((n_uint64_t)MAXBUF,bytes_left);
      read_bytes = bytes_to_copy;
      status = limeReaderReadData((void *)buf, &read_bytes, reader);
    
      if( status < 0 && status != LIME_EOR ) { 
	fprintf(stderr, "LIME read error occurred: status= %d\n", status);
	return EXIT_FAILURE;
      }
      if (read_bytes != bytes_to_copy) {
	fprintf(stderr, "Read error %lld bytes wanted,%lld read\n", 
		(unsigned long long)nbytes, (unsigned long long)read_bytes);
	return EXIT_FAILURE;
      }
    
      /* Print to stdout */
      
      fwrite(buf, bytes_to_copy, 1, fpout);
      bytes_left -= bytes_to_copy;
    }

    /* Quit at this record */
    break;
  }

  limeDestroyReader(reader);
  DCAP(fclose)(fp);
  DCAP(fclose)(fpout);

  return EXIT_SUCCESS;
}   
    
  
