/* Unpack a LIME formated file 
   Each LIME record payload becomes a file */
/* F. Maresca 12/22/04 based on code by Balint Joo */

/* Usage...

   lime_unpack <lime_file>

   Files are unpacked into a directory named after the lime_file with
   names that encode the message number, record number, and LIME type.

   lime_file.contents/msgnn.recnn.lime_type

*/

#include <lime_config.h>
#include <stdio.h>
#include <lime.h>
#include <lime_fixed_types.h>
#include <string.h>
#include <stdlib.h>
#define MAXFILENAME 512
#define MAXDIRNAME 512
#define MAXCOMLINE 512
#define MAXBUF 1048576

n_uint64_t mino(n_uint64_t i, n_uint64_t j){
  return i < j ? i : j;
}

int main(int argc, char *argv[])
{
  FILE *fp,*fp_dest;
  char filename[MAXFILENAME], dirname[MAXDIRNAME];
  char *limefile;
  char *lime_type;
  int rec, msg, status;
  char buf[MAXBUF];
  char com[MAXCOMLINE];
  LimeReader *reader;
  n_uint64_t nbytes, bytes_left, bytes_to_copy, read_bytes;
  size_t wrote_bytes;
  //size_t bytes_pad;
  int MB_flag/*, ME_flag*/;
  
  if( argc < 2 )
    {
      fprintf(stderr, "Usage: %s <lime_file>\n", argv[0]);
      return EXIT_FAILURE;
    }

  limefile = argv[1];
  
  /* Open LIME file for reading */
  fp = DCAPL(fopen)(limefile, "r");
  if(fp == (FILE *)NULL) 
    {
      fprintf(stderr,"Unable to open file %s for reading\n", filename);
      return EXIT_FAILURE;
    }
  
  /* Open LIME reader */
  reader = limeCreateReader(fp);
  if( reader == (LimeReader *)NULL ) 
    {
      fprintf(stderr, "Unable to open LimeReader\n");
      return EXIT_FAILURE;
    }
  
  /* Create a directory from the filename */
  
  if(strlen(limefile) + strlen(".contents") > MAXDIRNAME - 1){
    fprintf(stderr,"Not enough room for directory name\n");
    return EXIT_FAILURE;
  }
  snprintf(dirname,MAXDIRNAME,"%s.contents",limefile);

  if(strlen(dirname) + strlen("mkdir -p ") > MAXCOMLINE - 1){
    fprintf(stderr,"Not enough room for command line\n");
    return EXIT_FAILURE;
  }
  snprintf(com,MAXCOMLINE,"mkdir -p %s",dirname);

  if ( (status = system(com)) != 0)
    {
      fprintf(stderr,"Error: command %s returned %d\n",com,status);
      return EXIT_FAILURE;
    }
  
  /* Loop over LIME records */
  rec = 0;
  msg = 0;

  printf("   bytes file\n");

  while( (status = limeReaderNextRecord(reader)) != LIME_EOF )
    {
      if( status != LIME_SUCCESS ) 
	{
	  fprintf(stderr, "limeReaderNextRecord returned status = %d\n",
		  status);
	  return EXIT_FAILURE;
	}
      rec++;
      
      nbytes    = limeReaderBytes(reader);
      lime_type = limeReaderType(reader);
      //bytes_pad = limeReaderPadBytes(reader);
      MB_flag   = limeReaderMBFlag(reader);
      //ME_flag   = limeReaderMEFlag(reader);
      
      if (MB_flag == 1)
	{
	  rec = 1;
	  msg++;
	  printf("\n");
	}
      
      
      /* Create a file name for the payload */
      if(strlen(dirname) + strlen("/msgnn.recnn.") + strlen(lime_type) > 
	 MAXFILENAME - 1)
	{
	  fprintf(stderr,"Not enough space for filename\n");
	  return EXIT_FAILURE;
	}
	
      snprintf(filename,MAXFILENAME,"%s/msg%02d.rec%02d.%s",
	       dirname,msg,rec,lime_type);


      /* Announce file */
      printf("%8llu %s\n", (unsigned long long)nbytes, filename);
      
      /* Open the payload file for writing */
      fp_dest = DCAPL(fopen)(filename,"w");
      if(fp_dest == NULL){
	fprintf(stderr,"Can't open %s for writing\n",filename);
	return EXIT_FAILURE;
      }
      
      /* Buffered copy */

      bytes_left = nbytes;
      while(bytes_left > (n_uint64_t)0){
	bytes_to_copy = mino((n_uint64_t)MAXBUF,bytes_left);
	read_bytes = bytes_to_copy;

	/* Read from the LIME file */
	status = limeReaderReadData((void *)buf, &bytes_to_copy, reader);
      
	if( status < 0 && status != LIME_EOR ) { 
	  fprintf(stderr, "LIME read error occurred: status= %d", status);
	  return EXIT_FAILURE;
	}
	if (read_bytes != bytes_to_copy) {
	  fprintf(stderr, "Read error %llu bytes wanted,%llu read\n", 
		  (unsigned long long)nbytes, 
		  (unsigned long long)read_bytes);
	  return EXIT_FAILURE;
	}
    
	/* Write to the payload file */
	wrote_bytes = fwrite(buf,1,bytes_to_copy,fp_dest);
	if((n_uint64_t)wrote_bytes != bytes_to_copy){
	  fprintf(stderr,"Error writing %s.  Wrote %llu bytes but wanted %llu\n",
		  filename,(unsigned long long)wrote_bytes, 
		  (unsigned long long)bytes_to_copy);
	  return EXIT_FAILURE;
	}
	
	bytes_left -= bytes_to_copy;
      }      
      
      DCAP(fclose)(fp_dest);
    }
  
  limeDestroyReader(reader);
  fclose(fp);
  
  
  return EXIT_SUCCESS;
}

