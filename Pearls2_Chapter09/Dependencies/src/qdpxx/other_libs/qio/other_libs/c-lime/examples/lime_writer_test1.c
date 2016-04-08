/* Create a LIME test file */
/* Balint Joo 2003 */
/* 12/31/04 C. DeTar added more records to file */
/* 4/17/05  C. DeTar changed lyrics and added seek test */

#include <lime_config.h>
#include <stdio.h>
#include <lime.h>
#include <string.h>
#include <stdlib.h>

int write_rec(LimeWriter *writer, int MB_flag, int ME_flag, int shuffle,
	      char message[], char lime_type[]){
  off_t totbytes = strlen(message);
  off_t seek;
  n_uint64_t bytes;
  LimeRecordHeader *h;
  int status=EXIT_SUCCESS;
  char *bufstart;

  /* Write record header */
  fprintf(stderr, "Creating Header\n");
  h = limeCreateHeader(MB_flag, ME_flag, lime_type, totbytes);
  
  fprintf(stderr, "Writing Header\n");
  status = limeWriteRecordHeader( h, writer );

  if( status < 0 ) { 
    fprintf(stderr, "LIME write header error %d\n", status);
    return EXIT_FAILURE;
  }

  limeDestroyHeader(h);

  /* Write the record in pieces just to test multiple calls */

  bufstart = message;
  bytes = totbytes/2;

  if(!shuffle)
    {
      fprintf(stderr, "Writing first part of data\n"); fflush(stderr);
      status = limeWriteRecordData(bufstart, &bytes, writer);
      
      if( status != LIME_SUCCESS ) { 
	fprintf(stderr, "LIME write error %d\n", status);
	return EXIT_FAILURE;
      }
      fprintf(stderr, "Wrote %llu bytes\n", (unsigned long long)bytes);

      bufstart += bytes;
      bytes = totbytes - bytes;
      fprintf(stderr, "Writing second part of data\n"); fflush(stderr);
      status = limeWriteRecordData(bufstart, &bytes, writer);

      if( status != LIME_SUCCESS ) { 
	fprintf(stderr, "LIME write error %d\n", status);
	return EXIT_FAILURE;
      }
      fprintf(stderr, "Wrote %llu bytes\n", (unsigned long long)bytes);
    }
  else
    {
      seek = strlen(message)-bytes;
      bufstart += seek;

      fprintf(stderr, "Seeking to second part of record\n"); fflush(stderr);
      status = limeWriterSeek(writer, seek, SEEK_SET);
      if( status != LIME_SUCCESS ) { 
	fprintf(stderr, "LIME seek error %d\n", status);
	return EXIT_FAILURE;
      }

      fprintf(stderr, "Writing second part of data\n"); fflush(stderr);
      status = limeWriteRecordData(bufstart, &bytes, writer);
      if( status != LIME_SUCCESS ) { 
	fprintf(stderr, "LIME write error %d\n", status);
	return EXIT_FAILURE;
      }
      fprintf(stderr, "Wrote %llu bytes\n", (unsigned long long)bytes);

      bufstart -= seek;
      bytes = seek;

      fprintf(stderr, "Seeking to first part of record\n"); fflush(stderr);
      status = limeWriterSeek(writer, -strlen(message), SEEK_CUR);
      if( status != LIME_SUCCESS ) { 
	fprintf(stderr, "LIME seek error %d\n", status);
	return EXIT_FAILURE;
      }

      fprintf(stderr, "Writing first part of data\n"); fflush(stderr);
      status = limeWriteRecordData(bufstart, &bytes, writer);
      if( status != LIME_SUCCESS ) { 
	fprintf(stderr, "LIME write error %d\n", status);
	return EXIT_FAILURE;
      }
      fprintf(stderr, "Wrote %llu bytes\n", (unsigned long long)bytes);

      status = limeWriterCloseRecord(writer);
    }
  return status;
}

int main(int argc, char *argv[]) 
{

  char lime_file[] = "lime_file_test";

  LimeWriter *writer;
  FILE *fp;

  /* Open the file for the LimeWriter */
  fprintf(stderr, "Opening file %s\n", lime_file);
  fp = DCAPL(fopen)(lime_file, "w");
  if(fp == (FILE *)NULL) { 
    fprintf(stderr, "Unable to open %s\n", lime_file);
    return EXIT_FAILURE;
  }

  /* Set up the LimeWriter */
  fprintf(stderr, "Creating Writer\n");
  writer = limeCreateWriter(fp);
  if( writer == (LimeWriter *)NULL ) { 
    fprintf(stderr, "Unable to initialise LIME\n");
    return EXIT_FAILURE;
  }

  /* Write some messages */

  write_rec(writer,1,0,0,"Doctor! Ain't there nothin' I can take, I say", "lime-test-text1");
  write_rec(writer,0,1,0,"Doctor! To relieve this bellyache, I say", "lime-test-text1");
  write_rec(writer,1,0,1,"You put the lime in the coconut", "lime-test-text1");
  write_rec(writer,0,1,1,"drink 'em both together", "lime-test-text2");
  write_rec(writer,1,1,0,"Harry Nilsson, 1971", "lime-test-text3");

  limeDestroyWriter(writer);
  DCAP(fclose)(fp);

  return EXIT_SUCCESS;
}
