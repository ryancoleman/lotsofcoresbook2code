/* Pack files into LIME format.  
   Each file consists of the payload for one LIME record */
/* F. Maresca 12/22/04 with some code by Balint Joo */

/* Usage...

   lime_pack <filelist_file> <lime_file>

   Files to be packed are listed in a file list file.
   Files are packed in the order listed.
   Each line of the file list file has the format

     <file_name>  <lime_type>

   A blank line indicates a message break.  Thus

      file1 type1
      file2 type2
      
      file3 type3

   generates two LIME messages, the first containing two records with
   payloads file1 and file2 and the second consisting of a single
   record with payload file3.

   According to LIME format, file names are not retained in the LIME
   file.  Only the LIME type and order of appearance identifies them.

*/

/* These defines probably need to be config'ed */
/* Needed by Gnu C for ftello and fseeko */
/* #define _LARGEFILE_SOURCE */
/* #define _FILE_OFFSET_BITS 64 */

#include <lime_config.h>
#include <lime.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "lime_fseeko.h"

#define MAXLINE 512
#define MAXFILENAME MAXLINE
#define MAXLIMETYPE MAXLINE
#define MAXBUF 1048576

typedef struct {
  int eof;
  int ss_count;
  char filename[MAXFILENAME];
  char lime_type[MAXLIMETYPE];
} List_Line;


/* Discover how many bytes there are in the file */
off_t file_size(FILE *fp)
{
  off_t oldpos = ftello(fp);
  off_t length;
  
  if (DCAPL(fseeko)(fp, 0L,SEEK_END) == -1)
    return -1;
  
  length = DCAPL(ftello)(fp);
  
  return ( DCAPL(fseeko)(fp,oldpos,SEEK_SET) == -1 ) ? -1 : length;
  
}

void read_list_line(List_Line *list, FILE *fp_list)
{
  char *status;
  char line[MAXLINE];

  /* Read ahead one line */
  status         = fgets(line, sizeof(line), fp_list);
  list->eof      = (status == NULL);
  list->ss_count = sscanf(line,"%s %s", 
			  list->filename, list->lime_type);
}


int write_hdr(n_uint64_t bytes, char *type, 
	      int MB_flag, int ME_flag, LimeWriter *dg)
{
  LimeRecordHeader *h;
  int status;

  h = limeCreateHeader(MB_flag,ME_flag,type,bytes);
  status = limeWriteRecordHeader( h, dg);
  
  if( status < 0 ) 
    {
      fprintf(stderr, "LIME write header error %d\n", status);
      return EXIT_FAILURE;
    }
  
  limeDestroyHeader(h);
  
  return EXIT_SUCCESS;
}

int write_buf(char *buf, n_uint64_t bytes, LimeWriter *dg)
{
  int status;
  
  /* Write the record */
  status = limeWriteRecordData(buf, &bytes, dg);
  
  if( status < 0 ) 
    {
      fprintf(stderr, "LIME write error %d\n", status);
      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}

int min(int i, int j){
  return i < j ? i : j;
}

int main(int argc, char *argv[])
{
  FILE *fp_list,*fp_dest,*fp_src;
  char *listfile, *limefile;
  List_Line curr, next;
  char buf[MAXBUF];
  n_uint64_t bytes, bytes_left, bytes_to_copy;
  LimeWriter *dg;
  int MB_flag, ME_flag;
  int msg, rec;
  int status;
  
  if( argc != 3 ) 
    {
      fprintf(stderr, "Usage: %s <filelist_file> <lime_file>\n", argv[0]);
      return EXIT_FAILURE;
    }
  
  listfile = argv[1];
  limefile = argv[2];
  
  /* Open the file that contains the list of files to pack */
  /* fprintf(stderr,"Open file %s for reading\n", listfile); */
  fp_list = fopen(argv[1], "r");
  if(fp_list == NULL) 
    {
      fprintf(stderr,"Unable to open file %s for reading\n", 
	      listfile);
      return EXIT_FAILURE;
    }
  
  /* Open the LIME file for writing */
  fp_dest = DCAPL(fopen)(limefile, "w");
  if(fp_dest == (FILE *)NULL) 
    {
      fprintf(stderr, "Unable to open file %s for writing\n", limefile);
      return EXIT_FAILURE;
    }
  
  /* Set up the LimeWriter */
  /* fprintf(stderr, "Creating Writer\n"); */
  dg = limeCreateWriter(fp_dest);
  if( dg == (LimeWriter *)NULL ) 
    {
      fprintf(stderr, "Unable to initialise LIME\n");
      return EXIT_FAILURE;
    }

  /* Initialize message begin/end flags */

  MB_flag = 1;
  ME_flag = 0;  /* Tentative */
  msg = 0;
  rec = 0;

  /* Loop over entries in the file list file */

  /* Read ahead one line */
  read_list_line(&next, fp_list);
  
  printf("msg rec   bytes type\n");

  while ( 1 )
    {
      /* Copy next line values to current line */
      curr = next;
      if (curr.eof) break;

      /* Read ahead one line in file list */

      read_list_line(&next, fp_list);
      
      if (curr.ss_count <= 0)
	{
	  /* printf("New message\n"); */
	  
	  /* Reinitialize message begin/end flags */
	  MB_flag = 1;
	  ME_flag = 0;  /* Tentative */
	  continue;     /* Blank line: No file to process */
	}
      
      else if (curr.ss_count == 1)
	{
	  printf("Error: a file list line has only one entry :\n %s\n", 
		 curr.filename);
	  break;
	}
      
      else if (curr.ss_count == 2 ) 
	{
	  if (MB_flag == 1)
	    {
	      rec = 0;
	      msg++;
	    }
	  
	  rec++;
	  
	  /* Open file and get its length */
	  fp_src = DCAPL(fopen)(curr.filename,"r");
	  if(fp_src == (FILE *)NULL) 
	    {
	      fprintf(stderr, "Unable to open %s\n",curr.filename);
	      break;
	    }
	  else
	    if ((bytes = file_size(fp_src)) == -1)
	      fprintf(stderr,"Can't compute length of file %s\n", 
		      curr.filename);

	  /* Reset stream pointer to beginning of file */
	  fseeko(fp_src,0L,SEEK_SET);
	  
	  
	  /* End of message?  Signaled by blank line or EOF on file list */
	  if ( next.ss_count <= 0 || next.eof ) ME_flag = 1;
	    
	  /* Announce file */
	  printf("%2d %3d %8ld %s\n\t\t%s\n",
		 msg, rec, (long int)bytes, curr.lime_type, curr.filename);

	  /* Write header */
	  status = write_hdr(bytes,curr.lime_type,MB_flag,ME_flag,dg);
	  if (status != 0)break;

	  /* Buffered copy */
	  bytes_left = bytes;
	  while(bytes_left > 0){
	    bytes_to_copy = min(MAXBUF,bytes_left);
	    if( bytes_to_copy != fread(buf,1,bytes_to_copy,fp_src))
	      {
		fprintf(stderr, "Error reading %s\n", curr.filename);
		return EXIT_FAILURE;
	      }
	  
	    status = write_buf(buf,bytes_to_copy,dg);
	    if (status != 0)break;

	    bytes_left -= bytes_to_copy;
	  }

	  /* The next record is not the begining of a message */
	  MB_flag = 0;

	  DCAP(fclose)(fp_src);
	}
    }
  
  
  limeDestroyWriter(dg);
  DCAP(fclose)(fp_dest);
  fclose(fp_list);
  return EXIT_SUCCESS;
}
