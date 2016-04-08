#include <lime_config.h>
#include <lime.h>
#include <stdio.h>
#include <stdlib.h>

#include <lime_utils.h>
#include <lime_binary_header.h>

/* If we don't have fseeko, define it */
#include <lime_fseeko.h>

//#undef LIME_DEBUG

/* Forward declarations for internal routines */
int skipReaderBytes(LimeReader *r, off_t bytes_to_skip);
int readAndParseHeader(LimeReader *r);

LimeReader* limeCreateReader(FILE *fp)
{
  LimeReader *ret_val;

  ret_val = (LimeReader *)malloc(sizeof(LimeReader));
  if( ret_val == (LimeReader *)NULL ) { 
    return NULL;
  }

  ret_val->first_read = 0;
  ret_val->is_last = 0;
  ret_val->header_nextP = 1;
  ret_val->fp = fp;
  ret_val->curr_header = (LimeRecordHeader *)NULL;
  ret_val->bytes_left = 0;     /* Obsolete value */
  ret_val->bytes_total = 0;
  ret_val->rec_ptr = 0;
  ret_val->rec_start = 0;
  ret_val->bytes_pad = 0;
  return ret_val;
}

/* Move reader to a new position in the file.  We assume the new
   position marks the beginning of a LIME record, so we reset the
   reader state accordingly.  For random access inside an open LIME
   record, see limeReaderSeek. */
   
int limeSetReaderPointer(LimeReader *r, off_t offset){
  int status;

  /* Position file as requested */
  status = DCAPL(fseeko)(r->fp, offset, SEEK_SET);  
  if(status < 0)return LIME_ERR_SEEK;
  r->first_read = 0;
  r->is_last = 0;
  r->header_nextP = 1;
  r->curr_header = (LimeRecordHeader *)NULL;
  r->bytes_left = 0;     /* Obsolete value */
  r->bytes_total = 0;
  r->rec_ptr = 0;
  r->rec_start = 0;
  r->bytes_pad = 0;
  
  return LIME_SUCCESS;
}

/* Report the file pointer position of the beginning of the next LIME
   record or return the current position if we haven't read a header,
   yet. */

off_t limeGetReaderPointer(LimeReader *r){
  /* If we haven't read anything yet, we don't know anything. */
  if(r->first_read == 0)
    return DCAPL(ftello)(r->fp);
  else{
    return r->rec_start + r->bytes_total + r->bytes_pad;
  }
}

void limeDestroyReader(LimeReader *r)
{
  if (r != (LimeReader *)NULL) { 

    if ( r->curr_header != (LimeRecordHeader *)NULL ) { 
      limeDestroyHeader(r->curr_header);
      r->curr_header = (LimeRecordHeader *)NULL;
    }
    free(r);

  }
}

int limeReaderNextRecord(LimeReader *r)
{
  int status;
  char myname[] = "limeReaderNextRecord";

  if( r == (LimeReader *)NULL ) { 
    printf("%s LIME_ERR_PARAM\n",myname);
    return LIME_ERR_PARAM;
  }
  
  if( r->first_read == 0 ) { 
    /* We haven't yet read the first record */
    /* Call an auxiliary function to read and parse the header 
       -- this call may allocate memory */
    status = readAndParseHeader(r);

    /* We allow the caller to handle the error condition. */
    /* An EOF condition could be normal here. */
    if ( status < 0 ) { 
      if( status != LIME_EOF )
	printf("%s returning %d\n",myname,status);
      return status;
    }
    r->first_read = 1;
  }
  else { 
    /* We have read the first record -- so presumably we have 
       already got a header in the reader.*/

    /* In this case what we must do is skip to the end of the current
       record */
    status = skipReaderBytes(r, r->bytes_total - r->rec_ptr);
    if ( status < 0 ) { 
      if( status != LIME_EOF )
	printf("%s returns %d\n",myname,status);
      return status;
    }
    
    /* Right we have now skipped to the end of the previous 
       record and we believe it not to be the last record 
       so we can now safely destroy the current record header 
       and try to read the next header */
    status = readAndParseHeader(r);
    if ( status < 0 ){
      if( status != LIME_EOF )
	printf("%s returns %d\n",myname,status);
      return status;
    }
  }

  if(r->curr_header == NULL) {
    printf("%s No header info on this node\n",myname);
    return LIME_ERR_PARAM;
  }

  r->is_last = r->curr_header->ME_flag;
  r->bytes_left = r->curr_header->data_length;
  r->bytes_total = r->curr_header->data_length;
  r->rec_ptr = 0;
  r->rec_start = DCAPL(ftello)(r->fp);
  r->bytes_pad = lime_padding(r->bytes_total);

  return status;
}


int limeReaderReadData(void *dest, n_uint64_t *nbytes, LimeReader *r)
{
  n_uint64_t bytes_to_read;
  n_uint64_t bytes_read;
  char myname[] = "limeReaderReadData";

  /* Check if we are at the end of the record */
  if( r->rec_ptr == r->bytes_total ) {
    return LIME_EOR;
  }

  /* If we are not at the end then read how much we have to */
  if( *nbytes > 0 ) { 
    if( *nbytes + r->rec_ptr > r->bytes_total ) {
      bytes_to_read = r->bytes_total - r->rec_ptr;
    }
    else { 
      bytes_to_read = *nbytes;
    }

    if((size_t)bytes_to_read != bytes_to_read){
      printf("%s Can't read %llu bytes\n",
	     myname,(unsigned long long)bytes_to_read);
      return LIME_ERR_READ;
    }

    /* Actually read */
    bytes_read = DCAP(fread)(dest, sizeof(unsigned char), 
		       (size_t)bytes_to_read, r->fp);
    *nbytes = bytes_read;
    
    if( bytes_read != bytes_to_read ){
      printf("%s tried to read %llu bytes but got %llu bytes\n",
	     myname,
	     (unsigned long long)bytes_to_read,
	     (unsigned long long)bytes_read);
      if(feof(r->fp))
	printf("Unexpected EOF encountered\n");
      return LIME_ERR_READ;
    }

    r->bytes_left -= bytes_read;
    r->rec_ptr += bytes_read;
  }

  return LIME_SUCCESS;
}

/* Advance to beginning of the next record */
int limeReaderCloseRecord(LimeReader *r)
{
  int status;
  n_uint64_t offset;
  char myname[] = "limeReaderCloseRecord";

  /* Advance to the beginning of the next record */
  offset = r->bytes_total - r->rec_ptr;

  if((off_t)offset != offset){
    printf("%s: can't skip %llu bytes\n",myname,(long long unsigned)offset);
    return LIME_ERR_SEEK;
  }

  status = skipReaderBytes(r, (off_t)offset);

  return status;
}
      
/* Skip bytes within the current record.  Positive and negative
   offsets are allowed.  If the skip takes us to the end of the data,
   skip any padding as well.  If the skip takes us out of bounds of
   the data in the current record, skip to nearest boundary of the
   record and return an error. */

int skipReaderBytes(LimeReader *r, off_t bytes_to_skip)
{
  int status;
  n_uint64_t new_rec_ptr;
  n_uint64_t offset;
  char myname[] = "skipReaderBytes";


  new_rec_ptr = r->rec_ptr + bytes_to_skip;

  /* Prevent skip past the end of the data */
  if( new_rec_ptr > r->bytes_total ){
    new_rec_ptr = r->bytes_total;
    printf("%s: Seeking past end of data\n",myname);fflush(stdout);
    status = LIME_ERR_SEEK;
  }

  /* Prevent skips before the beginning of the data */
  /* In this case set the new pointer to the beginning of the record */
  if(new_rec_ptr < 0){
    new_rec_ptr = 0;
    printf("%s: Seeking before beginning end of data\n",myname);fflush(stdout);
    status = LIME_ERR_SEEK;
  }

  /* Seek */
  /* If there will be no bytes left, include padding in the seek */
  if(new_rec_ptr == r->bytes_total)
    offset = r->rec_start + new_rec_ptr + r->bytes_pad;
  else{
    offset = r->rec_start + new_rec_ptr;
  }

  /* Guard against insufficient integer size */
  if((off_t)offset != offset){
    printf("%s: fseeko can't seek to %llu. off_t too small.\n",
	   myname,(unsigned long long)offset);
    return LIME_ERR_SEEK;
  }

  status = DCAPL(fseeko)(r->fp, (off_t)offset , SEEK_SET);

  if(status < 0){
    printf("%s: fseek returned %d\n",myname,status);fflush(stdout);
    return LIME_ERR_SEEK;
  }

  /* Update the reader state */
  r->rec_ptr = new_rec_ptr;
  
  return LIME_SUCCESS;
}

int limeReaderSeek(LimeReader *r, off_t offset, int whence)
{
  int status;
  char myname[] = "limeReaderSeek";

  if(whence == SEEK_CUR){
    status = skipReaderBytes(r, offset);
  }
  else if(whence == SEEK_SET){
    status = skipReaderBytes(r, offset - r->rec_ptr);
  }
  else{
    fprintf(stderr, "%s code %x not implemented yet\n",myname,whence);  
    status = LIME_ERR_SEEK;
  }
  return status;
}


/* Accessors for header information */

/* Return MB flag in current header */
int limeReaderMBFlag(LimeReader *r){
  if(r == NULL)return -1;
  if(r->curr_header == NULL)return -1;
  return r->curr_header->MB_flag;
}

/* Return ME flag in current header */
int limeReaderMEFlag(LimeReader *r){
  if(r == NULL)return -1;
  if(r->curr_header == NULL)return -1;
  return r->curr_header->ME_flag;
}

/* Return pointer to LIME type string in current header */
char *limeReaderType(LimeReader *r){
  if(r == NULL)return NULL;
  if(r->curr_header == NULL)return NULL;
  return r->curr_header->type;
}

/* Return number of total bytes in current record */
n_uint64_t limeReaderBytes(LimeReader *r){
  if(r == NULL)return 0;
  if(r->curr_header == NULL)return 0;
  return r->bytes_total;
}

/* Return number of padding bytes in current record */
size_t limeReaderPadBytes(LimeReader *r){
  if(r == NULL)return 0;
  if(r->curr_header == NULL)return 0;
  return r->bytes_pad;
}

/* Entrance assumption to this function is that:
   i) The stream pointer is pointing to the beginning of 
      a record header.

   ii) There is a header about to come 
*/
int readAndParseHeader(LimeReader *r)
{

  unsigned int i_version;
  int i_MB, i_ME;
  n_uint32_t i_magic_no;
  n_uint64_t  i_data_length;
  unsigned char *typebuf;
  int status;
  //int msg_begin = 1;
  char myname[] = "lime::readAndParseHeader";

  /* Destroy any old header structure kicking about */
  if( r->curr_header != (LimeRecordHeader *)NULL ) { 
    //msg_begin = r->is_last;   /* Remember whether we finished a message */
    limeDestroyHeader(r->curr_header);
    r->curr_header = (LimeRecordHeader *)NULL;
  }

  /* Read the entire header */

  status = DCAP(fread)((void *)lime_header.int64, 
		 sizeof(n_uint64_t), MAX_HDR64, r->fp);
  if( status != MAX_HDR64 ) {
    if( DCAP(feof)( r->fp ) ) return LIME_EOF;
    fprintf(stderr,"%s read %d but wanted %d\n",myname,status,MAX_HDR64);
    return LIME_ERR_READ;
  }

  /* Check magic number */

  i_magic_no = big_endian_long(*lime_hdr_magic_no);
  if(i_magic_no != LIME_MAGIC_NO){
    fprintf(stderr,"%s: wrong magic number: read %x but wanted %x\n",
	    myname,i_magic_no,LIME_MAGIC_NO);
    return LIME_ERR_READ;
  }

#ifdef LIME_DEBUG
  fprintf(stderr, "%s Magic number OK: %d\n ", myname, i_magic_no);
  fflush(stderr);
#endif

  /* LIME version number */

  i_version = big_endian_short(*lime_hdr_version);

#ifdef LIME_DEBUG
  fprintf(stderr, "%s Input Version: %d\n ", myname,(int)i_version);
  fflush(stderr);
#endif

  /* MB flag */

  if ( *lime_hdr_mbme & MB_MASK ) { 
    i_MB = 1; 
  } 
  else {
    i_MB = 0;
  }

#ifdef LIME_DEBUG
  fprintf(stderr, "%s MB Flag: %d\n ", myname, (int)i_MB);
  fflush(stderr);
#endif
  
  /* ME Flag */

  if( *lime_hdr_mbme & ME_MASK ) {
    i_ME = 1;
  }
  else { 
    i_ME = 0;
  }

#ifdef LIME_DEBUG
  fprintf(stderr, "%s ME Flag: %d\n ", myname, (int)i_ME);
  fflush(stderr);
#endif

  /* Data length */

  i_data_length = big_endian_long_long(*lime_hdr_data_len);

#ifdef LIME_DEBUG
  fprintf(stderr, "%s Data Length: %llu\n ", myname, (unsigned long long)i_data_length);
  fflush(stderr);  
#endif

  /* Record type. */
  typebuf = (unsigned char *)lime_hdr_rec_type;

  /* Sanity Checking */
  /* Check Version */
  if( i_version != LIME_VERSION ) { 
    fprintf(stderr, "%s Unknown Lime Version\n",myname);
    exit(EXIT_FAILURE);
  }

#ifdef LIME_DEBUG
  printf("%s: type %s\n",myname,typebuf);
#endif

  /* If we are the first packet we MUST have MB flag set */
  /**  suppressed to allow seeking withing file without penalty
  if( msg_begin != i_MB ) { 
    fprintf(stderr, "%s MB Flag incorrect: last ME = %d but current MB=%d \n",
	    myname, msg_begin, i_MB);
    exit(EXIT_FAILURE);
  }
  **/

  r->curr_header = limeCreateHeader(i_MB, i_ME,
				    (char*)typebuf,
				    i_data_length);


  if (r->curr_header == (LimeRecordHeader *)NULL ) { 
    fprintf(stderr, "%s ERROR; Couldn't create header\n",myname);
    exit(EXIT_FAILURE);
  }

  return LIME_SUCCESS;
}

/* Manipulator to set the reader to a prescribed state.  We use this
   functionality to synchronize multinode reading from the same file.
   To synchronize, have the master node call limeCreateReader. Have
   the master node broadcast the resulting LimeReader structure to the
   secondary nodes.  Have each secondary node call rdest =
   limeCreateReader and then call this procedure with rsrc, the
   broadcast master node's reader. */

int limeReaderSetState(LimeReader *rdest, LimeReader *rsrc )
{
  int status;
  char myname[] = "limeReaderSetState";

  /* Set rdest reader state from rsrc */
  /* We do not copy the file pointer member fp or the curr_header member */
  rdest->first_read   = rsrc->first_read   ;
  rdest->is_last      = rsrc->is_last      ;
  rdest->header_nextP = rsrc->header_nextP ;
  rdest->bytes_total  = rsrc->bytes_total  ;
  rdest->bytes_left   = rsrc->bytes_left   ;
  rdest->rec_ptr      = rsrc->rec_ptr      ;
  rdest->rec_start    = rsrc->rec_start    ;
  rdest->bytes_pad    = rsrc->bytes_pad    ;

  /* Now make the system state agree with the reader state */
  status = DCAPL(fseeko)(rdest->fp, rdest->rec_start + rdest->rec_ptr, 
			 SEEK_SET);
  if(status < 0){
    printf("%s: fseek returned %d\n",myname,status);fflush(stdout);
    return LIME_ERR_SEEK;
  }
  return LIME_SUCCESS;
}
