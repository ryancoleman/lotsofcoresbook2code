#include <qio_config.h>
#include <lrl.h>
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1 // for fdopen in stdio
#endif
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
//#include <qmp.h>

/** 
 * Open a file for reading 
 *
 * \param filename   file for writing  ( Read )
 *
 * \return null if failure
 */
LRL_FileReader *LRL_open_read_file(const char *filename)
{
  LRL_FileReader *fr = NULL;
  /* Open and check for a readable file */
  FILE *fpt = DCAPL(fopen)(filename,"r");
  if(fpt != NULL) {
    /* Set up LRL_FileReader structure */
    fr = (LRL_FileReader *) malloc(sizeof(LRL_FileReader));
    if(fr != NULL) {
      fr->file = fpt;
      fr->dr = limeCreateReader(fr->file);
      if(fr->dr == NULL) {
	free(fr);
	fr = NULL;
      }
    }
    if(fr == NULL) {
      DCAP(fclose)(fpt);
    }
  }
  return fr;
}

/**
 * Seek to a new file pointer position
 *
 * \param fr     File reader
 * \param offset Absolute seek position.
 */

int LRL_set_reader_pointer(LRL_FileReader *fr, off_t offset){
  int status;

  /* Position file as requested */
  status = limeSetReaderPointer(fr->dr, offset);
  if(status != LIME_SUCCESS)return LRL_ERR_SEEK;
  return LRL_SUCCESS;
}

/**
 * Get the current file pointer position
 *
 * \param fr     File reader
 */

off_t LRL_get_reader_pointer(LRL_FileReader *fr){
  return limeGetReaderPointer(fr->dr);  
}

/** 
 * Open a file for writing 
 *
 * \param filename   file for writing  ( Read )
 *
 * \return null if failure
 */
LRL_FileWriter *LRL_open_write_file(const char *filename, int mode)
{
  LRL_FileWriter *fw = (LRL_FileWriter *) malloc(sizeof(LRL_FileWriter));
  if(fw != NULL) {
    /* Open according to requested mode */
    if(mode == LRL_APPEND) {
      fw->file = DCAPL(fopen)(filename,"a");
    } else if(mode == LRL_NOTRUNC) {
      int fd = DCAP(open)(filename, O_WRONLY, 0666);
      fw->file = fdopen(fd, "w");
    } else {
      fw->file = DCAPL(fopen)(filename,"w");
    }
    if(fw->file == NULL) {
      printf("%s: failed to open %s for writing\n", __func__, filename);
      free(fw);
      fw = NULL;
    } else {
      fw->dg = limeCreateWriter(fw->file);
      if(fw->dg == NULL) {
	printf("%s: limeCreateWriter failed\n", __func__);
	DCAP(fclose)(fw->file);
	free(fw);
	fw = NULL;
      }
    }
  }
  return fw;
}

/** 
 * Open a record for reading
 *
 * \param fr         LRL file reader  ( Read )
 * \param rec_size   record size ( Modify )
 * \param tag        tag for the record header ( Read )
 *
 * \return null if failure and set status to error flag
 */
LRL_RecordReader *LRL_open_read_record(LRL_FileReader *fr,
				       uint64_t *rec_size, 
				       LIME_type *lime_type, int *status)
{
  LRL_RecordReader *rr;
  int lime_status;
  char myname[] = "LRL_open_read_record";

  if (fr == NULL)
    return NULL;

  rr = (LRL_RecordReader *)malloc(sizeof(LRL_RecordReader));
  if (rr == NULL){
    printf("%s: Can't malloc reader\n",myname);
    return NULL;
  }
  rr->fr = fr;
  
  /* Get next record header */
  lime_status = limeReaderNextRecord(rr->fr->dr);
  if (lime_status != LIME_SUCCESS){
    if(lime_status == LIME_EOF)*status = LRL_EOF;
    else{
      printf("%s lime error %d getting next record header\n",myname,
	     lime_status);
      *status = LRL_ERR_READ;
    }
    return NULL;
  }

  /* Extract record information */
  *rec_size = limeReaderBytes(rr->fr->dr);
  *lime_type = limeReaderType(rr->fr->dr);

  *status = LRL_SUCCESS;
  return rr;
}


/** 
 * Search for the next record with one of the specified LIME types and
 * open it for reading
 *
 * \param fr               LRL file reader  ( Read )
 * \param lime_type_list   list of LIME types ( Read )
 * \param ntypes           number in list ( Read )
 *                         If 0 accept any LIME type.
 * \param rec_size         record size ( Modify )
 * \param lime_type_found  LIME type found ( Modify )
 * \param status           LRL status ( Modify )
 *
 * \return null if failure and set status to error flag
 */
LRL_RecordReader *LRL_open_read_target_record(LRL_FileReader *fr,
	      LIME_type *lime_type_list, int ntypes, uint64_t *rec_size, 
	      LIME_type *lime_type_found, int *status)
{
  LRL_RecordReader *rr;
  int i, found;

  rr = LRL_open_read_record(fr, rec_size, lime_type_found, status);
  if(rr == NULL)return rr;

  found = 0;
  while(1){
    for(i = 0; i < ntypes; i++){
      if(strcmp(lime_type_list[i],*lime_type_found) == 0){
	found = 1;
	break;
      }
    }
    /* If the list of targets is empty, accept any record */
    if(ntypes == 0 || found)break;

    *status = LRL_close_read_record(rr);
    if(*status != LRL_SUCCESS)return NULL;
    rr = LRL_open_read_record(fr, rec_size, lime_type_found, status);
    if(rr == NULL)return rr;
  }

  *status = LIME_SUCCESS;
  return rr;
}


/** 
 * Create a record writer.  Don't write anything.
 *
 * \param fw         LRL file writer  ( Read )
 *
 * \return null if failure
 */
LRL_RecordWriter *LRL_create_record_writer(LRL_FileWriter *fw)
{
  LRL_RecordWriter *rw;
  char myname[] = "LRL_create_record_writer";
  
  if (fw == NULL)
    return NULL;

  rw = (LRL_RecordWriter *)malloc(sizeof(LRL_RecordWriter));
  if (rw == NULL){
    printf("%s: Can't malloc writer\n",myname);
    return NULL;
  }
  rw->fw = fw;

  return rw;
}

/** 
 * Write the record header
 *
 * \param rw         LRL record writer  ( Read )
 * \param rec_size   record size ( Modify )
 * \param tag        tag for the record header ( Read )
 *
 * \return LRL_ERR_WRITE if failure, LRL_SUCCESS if successful.
 */
int LRL_write_record_header(LRL_RecordWriter *rw, 
			    int msg_begin, int msg_end, 
			    uint64_t rec_size, 
			    LIME_type lime_type)
{
  LimeRecordHeader *h;
  int status;
  char myname[] = "LRL_write_record_header";

  /* Create and write record header */
  h = limeCreateHeader(msg_begin, msg_end, lime_type, rec_size);
  status = limeWriteRecordHeader(h, rw->fw->dg);

  if (status < 0)
  { 
    printf("%s: fatal error. LIME status is: %d\n", myname, status);
    return LRL_ERR_WRITE;
  }

  limeDestroyHeader(h);

  return LRL_SUCCESS;
}

/** 
 * Open a record for writing and write the record header
 *
 * \param fw         LRL file writer  ( Read )
 * \param rec_size   record size ( Modify )
 * \param tag        tag for the record header ( Read )
 *
 * \return null if failure
 */
LRL_RecordWriter *LRL_open_write_record(LRL_FileWriter *fw, 
					int msg_begin, int msg_end, 
					uint64_t rec_size, 
					LIME_type lime_type)
{
  LRL_RecordWriter *rw;
  int status;
  
  if (fw == NULL)
    return NULL;

  rw = LRL_create_record_writer(fw);
  if (rw == NULL)
    return NULL;

  /* Create and write record header */
  status = LRL_write_record_header(rw, msg_begin, msg_end, rec_size, 
				   lime_type);
  if (status != LRL_SUCCESS)
    return NULL;

  return rw;
}

/** 
 * Copy reader
 *
 * \param rr         LRL record reader
 * \param state_ptr  Returned opaque state
 * \param state_size Returned size
 *
 */

void LRL_get_reader_state(LRL_RecordReader *rr,
			  void **state_ptr, size_t *state_size)
{
  LimeReader *spt;
  char myname[] = "LRL_get_reader_state";

  /* For now the LIME reader structure defines the state.  It would be
     cleaner to define a new state structure and copy just the reader
     components we need */
  spt = (LimeReader *)malloc(sizeof(LimeReader));
  if (spt == NULL){
    printf("%s: Can't malloc reader state\n",myname);
    *state_ptr = NULL;
    *state_size = 0;
  }
  else{
    if(rr) *spt = *(rr->fr->dr);
    *state_ptr = (void *)spt;
    *state_size = sizeof(LimeReader);
  }
}


/** 
 * Copy writer
 *
 * \param rw         LRL record writer
 * \param state_ptr  Returned opaque state
 * \param state_size Returned size
 *
 * \return copy of writer
 */

void LRL_get_writer_state(LRL_RecordWriter *rw,
			  void **state_ptr, size_t *state_size)
{
  LimeWriter *spt;
  char myname[] = "LRL_get_writer_state";

  /* For now the LIME writer structure defines the state.  It would be
     cleaner to define a new state structure and copy just the reader
     components we need */
  spt = (LimeWriter *)malloc(sizeof(LimeWriter));
  if (spt == NULL){
    printf("%s: Can't malloc writer state\n",myname);
    *state_ptr = NULL;
    *state_size = 0;
  }
  else{
    if(rw != NULL)
      *spt = *(rw->fw->dg);
    *state_ptr = (void *)spt;
    *state_size = sizeof(LimeWriter);
  }
}

/** 
 * Set reader state from a given state
 *
 * \param rw         LRL record reader  ( Write )
 * \param state_ptr  Reader state ( Read )
 *
 * \return LRL status
 */
int LRL_set_reader_state(LRL_RecordReader *rr, void *state_ptr){
  LimeReader *rsrc = (LimeReader *)state_ptr;
  int status;

  if(rr) {
    /* Set the LIME reader state to the state specified by state_ptr */
    status = limeReaderSetState(rr->fr->dr, rsrc);
    if(status != LIME_SUCCESS)return LRL_ERR_SETSTATE;

    status = limeReaderSeek(rr->fr->dr, 0, SEEK_SET);
    if(status != LIME_SUCCESS)return LRL_ERR_SEEK;
  }

  return LRL_SUCCESS;
}

/** 
 * Set writer state from a given state
 *
 * \param rw         LRL record writer  ( Write )
 * \param state_ptr  Writer state ( Read )
 *
 * \return LRL status
 */
int LRL_set_writer_state(LRL_RecordWriter *rw, void *state_ptr){
  LimeWriter *wsrc = (LimeWriter *)state_ptr;
  int status;

  if(rw) {
    /* Set the LIME writer state to the state specified by state_ptr */
    status = limeWriterSetState(rw->fw->dg, wsrc);
    if(status != LIME_SUCCESS)return LRL_ERR_SETSTATE;

    status = limeWriterSeek(rw->fw->dg, 0, SEEK_SET);
    if(status != LIME_SUCCESS)return LRL_ERR_SEEK;
  }

  return LRL_SUCCESS;
}

/** 
 * Read bytes
 *
 * \param rr         LRL record reader  ( Read )
 * \param buf        buffer for reading ( Write )
 * \param nbytes     number of bytes to read ( Read )
 *
 * \return number of bytes read
 */
uint64_t LRL_read_bytes(LRL_RecordReader *rr, char *buf, 
		      uint64_t nbytes)
{
  int status;
  uint64_t nbyt = nbytes;
  char myname[] = "LRL_read_bytes";

  if (rr == NULL)
    return 0;

  //printf("node %i pos %i reading %i bytes\n", QMP_get_node_number(), rr->fr->dr->rec_ptr, nbytes);
  status = limeReaderReadData((void *)buf, &nbyt, rr->fr->dr);
  if( status != LIME_SUCCESS ) 
  { 
    printf("%s: LIME error %d has occurred\n", myname, status);
    exit(EXIT_FAILURE);
  }

  return nbyt;
}


/** 
 * Write bytes
 *
 * \param rw         LRL record writer  ( Read )
 * \param buf        buffer for writing ( Read )
 * \param nbytes     number of bytes to write ( Read )
 *
 * \return number of bytes written
 */
uint64_t LRL_write_bytes(LRL_RecordWriter *rw, char *buf, 
		       uint64_t nbytes)
{
  int status;
  uint64_t nbyt = nbytes;

  if (rw == NULL)
    return 0;

  status = limeWriteRecordData(buf, &nbyt, rw->fw->dg);

  if( status != LIME_SUCCESS ) 
  { 
    printf("LRL_write_bytes: some error has occurred. status is: %d\n", status);
    exit(EXIT_FAILURE);
  }

  return nbyt;
}


/* For seeking to offset bytes from the beginning of the record
   payload.  We are not allowed to go beyond the end of the
   payload. */

int LRL_seek_read_record(LRL_RecordReader *rr, off_t offset)
{
  int status;

  if (rr == NULL || rr->fr == NULL)return LRL_ERR_SEEK;
  status = limeReaderSeek(rr->fr->dr, offset, SEEK_SET);

  if( status != LIME_SUCCESS ) 
  { 
    printf("LRL_seek_read_record: some error has occurred. status is: %d\n", status);
    return LRL_ERR_SEEK;
  }
  return LRL_SUCCESS;
}

/* For seeking to offset bytes from the beginning of the record
   payload.  We are not allowed to go beyond the end of the
   payload. */

int LRL_seek_write_record(LRL_RecordWriter *rw, off_t offset)
{
  int status;
  char myname[] = "LRL_seek_write_record";

  if (rw == NULL){
    printf("%s: null record writer\n",myname);
    return LRL_ERR_SEEK;
  }
  if(rw->fw == NULL){
    printf("%s: null file writer\n",myname);
    return LRL_ERR_SEEK;
  }

  status = limeWriterSeek(rw->fw->dg, offset, SEEK_SET);

  if( status != LIME_SUCCESS ) 
  { 
    printf("%s: LIME error %d\n", myname,status);
    return LRL_ERR_SEEK;
  }
  return LRL_SUCCESS;
}

/* For skipping to the end of the current message */

int LRL_next_message(LRL_FileReader *fr)
{
  int status;
  int msg_end = 0;
  char myname[] = "LRL_next_message";

  if(fr == NULL)return LRL_ERR_SKIP;
  while(msg_end == 0){
    status = limeReaderNextRecord(fr->dr);
    msg_end = limeReaderMEFlag(fr->dr);
    if( status != LIME_SUCCESS ) 
      { 
	printf("%s: LIME error %d\n", myname,status);
	return LRL_ERR_SKIP;
      }
  }
  return LRL_SUCCESS;
}

/* For skipping to the beginning of the next record? */

int LRL_next_record(LRL_RecordReader *rr)
{
  int status;
  LRL_FileReader *fr;

  if(rr == NULL)return LRL_ERR_SKIP;
  fr = rr->fr;
  if(fr == NULL)return LRL_ERR_SKIP;
  status = limeReaderNextRecord(fr->dr);

  if( status != LIME_SUCCESS ) 
  { 
    printf("LRL_next_record: LIME error %d\n", status);
    return LRL_ERR_SKIP;
  }
  return LRL_SUCCESS;
}

/** 
 * Destroy the LRL reader state (copy)
 *
 * \param state_ptr         LRL reader state
 *
 */
void LRL_destroy_reader_state_copy(void *state_ptr){
  if(state_ptr == NULL)return;
  free(state_ptr);
}

/** 
 * Destroy the LRL writer state (copy)
 *
 * \param state_ptr         LRL writer state
 *
 */
void LRL_destroy_writer_state_copy(void *state_ptr){
  free(state_ptr);
}

/**
 *
 * \param rr     LRL record reader
 *
 * return 0 for success and 1 for failure
 */
int
LRL_close_read_record(LRL_RecordReader *rr)
{
  int status = LRL_SUCCESS;
  if(rr == NULL) {
    status = LRL_ERR_CLOSE;
  } else {
    if(limeReaderCloseRecord(rr->fr->dr) != LIME_SUCCESS)
      status = LRL_ERR_CLOSE;
    free(rr);
  }
  return status;
}

/**
 *
 * \param rw     LRL record writer
 *
 */
int
LRL_close_write_record(LRL_RecordWriter *rw)
{
  int status = LRL_SUCCESS;
  if(rw == NULL) {
    status = LRL_ERR_CLOSE;
  } else {
    if(limeWriterCloseRecord(rw->fw->dg) != LIME_SUCCESS)
      status = LRL_ERR_CLOSE;
    free(rw);
  }
  return status;
}

/**
 * Close a file for reading
 *
 * \param fr   file buffer  ( Modify )
 *
 * \return status
 */
int
LRL_close_read_file(LRL_FileReader *fr)
{
  int status = LRL_SUCCESS;
  if(fr != NULL) {
    limeDestroyReader(fr->dr);
    if(DCAP(fclose)(fr->file)!=0) status = LRL_ERR_CLOSE;
    free(fr);
  }
  return status;
}

/** 
 * Close a file for writing
 *
 * \param fw   file buffer  ( Modify )
 *
 * \return status
 */
int
LRL_close_write_file(LRL_FileWriter *fw)
{
  int status = LRL_SUCCESS;
  if(fw != NULL) {
    limeDestroyWriter(fw->dg);
    int fstatus = DCAP(fclose)(fw->file);
    if(fstatus!=0) {
      status = LRL_ERR_CLOSE;
    }
#ifdef JCO_DEBUG
    fprintf(stderr, "%s: fclose return: %i\n", __func__, fstatus);
#endif
    free(fw);
  }
  return status;
}
