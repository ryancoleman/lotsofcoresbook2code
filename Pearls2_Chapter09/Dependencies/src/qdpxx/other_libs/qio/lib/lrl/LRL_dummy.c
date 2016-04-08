/* LRL_dummy.c */
/* Dummy */
/* LIME ignored */

#include <lrl.h>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

LRL_FileReader *LRL_open_read_file(const char *filename, off_t seek){
  LRL_FileReader *fr;

  fr = (LRL_FileReader *)malloc(sizeof(LRL_FileReader));
  if(fr == NULL)return NULL;

  /*** Ignore mode for now ***/
  fr->file = fopen(filename,"r");
  if(fr->file == NULL)return NULL;

  fseeko(fr->file, seek, SEEK_SET);

  return fr;
}

LRL_FileWriter *LRL_open_write_file(const char *filename, int mode){
  LRL_FileWriter *fr;

  fr = (LRL_FileWriter *)malloc(sizeof(LRL_FileWriter));
  if(fr == NULL)return NULL;

  /*** Ignore mode for now ***/
  fr->file = fopen(filename,"w");
  if(fr->file == NULL)return NULL;

  return fr;
}

LRL_RecordReader *LRL_open_read_record(LRL_FileReader *fr, size_t *rec_size, 
				       LIME_tag tag){
  LRL_RecordReader *rr;

  if(fr == NULL)return NULL;

  rr = (LRL_RecordReader *)malloc(sizeof(LRL_RecordReader));
  if(rr == NULL)return NULL;
  rr->fr = fr;
  
  /* Read and get byte size of record */
  if(fread(rec_size, 1, sizeof(size_t), fr->file) != sizeof(size_t))
    return NULL;  

  return rr;
}

LRL_RecordWriter *LRL_open_write_record(LRL_FileWriter *fr, size_t rec_size, 
					LIME_tag tag){
  LRL_RecordWriter *rr;
  
  if(fr == NULL)return NULL;

  rr = (LRL_RecordWriter *)malloc(sizeof(LRL_RecordWriter));
  if(rr == NULL)return NULL;
  rr->fr = fr;

  /* Write byte size of record */
  if(fwrite(&rec_size, 1, sizeof(size_t), fr->file) != sizeof(size_t))
    return NULL;  

  return rr;
}

size_t LRL_write_bytes(LRL_RecordWriter *rr, char *buf, size_t nbytes){
  return fwrite(buf, 1, nbytes, rr->fr->file);
}

size_t LRL_read_bytes(LRL_RecordReader *rr, char *buf, size_t nbytes){
  return fread(buf, 1, nbytes, rr->fr->file);
}

int LRL_seek_write_record(LRL_FileWriter *fr, size_t offset){
  if (fr == NULL)return LRL_ERR_SEEK;
  if(fseeko(fr->file, offset, SEEK_CUR)!=0)return LRL_ERR_SEEK;
  return LRL_SUCCESS;
}

int LRL_seek_read_record(LRL_FileReader *fr, size_t offset){
  if (fr == NULL)return LRL_ERR_SEEK;
  if(fseeko(fr->file, offset, SEEK_CUR)!=0)return LRL_ERR_SEEK;
  return LRL_SUCCESS;
}

int LRL_next_record(LRL_FileReader *fr){
  printf("LRL_next_record: not implemented\n");
  return 1;
}

int LRL_close_read_record(LRL_RecordReader *rr){
  free(rr);
  return LRL_SUCCESS;
}

int LRL_close_write_record(LRL_RecordWriter *rr){
  free(rr);
  return LRL_SUCCESS;
}

off_t LRL_close_read_file(LRL_FileReader *fr){
  off_t offset;

  if(fr == NULL)return 0;
  offset = ftello(fr->file);
  fclose(fr->file);
  free(fr);
  return offset;
}

int LRL_close_write_file(LRL_FileWriter *fr){
  if(fr == NULL)return LRL_SUCCESS;
  fclose(fr->file);
  free(fr);
  return LRL_SUCCESS;
}

