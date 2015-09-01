#ifndef LRL_H
#define LRL_H

#include <qio_config.h>
#include <qio_stdint.h>
#include <stdio.h>
#include <sys/types.h>
#include <lime.h>

/* Return codes */
#define LRL_SUCCESS      ( 0)
#define LRL_EOF          (-1)
#define LRL_ERR_READ     (-2)
#define LRL_ERR_SEEK     (-3)
#define LRL_ERR_SKIP     (-4)
#define LRL_ERR_CLOSE    (-5)
#define LRL_ERR_WRITE    (-6)
#define LRL_ERR_SETSTATE (-7)

/* For writing, either append or truncate */
#define LRL_CREAT      0
#define LRL_TRUNC      1
#define LRL_APPEND     2
#define LRL_NOTRUNC    3

#ifdef __cplusplus
extern "C"
{
#endif

/* Dummy LIME tag for now */
#define MAX_LIME_TYPE_LEN 32
typedef char* LIME_type;

typedef struct {
  FILE *file;
  LimeWriter *dg;
} LRL_FileWriter;

typedef struct {
  LRL_FileWriter *fw;
} LRL_RecordWriter;

typedef struct {
  FILE *file;
  LimeReader *dr;
} LRL_FileReader;

typedef struct {
  LRL_FileReader *fr;
} LRL_RecordReader;

LRL_FileReader *LRL_open_read_file(const char *filename);
int LRL_set_reader_pointer(LRL_FileReader *, off_t offset);
off_t LRL_get_reader_pointer(LRL_FileReader *fr);
LRL_FileWriter *LRL_open_write_file(const char *filename, int mode);
LRL_RecordReader *LRL_open_read_record(LRL_FileReader *fr, uint64_t *rec_size, 
				       LIME_type *lime_type, int *status);
LRL_RecordReader *LRL_open_read_target_record(LRL_FileReader *fr,
	      LIME_type *lime_type_list, int ntypes, uint64_t *rec_size, 
	      LIME_type *lime_type_found, int *status);
LRL_RecordWriter *LRL_create_record_writer(LRL_FileWriter *fw);
int LRL_write_record_header(LRL_RecordWriter *rw, 
			    int msg_begin, int msg_end, 
			    uint64_t rec_size, 
			    LIME_type lime_type);
LRL_RecordWriter *LRL_open_write_record(LRL_FileWriter *fr, 
					int msg_begin, int msg_end, 
					uint64_t rec_size, 
					LIME_type lime_type);
void LRL_get_writer_state(LRL_RecordWriter *rw,
			  void **state_ptr, size_t *state_size);
void LRL_get_reader_state(LRL_RecordReader *rr,
			  void **state_ptr, size_t *state_size);
int LRL_set_reader_state(LRL_RecordReader *rr, void *state_ptr);
int LRL_set_writer_state(LRL_RecordWriter *rw, void *state_ptr);
uint64_t LRL_write_bytes(LRL_RecordWriter *rr, char *buf, 
		       uint64_t nbytes);
uint64_t LRL_read_bytes(LRL_RecordReader *rr, char *buf, 
		      uint64_t nbytes);
int LRL_seek_write_record(LRL_RecordWriter *rr, off_t offset);
int LRL_seek_read_record(LRL_RecordReader *rr, off_t offset);
void LRL_destroy_reader_state_copy(void *state_ptr);
void LRL_destroy_writer_state_copy(void *state_ptr);
int LRL_next_record(LRL_RecordReader *rr);
int LRL_next_message(LRL_FileReader *fr);
int LRL_close_read_record(LRL_RecordReader *rr);
int LRL_close_write_record(LRL_RecordWriter *rr);
int LRL_close_read_file(LRL_FileReader *fr);
int LRL_close_write_file(LRL_FileWriter *fr);

#ifdef __cplusplus
}
#endif

#endif /* LRL_H */
