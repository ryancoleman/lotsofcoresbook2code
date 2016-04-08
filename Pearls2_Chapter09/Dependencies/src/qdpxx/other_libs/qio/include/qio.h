#ifndef QIO_H
#define QIO_H

#include <qio_config.h>
#include <qioxml.h>
#include <qio_string.h>
#include <lrl.h>
#include <dml.h>
#include <lime.h>

#define QIO_UNKNOWN    DML_UNKNOWN
#define QIO_SINGLEFILE DML_SINGLEFILE 
#define QIO_MULTIFILE  DML_MULTIFILE  
#define QIO_PARTFILE   DML_PARTFILE

#define QIO_FIELD      DML_FIELD
#define QIO_GLOBAL     DML_GLOBAL
#define QIO_HYPER      DML_HYPER

#define QIO_SERIAL     DML_SERIAL
#define QIO_PARALLEL   DML_PARALLEL

/* ILDG style file */
#define QIO_ILDGNO     0
#define QIO_ILDGLAT    1

/* SciDAC native or ILDG alien file (reading) */
#define QIO_SCIDAC_NATIVE  1
#define QIO_ILDG_ALIEN     2

#define QIO_CREAT      LRL_CREAT
#define QIO_TRUNC      LRL_TRUNC
#define QIO_APPEND     LRL_APPEND

/* Return codes */

#define QIO_SUCCESS               (  0)
#define QIO_EOF                   ( -1)
#define QIO_ERR_BAD_WRITE_BYTES   ( -2)
#define QIO_ERR_OPEN_READ         ( -3)
#define QIO_ERR_OPEN_WRITE        ( -4)
#define QIO_ERR_BAD_READ_BYTES    ( -5)
#define QIO_ERR_ALLOC             ( -6)
#define QIO_ERR_CLOSE             ( -7)
#define QIO_ERR_INFO_MISSED       ( -8)
#define QIO_ERR_BAD_SITELIST      ( -9)
#define QIO_ERR_PRIVATE_FILE_INFO (-10)
#define QIO_ERR_PRIVATE_REC_INFO  (-11)
#define QIO_BAD_XML               (-12)
#define QIO_BAD_ARG               (-13)
#define QIO_CHECKSUM_MISMATCH     (-14)
#define QIO_ERR_FILE_INFO         (-15)
#define QIO_ERR_REC_INFO          (-16)
#define QIO_ERR_CHECKSUM_INFO     (-17)
#define QIO_ERR_SKIP              (-18)
#define QIO_ERR_BAD_TOTAL_BYTES   (-19)
#define QIO_ERR_BAD_GLOBAL_TYPE   (-20)
#define QIO_ERR_BAD_VOLFMT        (-21)
#define QIO_ERR_BAD_IONODE        (-22)
#define QIO_ERR_BAD_SEEK          (-23)
#define QIO_ERR_BAD_SUBSET        (-24)

/* LIME types for SciDAC records */

#define QIO_LIMETYPE_PRIVATE_FILE_XML   "scidac-private-file-xml"
#define QIO_LIMETYPE_SITELIST           "scidac-sitelist"
#define QIO_LIMETYPE_FILE_XML           "scidac-file-xml"
#define QIO_LIMETYPE_PRIVATE_RECORD_XML "scidac-private-record-xml"
#define QIO_LIMETYPE_RECORD_XML         "scidac-record-xml"
#define QIO_LIMETYPE_BINARY_DATA        "scidac-binary-data"

/* LIME types for ILDG compatibility */
#define QIO_LIMETYPE_ILDG_FORMAT        "ildg-format"
#define QIO_LIMETYPE_ILDG_BINARY_DATA   "ildg-binary-data"
#define QIO_LIMETYPE_ILDG_DATA_LFN      "ildg-data-lfn"

#ifdef __cplusplus
extern "C"
{
#endif
  
/* For collecting and passing layout information */
typedef struct {
  int (*node_number)(const int coords[]);
  int (*node_index)(const int coords[]);
  void (*get_coords)(int coords[], int node, int index);
  int (*num_sites)(int node);
  int *latsize;
  int latdim;
  size_t volume;
  size_t sites_on_node;
  int this_node;
  int number_of_nodes;
} QIO_Layout;

typedef struct {
  LRL_FileWriter *lrl_file_out;
  int volfmt;
  int serpar;
  int ildgstyle;
  QIO_String *ildgLFN;
  DML_Layout *layout;
  DML_SiteList *sites;
  DML_Checksum last_checksum;
  DML_RecordWriter *dml_record_out;
} QIO_Writer;

#define QIO_RECORD_INFO_PRIVATE_NEXT 0
#define QIO_RECORD_INFO_USER_NEXT 1
#define QIO_RECORD_ILDG_INFO_NEXT 2
#define QIO_RECORD_DATA_NEXT 3
#define QIO_RECORD_CHECKSUM_NEXT 4

typedef struct {
  LRL_FileReader *lrl_file_in;
  int volfmt;
  int serpar;
  int format;
  int ildgstyle;
  int ildg_precision;
  DML_Layout *layout;
  DML_SiteList *sites;
  int read_state;
  QIO_String *xml_record;
  QIO_String *ildgLFN;
  QIO_RecordInfo record_info;
  DML_Checksum last_checksum;
  DML_RecordReader *dml_record_in;
} QIO_Reader;

typedef struct {
  int serpar;
  int volfmt;
} QIO_Iflag;

typedef struct {
  int serpar;
  int mode;
  int ildgstyle;
  QIO_String *ildgLFN;
} QIO_Oflag;

/* Support for host file conversion */

#define QIO_SINGLE_PATH 0
#define QIO_MULTI_PATH 1

typedef struct {
  int number_io_nodes;
  int type;                             /* Is node_path specified? */
  DML_io_node_t my_io_node;             /* Mapping as on compute nodes */
  DML_master_io_node_t master_io_node;  /* As on compute nodes */
  int *io_node;                         /* Only if number_io_nodes !=
					 number_of_nodes */
  char **node_path;                     /* Only if type = QIO_MULTI_PATH */
} QIO_Filesystem;

/* Internal host file conversion utilities in QIO_host_utils.c */
QIO_Layout *QIO_create_ionode_layout(QIO_Layout *layout, QIO_Filesystem *fs);
void QIO_delete_ionode_layout(QIO_Layout* layout);
QIO_Layout *QIO_create_scalar_layout(QIO_Layout *layout, QIO_Filesystem *fs);
void QIO_delete_scalar_layout(QIO_Layout *layout);
int QIO_ionode_io_node(int node);
int QIO_get_io_node_rank(int node);
int QIO_ionode_to_scalar_index(int ionode_node, int ionode_index);
int QIO_scalar_to_ionode_index(int scalar_node, int scalar_index);

/* Verbosity */
int QIO_verbose(int level);
int QIO_verbosity(void);

/* Enumerate in order of increasing verbosity */
#define QIO_VERB_OFF    0
#define QIO_VERB_LOW    1
#define QIO_VERB_MED    2
#define QIO_VERB_REG    3
#define QIO_VERB_HIGH   4
#define QIO_VERB_DEBUG  5

/* HostAPI */
int QIO_single_to_part( const char filename[], QIO_Filesystem *fs,
			QIO_Layout *layout);
int QIO_part_to_single( const char filename[], int ildgstyle, 
			QIO_String *ildgLFN, 
			QIO_Filesystem *fs, QIO_Layout *layout);
char *QIO_set_filepath(QIO_Filesystem *fs, 
		       const char * const filename, int node);

/* MPP API */
QIO_Writer *QIO_open_write(QIO_String *xml_file, const char *filename, 
			   int volfmt, QIO_Layout *layout, 
			   QIO_Filesystem *fs, QIO_Oflag *oflag);
QIO_Reader *QIO_open_read(QIO_String *xml_file, const char *filename, 
			  QIO_Layout *layout, QIO_Filesystem *fs,
			  QIO_Iflag *iflag);
int QIO_get_reader_latdim(QIO_Reader *in);
int *QIO_get_reader_latsize(QIO_Reader *in);
uint32_t QIO_get_reader_last_checksuma(QIO_Reader *in);
uint32_t QIO_get_reader_last_checksumb(QIO_Reader *in);
int QIO_set_reader_pointer(QIO_Reader *qio_in, off_t offset);
off_t QIO_get_reader_pointer(QIO_Reader *qio_in);
char *QIO_get_ILDG_LFN(QIO_Reader *qio_in);
int QIO_get_ildgstyle(QIO_Reader *in);
int QIO_get_reader_volfmt(QIO_Reader *in);
int QIO_get_reader_format(QIO_Reader *in);
void QIO_set_record_info(QIO_Reader *in, QIO_RecordInfo *rec_info);

uint32_t QIO_get_writer_last_checksuma(QIO_Writer *out);
uint32_t QIO_get_writer_last_checksumb(QIO_Writer *out);
void QIO_reset_writer_ILDG_flags(QIO_Writer *out, QIO_Oflag *oflag);

int QIO_close_write(QIO_Writer *out);
int QIO_close_read(QIO_Reader *in);

int QIO_write(QIO_Writer *out, QIO_RecordInfo *record_info,
	      QIO_String *xml_record, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      size_t datum_size, int word_size, void *arg);
int QIO_read(QIO_Reader *in, QIO_RecordInfo *record_info,
	     QIO_String *xml_record, 
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     size_t datum_size, int word_size, void *arg);
int QIO_read_record_info(QIO_Reader *in, QIO_RecordInfo *record_info,
			 QIO_String *xml_record);
int QIO_read_record_data(QIO_Reader *in, 
		 void (*put)(char *buf, size_t index, int count, void *arg),
		 size_t datum_size, int word_size, void *arg);
int QIO_next_record(QIO_Reader *in);

LRL_RecordWriter *QIO_open_write_field(QIO_Writer *out, 
    int msg_begin, int msg_end, size_t datum_size,
    const LIME_type lime_type, int *do_output, int *status);
int QIO_init_read_field(QIO_Reader *in, size_t datum_size, 
			LIME_type *lime_type_list, int ntypes,
			DML_Checksum *checksum, LIME_type *lime_type);
int QIO_seek_read_field_datum(QIO_Reader *in, 
	      DML_SiteRank seeksite,
	      void (*put)(char *buf, size_t index, int count, void *arg),
	      int count, size_t datum_size, int word_size, void *arg);
int QIO_close_read_field(QIO_Reader *in, uint64_t *nbytes);

int QIO_init_write_field(QIO_Writer *out, int msg_begin, int msg_end,
	    size_t datum_size, DML_Checksum *checksum,
	    const LIME_type lime_type);
int QIO_seek_write_field_datum(QIO_Writer *out, 
	      DML_SiteRank seeksite,
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      int count, size_t datum_size, int word_size, void *arg);
int QIO_close_write_field(QIO_Writer *out, uint64_t *nbytes);


/* Internal utilities  */
double QIO_time(void);
void QIO_wait(double sec);
void QIO_suppress_global_broadcast(QIO_Reader *qio_in);

QIO_Reader *QIO_open_read_master(const char *filename, QIO_Layout *layout,
				 QIO_Iflag *iflag, int (*io_node)(int),
				 int (*master_io_node)(void));
int QIO_open_read_nonmaster(QIO_Reader *qio_in, const char *filename,
			    QIO_Iflag *iflag);
int QIO_read_check_sitelist(QIO_Reader *qio_in);
int QIO_read_user_file_xml(QIO_String *xml_file, QIO_Reader *qio_in);
QIO_Writer *QIO_generic_open_write(const char *filename, 
				   int volfmt, QIO_Layout *layout, 
				   QIO_Oflag *oflag, 
				   int (*io_node)(int), 
				   int (*master_io_node)(void));
int QIO_reader_insert_hypercube_data(QIO_Reader *in, 
				     QIO_RecordInfo *record_info);
int QIO_writer_insert_hypercube_data(QIO_Writer *out, 
				     QIO_RecordInfo *record_info);
int QIO_write_file_header(QIO_Writer* qio_out, QIO_String *xml_file);
int QIO_read_private_record_info(QIO_Reader *in, QIO_RecordInfo *record_info);
int QIO_read_user_record_info(QIO_Reader *in, QIO_String *xml_record);
int QIO_read_ILDG_LFN(QIO_Reader *in);
int QIO_generic_read_record_data(QIO_Reader *in, 
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     size_t datum_size, int word_size, void *arg,
 	     DML_Checksum *checksum, uint64_t *nbytes);
QIO_ChecksumInfo *QIO_read_checksum(QIO_Reader *in);
int QIO_compare_checksum(int this_node,
	 QIO_ChecksumInfo *checksum_info_expect, DML_Checksum *checksum);

int QIO_generic_write(QIO_Writer *out, QIO_RecordInfo *record_info, 
	      QIO_String *xml_record, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      size_t datum_size, int word_size, void *arg,
	      DML_Checksum *checksum, uint64_t *nbytes,
	      int *msg_begin, int *msg_end);
int QIO_write_record_info(QIO_Writer *out, QIO_RecordInfo *record_info, 
              size_t datum_size, int word_size,
	      QIO_String *xml_record, 
	      int *msg_begin, int *msg_end);
int QIO_write_record_data(QIO_Writer *out, QIO_RecordInfo *record_info, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      size_t datum_size, int word_size, void *arg,
	      DML_Checksum *checksum, uint64_t *nbytes,
	      int *msg_begin, int *msg_end);
int QIO_write_checksum(QIO_Writer *out, DML_Checksum *checksum);

char *QIO_filename_edit(const char *filename, int volfmt, int this_node);
int QIO_write_string(QIO_Writer *out, int msg_begin, int msg_end,
		     QIO_String *xml,
		     const LIME_type lime_type);
int QIO_read_string(QIO_Reader *in,
		    QIO_String *xml, LIME_type *lime_type);
LRL_RecordReader *QIO_read_record_type(QIO_Reader *in, LIME_type *lime_type,
			       uint64_t *expected_rec_size, int *status);
LRL_RecordReader *QIO_open_read_target_record(QIO_Reader *in, 
    LIME_type *lime_type_list, int ntypes, LIME_type *lime_type,
    uint64_t *expected_rec_size, int *status);
int QIO_read_string_data(QIO_Reader *in, LRL_RecordReader *lrl_record_in, 
			 QIO_String *xml,  uint64_t expected_rec_size);
int QIO_skip_data(LRL_RecordReader *lrl_record_in);

DML_SiteList *QIO_create_sitelist(DML_Layout *layout, int volfmt, int serpar);
int QIO_close_read_record(LRL_RecordReader *lrl_record_in);
int QIO_read_sitelist(QIO_Reader *in, LIME_type *lime_type);
int QIO_write_sitelist(QIO_Writer *out, int msg_begin, int msg_end, 
		       const LIME_type lime_type);
LRL_RecordReader *QIO_open_read_field(QIO_Reader *in, size_t datum_size, 
  	       LIME_type *lime_type_list, int ntypes,
               LIME_type *lime_type, int *status);
int QIO_read_field(QIO_Reader *in, 
	   void (*put)(char *buf, size_t index, int count, void *arg),
	   int count, size_t datum_size, int word_size, void *arg, 
	   DML_Checksum *checksum, uint64_t* nbytes,
	   LIME_type *lime_type);
int QIO_read_field_data(QIO_Reader *in, LRL_RecordReader *lrl_record_in,
	   void (*put)(char *buf, size_t index, int count, void *arg),
	   int count, size_t datum_size, int word_size, void *arg, 
 	   DML_Checksum *checksum, uint64_t* nbytes);
int QIO_write_field_data(QIO_Writer *out, LRL_RecordWriter *lrl_record_out,
	    void (*get)(char *buf, size_t index, int count, void *arg),
	    int count, size_t datum_size, int word_size, void *arg, 
	    DML_Checksum *checksum, uint64_t *nbytes);
int QIO_write_field(QIO_Writer *out, int msg_begin, int msg_end,
	    void (*get)(char *buf, size_t index, int count, void *arg),
	    int count, size_t datum_size, int word_size, void *arg, 
	    DML_Checksum *checksum, uint64_t *nbytes,
	    const LIME_type lime_type);
#ifdef __cplusplus
}
#endif

#endif /* QIO_H */
