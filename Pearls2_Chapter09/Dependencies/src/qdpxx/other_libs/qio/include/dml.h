#ifndef DML_H
#define DML_H

#include <qio_config.h>
#include <lrl.h>
#include <qio_stdint.h>
#include <stdlib.h>

/* File fragmentation */
#define DML_UNKNOWN   -1
#define DML_SINGLEFILE 0
#define DML_MULTIFILE  1
#define DML_PARTFILE   2

/* Data distribution */
#define DML_FIELD      0
#define DML_GLOBAL     1
#define DML_HYPER      2

/* The following choices are ORed to create the oflag or iflag */

/* Multinode access to the same file */
#define DML_SERIAL     0
#define DML_PARALLEL   1

/* Limits size of read and write buffers (bytes) 2^18 for now */
#ifndef QIO_DML_BUF_BYTES
#define DML_BUF_BYTES  262144
#else
#define DML_BUF_BYTES  QIO_DML_BUF_BYTES
#endif

/* Size of message buffers (bytes).  Doesn't usually need to be bigger
   than the x dimension of the local subvolume times the size of the
   data per site */
#define DML_TBUF_BYTES 65536

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct {
  uint32_t suma;
  uint32_t sumb;
} DML_Checksum;

/* Type for sitelist values */
/* Change only if creating new file format */
typedef uint32_t DML_SiteRank;

/* For collecting and passing layout information */
/* See qio.h for QIO_Layout */
typedef struct {
  /* Data distribution fixed for the entire file */
  int (*node_number)(const int coords[]);
  int (*node_index)(const int coords[]);
  void (*get_coords)(int coords[], int node, const int index);
  int (*num_sites)(int node);
  int *latsize;
  int latdim;
  size_t volume;
  size_t sites_on_node;
  int this_node;
  int number_of_nodes;
  int broadcast_globaldata;
  int discover_dims_mode;

  /* Subset parameters that can change with each record */
  int recordtype;
  int *hyperlower;
  int *hyperupper;
  size_t subsetvolume;

  /* I/O partitions */
  int (*ionode)(int node);
  int master_io_node;
} DML_Layout;

typedef struct {
  /* Constant for the file */
  DML_SiteRank *list;            /* List of sites assigned to file */
  int use_list;
  DML_SiteRank first, current_rank;
  size_t number_of_io_sites;
  int number_of_my_ionodes;

  /* Variable */
  size_t current_index;

  /* Constant for a record */
  int use_subset;
  int *subset_rank;     /* Rank order of sites in subset */
  int subset_io_sites;
} DML_SiteList;


/* For saving the state of DML_partition_out */
typedef struct {
  LRL_RecordWriter *lrl_rw; /* LRL record writer */
  char *outbuf;             /* Allocated output buffer */
  char *buf;                /* Current location in output buffer */
  int *coords;              /* Workspace for coordinates */
  DML_Checksum *checksum;   /* Running checksum for this node */
  int current_node;         /* Current output node */
  int my_io_node;           /* The node to which I write */
  uint64_t nbytes;          /* Running total of bytes for this node */
  size_t isite;             /* Count of sites processed */
  size_t buf_sites;         /* Current site in the output buffer */
  size_t max_buf_sites;     /* Size of output buffer in sites */
  size_t max_dest_sites;    /* Total sites written by this node */
} DML_RecordWriter;

/* For saving the state of DML_partition_in */
typedef struct {
  LRL_RecordReader *lrl_rr; /* LRL record reader */
  char *inbuf;              /* Allocated input buffer */
  int *coords;              /* Workspace for coordinates */
  DML_Checksum *checksum;   /* Running checksum for this node */
  int current_node;         /* Current input node */
  int my_io_node;           /* The node to which I write */
  DML_SiteRank rcv_coords;  /* Current input coordinate rank */
  uint64_t nbytes;          /* Running total of bytes for this node */
  size_t isite;             /* Count of sites processed */
  size_t buf_extract;       /* Number of sites in the input buffer */
  size_t buf_sites;         /* Current site in the input buffer */
  size_t max_buf_sites;     /* Size of input buffer in sites */
  size_t max_send_sites;    /* Total sites to be read by this node */
} DML_RecordReader;

uint64_t DML_stream_out(LRL_RecordWriter *lrl_record_out, int recordtype,
	   void (*get)(char *buf, size_t index, int count, void *arg),
           int count, size_t size, int word_size, void *arg, 
	   DML_Layout *layout, DML_SiteList *sites,
 	   int volfmt, int serpar, DML_Checksum *checksum);

size_t DML_stream_global_out(LRL_RecordWriter *lrl_record_out, 
			     void *buf, 
			     size_t size, int word_size, DML_Layout *layout,
			     DML_Checksum *checksum);

uint64_t DML_stream_in(LRL_RecordReader *lrl_record_in, int recordtype,
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     int count, size_t size, int word_size, void *arg, 
             DML_Layout *layout, DML_SiteList *sites, 
	     int volfmt, int serpar, int broadcast_global, 
             DML_Checksum *checksum);


/* DML internal utilities */

void DML_lex_init(int *dim, int coords[], int latdim, int lower[]);
int DML_lex_next(int *dim, int coords[], int latdim, int lower[], 
  	   int upper[]);
void DML_lex_coords(int coords[], const int latdim, const int latsize[], 
		    const DML_SiteRank rcv_coords);
DML_SiteRank DML_lex_rank(const int coords[], int latdim, int latsize[]);
int *DML_allocate_coords(int latdim, const char *myname, int this_node);
char *DML_allocate_msg(size_t size, char *myname, int this_node);
size_t DML_msg_sizeof(size_t size);
char *DML_msg_datum(char *msg, size_t size);
DML_SiteRank *DML_msg_rank(char *msg, size_t size);
DML_SiteList *DML_init_sitelist(int volfmt, int serpar, DML_Layout *layout);
void DML_free_sitelist(DML_SiteList *sites);
int DML_fill_sitelist(DML_SiteList *sites, int volfmt, int serpar,
		      DML_Layout *layout);
int DML_read_sitelist(DML_SiteList *sites, LRL_FileReader *lrl_file_in,
		      int volfmt, DML_Layout *layout,
		      LIME_type *lime_type);
int DML_compare_sitelists(DML_SiteRank *lista, DML_SiteRank *listb, size_t n);
int DML_insert_subset_data(DML_Layout *layout, int recordtype,
			   int *lower, int *upper, int n);
void DML_checksum_init(DML_Checksum *checksum);
void DML_checksum_accum(DML_Checksum *checksum, DML_SiteRank rank, 
			char *buf, size_t size);
void DML_checksum_combine(DML_Checksum *checksum);
void DML_checksum_peq(DML_Checksum *total, DML_Checksum *checksum);
int DML_create_subset_rank(DML_SiteList *sites, DML_Layout *layout,
			   int volfmt, int serpar);
void DML_destroy_subset_rank(DML_SiteList *sites);
void DML_global_xor(uint32_t *x);
int DML_big_endian(void);
void DML_byterevn(char *buf, size_t size, int word_size);
size_t DML_max_buf_sites(size_t size, int factor);
char *DML_allocate_buf(size_t size, size_t *max_buf_sites);
int DML_write_buf_seek(LRL_RecordWriter *lrl_record_out, 
		       DML_SiteRank seeksite, 
		       char *lbuf, size_t buf_sites, size_t size,
		       uint64_t *nbytes, char *myname, 
		       int this_node);
int DML_write_buf_current(LRL_RecordWriter *lrl_record_out, 
			  char *lbuf, size_t buf_sites, size_t size,
			  uint64_t *nbytes, char *myname, 
			  int this_node);
size_t DML_seek_read_buf(LRL_RecordReader *lrl_record_in, 
			 DML_SiteRank seeksite, size_t size,
			 char *lbuf, size_t *buf_extract, size_t buf_sites, 
			 size_t max_buf_sites, size_t isite, 
			 size_t max_send_sites, 
			 uint64_t *nbytes, char *myname, int this_node,
			 int *err);
size_t DML_read_buf_next(LRL_RecordReader *lrl_record_in, size_t size,
			 char *lbuf, size_t *buf_extract, size_t buf_sites, 
			 size_t max_buf_sites, size_t isite, 
			 size_t max_send_sites, 
			 uint64_t *nbytes, char *myname, int this_node,
			 int *err);
int DML_my_ionode(int volfmt, int serpar, DML_Layout *layout);
           DML_RecordWriter *DML_partition_open_out(
	   LRL_RecordWriter *lrl_record_out, size_t size, 
	   size_t set_buf_sites, DML_Layout *layout, DML_SiteList *sites,
	   int volfmt, int serpar, DML_Checksum *checksum);
int DML_partition_sitedata_out(DML_RecordWriter *dml_record_out,
	   void (*get)(char *buf, size_t index, int count, void *arg),
           DML_SiteRank snd_coords, int count, size_t size, int word_size, 
           void *arg, DML_Layout *layout, DML_SiteList *sites);
int DML_partition_allsitedata_out(DML_RecordWriter *dml_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
 	   DML_Layout *layout, DML_SiteList *sites);
uint64_t DML_partition_close_out(DML_RecordWriter *dml_record_out);
uint64_t DML_partition_out(LRL_RecordWriter *lrl_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
	   DML_Layout *layout, DML_SiteList *sites, int volfmt,
	   int serpar, DML_Checksum *checksum);
size_t DML_global_out(LRL_RecordWriter *lrl_record_out, 
	   void (*get)(char *buf, size_t index, int count, void *arg),
	   int count, size_t size, int word_size, void *arg, 
           DML_Layout *layout, int volfmt, 
	   DML_Checksum *checksum);
uint64_t DML_multifile_out(LRL_RecordWriter *lrl_record_out, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      int count, size_t size, int word_size, void *arg, 
 	      DML_Layout *layout, DML_Checksum *checksum);
uint64_t DML_multifile_in(LRL_RecordReader *lrl_record_in, 
	     DML_SiteRank sitelist[],
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     int count, size_t size, int word_size, void *arg, 
	     DML_Layout *layout, DML_Checksum *checksum);
int DML_synchronize_out(LRL_RecordWriter *lrl_record_out, DML_Layout *layout);
int DML_synchronize_in(LRL_RecordReader *lrl_record_in, DML_Layout *layout);
DML_RecordReader *DML_partition_open_in(LRL_RecordReader *lrl_record_in, 
	  size_t size, size_t set_buf_sites, DML_Layout *layout, 
	  DML_SiteList *sites, int volfmt, int serpar, DML_Checksum *checksum);
int DML_partition_sitedata_in(DML_RecordReader *dml_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  DML_SiteRank rcv_coords, int count, size_t size, int word_size, 
          void *arg, DML_Layout *layout, DML_SiteList *sites);
int DML_partition_allsitedata_in(DML_RecordReader *dml_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  int count, size_t size, int word_size, void *arg, 
	  DML_Layout *layout, DML_SiteList *sites, int volfmt,
	  DML_Checksum *checksum);
uint64_t DML_partition_close_in(DML_RecordReader *dml_record_in);
uint64_t DML_partition_in(LRL_RecordReader *lrl_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  int count, size_t size, int word_size, void *arg, 
	  DML_Layout *layout, DML_SiteList *sites, int volfmt,
	  int serpar, DML_Checksum *checksum);
size_t DML_global_in(LRL_RecordReader *lrl_record_in, 
	  void (*put)(char *buf, size_t index, int count, void *arg),
	  int count, size_t size, int word_size, void *arg, 
          DML_Layout* layout, int volfmt, int broadcast_global,
	  DML_Checksum *checksum);

void DML_broadcast_bytes(char *buf, size_t size, int this_node, int from_node);
void DML_sum_uint64_t(uint64_t *ipt);
void DML_sum_int(int *ipt);
int DML_send_bytes(char *buf, size_t size, int tonode);
int DML_get_bytes(char *buf, size_t size, int fromnode);
int DML_clear_to_send(char *buf, size_t size, int my_io_node, int tonode);
void DML_sync(void);

/* I/O layout */
typedef int (*DML_io_node_t)(const int);
typedef int (*DML_master_io_node_t)(void);

int DML_io_node(const int node);
int DML_master_io_node(void);

uint32_t DML_crc32(uint32_t crc, const unsigned char *buf, uint32_t len);

/* Hacks to be removed */
int DML_grid_route(char *buf, size_t size, int fromnode, int tonode);
int DML_route_bytes(char *buf, size_t size, int fromnode, int tonode);

#ifdef __cplusplus
}
#endif

#endif /* DML_H */

