/* QIO Utilities */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <sys/time.h>

static int QIO_verbosity_level = QIO_VERB_OFF;

double QIO_time (void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double s = (double) tv.tv_sec;
  double u = (double) tv.tv_usec;
  return s + u*1.e-6;
}

/* Wait for "sec" seconds */
void QIO_wait (double sec){
  double finish = QIO_time() + sec;
  while(QIO_time() < finish);
}

/* Set verbosity level. */
int QIO_verbose (int level)
{
  int old = QIO_verbosity_level;
  QIO_verbosity_level = level;
  return old;
}

/* Check verbosity level */
int QIO_verbosity(){
  return QIO_verbosity_level;
}

/*------------------------------------------------------------------*/

/* In case of multifile format we use a common file name stem and add
   a suffix that depends on the node number */

char *QIO_filename_edit(const char *filename, int volfmt, int this_node){

  /* Caller must clean up returned filename */
  int n = strlen(filename) + 12;
  char *newfilename = (char *)malloc(n);

  if(!newfilename){
    printf("QIO_filename_edit: Can't malloc newfilename\n");
    return NULL;
  }

  /* No change for singlefile format */
  if(volfmt == QIO_SINGLEFILE){
    strncpy(newfilename,filename,strlen(filename));
    newfilename[strlen(filename)] = '\0';
  }
  /* Add volume suffix for multifile and partfile formats */
  else{
    snprintf(newfilename,n,"%s.vol%04d",filename,this_node);
  }
  return newfilename;
}


/*------------------------------------------------------------------*/

/* Write an XML record */

int QIO_write_string(QIO_Writer *out, 
		     int msg_begin, int msg_end,
		     QIO_String *xml,
		     const LIME_type lime_type)
{
  LRL_RecordWriter *lrl_record_out;
  char *buf;
  uint64_t actual_rec_size, planned_rec_size;
  char myname[] = "QIO_write_string";

  buf = QIO_string_ptr(xml);
  planned_rec_size = strlen(buf)+1;  /* Include terminating null */

  lrl_record_out = LRL_open_write_record(out->lrl_file_out, 
					 msg_begin,
					 msg_end, planned_rec_size, 
					 lime_type);
  actual_rec_size = LRL_write_bytes(lrl_record_out, buf, planned_rec_size);
  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): wrote %llu bytes\n",myname,out->layout->this_node,
	   (unsigned long long)actual_rec_size);fflush(stdout);

  /* Check byte count */
  if(actual_rec_size != planned_rec_size){
    printf("%s(%d): bytes written %llu != planned rec_size %llu\n",
	   myname, out->layout->this_node, 
	   (unsigned long long)actual_rec_size, 
	   (unsigned long long)planned_rec_size);
    return QIO_ERR_BAD_WRITE_BYTES;
  }
  LRL_close_write_record(lrl_record_out);
  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): closed string record\n",myname,out->layout->this_node);
    fflush(stdout);
  }
  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/
/* Create list of sites expected for this node */

DML_SiteList *QIO_create_sitelist(DML_Layout *layout, int volfmt, int serpar){
  DML_SiteList *sites;
  char myname[] = "QIO_create_sitelist";

  /* Initialize sitelist structure */
  sites = DML_init_sitelist(volfmt, serpar, layout);
  if(sites == NULL){
    printf("%s(%d): Error creating the sitelist structure\n", 
	   myname,layout->this_node);fflush(stdout);
    return sites;
  }

  /* Populate the sitelist */
  /* unless we are reading in discovery mode */
  if(!layout->discover_dims_mode){
    if(DML_fill_sitelist(sites, volfmt, serpar, layout)){
      printf("%s(%d): Error building the site list\n", 
	     myname,layout->this_node);fflush(stdout);
      DML_free_sitelist(sites);
      return NULL;
    }
  }

  return sites;
}

/*------------------------------------------------------------------*/
/* Write list of sites (used with multifile and partitioned file formats) */
/* Returns number of bytes written */

int QIO_write_sitelist(QIO_Writer *out, int msg_begin, int msg_end, 
		       const LIME_type lime_type){
  LRL_RecordWriter *lrl_record_out;
  uint64_t nbytes;
  size_t rec_size;
  DML_SiteRank *outputlist;
  DML_SiteList *sites = out->sites;
  int volfmt = out->volfmt;
  int this_node = out->layout->this_node;
  char myname[] = "QIO_write_sitelist";

  /* Quit if we aren't writing the sitelist */
  /* Single file writes no sitelist.  Multifile always writes one.
     Partitioned file writes only if an I/O node */

  if(volfmt == QIO_SINGLEFILE)return 0;
  if(volfmt == QIO_PARTFILE)
    if(this_node != out->layout->ionode(this_node))return 0;

  /* Make a copy in case we have to byte reverse */
  rec_size = sites->number_of_io_sites * sizeof(DML_SiteRank);
  if(this_node == DML_master_io_node() && QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d) allocating %llu for output sitelist\n",myname,this_node,
	   (unsigned long long)rec_size);fflush(stdout);
  }

  outputlist = (DML_SiteRank *)malloc(rec_size);
  if(outputlist == NULL){
    printf("%s(%d) no room for output sitelist\n",myname,this_node);
    fflush(stdout);
    return QIO_ERR_ALLOC;
  }

  memcpy(outputlist, sites->list, rec_size);

  /* Byte reordering for entire sitelist */
  if (! DML_big_endian())
    DML_byterevn((char *)outputlist, rec_size, sizeof(DML_SiteRank));

  /* Write site list */
  lrl_record_out = LRL_open_write_record(out->lrl_file_out, msg_begin,
					 msg_end, rec_size, lime_type);
  nbytes = LRL_write_bytes(lrl_record_out, (char *)outputlist, rec_size);

  if(nbytes != rec_size){
    printf("%s(%d): Error writing site list. Wrote %llu bytes expected %lu\n", 
	   myname,out->layout->this_node,
	   (unsigned long long)nbytes,(unsigned long)rec_size);
    free(outputlist);
    return QIO_ERR_BAD_WRITE_BYTES;
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): wrote sitelist\n", myname, out->layout->this_node);
  
  /* Close record when done and clean up*/
  LRL_close_write_record(lrl_record_out);

  free(outputlist); 
  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/

LRL_RecordWriter *
QIO_open_write_field(QIO_Writer *out, 
		     int msg_begin, int msg_end, size_t datum_size,
		     const LIME_type lime_type, int *do_output, int *status)
{
  LRL_RecordWriter *lrl_record_out = NULL;
  uint64_t planned_rec_size;
  int this_node = out->layout->this_node;
  int recordtype = out->layout->recordtype;
  int volfmt = out->volfmt;
  int serpar = out->serpar;
  DML_SiteList *sites = out->sites;
  char myname[] = "QIO_open_write_field";

  /* Compute record size */
  if(recordtype == QIO_GLOBAL){
    /* Global data */
    planned_rec_size = datum_size;
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): global data: size %lu\n",myname,this_node,
	     (unsigned long)datum_size);
    }
  }
  else{
    /* Create list of sites in subset for output and count them */
    if(DML_create_subset_rank(out->sites, out->layout, volfmt, serpar) == 1){
      printf("%s(%d) No room for subset rank list\n",
	     myname,this_node);
      return NULL;
    }
    /* Field or subset data */
    if(out->serpar == QIO_SERIAL)
      /* Serial output.  Record size equals the size we write. */
      planned_rec_size = ((uint64_t)sites->subset_io_sites) * datum_size;
    else
      /* Parallel output.  Record size equals the total volume
	 NOTE: If we later decide to write partitions in parallel,
	 this has to be changed to the size for the partition. */
      planned_rec_size = ((uint64_t)out->layout->subsetvolume) * datum_size;
    
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): field data: sites %lu datum %lu\n",
	     myname,this_node,
	     (unsigned long)sites->subset_io_sites,
	     (unsigned long)datum_size);
    }
  }
  
  /* For global data only the master node opens and writes the record.
     Otherwise, all nodes process output, even though only some nodes
     actually write */
  *do_output = (recordtype == QIO_FIELD) || (recordtype == QIO_HYPER) ||
    (this_node == out->layout->master_io_node);
  
  /* Open record only if we have a file handle and are writing */
  /* Nodes that do not write to a file will have a NULL file handle */
  
  if(!out->lrl_file_out || !(*do_output)) {
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("%s(%d): skipping LRL_open_write_record\n",
	     myname,this_node);
  } else {
    /* For serial output, the io_nodes open their records */
    /* For parallel output, only the master node opens its record */
    if( ( out->serpar == DML_SERIAL && 
	  this_node == out->layout->ionode(this_node) ) ||
	( out->serpar == DML_PARALLEL &&
	  this_node == out->layout->master_io_node ) ) {
      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): calling LRL_open_write_record size %llu\n",
	       myname,this_node,(unsigned long long)planned_rec_size);
      lrl_record_out = 
	LRL_open_write_record(out->lrl_file_out, msg_begin, msg_end, 
			      planned_rec_size, lime_type);
      if(lrl_record_out == NULL){
	*status = QIO_ERR_OPEN_WRITE;
	return NULL;
      }
    }
    /* For field data other nodes just create the reader without
       writing a header */
    else if( recordtype == QIO_FIELD || recordtype == QIO_HYPER ) {
      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): calling LRL_create_record_header\n",
	       myname,this_node);
      lrl_record_out = LRL_create_record_writer(out->lrl_file_out);
    }
  }

  /* Then, if we are writing a field in parallel mode, we have to
     synchronize the writers */
  // all nodes must synchronize
  if(out->serpar == DML_PARALLEL && (recordtype == QIO_FIELD ||
				     recordtype == QIO_HYPER)) {
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("%s(%d): Doing parallel write sync\n",myname,this_node);
    if(DML_synchronize_out(lrl_record_out,out->layout) != 0) {
      printf("%s(%d): DML_synchronize returns error\n",
	     myname,this_node);
      *status = QIO_ERR_OPEN_WRITE;
      return NULL;
    }
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): finished\n",myname,this_node);

  *status = QIO_SUCCESS;
  return lrl_record_out;
}  


/*------------------------------------------------------------------*/

/* The next three procedures are available for the API and allow
   random access writing to the binary payload of a lattice field
   but not of global data.

   The first opens the record and initializes the data movement.
   The second is called for each site datum.
   The third closes the record. 

   To write the whole field at once, use QIO_write_field instead,
   which is accessed through QIO_generic_write.

*/

/* Initialize a binary field record for writing data site-by-site */
/* Opens the field and initializes data movement */

int QIO_init_write_field(QIO_Writer *out, int msg_begin, int msg_end,
	    size_t datum_size, DML_Checksum *checksum,
	    const LIME_type lime_type){
  
  LRL_RecordWriter *lrl_record_out;
  DML_RecordWriter *dml_record_out;
  int this_node = out->layout->this_node;
  int do_output;
  int status;
  char myname[] = "QIO_init_write_field";

  /* NOTE: we aren't currently returning do_output */
  lrl_record_out = QIO_open_write_field(out, msg_begin, msg_end, datum_size,
					lime_type, &do_output, &status);

  if(lrl_record_out == NULL)
    return status;

  /* Next we initialize the DML engine */

  dml_record_out = DML_partition_open_out(lrl_record_out,
	  datum_size, 1, out->layout, out->sites, out->volfmt, 
	  out->serpar, checksum);

  if(dml_record_out == NULL)
    {
      printf("%s(%d): Open record failed\n",myname,this_node);
      return QIO_ERR_OPEN_WRITE;
    }

  out->dml_record_out = dml_record_out;

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): finished\n",myname,this_node);
  return QIO_SUCCESS;
}  


/*------------------------------------------------------------------*/

/* Random access write.

   Currently used only in file format conversion.

   Write a single site's data to a location in the binary payload
   specified by a site rank parameter.  The record must first be
   initialized with QIO_init_write_field and closed with
   QIO_close_write_field

*/

int QIO_seek_write_field_datum(QIO_Writer *out, 
	      DML_SiteRank snd_coords,
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      int count, size_t datum_size, int word_size, void *arg)
{
  DML_RecordWriter *dml_record_out = out->dml_record_out;
  int this_node                    = out->layout->this_node;
  int status;
  char myname[] = "QIO_seek_write_site_data";

  status = DML_partition_sitedata_out(dml_record_out, get, snd_coords, 
	      count, datum_size, word_size, arg, out->layout, out->sites);

  if(status != QIO_SUCCESS){
    printf("%s(%d): Error writing site datum\n",myname,this_node);
    return status;
  }

  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/

int QIO_close_write_field(QIO_Writer *out, uint64_t *nbytes)
{
  DML_RecordWriter *dml_record_out = out->dml_record_out;
  LRL_RecordWriter *lrl_record_out = dml_record_out->lrl_rw;

  /* Copy most recent node checksum into writer */
  out->last_checksum = *(dml_record_out->checksum);
  *nbytes = DML_partition_close_out(dml_record_out);
  out->dml_record_out = NULL;

  /* Close record when done and clean up*/
  if(out->lrl_file_out)
    LRL_close_write_record(lrl_record_out);

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("QIO_close_write_field(%d): finished\n",out->layout->this_node);

  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/
/* Write the whole binary lattice field or global data at once */

int QIO_write_field_data(QIO_Writer *out, LRL_RecordWriter *lrl_record_out,
	    void (*get)(char *buf, size_t index, int count, void *arg),
	    int count, size_t datum_size, int word_size, void *arg, 
	    DML_Checksum *checksum, uint64_t *nbytes)
{
  int this_node = out->layout->this_node;
  int recordtype = out->layout->recordtype;
  char myname[] = "QIO_write_field_data";

  /* Initialize byte count and checksum */
  *nbytes = 0;
  DML_checksum_init(checksum);

  /* Write all bytes */

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): starting DML call\n", myname,this_node);

  /* Global data type. */
  if(recordtype == DML_GLOBAL){
    *nbytes = DML_global_out(lrl_record_out, get, count, datum_size, word_size,
			     arg, out->layout, out->volfmt, checksum);
  }

  /* Lattice field data type */
  else {
    *nbytes = DML_partition_out(lrl_record_out, get, count, datum_size, 
				word_size, 
				arg, out->layout, out->sites, out->volfmt, 
				out->serpar, checksum);
    if(out->serpar==QIO_PARALLEL) DML_sync();
  }

  /* Close record when done and clean up*/
  if(out->lrl_file_out)
    LRL_close_write_record(lrl_record_out);

  return QIO_SUCCESS;
}


/*------------------------------------------------------------------*/
/* Write the whole binary lattice field or global data at once */

int QIO_write_field(QIO_Writer *out, int msg_begin, int msg_end,
	    void (*get)(char *buf, size_t index, int count, void *arg),
	    int count, size_t datum_size, int word_size, void *arg, 
	    DML_Checksum *checksum, uint64_t *nbytes,
	    const LIME_type lime_type){
  
  LRL_RecordWriter *lrl_record_out = NULL;
  int recordtype = out->layout->recordtype;
  int do_output = 0;
  int status = QIO_SUCCESS;

  lrl_record_out = QIO_open_write_field(out, msg_begin, msg_end, 
		       datum_size, lime_type, &do_output, &status);

  if(status != QIO_SUCCESS)
    return status;

  /* Initialize byte count and checksum */
  *nbytes = 0;
  DML_checksum_init(checksum);
  status = QIO_SUCCESS;

  /* Write data and close LRL record writer */
  if(do_output)
    status = QIO_write_field_data(out, lrl_record_out, 
			 get, count, datum_size, word_size, arg,
			 checksum, nbytes);

  /* Clean up */
  if(recordtype != QIO_GLOBAL)
    DML_destroy_subset_rank(out->sites);

  return status;
}


/*------------------------------------------------------------------*/
/* Open a record and get the LIME type and expected record size */

LRL_RecordReader *QIO_read_record_type(QIO_Reader *in, LIME_type *lime_type,
			       uint64_t *expected_rec_size, int *status){
  LRL_RecordReader *lrl_record_in;
  int lrl_status;

  if(!in->lrl_file_in){
    *status = QIO_SUCCESS;
    return NULL;
  }

  /* Open record and find type and record size */
  lrl_record_in = LRL_open_read_record(in->lrl_file_in, expected_rec_size, 
				       lime_type, &lrl_status);
  if(!lrl_record_in){
    if(lrl_status == LRL_EOF)*status = QIO_EOF;
    else *status = QIO_ERR_OPEN_READ;
  }
  else
    *status = QIO_SUCCESS;
  
  return lrl_record_in;
}

/*------------------------------------------------------------------*/
/* Skip ahead to the record with one of the specified LIME types and
   read just the LIME header information */

LRL_RecordReader *QIO_open_read_target_record(QIO_Reader *in, 
    LIME_type *lime_type_list, int ntypes, LIME_type *lime_type,
    uint64_t *expected_rec_size, int *status){
  LRL_RecordReader *lrl_record_in;
  int lrl_status;

  /* Open record and find record size */
  if(!in->lrl_file_in){
    *status = QIO_SUCCESS;
    return NULL;
  }

  lrl_record_in = LRL_open_read_target_record(in->lrl_file_in,
		      lime_type_list, ntypes,
		      expected_rec_size, lime_type, &lrl_status);

  if(!lrl_record_in){
    if(lrl_status == LRL_EOF)*status = QIO_EOF;
    else *status = QIO_ERR_OPEN_READ;
  }
  else
    *status = QIO_SUCCESS;
  
  return lrl_record_in;
}

/*------------------------------------------------------------------*/
/* Read an XML record */

int QIO_read_string(QIO_Reader *in, QIO_String *xml, LIME_type *lime_type){
  LRL_RecordReader *lrl_record_in;
  uint64_t expected_rec_size;
  int status;

  /* Open record and find record size */
  if(!in->lrl_file_in)return QIO_SUCCESS;
  lrl_record_in = LRL_open_read_record(in->lrl_file_in, &expected_rec_size, 
				       lime_type, &status);
  if(!lrl_record_in){
    if(status == LRL_EOF)return QIO_EOF;
    else return QIO_ERR_OPEN_READ;
  }

  status = QIO_read_string_data(in, lrl_record_in, xml, expected_rec_size);

  return status;
}

/*------------------------------------------------------------------*/
/* Read string data from a previously opened record */

int QIO_read_string_data(QIO_Reader *in, LRL_RecordReader *lrl_record_in, 
			 QIO_String *xml,  uint64_t expected_rec_size){
  char *buf;
  size_t buf_size;
  uint64_t actual_rec_size;
  char myname[] = "QIO_read_string_data";

  buf_size = QIO_string_length(xml);   /* The size allocated for the string */
  buf      = QIO_string_ptr(xml);

  /* Realloc if necessary */
  if(expected_rec_size+1 > buf_size){
    QIO_string_realloc(xml,(size_t)expected_rec_size+1);  /* +1 for null termination */
  }
  /* Get this again: Guard against truncation in (size_t) conversion */
  buf_size = QIO_string_length(xml);
  buf      = QIO_string_ptr(xml);

  actual_rec_size = LRL_read_bytes(lrl_record_in, buf, buf_size-1);
  LRL_close_read_record(lrl_record_in);

  if(actual_rec_size != expected_rec_size){
    printf("%s(%d): bytes read %llu != expected rec_size %lu\n",
	   myname, in->layout->this_node, (unsigned long long)actual_rec_size, 
	   (unsigned long)expected_rec_size);
    return QIO_ERR_BAD_READ_BYTES;
  }

  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/
/* Skip the LIME record data and close the record reader */

int QIO_close_read_record(LRL_RecordReader *lrl_record_in){
  int status;

  status = LRL_close_read_record(lrl_record_in);

  if(status == LRL_SUCCESS)return QIO_SUCCESS;
  else return QIO_ERR_SKIP;
}

/*------------------------------------------------------------------*/
/* Read site list */

int
QIO_read_sitelist(QIO_Reader *in, LIME_type *lime_type)
{
  int this_node = in->layout->this_node;
  //int number_of_nodes = in->layout->number_of_nodes;
  int volfmt = in->volfmt;
  /* char myname[] = "QIO_read_sitelist"; */
  int status = QIO_SUCCESS;

  /* SINGLEFILE format has no sitelist */
  if(volfmt == QIO_SINGLEFILE) return QIO_SUCCESS;

  /* Only I/O nodes read and verify the sitelist */
  if((volfmt == QIO_MULTIFILE) || 
     ((volfmt == QIO_PARTFILE) 
      && (this_node == in->layout->ionode(this_node)))){
    /* Time release */
    /* double lapse = 1;
       QIO_wait(this_node*lapse); */
    status = DML_read_sitelist(in->sites, 
			       in->lrl_file_in, in->volfmt, 
			       in->layout, lime_type);
    /* QIO_wait((number_of_nodes - this_node)*lapse); */
  }

  return status;
}

/*------------------------------------------------------------------*/
/* Open field data record */

LRL_RecordReader *QIO_open_read_field(QIO_Reader *in, size_t datum_size, 
  	       LIME_type *lime_type_list, int ntypes,
               LIME_type *lime_type, int *status)
{
  LRL_RecordReader *lrl_record_in = NULL;
  DML_SiteList *sites = in->sites;
  uint64_t announced_rec_size, expected_rec_size;
  int this_node = in->layout->this_node;
  int recordtype = in->layout->recordtype;
  int volfmt = in->volfmt;
  int serpar = in->serpar;
  int do_open, do_read;
  int lrl_status;
  int open_fail, open_eof;
  char myname[] = "QIO_open_read_field";

  *status = QIO_SUCCESS;

  /* Will we read the data?
     For field or hypercube data all nodes with readers will read.
     For global data only the master node reads. */

  do_read = ( in->lrl_file_in && (recordtype == QIO_FIELD ||
				  recordtype == QIO_HYPER) ) || 
    (this_node == in->layout->master_io_node);

  /* Should we open the next record with one of the desired LIME types? */

  if(in->serpar == QIO_SERIAL)
    /* For serial I/O we open the record if we will actually read it. */
    do_open = do_read;
  else
    /* For parallel I/O some nodes will have a file reader.
       For field or hypercube data all nodes read.
       For global data only the master node reads but
       other nodes must also open and seek to maintain 
       record-level synchronization. */
    do_open = (in->lrl_file_in != NULL);
  
  open_fail = open_eof = 0;

  if(!do_open) {
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("%s(%d): skipping LRL_open_read_target_record\n",
	     myname,this_node);
  }
  else{
    lrl_record_in = LRL_open_read_target_record(in->lrl_file_in,
			lime_type_list, ntypes, 
			&announced_rec_size, lime_type, &lrl_status);
    
    /* An EOF condition is acceptable here */
    if(lrl_record_in == NULL){
      open_fail = 1;
      if(lrl_status == LRL_EOF){
	open_eof = 0;
      }
      else{
	open_eof = 1;
      }
    }
  }

  /* Poll all nodes for status of open */

  DML_sum_int(&open_fail); DML_sum_int(&open_eof);

  /* Bail out if open on any node returned an error */

  if(open_fail > 0){
    if(open_eof > 0)*status = QIO_ERR_OPEN_READ;
    else *status = QIO_EOF;
    return NULL;
  }

  /* If we are doing parallel I/O, we have to
     synchronize all the the readers to the master node reader. */

  if(in->serpar == DML_PARALLEL){
    if(DML_synchronize_in(lrl_record_in,in->layout) != 0){
      *status = QIO_ERR_OPEN_READ;
      return NULL;
    }
  }

  /* Create list of sites in subset for output and count them */
  if(recordtype != QIO_GLOBAL){
    if(DML_create_subset_rank(in->sites, in->layout, volfmt, serpar) == 1){
      printf("%s(%d) No room for subset rank list\n",
	     myname,this_node);
      return NULL;
    }
  }

  /* Now check the record size for consistency */

  if( do_read )
  {
    /* Check that the record size matches the expected size of the data */
    if(recordtype == QIO_GLOBAL)
      /* Global data */
      expected_rec_size = datum_size;
    else {
      /* Field data or hypercube data */
      if(in->serpar == QIO_SERIAL)
	/* Serial input. Record size equals the size we actually read */
	expected_rec_size = ((uint64_t)sites->subset_io_sites) * datum_size; 
      else
	/* Parallel input.  Record size equals the total volume
	   NOTE: If we later decide to read partitions in parallel,
	   this has to be changed to the size for the partition. */
	expected_rec_size = ((uint64_t)in->layout->subsetvolume) * datum_size;
    }
    if (announced_rec_size != expected_rec_size){
      printf("%s(%d): rec_size mismatch: found %llu expected %llu\n",
	     myname, this_node, (unsigned long long)announced_rec_size, 
	     (unsigned long long)expected_rec_size);
      open_fail = 1;
    }
  }

  /* Poll all nodes for result of size test */

  DML_sum_int(&open_fail);

  if(open_fail > 0){
    *status = QIO_ERR_BAD_READ_BYTES;
    return NULL;
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): finished\n",myname,this_node);

  return lrl_record_in;
}

/*------------------------------------------------------------------*/

/* The next three procedures are available for the API and allow
   random access reading from the binary payload of a lattice field
   but not of global data.

   The first opens the record and initializes the data movement.
   The second is called for each site datum.
   The third closes the record. 

   To read the whole field at once, use QIO_read_field instead,
   which is accessed through QIO_generic_read_record_data.
*/

/* Initialize a binary field record for reading data site-by-site */
/* Opens the field and initializes data movement */

/* Read binary data for a lattice field */

/* Start reading a field */

int QIO_init_read_field(QIO_Reader *in, size_t datum_size, 
			LIME_type *lime_type_list, int ntypes,
			DML_Checksum *checksum, LIME_type *lime_type)
{
  LRL_RecordReader *lrl_record_in = NULL;
  DML_RecordReader *dml_record_in;
  int status;
  int this_node = in->layout->this_node;
  char myname[] = "QIO_init_read_field";

  lrl_record_in = QIO_open_read_field(in, datum_size, 
              lime_type_list, ntypes, lime_type, &status);

  if(lrl_record_in == NULL){
    printf("%s(%d): QIO_open_read_field failed\n",myname,this_node);
    return QIO_ERR_OPEN_READ;
  }

  dml_record_in = DML_partition_open_in(lrl_record_in,
	  datum_size, 1, in->layout, in->sites, in->volfmt, in->serpar,
	  checksum);

  if(dml_record_in == NULL)
    {
      printf("%s(%d): Open record failed\n",myname,this_node);
      return QIO_ERR_OPEN_READ;
    }

  in->dml_record_in = dml_record_in;

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("%s(%d): finished\n",myname,this_node);
  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/

/* Random access read.

   Read a single site's datum from a location in the binary payload
   specified by a site rank parameter.  The record must first be
   initialized with QIO_init_read_field and closed with
   QIO_close_read_field

*/

int QIO_seek_read_field_datum(QIO_Reader *in, DML_SiteRank rcv_coords,
	      void (*put)(char *buf, size_t index, int count, void *arg),
	      int count, size_t datum_size, int word_size, void *arg)
{

  DML_RecordReader *dml_record_in = in->dml_record_in;
  int this_node                   = in->layout->this_node;
  int status;
  char myname[] = "QIO_seek_read_field_datum";

  status = DML_partition_sitedata_in(dml_record_in, put, rcv_coords, 
				     count, datum_size, word_size, arg, 
				     in->layout, in->sites);

  if(status != 0){
    printf("%s(%d): DML error %d reading site datum\n",myname,this_node,
	   status);
    return QIO_ERR_BAD_READ_BYTES;
  }

  in->read_state = QIO_RECORD_CHECKSUM_NEXT;
  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/

int QIO_close_read_field(QIO_Reader *in, uint64_t *nbytes)
{
  DML_RecordReader *dml_record_in = in->dml_record_in;
  LRL_RecordReader *lrl_record_in = dml_record_in->lrl_rr;

  /* Copy most recent node checksum into reader */
  in->last_checksum = *(dml_record_in->checksum);

  *nbytes = DML_partition_close_in(dml_record_in);
  in->dml_record_in = NULL;

  /* Close record when done and clean up*/
  if(in->lrl_file_in)
    LRL_close_read_record(lrl_record_in);

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("QIO_close_read_field(%d): finished\n",in->layout->this_node);
  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/
/* Read binary data for a previously opened lattice field.  
   Close record when done. */

int QIO_read_field_data(QIO_Reader *in, LRL_RecordReader *lrl_record_in,
	   void (*put)(char *buf, size_t index, int count, void *arg),
	   int count, size_t datum_size, int word_size, void *arg, 
 	   DML_Checksum *checksum, uint64_t* nbytes){

  int this_node = in->layout->this_node;
  int recordtype = in->layout->recordtype;
  char myname[] = "QIO_read_field_data";

  /* Initialize byte count and checksum */
  *nbytes = 0;
  DML_checksum_init(checksum);
  
  /* All nodes process input.  Compute checksum and byte count
     for node*/

  /* Global data */
  if(recordtype == QIO_GLOBAL){
    *nbytes = DML_global_in(lrl_record_in,
                         put, count, datum_size, word_size, arg, in->layout, 
			 in->volfmt, in->layout->broadcast_globaldata, 
			 checksum);
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): done with DML_global_in\n", myname,this_node);
    }
  }

  /* Field data */
  else{
    /* Partition I/O only.  Nodes are assigned to disjoint I/O partitions */
    *nbytes = DML_partition_in(lrl_record_in, 
		       put, count, datum_size, word_size, arg, in->layout, 
		       in->sites, in->volfmt, in->serpar, checksum);
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): done with DML_partition_in\n", myname,this_node);
    }
  }

  /* Close record when done and clean up*/

  if(recordtype != QIO_GLOBAL)
    DML_destroy_subset_rank(in->sites);

  if(lrl_record_in)
    LRL_close_read_record(lrl_record_in);

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): record closed\n", myname,this_node);
  }
    
  return QIO_SUCCESS;
}

/*------------------------------------------------------------------*/
/* Read binary data for a lattice field */

int QIO_read_field(QIO_Reader *in, 
	   void (*put)(char *buf, size_t index, int count, void *arg),
	   int count, size_t datum_size, int word_size, void *arg, 
	   DML_Checksum *checksum, uint64_t* nbytes,
	   LIME_type *lime_type){

  LRL_RecordReader *lrl_record_in = NULL;
  int status;

  lrl_record_in = QIO_open_read_field(in, datum_size, NULL, 0, 
				      lime_type, &status);

  status = QIO_read_field_data(in, lrl_record_in, 
			       put, count, datum_size, word_size, arg,
			       checksum, nbytes);
  return status;
}

