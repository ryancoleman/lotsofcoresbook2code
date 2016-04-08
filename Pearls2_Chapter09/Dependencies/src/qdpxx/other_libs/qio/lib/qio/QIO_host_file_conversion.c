/* Single processor code for converting between SciDAC SINGLEFILE and
   SciDAC PARTFILE format */

#include <qio_config.h>
#include <qio.h>
#include <dml.h>
#include <qio_string.h>
#include <qioxml.h>
#include <stdio.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

typedef struct
{
  char *data;
  size_t datum_size;
  int word_size;
  size_t volume;
} s_field;

static int QIO_init_scalar_field(s_field *field,int vol,size_t datum_size, 
			  int word_size)
{
  size_t bytes;
  field->datum_size  = datum_size;
  field->word_size   = word_size;
  field->volume      = vol;
  bytes = datum_size*vol;

  field->data = (char *) malloc(bytes);
  if(!field->data){
    printf("QIO_init_scalar_field: Can't malloc field data (%f MB)\n",
	   (float)bytes/1e6);
    return 1;
  }
  return 0;
}

static void QIO_free_scalar_field(s_field *field)
{
  free(field->data);
}

typedef struct
{
  s_field *field;
  int node;
  int master_io_node;
} get_put_arg;

typedef struct
{
  get_put_arg *arg;
  int node;
  int master_io_node;
  QIO_Reader *reader;
} read_seek_arg;

typedef struct
{
  get_put_arg *arg;
  int node;
  int master_io_node;
  QIO_Writer *writer;
} write_seek_arg;

/* Populate a pass-through structure for factory functions */
static void QIO_init_get_put_arg(get_put_arg *arg, s_field *field, int node,
				 int master_io_node)
{
  arg->field = field;
  arg->node = node;
  arg->master_io_node = master_io_node;
}

static void QIO_init_read_seek_arg(read_seek_arg *arg_seek, get_put_arg *arg,
				   QIO_Reader *reader, int node,
				   int master_io_node)
{
  arg_seek->arg            = arg;
  arg_seek->node           = node;
  arg_seek->master_io_node = master_io_node;
  arg_seek->reader         = reader;
}

static void QIO_init_write_seek_arg(write_seek_arg *arg_seek, get_put_arg *arg,
				    QIO_Writer *writer, int node, 
				    int master_io_node)
{
  arg_seek->arg            = arg;
  arg_seek->node           = node;
  arg_seek->master_io_node = master_io_node;
  arg_seek->writer         = writer;
}

/* Convert precision code to bytes */
int QIO_bytes_of_word(char *type)
{
  int value=0;
  switch(*type)
    {
    case 'I':
      value = 4;
      break;
      
    case 'F':
      value = 4;
      break;
      
    case 'D':
      value = 8;
      break;
      
    case 'S':
      value = 4;
      break;
    }
  
  return value;
}


/* my_io_node function for host should be called only for node 0 */
int QIO_host_my_io_node(int node)
{
  return node;
}

/* master I/O node for host */
int QIO_host_master_io_node( void )
{
  return 0;
}

/* Factory functions for readers and writers */

/* Copy a chunk of data of length "datum_size" from the input buffer
   to the field */
void QIO_scalar_put( char *s1 , size_t scalar_index, int count, void *s2 )
{
  get_put_arg *arg = (get_put_arg *)s2;
  s_field *field = arg->field;
  size_t datum_size = field->datum_size;
  /* Since we are processing data one site at a time, our "field" holds
     data for only one site and we ignore the scalar index */
  char *dest = field->data;

  memcpy(dest,s1,datum_size);
}


/* Copy a chunk of data of length "datum_size" from the input buffer
   to the field */
void QIO_scalar_put_global( char *s1 , size_t scalar_index, 
			    int count, void *s2 )
{
  get_put_arg *arg = (get_put_arg *)s2;
  s_field *field = arg->field;

  size_t datum_size = field->datum_size;
  /* For global data we ignore the scalar_index */
  char *dest        = field->data;

  memcpy(dest,s1,datum_size);
}


/* Seek and read one site's worth of data from the input file and copy
   it to the output buffer */
void QIO_scalar_get( char *s1, size_t ionode_index, int count, void *s2 )
{
  read_seek_arg *arg_seek = (read_seek_arg *)s2;
  get_put_arg *arg   = arg_seek->arg;    /* The arg for input */
  int ionode_node    = arg_seek->node; 
  s_field *field_in  = arg->field;
  QIO_Reader *infile = arg_seek->reader;
  size_t datum_size = field_in->datum_size;
  char *src         = field_in->data;
  int word_size     = field_in->word_size;
  int scalar_index;
  int status;

  /* Convert site rank ionode_index to scalar_index */
  scalar_index = QIO_ionode_to_scalar_index(ionode_node,ionode_index);

  /* Read the field at location "scalar_index" and put it into
     field_in using the QIO_scalar_put factory function */
  status = 
    QIO_seek_read_field_datum(infile, scalar_index, QIO_scalar_put,
			      count, datum_size, word_size, 
			      (void *)arg);
  if(status != QIO_SUCCESS){
    printf("QIO_scalar_get: QIO_seek_read_field_datum returned %d\n", status);
  }
  else{
    memcpy(s1,src,datum_size);
  }
}

/* Copy the global data of length "datum_size" from the field to the
   output buffer */
void QIO_scalar_get_global( char *s1 , size_t ionode_index, 
			    int count, void *s2 )
{
  get_put_arg *arg = (get_put_arg *)s2;
  s_field *field     = arg->field;
  int node           = arg->node;
  int master_io_node = arg->master_io_node;
  size_t datum_size = field->datum_size;
  char *src         = field->data;

  /* Copy buffer only for the master node and ignore the ionode_index */
  if(node == master_io_node)
    memcpy(s1,src,datum_size);
}

void QIO_part_get( char *s1 , size_t scalar_index, int count, void *s2 )
{
  get_put_arg *arg = (get_put_arg *)s2;
  s_field *field = arg->field;
  size_t datum_size = field->datum_size;
  /* Since we are processing data one site at a time, our "field" holds
     data for only one site and we ignore the scalar index */
  char *src = field->data;;

  memcpy(s1,src,datum_size);
}


/* Copy a chunk of data of length "datum_size" from the input buffer
   and seek and write it to the single file */
void QIO_part_put( char *s1 , size_t ionode_index, int count, void *s2 )
{
  write_seek_arg *arg_seek = (write_seek_arg *)s2;
  get_put_arg *arg         = arg_seek->arg;
  int ionode_node          = arg_seek->node;
  QIO_Writer *outfile      = arg_seek->writer;
  s_field *field_in        = arg->field;
  size_t datum_size        = field_in->datum_size;
  char *dest               = field_in->data;
  int word_size            = field_in->word_size;
  int scalar_index;
  int status;

  /* Copy the input buffer to field_in */
  memcpy(dest,s1,datum_size);

  /* Convert ionode_index to scalar_index */
  scalar_index = QIO_ionode_to_scalar_index(ionode_node,ionode_index);

  /* Write the field to the host single file at the location
     "scalar_index", getting it from field_in using the QIO_part_get
     factory function */
  status =
    QIO_seek_write_field_datum(outfile, scalar_index, QIO_part_get,
			       count, datum_size, word_size,
			       (void *)arg);
  if(status != QIO_SUCCESS){
    printf("QIO_part_put: QIO_seek_write_field_datum returned %d\n", status);
  }
}


/* Copy a chunk of global data of length "datum_size" from the input buffer
   to the field */
void QIO_part_put_global( char *s1 , size_t ionode_index, int count, void *s2 )
{
  get_put_arg *arg = (get_put_arg *)s2;
  s_field *field = arg->field;
  int node = arg->node;
  int master_io_node = arg->master_io_node;
  size_t datum_size = field->datum_size;
  char *dest = field->data;

  /* Copy buffer only for the master node and ignore the ionode_index */
  if(node == master_io_node)
    memcpy(dest,s1,datum_size);
}


int QIO_set_this_node(QIO_Filesystem *fs, const QIO_Layout *layout, int node)
{
  if ( fs->number_io_nodes < layout->number_of_nodes)
    return fs->io_node[node];
  else 
    return node;
}

/* Append the file name to a possibly node-dependent directory path */
char *QIO_set_filepath(QIO_Filesystem *fs, 
		  const char * const filename, int node)
{
  char *path, *newfilename=NULL;
  int fnlength = strlen(filename);
  int drlength;
  
  if (fs->type == QIO_MULTI_PATH)
    {
      path = fs->node_path[node];
      drlength = strlen(path);
      newfilename = (char *) malloc(fnlength+drlength+2);
      if(!newfilename){
	printf("QIO_set_filepath: Can't malloc newfilename\n");
	return NULL;
      }
      newfilename[0] = '\0';
      if(drlength > 0){
	strncpy(newfilename, path, drlength+1);
	strncat(newfilename, "/", 2);
      }
      strncat(newfilename, filename, fnlength+1);
    }
  else if (fs->type == QIO_SINGLE_PATH)
    {
      newfilename = (char *) malloc(fnlength+1);
      strncpy(newfilename,filename,fnlength+1);
    }
  
  return newfilename;
}


/* Open a partition file for reading and read the file header and sitelist */

static QIO_Reader *QIO_open_read_partfile(int io_node_rank, QIO_Iflag *iflag,
					  const char *filename,
					  QIO_Layout *ionode_layout,
					  const QIO_Layout *layout,
					  QIO_Filesystem *fs)
{
  char *newfilename;
  int volfmt;
  QIO_Reader *infile;
  int status;
  char myname[] = "QIO_open_read_partfile";

  /* Pretend we are the ionode */
  ionode_layout->this_node = 
    QIO_set_this_node(fs,layout,io_node_rank);
  
  /* Set output path according to MULTI/SINGLE PATH flag */
  newfilename = QIO_set_filepath(fs,filename,io_node_rank);
  
  /* Open master ionode file to read */
  infile = QIO_open_read_master(newfilename,ionode_layout,
				iflag,fs->my_io_node,fs->master_io_node);
  if(infile == NULL)return NULL;

  /* Check the volume format */
  volfmt = QIO_get_reader_volfmt(infile);

  if (volfmt != QIO_PARTFILE){
    printf("%s(%d) File %s volume format must be PARTFILE.  Found %d\n",
	   myname, io_node_rank, newfilename, QIO_get_reader_volfmt(infile));
    return NULL;
  }

  /* Open nonmaster ionode file to read */
  status = QIO_open_read_nonmaster(infile, newfilename, iflag);
  if(status != QIO_SUCCESS)return NULL;
  
  /* Read site list from master ionode file */
  status = QIO_read_check_sitelist(infile);
  if(status != QIO_SUCCESS)return NULL;

  free(newfilename);
  return infile;
}


/* Open a partition file for writing */

static QIO_Writer *QIO_open_write_partfile(int io_node_rank, QIO_Oflag *oflag,
					   int volfmt, const char *filename,
					   QIO_Layout *ionode_layout,
					   const QIO_Layout *layout,
					   QIO_Filesystem *fs)
{
  char *newfilename;
  QIO_Writer *outfile;

  /* Pretend we are the ionode */
  ionode_layout->this_node = QIO_set_this_node(fs,layout,io_node_rank);
  
  /* Set output path according to MULTI/SINGLE PATH flag */
  newfilename = QIO_set_filepath(fs,filename,io_node_rank);
  
  /* Open to write by appending or with truncation if the file exists */
  outfile = QIO_generic_open_write(newfilename,volfmt,
				   ionode_layout,oflag,
				   fs->my_io_node,fs->master_io_node);
  return outfile;
}


/********************************************************************/
/*  Convert SINGLEFILE to PARTFILE                                  */
/********************************************************************/

int QIO_single_to_part( const char filename[], QIO_Filesystem *fs,
			QIO_Layout *layout)
{
  QIO_Layout *scalar_layout, *ionode_layout;
  QIO_String *xml_file_in, *xml_record_in;
  QIO_String *xml_file_out, *xml_record_out;
  QIO_Reader *infile;
  QIO_Writer *outfile;
  QIO_RecordInfo rec_info;
  QIO_Oflag oflag;
  DML_Checksum checksum, checksum_out, checksum_in;
  uint64_t nbytes_out, totnbytes_out, nbytes_in, totnbytes_in;
  int *msg_begin, *msg_end;
  int i,status,master_io_node_rank;
  int number_io_nodes = fs->number_io_nodes;
  int master_io_node = fs->master_io_node();
  uint64_t total_bytes;
  size_t datum_size;
  int typesize,datacount,recordtype,word_size,volfmt;
  int ntypes = 2;
  LIME_type lime_type_list[2] = {
    QIO_LIMETYPE_BINARY_DATA,
    QIO_LIMETYPE_ILDG_BINARY_DATA
  };
  LIME_type lime_type = NULL;
  s_field field_in;
  get_put_arg arg;
  read_seek_arg arg_seek;
  QIO_ChecksumInfo *checksum_info_expect;
  char myname[] = "QIO_single_to_part";
 
  /* Default values */
  oflag.mode = QIO_TRUNC;
  oflag.serpar = QIO_SERIAL;
  oflag.ildgstyle = QIO_ILDGNO;
  oflag.ildgLFN = NULL;
  
  if(number_io_nodes <= 1){
   printf("%s: No conversion since number_io_nodes %d <= 1\n",
	  myname,number_io_nodes);
   return 1;
  }
 
  /* Create the file XML */
  xml_file_in = QIO_string_create();
  
  /* Create the MPP io_node layout structure */
  ionode_layout = QIO_create_ionode_layout(layout, fs);
  
  /* Create the scalar host layout structure */
  scalar_layout = QIO_create_scalar_layout(layout, fs);
  
  /* Which entry in the table is the MPP master ionode? */
  master_io_node_rank = QIO_get_io_node_rank(master_io_node);
  if(master_io_node_rank < 0){
    printf("%s: Bad Filesystem structure.  Master node %d is not an I/O node\n",
	   myname,master_io_node);
    return QIO_ERR_BAD_IONODE;
  }

  /* Open the input master file for reading */
  infile = QIO_open_read_master(filename, scalar_layout, NULL,
				QIO_host_my_io_node,
				QIO_host_master_io_node);
  if(infile == NULL)return QIO_ERR_OPEN_READ;
  
  if (QIO_get_reader_volfmt(infile) != QIO_SINGLEFILE)
    {
      printf("%s: File %s format %d is not SINGLEFILE\n",
	     myname, filename, QIO_get_reader_volfmt(infile));
      return QIO_ERR_BAD_VOLFMT;
    }
  
  /* Read site list */
  status = QIO_read_check_sitelist(infile);
  if(status != QIO_SUCCESS)return status;

  /* Read user file XML from input master */
  status = QIO_read_user_file_xml(xml_file_in, infile);
  if(status != QIO_SUCCESS)return status;
  
  /* Copy user file XML for eventual output */
  xml_file_out = QIO_string_create();
  QIO_string_copy(xml_file_out, xml_file_in);
  
  /* Set the output volfmt */
  volfmt = QIO_PARTFILE;
  
  /* Make space for message flags */
  msg_begin = (int *)malloc(sizeof(int)*number_io_nodes);
  msg_end   = (int *)malloc(sizeof(int)*number_io_nodes);
  if(!msg_begin || !msg_end){
    printf("%s: Can't malloc msg_begin, msg_end\n", myname);
    return QIO_ERR_ALLOC;
  }
  
  /* Open all partition files and write the file header information
     and site list.  Then close for now */
  for (i=0; i < number_io_nodes; i++)
    {
      /* Open the partition file for writing */
      oflag.mode = QIO_TRUNC;
      outfile = QIO_open_write_partfile(i, &oflag, volfmt, filename,
					ionode_layout, layout, fs);
      if(outfile == NULL)return QIO_ERR_OPEN_WRITE;
	
      /* Write the appropriate file header including site list */
      status = QIO_write_file_header(outfile, xml_file_out);
      if(status != QIO_SUCCESS){
	printf("%s: Can't write file header to %s part %i\n", myname,
	       filename, i);
	return status;
      }

      /* Close the file for now */
      QIO_close_write(outfile);
    }
  
  /***** iterate on field/hypercube/globaldata records up to EOF *******/
  /* We assume the canonical order of SciDAC records */
  while (1)
    {
      
      /* Initialize message flags for output records */
      for(i = 0; i < number_io_nodes; i++){
	msg_begin[i] = 1; msg_end[i] = 0;
      }

      /* Read the private record info or if EOF, quit the loop */
      status = QIO_read_private_record_info(infile, &rec_info);
      if (status==QIO_EOF)break;
      if (status!=QIO_SUCCESS) return status;

      /* Parse the record info */
      datacount  = QIO_get_datacount(&rec_info);
      typesize   = QIO_get_typesize(&rec_info);
      word_size  = QIO_bytes_of_word(QIO_get_precision(&rec_info));
      recordtype = QIO_get_recordtype(&rec_info);
      datum_size = typesize*datacount;
      
      /* Create and read the input user record XML */
      xml_record_in = QIO_string_create();
      status = QIO_read_user_record_info(infile, xml_record_in);
      if (status!=QIO_SUCCESS) return status;
      
      /* Read the ILDG LFN if present */
      status = QIO_read_ILDG_LFN(infile);
      if (status!=QIO_SUCCESS) return status;

      /* If ILDG style, switch output to ILDG style */
      oflag.ildgstyle = QIO_get_ildgstyle(infile);

      /* Copy the ILDG LFN if we have it */
      if(QIO_get_ILDG_LFN(infile) != NULL)
	if(strlen(QIO_get_ILDG_LFN(infile)) > 0)
	  {
	    oflag.ildgLFN = QIO_string_create();
	    QIO_string_set(oflag.ildgLFN, QIO_get_ILDG_LFN(infile));
	    QIO_reset_writer_ILDG_flags(outfile, &oflag);
	  }

      /* Reopen the master ionode file */
      oflag.mode = QIO_APPEND;
      outfile = QIO_open_write_partfile(master_io_node, &oflag, volfmt, 
					filename, ionode_layout, layout, fs);
      if(outfile == NULL)return QIO_ERR_OPEN_WRITE;
      
      /* Copy the user record XML */
      xml_record_out = QIO_string_create();
      QIO_string_copy(xml_record_out, xml_record_in);
      
      /* Write the record header to the master_io_node file */
      /* This includes the private record XML and user record XML */
      status = QIO_write_record_info(outfile, &rec_info, datum_size, word_size,
				     xml_record_out, 
				     &msg_begin[master_io_node_rank], 
				     &msg_end[master_io_node_rank]);
      if(status != QIO_SUCCESS)return status;

      /* Next we process the payload.  Processing depends on whether
	 we have global data or a field */

      if( recordtype == QIO_GLOBAL){

	/* Read global data in its entirety and write it all to the
	   the master_io_node file */

	/* Allocate space for the global data */
	status = QIO_init_scalar_field(&field_in,1,datum_size,word_size);
	if(status != QIO_SUCCESS)return status;

	/* Prepare to read */
	QIO_init_get_put_arg(&arg, &field_in, QIO_host_my_io_node(0),
			     QIO_host_master_io_node());
	/* Read the data from the host file */
	QIO_suppress_global_broadcast(infile);  /* Scalar operation here */
	status = 
	  QIO_generic_read_record_data(infile,QIO_scalar_put_global,datum_size,
				       word_size,&arg, &checksum_in, 
				       &nbytes_in);
	if(status != QIO_SUCCESS)return status;

	/* Expected output byte count */
	total_bytes = datum_size;
	totnbytes_in = nbytes_in;

	/* Prepare to write */
	QIO_init_get_put_arg(&arg, &field_in, ionode_layout->this_node,
			     master_io_node);

	/* Write the global data */
	status = QIO_write_record_data(outfile, &rec_info,
				       QIO_scalar_get_global, 
				       datum_size, word_size, 
				       &arg, &checksum_out, &nbytes_out,
				       &msg_begin[master_io_node_rank], 
				       &msg_end[master_io_node_rank]);
	if(status != QIO_SUCCESS)return status;
	totnbytes_out = nbytes_out;

      } /* global data */

      else{

	/* Write the field or hypercube data. */
	
	/* Prepare the input file for reading the site data via random
	   access */
	
	status = QIO_init_read_field(infile, datum_size, 
				     lime_type_list, ntypes, &checksum_in, 
				     &lime_type);
	if(status != QIO_SUCCESS)return status;
	
	/* Allocate space for the field datum for one site */
	status = QIO_init_scalar_field(&field_in,1,datum_size,word_size);
	if(status != QIO_SUCCESS)return status;
	
	/* Prepare to read */
	QIO_init_get_put_arg(&arg, &field_in, 
			     QIO_host_my_io_node(0),
			     QIO_host_master_io_node());
	
	/* Expected total for the entire field */
	total_bytes = ((uint64_t)infile->layout->subsetvolume) * datum_size;
	totnbytes_out = 0;
	DML_checksum_init(&checksum_out);
	
	/*  Cycle through all the partition files, copying one site at
	    a time */
	
	for(i = 0; i < number_io_nodes; i++){
	  
	  /* Reopen the partition file for appending */
	  oflag.mode = QIO_APPEND;
	  outfile = QIO_open_write_partfile(i, &oflag, 
					    volfmt, filename, 
					    ionode_layout, layout, fs);
	  if(outfile == NULL)return QIO_ERR_OPEN_WRITE;

	  /* Prepare part file output */
	  QIO_init_read_seek_arg(&arg_seek, &arg, infile,
				 ionode_layout->this_node,
				 master_io_node);
	  
	  /* Copy hypercube data from record_info structure to writer */
	  status = QIO_writer_insert_hypercube_data(outfile, &rec_info);
	  
	  if(status != QIO_SUCCESS)return status;
	  
	  /* Write the data.  The factory function QIO_scalar_get
	     seeks and reads from the input file */
	  status = 
	    QIO_write_record_data(outfile, &rec_info, QIO_scalar_get, 
				  datum_size, word_size, &arg_seek, 
				  &checksum, &nbytes_out, 
				  &msg_begin[i], &msg_end[i]);
	  if(status != QIO_SUCCESS)return status;
	  
	  /* Add partial byte count to total output bytes */
	  totnbytes_out += nbytes_out;
	  /* Add partial checksum to total */
	  DML_checksum_peq(&checksum_out, &checksum);

	  /* Close the file for now */
	  QIO_close_write(outfile);
	}
	
	/* Close the input field. (File remains open) */
	status = QIO_close_read_field(infile, &totnbytes_in);
	if(status != QIO_SUCCESS)return status;

	/* Reopen the master_io_node file for writing the checksum */
	oflag.mode = QIO_APPEND;
	outfile = QIO_open_write_partfile(master_io_node, &oflag, volfmt, 
					  filename, ionode_layout, layout, fs);
	if(outfile == NULL)return QIO_ERR_OPEN_WRITE;

      } /* field or hypercube data */
	
      /* Check that input byte count matches total expected */
      if(total_bytes != totnbytes_in){
	printf("%s: Input byte count %llu does not match expected %llu\n",
	       myname,
	       (unsigned long long)totnbytes_in, 
	       (unsigned long long)total_bytes);
	return QIO_ERR_BAD_READ_BYTES;
      }
      
      /* Check that input byte and output byte counts match */
      if(totnbytes_out != totnbytes_in){
	printf("%s: Input byte count %llu does not match output %llu\n",
	       myname,
	       (unsigned long long)totnbytes_in, 
	       (unsigned long long)totnbytes_out);
	return QIO_ERR_BAD_READ_BYTES;
      }
      
      /* Compare checksums */
      if(checksum_in.suma != checksum_out.suma ||
	 checksum_in.sumb != checksum_out.sumb){
	printf("%s Input checksum %0x %0x != output checksum %0x %0x.\n",
	       myname,
	       checksum_in.suma, checksum_in.sumb,
	       checksum_out.suma, checksum_out.sumb);
	return QIO_CHECKSUM_MISMATCH;
      }
      
      /* Read the checksum record and compare with what we got */
      
      checksum_info_expect = QIO_read_checksum(infile);
      if(checksum_info_expect == NULL)return QIO_ERR_CHECKSUM_INFO;
      
      if(QIO_get_reader_format(infile) == QIO_SCIDAC_NATIVE){
	/* Only native SciDAC files have a checksum to compare */
	status = QIO_compare_checksum(0, checksum_info_expect, &checksum_in);
	if(status != QIO_SUCCESS){
	  printf("%s Input data checksum does not match input file checksum\n",
		 myname);
	  return status;
	}
      }
      QIO_destroy_checksum_info(checksum_info_expect);
      
      /* Write checksum record */
      status = QIO_write_checksum(outfile, &checksum_out);
      
      if(QIO_verbosity() >= QIO_VERB_LOW){
	if(recordtype == QIO_GLOBAL)printf("%s: Global data\n",myname);
	else printf("%s: Field data\n",myname);
	printf("  %s\n  Datatype %s\n  precision %s colors %d spins %d count %d\n",
	       QIO_string_ptr(xml_record_in),
	       QIO_get_datatype(&rec_info),
	       QIO_get_precision(&rec_info),
	       QIO_get_colors(&rec_info),
	       QIO_get_spins(&rec_info),
	       QIO_get_datacount(&rec_info));
	
	printf("  Checksums %0x %0x\n",
	       checksum_out.suma, checksum_out.sumb);
      }
      
      QIO_free_scalar_field(&field_in);
      
      /* Close the master_io_node file for now */
      QIO_close_write(outfile);
      
      QIO_string_destroy(xml_record_in);
      QIO_string_destroy(xml_record_out);

      /* If this is a non-native ILDG file, only one data field (the
	 lattice) is permitted, so we bail out here */
      if(QIO_get_reader_format(infile) == QIO_ILDG_ALIEN)break;
    }
  /************* end iteration on records *********/
  
  QIO_string_destroy(xml_file_in);
  
  QIO_delete_scalar_layout(scalar_layout);
  QIO_delete_ionode_layout(ionode_layout);
  
  return QIO_SUCCESS;
}

/********************************************************************/
/*  Convert PARTFILE to SINGLEFILE                                  */
/********************************************************************/

int QIO_part_to_single( const char filename[], int ildgstyle, 
			QIO_String *ildgLFN_override,
			QIO_Filesystem *fs, QIO_Layout *layout)
{
  QIO_Layout *scalar_layout, *ionode_layout;
  QIO_String *xml_file_in, *xml_record_in;
  QIO_String *xml_file_out, *xml_record_out;
  QIO_Reader *infile;
  QIO_Writer *outfile;
  QIO_RecordInfo rec_info, rec_info_in;
  QIO_Iflag iflag;
  QIO_Oflag oflag;
  DML_Checksum checksum_out, checksum_in, checksum;
  QIO_ChecksumInfo *checksum_info_expect=NULL, *checksum_info_tmp;
  int read_format = QIO_SCIDAC_NATIVE;
  uint64_t nbytes_in, nbytes_out, totnbytes_out, totnbytes_in,
    total_bytes;
  int msg_begin, msg_end;
  int i,status,master_io_node_rank;
  int number_io_nodes = fs->number_io_nodes;
  int master_io_node = fs->master_io_node();
  size_t datum_size;
  int typesize,datacount,recordtype,word_size;
  int ntypes = 2;
  LIME_type lime_type_list[2] = {
    QIO_LIMETYPE_BINARY_DATA,
    QIO_LIMETYPE_ILDG_BINARY_DATA
  };
  LIME_type lime_type_out = NULL, lime_type_in;
  off_t *offset;
  s_field field_in;
  get_put_arg arg;
  write_seek_arg arg_seek;
  FILE *check;
  char myname[] = "QIO_part_to_single";

  /* Default values */
  iflag.serpar = QIO_SERIAL;
  iflag.volfmt = QIO_PARTFILE;

  oflag.serpar = QIO_SERIAL;
  oflag.mode = QIO_TRUNC;
  oflag.ildgstyle = ildgstyle;
  oflag.ildgLFN = NULL;

  /* Sanity checks */

  /* No conversion if the number of nodes is not greater than 1 */
  if(number_io_nodes <= 1){
   printf("%s: No conversion since number_io_nodes %d <= 1\n",
	  myname,number_io_nodes);
   return 1;
  }

  /* The single file target must not exist.  Otherwise, because it is
     designed to autodetect the file format, QIO_open_read can't tell
     whether to open it or the input partition file. */
  
  check = DCAPL(fopen)(filename,"r");
  if(check){
    printf("%s: No conversion since the file %s already exists\n",
	   myname, filename);
    DCAP(fclose)(check);
    return 1;
  }

  /* Allocate space for input file offsets */
  /* We use these to mark our place in each file so we can reopen them
     and continue reading from where we left off */
  offset = (off_t *)malloc(number_io_nodes * sizeof(off_t));
  if(offset == NULL){
    printf("%s No space for offsets\n",myname);
    return QIO_ERR_ALLOC;
  }
  /* Indicates we haven't determined the offset, yet */
  for(i = 0; i < number_io_nodes; i++)
    offset[i] = -1;
  
  /* Create the MPP ionode layout structure */
  ionode_layout = QIO_create_ionode_layout(layout, fs);
  if(ionode_layout == NULL)return QIO_ERR_ALLOC;

  /* Create the scalar layout structure */
  scalar_layout = QIO_create_scalar_layout(layout, fs);
  if(scalar_layout == NULL)return QIO_ERR_ALLOC;

  /* Which entry in the table is the MPP master ionode? */
  master_io_node_rank = QIO_get_io_node_rank(master_io_node);
  if(master_io_node_rank < 0){
    printf("%s: Bad Filesystem structure.  Master node %d is not an I/O node\n",
           myname, master_io_node);
    return QIO_ERR_BAD_IONODE;
  }

  /* Open the master ionode file, read private file info and sitelist */
  infile = QIO_open_read_partfile(master_io_node_rank, &iflag, filename, 
				  ionode_layout, layout, fs);
  if(infile == NULL)return QIO_ERR_OPEN_READ;

  /* Read user file XML from input master */
  xml_file_in = QIO_string_create();
  status = QIO_read_user_file_xml(xml_file_in, infile);
  if(status != QIO_SUCCESS)return status;
  
  /* Close the master ionode file temporarily, saving the current
     position, which should be just before the first data message */
  offset[master_io_node_rank] = QIO_get_reader_pointer(infile);
  status = QIO_close_read(infile);
  if(status != QIO_SUCCESS)return status;

  /* Copy user file XML for eventual output */
  xml_file_out = QIO_string_create();
  QIO_string_copy(xml_file_out, xml_file_in);
  
  /* Copy user file XML for output */
  xml_file_out = QIO_string_create();
  QIO_string_copy(xml_file_out, xml_file_in);

  /* Open the host single file and write the header, including user file
     XML */

  outfile =  QIO_generic_open_write(filename,QIO_SINGLEFILE,
				    scalar_layout, &oflag,
				    QIO_host_my_io_node,
				    QIO_host_master_io_node);
  
  if(outfile == NULL)return QIO_ERR_OPEN_WRITE;

  status = QIO_write_file_header(outfile, xml_file_out);
  if(status != QIO_SUCCESS)return status;

  /***** iterate on records up to EOF ***********/
  
  while (1)
    {
      
      /* Reopen the master ionode file, read private file info and sitelist */
      infile = QIO_open_read_partfile(master_io_node_rank, &iflag, filename,
				      ionode_layout, layout, fs);
      if(infile == NULL)return QIO_ERR_OPEN_READ;
      
      /* Skip to where we left off */
      status = QIO_set_reader_pointer(infile,offset[master_io_node_rank]);
      if(status != QIO_SUCCESS)return status;

      /* Read the next record info from the master ionode file.
	 Stop when we reach the end of file. */
      status = QIO_read_private_record_info(infile, &rec_info_in);
      if (status==QIO_EOF) break;
      if (status!=QIO_SUCCESS) return status;
      
      /* Collect record format data */
      datacount   = QIO_get_datacount(&rec_info_in);
      typesize    = QIO_get_typesize(&rec_info_in);
      word_size   = QIO_bytes_of_word(QIO_get_precision(&rec_info_in));
      recordtype  = QIO_get_recordtype(&rec_info_in);
      datum_size  = typesize*datacount;

      /* Read the user record info from the master ionode file. */
      xml_record_in = QIO_string_create();
      status = QIO_read_user_record_info(infile, xml_record_in);
      if(status != QIO_SUCCESS)return status;

      /* Set the output record XML */
      xml_record_out = QIO_string_create();
      QIO_string_copy(xml_record_out, xml_record_in);
      
      /* Read the ILDG LFN from the master ionode file. */
      status = QIO_read_ILDG_LFN(infile);
      if(status != QIO_SUCCESS)return status;
      
      /* If we want ILDG format output and we have specified an output
	 LFN, then use the specified output LFN regardless of the
	 input LFN. */
      if( ildgstyle == QIO_ILDGLAT && 
	  ildgLFN_override != NULL)
	if(QIO_string_length(ildgLFN_override) > 0){
	  oflag.ildgLFN = QIO_string_create();
	  QIO_string_copy(oflag.ildgLFN, ildgLFN_override);
	  QIO_reset_writer_ILDG_flags(outfile, &oflag);
	}
      
      /* If we want ILDG format output and we have an input LFN and 
	 we don't have an output LFN, yet,
         then copy the input LFN to the output LFN */
      if(ildgstyle == QIO_ILDGLAT && 
	 QIO_get_ILDG_LFN(infile) != NULL &&
	 oflag.ildgLFN == NULL)
	if(strlen(QIO_get_ILDG_LFN(infile)) > 0){
	  oflag.ildgLFN = QIO_string_create();
	  QIO_string_set(oflag.ildgLFN, QIO_get_ILDG_LFN(infile));
	}

      /* If we now have an output LFN, give it to the writer */
      if(oflag.ildgLFN != NULL)
	QIO_reset_writer_ILDG_flags(outfile, &oflag);

      /* Set the output record XML and write the private and user
	 record XML to the host single file */
      status = QIO_write_record_info(outfile, &rec_info_in, 
				     datum_size, word_size,
				     xml_record_out,
				     &msg_begin, &msg_end);
      if(status != QIO_SUCCESS)return status;

      /* Process the record according to type: global data or field data */

      if( recordtype == QIO_GLOBAL)
	{
	  /* Global data.  Read the data in its entirety and write it */

	  /* Allocate space for the global data */
	  status = QIO_init_scalar_field(&field_in,1,datum_size,word_size);
	  if(status != QIO_SUCCESS)return status;

	  /* Read the data */
	  QIO_init_get_put_arg(&arg, &field_in, ionode_layout->this_node,
			       master_io_node);
	  QIO_suppress_global_broadcast(infile);  /* Scalar operation here */
	  status = 
	    QIO_generic_read_record_data(infile, QIO_part_put_global,
					 datum_size, word_size,
					 &arg,&checksum_in,&nbytes_in);
	  if(status != QIO_SUCCESS)return status;

	  /* Read checksum data */
	  checksum_info_expect = QIO_read_checksum(infile);
	  if(checksum_info_expect == NULL)return QIO_ERR_CHECKSUM_INFO;
	  
	  /* Expected output byte count */
	  total_bytes = datum_size;
	  totnbytes_in = nbytes_in;
	  
	  /* Write the global data to the host single file */
	  QIO_init_get_put_arg(&arg, &field_in, QIO_host_my_io_node(0),
			       QIO_host_master_io_node());

	  status = QIO_write_record_data(outfile, &rec_info_in,
					 QIO_scalar_get_global, 
					 datum_size, word_size, 
					 &arg, &checksum_out, &nbytes_out,
					 &msg_begin, &msg_end);
	  if(status != QIO_SUCCESS)return status;
	  totnbytes_out = nbytes_out;

	  /* Close the master ionode file temporarily, saving the current
	     position, which should be just before the next data message */
	  offset[master_io_node_rank] = QIO_get_reader_pointer(infile);
	  status = QIO_close_read(infile);
	  if(status != QIO_SUCCESS)return status;
	}
      else
	{
	  /* Write the field or hypercube data */

	  /* We need the LIME type for the record data. So we peek at
	     the header for LIME record for the binary data in the
	     master ionode file.  Also initialize checksum_in.
	     We may need to search ahead for this record, since
	     there may be an intervening ILDG format and LFN record. */
	  
	  status = QIO_init_read_field(infile, datum_size, 
				       lime_type_list, ntypes,
				       &checksum_in, &lime_type_in);
	  if(status != QIO_SUCCESS)return status;

	  /* Expected byte count */
	  total_bytes = ((uint64_t)infile->layout->subsetvolume) * datum_size;

	  /* Copy LIME type */
	  lime_type_out = (char *)malloc(strlen(lime_type_in)+1);
	  strncpy(lime_type_out,lime_type_in,strlen(lime_type_in)+1);

	  /* Now close the master ionode file.  We will reread the
	     private and user file xml later */
	  status = QIO_close_read(infile);
	  if(status != QIO_SUCCESS)return status;

	  /* Prepare to write the output field to the host single file */
	  status = QIO_init_write_field(outfile, msg_begin, msg_end,
					datum_size, &checksum_out,
					lime_type_out);
	  if(status != QIO_SUCCESS)return status;
	  free(lime_type_out);
	  lime_type_out = NULL;

	  /* Input byte counting and checksums */
	  totnbytes_in = 0;
	  DML_checksum_init(&checksum_in);

	  /* Create space for holding the binary data for one site */
	  status = QIO_init_scalar_field(&field_in,1,datum_size,word_size);
	  if(status != QIO_SUCCESS)return status;
      
	  /* Prepare to write */
	  QIO_init_get_put_arg(&arg, &field_in, 
			       QIO_host_my_io_node(0),
			       QIO_host_master_io_node());
	
	  /* Cycle through all the partition files, reading and
	     copying one site at a time.  The factory "put" function
	     writes the site data to the large file. */
	  
	  for(i = 0; i < number_io_nodes; i++){
	  
	    /* Open the ionode file for reading */
	    infile = QIO_open_read_partfile(i, &iflag, filename, 
					    ionode_layout, layout, fs);
	    if(infile == NULL)return QIO_ERR_OPEN_READ;
	    
	    /* Position the file for the next data message */
	    if(offset[i] >= 0)
	      status = QIO_set_reader_pointer(infile,offset[i]);
	    if(status != QIO_SUCCESS)return status;

	    /* Reread the record info or set state */
	    status = QIO_read_private_record_info(infile, &rec_info);
	    if(status != QIO_SUCCESS)return status;

	    /* The nonmaster node files don't have any record info, so
               we copy the master ionode record info into the reader.
	    */
	    QIO_set_record_info(infile, &rec_info_in);
	    
	    /* Copy hypercube data into reader (if any) */
	    status = QIO_reader_insert_hypercube_data(infile, &rec_info_in);
	    if(status != QIO_SUCCESS)return status;

	    /* Reread the user record info or set state */
	    status = QIO_read_user_record_info(infile, xml_record_in);
	    if(status != QIO_SUCCESS)return status;

	    /* Reread the ILDG LFN from the master ionode file. */
	    status = QIO_read_ILDG_LFN(infile);
	    if(status != QIO_SUCCESS)return status;
      
	    /* Prepare host single file output */
	    QIO_init_write_seek_arg(&arg_seek, &arg, outfile,
				    ionode_layout->this_node,
				    master_io_node);
	    
	    /* Read the record data, writing it to the host single
	       file via the factory function QIO_part_put */
	    status = 
	      QIO_generic_read_record_data(infile,QIO_part_put,datum_size,
					   word_size,&arg_seek,
					   &checksum, &nbytes_in);
	    if(status != QIO_SUCCESS)return status;

	    /* Add partial checksum to total */
	    DML_checksum_peq(&checksum_in, &checksum);

	    /* Read checksum data. Only the master ionode has it, so
	       grab it when we process the master ionode file */
	    checksum_info_tmp = QIO_read_checksum(infile);
	    if(i == master_io_node_rank){
	      if(checksum_info_tmp == NULL)return QIO_ERR_CHECKSUM_INFO;
	      checksum_info_expect = checksum_info_tmp;
	      /* Save read format for later */
	      read_format = QIO_get_reader_format(infile);
	    }
	    
	    /* Add partial byte count to total output bytes */
	    totnbytes_in += nbytes_in;

	    /* Close the file temporarily */
	    offset[i] = QIO_get_reader_pointer(infile);
	    status = QIO_close_read(infile);
	    if(status != QIO_SUCCESS)return status;
	  }

	  /* Close the output field (file remains open) */
	  status = QIO_close_write_field(outfile, &totnbytes_out);
	  if(status != QIO_SUCCESS)return QIO_ERR_CLOSE;

	} /* else field */
      
      /* Compare output byte count with expected total record size */
      if(totnbytes_out != total_bytes)
	{
	  printf("%s: bytes written %llu != expected rec_size %llu\n",
		 myname, (unsigned long long)totnbytes_out, 
		 (unsigned long long)total_bytes);
	  return QIO_ERR_BAD_WRITE_BYTES;
	}

      /* Compare input and output byte counts */
      if(totnbytes_in != totnbytes_out){
	  printf("%s: bytes written %llu != bytes read %llu\n",
		 myname, (unsigned long long)totnbytes_in, 
		 (unsigned long long)totnbytes_out);
	  return QIO_ERR_BAD_WRITE_BYTES;
      }
      
      /* Compare input and output checksums */
      if(checksum_in.suma != checksum_out.suma ||
	 checksum_in.sumb != checksum_out.sumb){
	printf("%s Input checksum %0x %0x != output checksum %0x %0x.\n",
	       myname,
	       checksum_in.suma, checksum_in.sumb,
	       checksum_out.suma, checksum_out.sumb);
	return QIO_CHECKSUM_MISMATCH;
      }
      
      /* Compare the computed input checksums with checksum on the
	 input file.  Only for SciDAC native files! */
      if(read_format == QIO_SCIDAC_NATIVE){
	status = QIO_compare_checksum(master_io_node, 
				      checksum_info_expect, &checksum_in);
	if (status != QIO_SUCCESS) return status;
      }
      
      /* Write checksum record to host single file */
      status = QIO_write_checksum(outfile, &checksum_out);
      if(status != QIO_SUCCESS)return status;      

      if(QIO_verbosity() >= QIO_VERB_LOW){
	if(recordtype == QIO_GLOBAL)printf("Global data\n");
	else printf("Field data\n");
	printf("  %s\n  Datatype %s\n  precision %s colors %d spins %d count %d\n",
	       QIO_string_ptr(xml_record_out),
	       QIO_get_datatype(&rec_info),
	       QIO_get_precision(&rec_info),
	       QIO_get_colors(&rec_info),
	       QIO_get_spins(&rec_info),
	       QIO_get_datacount(&rec_info));
	
	printf("  Checksums %0x %0x\n",
	       checksum_out.suma, checksum_out.sumb);
      }
      
      QIO_free_scalar_field(&field_in);
      
      QIO_string_destroy(xml_record_in);
      QIO_string_destroy(xml_record_out);
     
      QIO_destroy_checksum_info(checksum_info_expect);
    }
  /************* end iteration on records *********/

  /* Close the master_io_node file */
  QIO_close_write(outfile);
  
  QIO_string_destroy(xml_file_in);
  
  QIO_delete_scalar_layout(scalar_layout);
  QIO_delete_ionode_layout(ionode_layout);
  
  return QIO_SUCCESS;
}
