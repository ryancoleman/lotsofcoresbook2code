/* QIO_open_write.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <qioxml.h>
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif


/* Accessors for writer information */
uint32_t QIO_get_writer_last_checksuma(QIO_Writer *out){
  return out->last_checksum.suma;
}

uint32_t QIO_get_writer_last_checksumb(QIO_Writer *out){
  return out->last_checksum.sumb;
}

void QIO_reset_writer_ILDG_flags(QIO_Writer *out, QIO_Oflag *oflag){
  if(oflag == NULL)return;
  out->ildgstyle = oflag->ildgstyle;
  out->ildgLFN = oflag->ildgLFN;
}

/* Opens a file for writing, whether in MPP mode or on host */
/* Constructs the site list if multifile format */

QIO_Writer *QIO_generic_open_write(const char *filename, 
				   int volfmt, QIO_Layout *layout, 
				   QIO_Oflag *oflag, 
				   DML_io_node_t io_node, 
				   DML_master_io_node_t master_io_node)
{
  QIO_Writer *qio_out = NULL;
  LRL_FileWriter *lrl_file_out = NULL;
  DML_Layout *dml_layout;
  int *latsize, *upper, *lower;
  int latdim = layout->latdim;
  int this_node = layout->this_node;
  int i;
  int mode;
  int serpar;
  char *newfilename;
  char myname[] = "QIO_generic_open_write";

  /* Make a local copy of lattice size */
  latsize = (int *)malloc(sizeof(int)*latdim);
  for(i=0; i < latdim; ++i)
    latsize[i] = layout->latsize[i];

  /* Set up hypercube subset dimensions in case they are needed */
  upper = (int *)malloc(sizeof(int)*latdim);
  lower = (int *)malloc(sizeof(int)*latdim);

  /* Construct the layout data from the QIO_Layout structure*/
  dml_layout = (DML_Layout *)malloc(sizeof(DML_Layout));
  if (!layout){
    printf("%s(%d): can't malloc dml_layout\n",myname,this_node);
    return NULL;
  }
  
  /* Force single file format if there is only one node */
  if(layout->number_of_nodes==1) volfmt = QIO_SINGLEFILE;
  
  dml_layout->node_number          = layout->node_number;
  dml_layout->node_index           = layout->node_index;
  dml_layout->get_coords           = layout->get_coords;
  dml_layout->num_sites            = layout->num_sites;
  dml_layout->latsize              = latsize;
  dml_layout->latdim               = layout->latdim;
  dml_layout->volume               = layout->volume;
  dml_layout->sites_on_node        = layout->sites_on_node;
  dml_layout->this_node            = layout->this_node;
  dml_layout->number_of_nodes      = layout->number_of_nodes;
  dml_layout->broadcast_globaldata = 0; /* Not applicable for writing */
  dml_layout->discover_dims_mode   = 0;
  				   
  dml_layout->hyperlower           = lower;
  dml_layout->hyperupper           = upper;
  dml_layout->subsetvolume         = layout->volume;
				   
  dml_layout->ionode               = io_node;
  dml_layout->master_io_node       = master_io_node();
  
  /* Construct the writer handle */
  qio_out = (QIO_Writer *)malloc(sizeof(QIO_Writer));
  if(qio_out == NULL) return NULL;
  qio_out->lrl_file_out   = NULL;
  qio_out->volfmt         = volfmt;
  qio_out->layout         = dml_layout;
  qio_out->dml_record_out = NULL;
  DML_checksum_init(&(qio_out->last_checksum));

  /* Unpack the QIO_Oflag parameter */
  if(oflag == NULL){
    /* For backward compatibility.  Default values */
    mode = QIO_TRUNC;
    qio_out->serpar = QIO_SERIAL;
    qio_out->ildgstyle = QIO_ILDGLAT;
    qio_out->ildgLFN = NULL;
  }
  else {
    mode = oflag->mode;
    qio_out->serpar = oflag->serpar;
    qio_out->ildgstyle = oflag->ildgstyle;
    /* This should be a deep rather than pointer copy I think */
    /* qio_out->ildgLFN = oflag->ildgLFN; */
    if( oflag->ildgLFN != NULL ) { 
      qio_out->ildgLFN=QIO_string_create();
      QIO_string_copy(qio_out->ildgLFN, oflag->ildgLFN);
    }
    else { 
      qio_out->ildgLFN = NULL ; /* NO user supplied LFN */
    }

  }
  serpar = qio_out->serpar;

  /* For now parallel I/O is supported only for SINGLEFILE */
  if(qio_out->volfmt != QIO_SINGLEFILE && 
     qio_out->serpar == QIO_PARALLEL){
    serpar = qio_out->serpar = QIO_SERIAL;
    if(QIO_verbosity() >= QIO_VERB_REG){
      printf("%s(%d): changed mode from QIO_PARALLEL to QIO_SERIAL \n",
	     myname,this_node);
    }
  }
  
  /*****************************/
  /* Open the file for writing */
  /*****************************/
  /* Which node does this depends on the user request. */
  /* If parallel write, ionodes open.
     If multifile, all nodes open.
     If writing by partitions, the partition I/O node does.
     In all cases, the master I/O node opens the file. */
  
  if( (qio_out->volfmt == QIO_MULTIFILE)
      || ((qio_out->volfmt == QIO_PARTFILE) 
	  && (dml_layout->ionode(this_node) == this_node))
      || ((serpar == QIO_PARALLEL)
	  && (dml_layout->ionode(this_node) == this_node))
      || (this_node == dml_layout->master_io_node) ) {
    /* Modifies filename for non master nodes */
    newfilename = QIO_filename_edit(filename, volfmt, dml_layout->this_node);
    if(qio_out->volfmt==QIO_SINGLEFILE) {
      if(this_node == dml_layout->master_io_node) {
	// optimization: one node creates file (if necessary) and closes
	lrl_file_out = LRL_open_write_file(newfilename, mode);
	LRL_close_write_file(lrl_file_out);
      }
      // the initial open will create or truncate the file if needed
      // so now everyone can open it without extra creation or truncation
      mode = LRL_NOTRUNC;
    }
    DML_sync();
    lrl_file_out = LRL_open_write_file(newfilename, mode);
    if(lrl_file_out == NULL) {
      printf("%s(%d): failed to open file for writing\n",myname,this_node);
      return NULL;
    }
    qio_out->lrl_file_out = lrl_file_out;
    free(newfilename);
  } else {
    DML_sync();
  }
  DML_sync();

  if(this_node == dml_layout->master_io_node && 
     QIO_verbosity() >= QIO_VERB_MED)
    {
      if(qio_out->volfmt == QIO_SINGLEFILE)
	printf("%s(%d): Opened %s for writing in singlefile format\n",
	       myname,this_node,filename);
      else if(qio_out->volfmt == QIO_MULTIFILE)
	printf("%s(%d): Opened %s for writing in multifile format\n",
	       myname,this_node,filename);
      else if(qio_out->volfmt == QIO_PARTFILE)
	printf("%s(%d): Opened %s for writing in partfile format\n",
	       myname,this_node,filename);
    }

  /* Determine sites to be written and create site list if needed */
  
  qio_out->sites = QIO_create_sitelist(qio_out->layout,qio_out->volfmt,
				       qio_out->serpar);
  if(qio_out->sites == NULL){
    printf("%s(%d): error creating sitelist\n",
	   myname,this_node);
    return NULL;
  }
  
  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): sitelist structure created \n",
	   myname,this_node);
    printf("%s(%d): I/O for %lu sites \n",
	   myname,this_node,(unsigned long)qio_out->sites->number_of_io_sites);
  }

  return qio_out;
}

/* Writes the private file XML record */
/* Writes the site list if multifile format */
/* Writes the user file XML record */

int QIO_write_file_header(QIO_Writer* qio_out, QIO_String *xml_file)
{
  QIO_String *xml_file_private;
  DML_Layout *dml_layout = qio_out->layout;
  int volfmt             = qio_out->volfmt;
  int serpar             = qio_out->serpar;
  int *latsize       = dml_layout->latsize;
  int latdim         = dml_layout->latdim;
  int this_node      = dml_layout->this_node;
  int master_io_node = dml_layout->master_io_node;

  QIO_FileInfo *file_info;
  int msg_begin, msg_end;
  char myname[] = "QIO_write_file_header";

  
  /****************************************/
  /* Load the private file info structure */
  /****************************************/
  
  file_info = QIO_create_file_info(latdim, latsize, volfmt);
  if(!file_info){
    printf("%s(%d): Can't create file info structure\n",myname,this_node);
    return QIO_ERR_FILE_INFO;
  }

  /*******************************/
  /* Encode the private file XML */
  /*******************************/

  xml_file_private = QIO_string_create();
  QIO_string_realloc(xml_file_private,QIO_STRINGALLOC);
  QIO_encode_file_info(xml_file_private, file_info);
  QIO_destroy_file_info(file_info);
  
  /*******************************************/
  /* Write the file header as a LIME message */
  /*******************************************/

  /* The master I/O node writes the full header consisting of
     (1) private file XML
     (2) site list (only if more than one volume)
     (3) user file XML
     If there are other I/O nodes, they write only
     (1) site list
  */
     
  /* First and last records in a message are flagged */
  msg_begin = 1; msg_end = 0;
  
  /* Master node writes the private file XML record */
  if(this_node == master_io_node){
    if(QIO_write_string(qio_out, msg_begin, msg_end, 
			xml_file_private, 
			(LIME_type)QIO_LIMETYPE_PRIVATE_FILE_XML)){
      printf("%s(%d): error writing private file XML\n",
	     myname,this_node);
      return QIO_ERR_PRIVATE_FILE_INFO;
    }
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): private file XML = %s\n",
	     myname,this_node,QIO_string_ptr(xml_file_private));
    }
    msg_begin = 0;
  }
  QIO_string_destroy(xml_file_private);

  /* For parallel I/O all nodes pretend to have written the XML */
  if(serpar == QIO_PARALLEL)
    msg_begin = 0;
  
  /******************************************/
  /* Write the sitelist if needed */
  /******************************************/

  /* Next record is last in message for all but master I/O node */
  /* For parallel I/O all nodes pretend to continue writing */
  if (this_node != master_io_node && serpar == QIO_SERIAL)msg_end = 1;
  
  if(QIO_write_sitelist(qio_out, msg_begin, msg_end, 
			(LIME_type)QIO_LIMETYPE_SITELIST)){
    printf("%s(%d): error writing the site list\n", myname,this_node);
    return QIO_ERR_BAD_SITELIST;
  }
  
  msg_begin = 0;
  
  /* Not really necessary */
  if(serpar == QIO_PARALLEL)
    msg_end = 1;
  
  /* Master node writes the user file XML record */
  /* For parallel output the other nodes pretend to write */
  if(this_node == master_io_node){
    msg_end = 1;
    if(QIO_write_string(qio_out, msg_begin, msg_end, xml_file, 
			(LIME_type)QIO_LIMETYPE_FILE_XML)){
      printf("%s(%d): error writing the user file XML\n",
	     myname,this_node);
      return QIO_ERR_FILE_INFO;
    }
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("%s(%d): wrote user file XML  = \"%s\"\n",
	     myname,this_node, QIO_string_ptr(xml_file));
  }

  return QIO_SUCCESS;
}
 
/* Opens a file for writing, whether in MPP mode or on host */
/* Writes the private file XML record */
/* Writes the site list if multifile format */
/* Writes the user file XML record */

QIO_Writer *QIO_open_write(QIO_String *xml_file, const char *filename, 
			   int volfmt, QIO_Layout *layout, 
			   QIO_Filesystem *fs, QIO_Oflag *oflag)
{
  QIO_Writer *qio_out;
  int status;
  int this_node = layout->this_node;
  DML_io_node_t my_io_node;
  DML_master_io_node_t master_io_node;

  /* Assign default behavior for io_node functions if needed */
  if(fs == NULL){
    my_io_node = DML_io_node;
    master_io_node = DML_master_io_node;
  }
  else{
    if(fs->my_io_node == NULL)
      my_io_node = DML_io_node;
    else
      my_io_node = fs->my_io_node;
    if(fs->master_io_node == NULL)
      master_io_node = DML_master_io_node;
    else
      master_io_node = fs->master_io_node;
  }

  qio_out = QIO_generic_open_write(filename, volfmt, layout, oflag,
				   my_io_node, master_io_node);

#if 0
  /* Prevent premature file truncation in parallel writes */
  /* Note, the test will cause a hang if the oflag->serpar value is
     not the same for all nodes */
  if(qio_out->serpar == QIO_PARALLEL) DML_sync();
#endif
  int fail = (qio_out==NULL);
  DML_sum_int(&fail);
  if(fail) {
    if(this_node == master_io_node()) {
      fprintf(stderr, "%s(%d): %i ranks failed to open file %s\n",
	      __func__, this_node, fail, filename);
    }
    // should free qio_out if not NULL
    return NULL;
  }

  status = QIO_write_file_header(qio_out, xml_file);
  fail = (status != QIO_SUCCESS);
  DML_sum_int(&fail);
  if(fail) {
    if(this_node == master_io_node()) {
      fprintf(stderr, "%s(%d): %i ranks failed to write header to file %s\n",
	      __func__, this_node, fail, filename);
    }
    // should free qio_out
    return NULL;
  }

  return qio_out;
}
