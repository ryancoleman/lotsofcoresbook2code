/* QIO_open_read.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <qioxml.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

/*****************************************/
/* Build QIO_Reader and open master file */
/*****************************************/

QIO_Reader *
QIO_create_reader(const char *filename, 
		  QIO_Layout *layout, QIO_Iflag *iflag,
		  int (*io_node)(int), int (*master_io_node)())
{
  QIO_Reader *qio_in;
  LRL_FileReader *lrl_file_in = NULL;
  DML_Layout *dml_layout;
  int i;
  int *latsize, *upper, *lower;
  int latdim = layout->latdim;
  int this_node = layout->this_node;
  int master_ionode = master_io_node();
  char *newfilename;
  char myname[] = "QIO_create_reader";

  int serpar=QIO_SERIAL, volfmt=QIO_UNKNOWN;
  if(iflag != NULL) {
    serpar = iflag->serpar;
    volfmt = iflag->volfmt;
  }

  /* First, only the global master node opens the file, regardless of
     whether it will be read by all nodes */

  if(this_node == master_ionode) {
    /* If the format is SINGLEFILE the file name is as specified by
       "filename".  If it is PARTFILE or MULTIFILE, the name is formed
       by adding a volume number extension to the "filename" string.
       If the volume format is unknown or single file, we try
       "filename" first. */
    if( volfmt == QIO_UNKNOWN || 
        volfmt == QIO_SINGLEFILE ) {
      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): Calling LRL_open_read_file %s\n",
	       myname,this_node,filename);
      lrl_file_in = LRL_open_read_file(filename);
      /* If the open succeeded with just "filename" the format
	 must be SINGLEFILE */
      if(lrl_file_in != NULL) volfmt = QIO_SINGLEFILE;
      else if( volfmt == QIO_UNKNOWN || 
	       volfmt == QIO_SINGLEFILE ) {
	printf("%s(%d): Can't open %s in singlefile format.\n Trying partfile or multifile.\n",
	       myname,this_node,filename);
      }
    }
    if(lrl_file_in == NULL) {
      /* If plain filename fails, try filename with volume number suffix */
      newfilename = QIO_filename_edit(filename, QIO_PARTFILE, this_node);
      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): Calling LRL_open_read_file %s\n",
	       myname,this_node,newfilename);
      lrl_file_in = LRL_open_read_file(newfilename);
      if(lrl_file_in == NULL) {
	printf("%s(%d): Can't open %s\n",myname,this_node,newfilename);
      }
      /* If we split an alien ILDG file, we currently put the
	 resulting PARTFILE volume in SciDAC format.  But if we ever
	 do split into an alien format, alien ILDG files with volume
	 extensions could only be PARTFILE, not MULTIFILE.  For native
	 SciDAC files, we will get the correct volume format later
	 when we decode the private file XML */
      if(volfmt == QIO_UNKNOWN)
	volfmt = QIO_PARTFILE;
      free(newfilename);
    }
  }
  int success = (lrl_file_in != NULL);
  DML_broadcast_bytes((char *)&success, sizeof(int), this_node, master_ionode);
  if(!success) return NULL;

  /* Make a local copy of lattice size */
  if(latdim != 0){
    latsize = (int *)malloc(sizeof(int)*latdim);
    for(i=0; i < latdim; ++i)
      latsize[i] = layout->latsize[i];
  }
  else
    latsize = NULL;

  /* Set up hypercube subset dimensions in case they are needed */
  if(latdim != 0){
    upper = (int *)malloc(sizeof(int)*latdim);
    lower = (int *)malloc(sizeof(int)*latdim);
  }  else {
    upper = NULL;
    lower = NULL;
  }

  /* Construct the layout data from the QIO_Layout structure*/
  dml_layout = (DML_Layout *)malloc(sizeof(DML_Layout));
  if (dml_layout == NULL || layout == NULL)
    return NULL;

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
  dml_layout->broadcast_globaldata = 1;
  dml_layout->discover_dims_mode   = (layout->latdim == 0);

  dml_layout->hyperlower           = lower;
  dml_layout->hyperupper           = upper;
  dml_layout->subsetvolume         = layout->volume;

  dml_layout->ionode               = io_node;
  dml_layout->master_io_node       = master_ionode;

  /* Construct the reader handle */
  qio_in = (QIO_Reader *)malloc(sizeof(QIO_Reader));
  if(qio_in == NULL){
    printf("%s(%d): Can't malloc QIO_Reader\n",myname,this_node);
    return NULL;
  }
  qio_in->lrl_file_in = NULL;
  qio_in->layout      = dml_layout;
  qio_in->sites       = NULL;
  qio_in->read_state  = QIO_RECORD_INFO_PRIVATE_NEXT;
  qio_in->xml_record  = NULL;
  qio_in->ildgLFN     = QIO_string_create();
  DML_checksum_init(&(qio_in->last_checksum));

  qio_in->serpar = serpar;
  qio_in->volfmt = volfmt;

  /* Force single file format if there is only one node */
  if(layout->number_of_nodes==1) qio_in->volfmt = QIO_SINGLEFILE;

  /* We start by assuming this is neither a foreign ILDG file
     nor a SciDAC ILDG file */
  qio_in->format = QIO_SCIDAC_NATIVE;
  qio_in->ildgstyle = QIO_ILDGNO;
  qio_in->ildg_precision = 0;

  qio_in->lrl_file_in = lrl_file_in;

  return qio_in;
}

/* Used during file conversion on a scalar machine */
void QIO_suppress_global_broadcast(QIO_Reader *qio_in){
  qio_in->layout->broadcast_globaldata = 0;
}

/*******************************************/
/* Read private file info from master file */
/*******************************************/

QIO_FileInfo *QIO_read_private_file_info(QIO_Reader *qio_in)
{
  QIO_String *xml_file_private;
  QIO_FileInfo *file_info_found = NULL;
  QIO_ILDGFormatInfo *ildg_info = NULL;
  DML_Layout *dml_layout = qio_in->layout;
  int this_node = dml_layout->this_node;
  LRL_RecordReader *lrl_record_in;
  LIME_type lime_type = NULL;
  uint64_t expected_rec_size;
  int status;
  int ildg_dims[4];
  int ntypes = 1;
  LIME_type lime_type_list[1] = {
    QIO_LIMETYPE_ILDG_FORMAT
  };
  char myname[] = "QIO_read_private_file_info";

  /* Master node reads and decodes the private file XML record */
  if(this_node == dml_layout->master_io_node){
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): Reading xml_file_private\n",myname,this_node);fflush(stdout);
    }

    xml_file_private = QIO_string_create();
    QIO_string_realloc(xml_file_private,QIO_STRINGALLOC);

    /* Get the first record type */
    lrl_record_in = 
      QIO_read_record_type(qio_in, &lime_type, &expected_rec_size, &status);
    if(lrl_record_in == NULL){
      printf("%s(%d): error %d finding private file XML\n",
	     myname,this_node,status);
      return NULL;
    }

    /* If SciDAC format, we expect the private file info type */
    if(strcmp(lime_type,QIO_LIMETYPE_PRIVATE_FILE_XML) == 0){
      /* Read the file info string */
      qio_in->format = QIO_SCIDAC_NATIVE;

      status = QIO_read_string_data(qio_in, lrl_record_in, xml_file_private, 
				    expected_rec_size);
      if(status != QIO_SUCCESS){
	printf("%s(%d): error %d reading private file XML\n",
	       myname,this_node,status);
	return NULL;
      }
      
      if(QIO_verbosity() >= QIO_VERB_DEBUG){
	printf("%s(%d): private file XML = %s\n",myname,this_node,
	       QIO_string_ptr(xml_file_private));
      }
      /* Decode the file info */
      file_info_found = QIO_create_file_info(0,NULL,0);
      QIO_decode_file_info(file_info_found, xml_file_private);
      QIO_string_destroy(xml_file_private);
    }
    
    else{
      qio_in->format = QIO_ILDG_ALIEN;
      if(QIO_verbosity() >= QIO_VERB_REG)
	printf("%s(%d): Non-SciDAC file.  Attempting to read as an alien ILDG lattice file\n",
	       myname,this_node);
      
      if(strcmp(lime_type,(LIME_type)QIO_LIMETYPE_ILDG_FORMAT)!=0){
	/* If the first record is not the ILDG format record,
	   keep reading until we find it */
	/* Close the LIME record we just read and don't want */
	status = QIO_close_read_record(lrl_record_in);
	if(status != QIO_SUCCESS){
	  printf("%s(%d): error %d closing alien record\n",
		 myname,this_node,status);
	  return NULL;
	}
	/* Scan for the ILDG format record */
	lrl_record_in = QIO_open_read_target_record(qio_in, 
	    lime_type_list, ntypes, &lime_type, &expected_rec_size, &status);
	
	/* One cause of failure is hitting an EOF */
	if(lrl_record_in == NULL){
	  printf("%s(%d): error %d finding ILDG format record\n",
		 myname,this_node,status);
	  return NULL;
	}
      }

      /* Found the ILDG format record.  Now read the XML string. */
      qio_in->ildgstyle = QIO_ILDGLAT;

      status = QIO_read_string_data(qio_in, lrl_record_in, xml_file_private, 
				    expected_rec_size);
      if(status != QIO_SUCCESS){
	printf("%s(%d): error %d reading ILDG format record\n",
	       myname,this_node,status);
	return NULL;
      }
      
      if(QIO_verbosity() >= QIO_VERB_DEBUG){
	printf("%s(%d): ILDG format XML = %s\n",myname,this_node,
	       QIO_string_ptr(xml_file_private));
      }
      /* Decode the ILDG format XML record */
      ildg_info = QIO_create_ildg_format_info(0, NULL);
      QIO_decode_ILDG_format_info(ildg_info, xml_file_private);
      QIO_string_destroy(xml_file_private);
      
      /* Extract the ILDG format information */
      /* We don't look at the field type for now, since only lattice
	 data is currently supported in the ILDG */
      qio_in->ildg_precision = QIO_get_ildgformat_precision(ildg_info);
      ildg_dims[0] = QIO_get_ildgformat_lx(ildg_info);
      ildg_dims[1] = QIO_get_ildgformat_ly(ildg_info);
      ildg_dims[2] = QIO_get_ildgformat_lz(ildg_info);
      ildg_dims[3] = QIO_get_ildgformat_lt(ildg_info);
      
      /* Create a SciDAC file info structure 
	 and populate it with the ILDG format data */
      
      file_info_found = QIO_create_file_info(4, ildg_dims, qio_in->volfmt);
    }
  }    
  return file_info_found;
}

/***********************************/
/* Accessors for header information */
/***********************************/
int QIO_get_reader_latdim(QIO_Reader *in){
  return in->layout->latdim;
}

int *QIO_get_reader_latsize(QIO_Reader *in){
  return in->layout->latsize;
}

uint32_t QIO_get_reader_last_checksuma(QIO_Reader *in){
  return in->last_checksum.suma;
}

uint32_t QIO_get_reader_last_checksumb(QIO_Reader *in){
  return in->last_checksum.sumb;
}

char *QIO_get_ILDG_LFN(QIO_Reader *in){
  return QIO_string_ptr(in->ildgLFN);
}

int QIO_get_ildgstyle(QIO_Reader *in){
  return in->ildgstyle;
}

int QIO_get_reader_volfmt(QIO_Reader *in){
  return in->volfmt;
}

int QIO_get_reader_format(QIO_Reader *in){
  return in->format;
}

void QIO_set_record_info(QIO_Reader *in, QIO_RecordInfo *rec_info){
  memcpy((char *)&(in->record_info),(char *)rec_info,sizeof(QIO_RecordInfo));
}

/*****************************************/
/* In case we are in discovery mode,
   set the QIO_Reader lattice dimensions */
/*****************************************/

int QIO_set_latdim(QIO_Reader *qio_in, int latdim, int *latsize)
{

  DML_Layout *dml_layout = qio_in->layout;
  int i;
  int this_node = dml_layout->this_node;
  int latdim_r = dml_layout->latdim;
  char myname[] = "QIO_set_latdim";

  /* Make space for lattice dimensions */
  if(latdim != latdim_r){
    /* Set correct space for given dimensions */
    dml_layout->latsize = (int *)realloc(dml_layout->latsize, 
					 latdim*sizeof(int));
    if(!dml_layout->latsize){
      printf("%s(%d): Can't malloc space for latsize dim %d\n",
	     myname,this_node,latdim);
      return QIO_ERR_ALLOC;
    }
    dml_layout->latdim = latdim;
  }

  for(i = 0; i < latdim; i++)
    dml_layout->latsize[i] = latsize[i];

  return QIO_SUCCESS;
}

/*****************************************/
/* In case we are reading intermittently,
   set the file pointer                  */
/*****************************************/

int QIO_set_reader_pointer(QIO_Reader *qio_in, off_t offset)
{
  int status;

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("QIO_set_reader_pointer(%d): Setting to %llu\n",
	   qio_in->layout->this_node, (unsigned long long)offset);
  status = LRL_set_reader_pointer(qio_in->lrl_file_in, offset);
  if(status != LRL_SUCCESS)
    return QIO_ERR_BAD_SEEK;
  return QIO_SUCCESS;
}

off_t QIO_get_reader_pointer(QIO_Reader *qio_in)
{
  off_t offset;

  offset = LRL_get_reader_pointer(qio_in->lrl_file_in);

  if(QIO_verbosity() >= QIO_VERB_DEBUG)
    printf("QIO_get_reader_pointer(%d): Got %llu\n",
	   qio_in->layout->this_node, (unsigned long long)offset);
  return offset;
}


/****************************************************************/
/* Compares lattice dimension and lattice size in file with the
   user-specified dimension and size */
/****************************************************************/

int QIO_check_file_info(DML_Layout *dml_layout, QIO_FileInfo *file_info_found)
{
  QIO_FileInfo *file_info_expect;
  int this_node = dml_layout->this_node;
  int status;
  char myname[] = "QIO_check_file_info";

  /* Create structure with what we expect */
  /* (We discover and respond to the volume format parameter.
     A zero here implies a matching value is not required.) */
  file_info_expect = 
    QIO_create_file_info(dml_layout->latdim, dml_layout->latsize, 0);
  
  /* Compare what we found and what we expected */
  if((status = 
      QIO_compare_file_info(file_info_found,file_info_expect,
			    myname,this_node))!=QIO_SUCCESS){
    printf("%s(%d): Error %d checking file info\n",myname,this_node,status);
    return QIO_ERR_FILE_INFO;
  }
  QIO_destroy_file_info(file_info_expect);
  return QIO_SUCCESS;
}

/************************************************************************/
/* All nodes create the reader.  Master I/O node opens the file and
   reads the private file XML */
/************************************************************************/
/* Other nodes just create the reader */

QIO_Reader *
QIO_open_read_master(const char *filename, 
		     QIO_Layout *layout, QIO_Iflag *iflag,
		     int (*io_node)(int), int (*master_io_node)())
{
  QIO_Reader *qio_in;
  DML_Layout *dml_layout;
  QIO_FileInfo *file_info_found;
  int this_node = layout->this_node;
  int status;
  char myname[] = "QIO_open_read_master";

  /* First, only the global master node opens the file, regardless of
     whether it will be read by all nodes */

  qio_in = QIO_create_reader(filename, layout, iflag, io_node, master_io_node);
  if(!qio_in)return NULL;

  dml_layout = qio_in->layout;

  if(this_node == dml_layout->master_io_node && !qio_in->lrl_file_in)
    return NULL;

  /* Master node reads and decodes the private file XML record */
  /* For parallel input the other nodes pretend to read */

  if(this_node == dml_layout->master_io_node) {

    /* Read the private file info from the master file */
    file_info_found = QIO_read_private_file_info(qio_in);
    if(file_info_found == NULL)return NULL;

    /* Compare what we found with what we expected */
    /* If we read in discovery mode, accept whatever the file gives
       us */

    if(!qio_in->layout->discover_dims_mode){
      status = QIO_check_file_info(dml_layout, file_info_found);
      if(status != QIO_SUCCESS)return NULL;
    }

    /* Get the volume format from the file info */
    qio_in->volfmt = QIO_get_volfmt(file_info_found);
    if(QIO_verbosity() >= QIO_VERB_DEBUG)
      printf("%s(%d): set volfmt to %d\n",myname,this_node,qio_in->volfmt);

    /* Get lattice dimensions from the file info
       if we are in discovery mode */
    if(qio_in->layout->discover_dims_mode){
      status = QIO_set_latdim(qio_in, QIO_get_spacetime(file_info_found),
			      QIO_get_dims(file_info_found));
      if(status != QIO_SUCCESS)return NULL;
    }

   QIO_destroy_file_info(file_info_found);
  }

  /* Return what we have read so far */
  return qio_in;
}

int QIO_broadcast_file_reader_info(QIO_Reader *qio_in, int discover_dims)
{
  DML_Layout *dml_layout = qio_in->layout;
  int this_node = dml_layout->this_node;
  int discover = discover_dims;
  char myname[] = "QIO_broadcast_file_reader_info";

  /* Master I/O node broadcasts the volume format to all the nodes, */
  /* inserting the value in the qio_in structure */
  DML_broadcast_bytes((char *)(&qio_in->volfmt), sizeof(int),
		      this_node, DML_master_io_node());
  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): volume format info was broadcast\n",
	   myname,this_node);fflush(stdout);
  }

  /* In discover dimensions mode the master broadcasts the lattice
     dimensions in its reader structure to all the nodes */

  /* Master node calls for discovery */
  DML_broadcast_bytes((char *)(&discover), sizeof(int),
		      this_node, DML_master_io_node());

  if(discover){
    /* Set the lattice dimension from the master */
    DML_broadcast_bytes((char *)(&dml_layout->latdim),
			sizeof(int), this_node, DML_master_io_node());
    /* Adjust space for the lattice dimensions */
    dml_layout->latsize = (int *)realloc(dml_layout->latsize,
					 sizeof(int)*dml_layout->latdim);
    if(!dml_layout->latsize){
      printf("%s(%d): Can't realloc latsize\n",myname,this_node);
      return QIO_ERR_ALLOC;
    }
    /* Broadcast the lattice dimensions */
    DML_broadcast_bytes((char *)(dml_layout->latsize),
			sizeof(int)*dml_layout->latdim, 
			this_node, DML_master_io_node());
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): lattice dimension info was broadcast\n",
	     myname,this_node);fflush(stdout);
    }
  }

  return QIO_SUCCESS;
}



/*********************************************************************/
/* Other nodes open their file (if needed)                           */
/*********************************************************************/

int
QIO_open_read_nonmaster(QIO_Reader *qio_in, const char *filename,
			QIO_Iflag *iflag)
{
  DML_Layout *dml_layout = qio_in->layout;
  int this_node = dml_layout->this_node;
  LRL_FileReader *lrl_file_in = NULL;
  char *newfilename;
  char myname[] = "QIO_open_read_nonmaster";

  /* Open any additional file handles as needed */
  /* Read the sitelist if needed */

  /* Parallel I/O is supported only for SINGLEFILE. */
  if(qio_in->volfmt != QIO_SINGLEFILE && 
     qio_in->serpar == QIO_PARALLEL){
    qio_in->serpar = QIO_SERIAL;
    if(QIO_verbosity() >= QIO_VERB_REG){
      printf("%s(%d): changed mode from QIO_PARALLEL to QIO_SERIAL \n",
	     myname,this_node);
    }
  }
  
  if(qio_in->volfmt == QIO_SINGLEFILE) {
    /* One file for all nodes */
    /* If parallel read, all io nodes (except master) open the file */
    if( qio_in->serpar == QIO_PARALLEL &&
	this_node != dml_layout->master_io_node &&
	this_node == dml_layout->ionode(this_node) ) {
      lrl_file_in = LRL_open_read_file(filename);
      if(lrl_file_in == NULL){
	printf("%s(%d): Can't open %s\n",myname,this_node,filename);
	return QIO_ERR_OPEN_READ;
      }
      qio_in->lrl_file_in = lrl_file_in; 
    }

    if(this_node == dml_layout->master_io_node){
      if(QIO_verbosity() >= QIO_VERB_MED)
	printf("%s(%d): Opened %s for reading in singlefile mode\n",
	       myname,this_node,filename);fflush(stdout);
    }
  }
  else if(qio_in->volfmt == QIO_PARTFILE) {
    /* One file per machine partition in lexicographic order */
    if(this_node == dml_layout->master_io_node){
      if(QIO_verbosity() >= QIO_VERB_MED)
	printf("%s(%d): Opened %s for reading in partfile mode\n",
	       myname,this_node,filename);fflush(stdout);
    }

    /* All the partition I/O nodes open their files.  */
    if(this_node == dml_layout->ionode(this_node)){
      /* (The global master has already opened its file) */
      if(this_node != dml_layout->master_io_node){
	/* Construct the file name based on the partition I/O node number */
	newfilename = QIO_filename_edit(filename, qio_in->volfmt, this_node);
	/* Open the file */
	if(QIO_verbosity() >= QIO_VERB_DEBUG)
	  printf("%s(%d): Calling LRL_open_read_file %s\n",
		 myname,this_node,newfilename);
	lrl_file_in = LRL_open_read_file(newfilename);
	if(QIO_verbosity() >= QIO_VERB_DEBUG)
	  printf("%s(%d): LRL_open_read_file returns %p\n",myname,this_node,
		 lrl_file_in);
	if(lrl_file_in == NULL){
	  printf("%s(%d): Can't open %s for reading\n",myname,this_node,
		 newfilename);
	  free(newfilename);
	  return QIO_ERR_OPEN_READ;
	}
	free(newfilename);
	qio_in->lrl_file_in = lrl_file_in; 
      }
    }
  }
  else if(qio_in->volfmt == QIO_MULTIFILE) {
    /* One file per node */
    if(this_node == dml_layout->master_io_node){
      if(QIO_verbosity() >= QIO_VERB_MED)
	printf("%s(%d): Opened %s for reading in multifile mode\n",
	       myname,this_node,filename);fflush(stdout);
    }
    /* The non-master nodes open their files */
    if(this_node != dml_layout->master_io_node){
      /* Edit file name */
      newfilename = QIO_filename_edit(filename, qio_in->volfmt, this_node);
      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): Calling LRL_open_read_file %s\n",
	       myname,this_node,newfilename);
      lrl_file_in = LRL_open_read_file(newfilename);
      if(lrl_file_in == NULL){
	printf("%s(%d): Can't open %s for reading\n",myname,this_node,
	       newfilename);
	free(newfilename);
	return QIO_ERR_OPEN_READ;
      }
      free(newfilename);
      qio_in->lrl_file_in = lrl_file_in; 
    }
  }
  else {
    printf("%s(%d): bad volfmt parameter %d\n",myname,this_node,
	   qio_in->volfmt);
    return QIO_ERR_FILE_INFO;
  }

  return QIO_SUCCESS;
}

/****************************************/
/* Nodes read their site lists (if any) */
/****************************************/

int QIO_read_check_sitelist(QIO_Reader *qio_in){

  int this_node = qio_in->layout->this_node;
  LIME_type lime_type = NULL;
  int status;
  char myname[] = "QIO_read_check_sitelist";

  /* Create the expected sitelist.  Input must agree exactly. */

  qio_in->sites = QIO_create_sitelist(qio_in->layout, qio_in->volfmt, 
				      qio_in->serpar);
  if(qio_in->sites == NULL){
    printf("%s(%d): error creating sitelist\n",
	   myname,this_node);
    return QIO_ERR_BAD_SITELIST;
  }
  
  /* We don't expect to find a sitelist in an alien ILDG file */
  if(qio_in->format == QIO_ILDG_ALIEN)return QIO_SUCCESS;

  if((status = QIO_read_sitelist(qio_in, &lime_type))!= QIO_SUCCESS){
    printf("%s(%d): Error %d reading site list\n",myname,this_node,status);
    return QIO_ERR_BAD_SITELIST;
  }
  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Sitelist was read\n",myname,this_node);fflush(stdout);
  }
  
  return QIO_SUCCESS;
}

/*************************************************************/
/* Master node reads and broadcasts the user file XML record */
/*************************************************************/

int QIO_read_user_file_xml(QIO_String *xml_file, QIO_Reader *qio_in){

  /* Calling program must allocate *xml_file */

  int this_node = qio_in->layout->this_node;
  LIME_type lime_type = NULL;
  int status;
  char myname[] = "QIO_read_user_file_xml";

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Reading user file XML\n",myname,this_node);fflush(stdout);
  }
  
  /* Assumes the xml_file structure was created by caller */
  if(this_node == qio_in->layout->master_io_node){
    if(qio_in->format == QIO_SCIDAC_NATIVE){
      if((status = 
	  QIO_read_string(qio_in, xml_file, &lime_type))!= QIO_SUCCESS){
	printf("%s(%d): error %d reading user file XML\n",
	       myname,this_node,status);
	return QIO_ERR_PRIVATE_FILE_INFO;
      }
    }
    else{
      QIO_string_set(xml_file,QIO_NONSCIDAC_FILE);
    }
    
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): user file XML = \"%s\"\n",
	     myname,this_node,QIO_string_ptr(xml_file));
    }
  }
  
  return QIO_SUCCESS;
}

/**************************************/
/* Main entry point for compute nodes */
/**************************************/

QIO_Reader *
QIO_open_read(QIO_String *xml_file, const char *filename, 
	      QIO_Layout *layout, QIO_Filesystem *fs,
	      QIO_Iflag *iflag)
{
  QIO_Reader *qio_in;
  DML_Layout *dml_layout;
  char myname[] = "QIO_open_read";
  int this_node = layout->this_node;
  int status;
  int length;
  DML_io_node_t my_io_node;
  DML_master_io_node_t master_io_node;

  /* Assign default behavior for io_node functions if needed */
  if(fs == NULL) {
    my_io_node = DML_io_node;
    master_io_node = DML_master_io_node;
  } else {
    if(fs->my_io_node == NULL)
      my_io_node = DML_io_node;
    else
      my_io_node = fs->my_io_node;
    if(fs->master_io_node == NULL)
      master_io_node = DML_master_io_node;
    else
      master_io_node = fs->master_io_node;
  }

  /* On the compute nodes, we use DML calls to specify the I/O nodes
     and the master I/O node */
  qio_in = QIO_open_read_master(filename, layout, iflag,
				my_io_node, master_io_node);
  if(qio_in == NULL) return NULL;

  /* Master I/O node broadcasts the volume format to all the nodes, */
  /* inserting the value in the qio_in structure */
  /* In discovery mode the master also broadcasts the lattice
     dimensions to all the nodes */
  status = QIO_broadcast_file_reader_info(qio_in, 
					  qio_in->layout->discover_dims_mode);

  /* Read the rest of the header */
  status = QIO_open_read_nonmaster(qio_in, filename, iflag);
  if(status != QIO_SUCCESS)return NULL;

  status = QIO_read_check_sitelist(qio_in);
  if(status != QIO_SUCCESS)return NULL;

  status = QIO_read_user_file_xml(xml_file, qio_in);
  if(status != QIO_SUCCESS)return NULL;

  /* Broadcast the user file XML to all nodes */

  dml_layout = qio_in->layout;
  if(this_node == dml_layout->master_io_node)
    length = QIO_string_length(xml_file);

  /* First broadcast length */
  DML_broadcast_bytes((char *)&length,sizeof(int),
		      this_node, dml_layout->master_io_node);
  
  /* Receiving nodes resize their strings */
  if(this_node != dml_layout->master_io_node){
    QIO_string_realloc(xml_file,length);
  }

  /* Then broadcast the string itself */
  DML_broadcast_bytes(QIO_string_ptr(xml_file),length,
		      this_node, dml_layout->master_io_node);
  
  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Done with user file XML\n",myname,this_node);
    fflush(stdout);
  }

  return qio_in;
}
