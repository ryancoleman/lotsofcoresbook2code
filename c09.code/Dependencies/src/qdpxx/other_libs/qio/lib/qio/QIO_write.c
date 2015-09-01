/* QIO_write.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <stdio.h>
#include <string.h>

/* Update the record type and hypercube information */
int QIO_writer_insert_hypercube_data(QIO_Writer *out, 
				     QIO_RecordInfo *record_info)
{

  int status;

  /* Copy subset data from record_info structure to layout structure
     and check it */

  status = DML_insert_subset_data(out->layout, 
				  QIO_get_recordtype(record_info),
				  QIO_get_hyperlower(record_info),
				  QIO_get_hyperupper(record_info),
				  QIO_get_hyper_spacetime(record_info));
  if(status != 0) return QIO_ERR_BAD_SUBSET;

  return QIO_SUCCESS;
}

/* Write the record metadata for a lattice field. */

/* Handles the write operation on the compute nodes as well as the host */
int QIO_write_record_info(QIO_Writer *out, QIO_RecordInfo *record_info, 
              size_t datum_size, int word_size,
	      QIO_String *xml_record, 
	      int *msg_begin, int *msg_end){

  QIO_String *xml_record_private;
  QIO_String *xml_ildg_format;
  QIO_ILDGFormatInfo *ildg_info;
  int this_node = out->layout->this_node;
  int master_io_node = out->layout->master_io_node;
  int ildg_precision;
  int status;
  int count = QIO_get_datacount(record_info);
  char myname[] = "QIO_write_record_info";

  /* Require consistency between the byte count specified in the
     private record metadata and the byte count per site to be written */
  if(datum_size != QIO_get_typesize(record_info) * count)
    {
      printf("%s(%d): bytes per site mismatch %lu != %d * %d\n",
	     myname,this_node,(unsigned long)datum_size,
	     QIO_get_typesize(record_info),
	     QIO_get_datacount(record_info));
      return QIO_ERR_BAD_WRITE_BYTES;
    }

  /* Copy hypercube data from record_info structure to writer */
  status = QIO_writer_insert_hypercube_data(out, record_info);
  if(status != QIO_SUCCESS)return status;

  /* A message consists of the XML, binary payload, and checksums */
  /* An ILDG lattice message includes an ILDG metadata record */
  /* First and last records in a message are flagged */
  *msg_begin = 1; *msg_end = 0;

  /* Create private record XML */
  xml_record_private = QIO_string_create();
  QIO_encode_record_info(xml_record_private, record_info);

  /* Master node writes the private record XML record */
  if(this_node == master_io_node){
    if ((status = 
	 QIO_write_string(out, *msg_begin, *msg_end, 
			  xml_record_private, 
			  (LIME_type)QIO_LIMETYPE_PRIVATE_RECORD_XML))
	!= QIO_SUCCESS){
      printf("%s(%d): Error writing private record XML\n",
	     myname,this_node);
      return status;
    }

    if(QIO_verbosity() >= QIO_VERB_REG){
      printf("%s(%d): private record XML = \"%s\"\n",
	     myname,this_node,QIO_string_ptr(xml_record_private));
    }
    *msg_begin = 0;
  }
  /* In singlefile parallel mode all nodes pretend they also wrote the
     private record XML */
  if(out->serpar == QIO_PARALLEL && out->volfmt == QIO_SINGLEFILE)
    *msg_begin = 0;

  QIO_string_destroy(xml_record_private);

  /* Master node writes the user's record metadata */
  if(this_node == master_io_node){
    if ((status = 
	 QIO_write_string(out, *msg_begin, *msg_end, xml_record, 
			  (LIME_type)QIO_LIMETYPE_RECORD_XML))
	!= QIO_SUCCESS){
      printf("%s(%d): Error writing user record XML\n",myname,this_node);
      return status;
    }
    if(QIO_verbosity() >= QIO_VERB_DEBUG){
      printf("%s(%d): user record XML = \"%s\"\n",
	     myname,this_node,QIO_string_ptr(xml_record));
    }
    *msg_begin = 0;
  }
  /* In singlefile parallel mode all nodes pretend they also wrote the
     user record XML */
  if(out->serpar == QIO_PARALLEL && out->volfmt == QIO_SINGLEFILE)
    *msg_begin = 0;

  /* In case of an ILDG lattice record, create the ILDG lattice metadata */
  if(this_node == master_io_node && out->ildgstyle == QIO_ILDGLAT){

    /* Check for a valid ILDG datatype, count, dimension and get precision */
    ildg_precision = 0;
    
    /* The ILDG lattice must have the correct data type */
    if(strcmp(QIO_get_datatype(record_info),"QDP_F3_ColorMatrix") == 0 ||
       strcmp(QIO_get_datatype(record_info),"USQCD_F3_ColorMatrix") == 0)
      ildg_precision = 32;

    else if(strcmp(QIO_get_datatype(record_info),"QDP_D3_ColorMatrix") == 0 ||
	    strcmp(QIO_get_datatype(record_info),"USQCD_D3_ColorMatrix") == 0) 
      ildg_precision = 64;

    /* There must be four color matrices per site */
    if(QIO_get_datacount(record_info) != 4)
      ildg_precision = 0;

    /* Only four-dimensional lattices are supported */
    if(out->layout->latdim != 4)
      ildg_precision = 0;

    /* A zero value cancels ILDG */
    if(ildg_precision == 0){
      /* Other datatypes not supported. */
      out->ildgstyle = QIO_ILDGNO;
      if(QIO_verbosity() >= QIO_VERB_LOW){
	printf("%s(%d): ILDG format not supported for datatype %s and datacount %d and dimension %d.  ILDG format cancelled.\n",
	       myname,this_node,QIO_get_datatype(record_info),
	       QIO_get_datacount(record_info),out->layout->latdim);
      }
    }
    
    else{
      /* Create data structure with ILDG information */
      ildg_info = QIO_create_ildg_format_info(ildg_precision, 
					      out->layout->latsize);
      /* Convert data structure to XML */
      xml_ildg_format = QIO_string_create();
      QIO_encode_ILDG_format_info(xml_ildg_format, ildg_info);
      QIO_destroy_ildg_format_info(ildg_info);
      
      /* Write the ildg-format record */
      if ((status = 
	   QIO_write_string(out, *msg_begin, *msg_end, xml_ildg_format, 
			    (LIME_type)QIO_LIMETYPE_ILDG_FORMAT))
	  != QIO_SUCCESS){
	printf("%s(%d): Error writing user record XML\n",myname,this_node);
	return status;
      }
      QIO_string_destroy(xml_ildg_format);

      /* Write the ildg-data-lfn record if the LFN is known */
      if(out->ildgLFN != NULL){
	if ((status = 
	     QIO_write_string(out, *msg_begin, *msg_end, 
			      out->ildgLFN, 
			      (LIME_type)QIO_LIMETYPE_ILDG_DATA_LFN))
	    != QIO_SUCCESS){
	  printf("%s(%d): Error writing ILDG LFN record\n",
		 myname,this_node);
	  return status;
	}
	if(QIO_verbosity() >= QIO_VERB_REG){
	  printf("%s(%d): LFN = \"%s\"\n",
		 myname,this_node,QIO_string_ptr(out->ildgLFN));
	}
      }
    }
  }

  return QIO_SUCCESS;
}

/* Write the binary payload for a lattice field, but not the checksum */

/* Handles the write operation on the compute nodes as well as the host */
int QIO_write_record_data(QIO_Writer *out, QIO_RecordInfo *record_info, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      size_t datum_size, int word_size, void *arg,
	      DML_Checksum *checksum, uint64_t *nbytes,
	      int *msg_begin, int *msg_end){

  int this_node = out->layout->this_node;
  int master_io_node = out->layout->master_io_node;
  int status;
  int count = QIO_get_datacount(record_info);
  char scidac_type[] = QIO_LIMETYPE_BINARY_DATA;
  char ildg_type[] = QIO_LIMETYPE_ILDG_BINARY_DATA;
  LIME_type lime_type;
  char myname[] = "QIO_write_record_data";

  /* Next one is last record in message for all but master node */
  /* But if we are writing in parallel mode all nodes think they
     are the master */
  if (this_node != master_io_node && out->serpar == QIO_SERIAL)
    *msg_end = 1;
  
  /* Set LIME type for the data record.  Depends whether we are creating
     an ILDG compatible file */
  if(out->ildgstyle == QIO_ILDGLAT)lime_type = ildg_type;
  else lime_type = scidac_type;

  status = 
    QIO_write_field(out, *msg_begin, *msg_end, 
		    get, count, datum_size, word_size, arg, 
		    checksum, nbytes, lime_type);
  if(status != QIO_SUCCESS){
    printf("%s(%d): Error writing field data\n",myname,this_node);
    return status;
  }
  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): wrote field\n",myname,this_node);fflush(stdout);
  }

  /* Copy most recent node checksum into writer */
  memcpy(&(out->last_checksum),checksum,sizeof(DML_Checksum));

  return QIO_SUCCESS;
}



/* Write records for a lattice field.  Writes metadata, binary payload,
   but NOT checksum */

/* Handles the write operation on the compute nodes or the host in
   case of file conversion */
int QIO_generic_write(QIO_Writer *out, QIO_RecordInfo *record_info, 
	      QIO_String *xml_record, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      size_t datum_size, int word_size, void *arg,
	      DML_Checksum *checksum, uint64_t *nbytes,
	      int *msg_begin, int *msg_end){

  int status;

  status = QIO_write_record_info(out, record_info, datum_size, word_size, 
				 xml_record, msg_begin, msg_end);
  if(status != QIO_SUCCESS)
    return status;

  status = QIO_write_record_data(out, record_info, get, datum_size, 
				     word_size, arg, checksum, nbytes, 
				     msg_begin, msg_end);
  return status;
}

/* Write the checksum record */

int QIO_write_checksum(QIO_Writer *out, DML_Checksum *checksum)
{
  QIO_ChecksumInfo *checksum_info;
  int this_node = out->layout->this_node;
  int master_io_node = out->layout->master_io_node;
  int status;
  int msg_end = 1;   /* The last record */
  int msg_begin = 0; /* Always a preceding record */
  QIO_String *xml_checksum;
  char myname[] = "QIO_write_checksum";

  /* Master node encodes and writes the checksum and announces result */
  if(this_node == master_io_node){
    checksum_info = QIO_create_checksum_info(checksum->suma,checksum->sumb);
    xml_checksum = QIO_string_create();
    QIO_encode_checksum_info(xml_checksum, checksum_info);

    if ((status = 
	 QIO_write_string(out, msg_begin, msg_end, xml_checksum,
			  (LIME_type)"scidac-checksum"))
	!= QIO_SUCCESS) {
      printf("%s(%d): Error writing checksum\n",myname,this_node);
      return status;
    }

    QIO_string_destroy(xml_checksum);
    QIO_destroy_checksum_info(checksum_info);
  }
  return QIO_SUCCESS;
}

/* Write a lattice field on the compute nodes. Includes checksum
   record */

int QIO_write(QIO_Writer *out, QIO_RecordInfo *record_info, 
	      QIO_String *xml_record, 
	      void (*get)(char *buf, size_t index, int count, void *arg),
	      size_t datum_size, int word_size, void *arg){

  DML_Checksum checksum;
  uint64_t nbytes;
  int this_node = out->layout->this_node;
  int msg_begin, msg_end;
  int status;
  int recordtype;
  uint64_t total_bytes;
  size_t volume;
  char myname[] = "QIO_write";

  status = QIO_generic_write(out, record_info, xml_record, get, datum_size, 
			     word_size, arg, &checksum, &nbytes, 
			     &msg_begin, &msg_end);

 if(status != QIO_SUCCESS)return status;

  recordtype = out->layout->recordtype;

  /* Combine checksums over all nodes */
  DML_checksum_combine(&checksum);

  /* Copy most recent combined checksum into writer */
  memcpy(&(out->last_checksum),&checksum,sizeof(DML_Checksum));

  /* Sum the bytes written by all nodes */
  DML_sum_uint64_t(&nbytes);

  /* Compute and compare byte count with expected record size */
  volume = out->layout->subsetvolume;
  if(recordtype == QIO_GLOBAL)
    total_bytes = datum_size;
  else 
    total_bytes = ((uint64_t)volume) * datum_size;

  /* 
   * Global data may only live on the master_io_node, so only return an
   * error if this node expects data to have been written
   */
  if (! (recordtype == QIO_GLOBAL && total_bytes == 0))
  {
    if(nbytes != total_bytes)
    {
      printf("%s(%d): bytes written %llu != expected rec_size %llu\n",
	     myname, this_node, (unsigned long long)nbytes, 
	     (unsigned long long)total_bytes);
      return QIO_ERR_BAD_WRITE_BYTES;
    }
  }

  /* Write the checksum */
  status = QIO_write_checksum(out, &out->last_checksum);

  /* Some useful information */
  if(QIO_verbosity() >= QIO_VERB_REG) {
    printf("%s(%d): Wrote field. datatype %s recordtype %d \n              precision %s colors %d spins %d count %d\n",
	   myname,this_node,
	   QIO_get_datatype(record_info),
	   QIO_get_recordtype(record_info),
	   QIO_get_precision(record_info),
	   QIO_get_colors(record_info),
	   QIO_get_spins(record_info),
	   QIO_get_datacount(record_info));
    
    printf("%s(%d): checksums %0x %0x\n",
	   myname,this_node,
	   QIO_get_writer_last_checksuma(out),
	   QIO_get_writer_last_checksumb(out));
  }

  return status;
}
