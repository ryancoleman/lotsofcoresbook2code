/* QIO_read_record_info.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <qioxml.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

/* Update the record type and hypercube information */
int QIO_reader_insert_hypercube_data(QIO_Reader *in, 
				     QIO_RecordInfo *record_info)
{

  int status;

  /* Copy subset data from record_info structure to layout structure
     and check it */

  status = DML_insert_subset_data(in->layout, 
				  QIO_get_recordtype(record_info),
				  QIO_get_hyperlower(record_info),
				  QIO_get_hyperupper(record_info),
				  QIO_get_hyper_spacetime(record_info));
  if(status != 0) return QIO_ERR_BAD_SUBSET;

  return QIO_SUCCESS;
}

/* Read user record XML */
/* Can be called separately from QIO_read for discovering the record
   contents without reading the field itself */

int QIO_read_private_record_info(QIO_Reader *in, QIO_RecordInfo *record_info)
{
  /* Caller must allocate *record_info */

  QIO_String *xml_record_private;
  int this_node = in->layout->this_node;
  int status;
  QIO_RecordInfo *alien_info;
  char myname[] = "QIO_read_private_record_info";
  LIME_type lime_type=NULL;
  
  /* Read private record XML if not already done */
  if(in->read_state == QIO_RECORD_INFO_PRIVATE_NEXT){

    /* Allocate space for internal copy of user record XML */
    /* This may already be in use, if so free it */
    if( in->xml_record != NULL ) { 
      QIO_string_destroy(in->xml_record);
    }

    in->xml_record = QIO_string_create();
    
    /* Master node reads and interprets the private record XML record */
    if(this_node == in->layout->master_io_node){
      if(in->format == QIO_SCIDAC_NATIVE){
	/* Initialize private record XML - will be allocated by read_string */
	xml_record_private = QIO_string_create();
	
	/* Read private record XML */
	if((status=QIO_read_string(in, xml_record_private, &lime_type ))
	   != QIO_SUCCESS)return status;
	
	if(QIO_verbosity() >= QIO_VERB_DEBUG){
	  printf("%s(%d): private XML = \"%s\"\n",myname,this_node,
		 QIO_string_ptr(xml_record_private));
	}
	
	/* Decode the private record XML */
	status = QIO_decode_record_info(&(in->record_info), xml_record_private);
	if(status != 0)return QIO_ERR_PRIVATE_REC_INFO;

	/* Free storage */
	QIO_string_destroy(xml_record_private);
      }
      /* For alien ILDG formats we create the members of the record
	 info structure, using the ildg_precision value we have read from
	 the ILDG format record */
      else{
	if(in->ildg_precision == 32)
	  /* Single precision SU(3) matrix */
	  alien_info = 
	    QIO_create_record_info(QIO_FIELD, NULL, NULL, 0,
				   "USQCD_F3_ColorMatrix",
				   "F", 3, 0, 72, 4);
	else
	  alien_info = 
	    QIO_create_record_info(QIO_FIELD, NULL, NULL, 0,
				   "USQCD_D3_ColorMatrix",
				   "D", 3, 0, 144, 4);
	memcpy(&(in->record_info), alien_info, sizeof(QIO_RecordInfo));
	QIO_destroy_record_info(alien_info);
      }
    }
    /* Set state in case record is reread */
    in->read_state = QIO_RECORD_INFO_USER_NEXT;
  }

  // broadcast just so everyone has it (possibly not needed)
  DML_broadcast_bytes((char *)&(in->record_info), sizeof(QIO_RecordInfo),
		      this_node, in->layout->master_io_node);

  /* Copy record info on all calls */
  /*  *record_info = in->record_info; */
  memcpy(record_info, &(in->record_info), sizeof(QIO_RecordInfo));
  
  /* Copy hypercube data into reader */
  status = QIO_reader_insert_hypercube_data(in, record_info);
  if(status != QIO_SUCCESS)return status;

  return QIO_SUCCESS;
}

int QIO_read_user_record_info(QIO_Reader *in, QIO_String *xml_record){
  /* Caller must allocate *xml_record */

  int this_node = in->layout->this_node;
  int status;
  char myname[] = "QIO_read_user_record_info";
  LIME_type lime_type=NULL;
  
  /* Read user record XML if not already done */
  if(in->read_state == QIO_RECORD_INFO_USER_NEXT){
    /* We don't expect alien ILDG files to have this record */
    if(in->format == QIO_SCIDAC_NATIVE){
      /* Master node reads the user XML record */
      if(this_node == in->layout->master_io_node){
	if((status=QIO_read_string(in, in->xml_record, &lime_type))
	   != QIO_SUCCESS){
	  printf("%s(%d): Error reading user record XML\n",myname,this_node);
	  return status;
	}
      }
    }
    else{
      QIO_string_set(in->xml_record,QIO_NONSCIDAC_RECORD);
    }
    
    if(QIO_verbosity() >= QIO_VERB_DEBUG && 
       QIO_string_ptr(in->xml_record) != NULL){
      printf("%s(%d): user XML = \"%s\"\n",myname,this_node,
	     QIO_string_ptr(in->xml_record));
    }
    
    /* Set state in case record is being reread */
    in->read_state = QIO_RECORD_ILDG_INFO_NEXT;
  }

  /* Copy user record info (for this node) */
  if(xml_record != NULL)
    QIO_string_copy(xml_record,in->xml_record);

  return QIO_SUCCESS;
}

    
/* Look for and read the ILDG format record and the ILDG
   logical file name record, if present.  The result is placed in
   the string member of the reader. */

int QIO_read_ILDG_LFN(QIO_Reader *in){

  int this_node = in->layout->this_node;
  int status;
  off_t offset;
  uint64_t rec_size;
  LRL_RecordReader *lrl_record_in = NULL;
  char myname[] = "QIO_read_ILDG_LFN";
  LIME_type lime_type=NULL;
  
  if(in->read_state == QIO_RECORD_ILDG_INFO_NEXT){
    /* Only the master I/O node reads in any event */
    /* At present we don't try to extract the LFN from a nonnative file */
    if(this_node == in->layout->master_io_node &&
       in->format == QIO_SCIDAC_NATIVE){

      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): looking for ILDG LFN \n",myname,this_node);
      
      /* Mark the current reader pointer and read the next record
	 header */
      offset = QIO_get_reader_pointer(in);
      lrl_record_in = LRL_open_read_record(in->lrl_file_in, &rec_size, 
					   &lime_type, &status);
      if(status != QIO_SUCCESS)return status;
      /* We don't actually read the ILDG format record */
      LRL_close_read_record(lrl_record_in);
      
      if(QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): found LIME type %s\n",myname,this_node,lime_type);
      
      /* If this is the ILDG format record, indicate we are reading an
	 ILDG format file and read on */
      if(strcmp(lime_type,QIO_LIMETYPE_ILDG_FORMAT) == 0){
	in->ildgstyle = QIO_ILDGLAT;
	
	/* Update the current reader pointer */
	offset = QIO_get_reader_pointer(in);
	
	/* Look at next record */
	lrl_record_in = LRL_open_read_record(in->lrl_file_in, &rec_size, 
					     &lime_type, &status);
	if(status != QIO_SUCCESS)return status;
	/* We will reread the record later */
	LRL_close_read_record(lrl_record_in);
	
	if(QIO_verbosity() >= QIO_VERB_DEBUG)
	  printf("%s(%d): found LIME type %s\n",myname,this_node,lime_type);
	
	/* If this is the ILDG LFN, grab it if requested. */
	if(strcmp(lime_type,QIO_LIMETYPE_ILDG_DATA_LFN) == 0 &&
	   in->ildgLFN != NULL) {
	  /* Back up to reread */
	  status = QIO_set_reader_pointer(in,offset);
	  if(status != QIO_SUCCESS)return status;
	  
	  status=QIO_read_string(in, in->ildgLFN, &lime_type);
	  if(status != QIO_SUCCESS){
	    printf("%s(%d): Error reading ILDG LFN\n",myname,this_node);
	    return status;
	  }
	  
	  /* Update the current reader pointer */
	  offset = QIO_get_reader_pointer(in);
	  
	  if(QIO_verbosity() >= QIO_VERB_DEBUG && 
	     QIO_string_ptr(in->ildgLFN) != NULL){
	    printf("%s(%d): ILDG LFN = \"%s\"\n",myname,this_node,
		   QIO_string_ptr(in->ildgLFN));
	  }
	}
      }
      
      /* Restore reader position */
      status = QIO_set_reader_pointer(in,offset);
    }
    else { 
      if(this_node == in->layout->master_io_node && 
	 QIO_verbosity() >= QIO_VERB_DEBUG)
	printf("%s(%d): Currently not grokking ILDG LFN from non SciDAC produced files\n",myname,this_node);     
    }
    
    /* Set state in case record is being reread */
    in->read_state = QIO_RECORD_DATA_NEXT;
  }

  return QIO_SUCCESS;
}

    
/* Read private and user record XML on compute node(s) */
/* Can be called separately from QIO_read for discovering the record
   contents without reading the field itself */

int QIO_read_record_info(QIO_Reader *in, QIO_RecordInfo *record_info,
			 QIO_String *xml_record){
  /* Caller must allocate *xml_record and *record_info */

  int this_node = in->layout->this_node;
  int length;
  int status;
  char myname[] = "QIO_read_record_info";
  
  /* Read private record XML if not already done */

  status = QIO_read_private_record_info(in, record_info);
  if(status != QIO_SUCCESS)return status;

  /* Broadcast the private record data to all nodes */
  DML_broadcast_bytes((char *)&(in->record_info), 
		      sizeof(QIO_RecordInfo), this_node, 
		      in->layout->master_io_node);
  
  /* Return the record info to caller */
  memcpy(record_info, &(in->record_info), sizeof(QIO_RecordInfo));

  /* Add record type and subset data to layout structure */
  QIO_reader_insert_hypercube_data(in, record_info);

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Done broadcasting private record info\n",
	   myname,this_node);
  }

  /* Read user record XML if not already done */

  status = QIO_read_user_record_info(in, xml_record);
  if(status != QIO_SUCCESS)return status;
  
  /* Broadcast the user xml record to all nodes */
  /* First broadcast length */

  length = QIO_string_length(in->xml_record);
  DML_broadcast_bytes((char *)&length, sizeof(int), this_node, 
		      in->layout->master_io_node);

  /* Receiving nodes resize their strings */
  /* if(this_node != in->layout->master_io_node){ */

  /* QIO_string_realloc is supposedly non-destructive. Can do it on
     all nodes */
  QIO_string_realloc(in->xml_record,length);

  DML_broadcast_bytes(QIO_string_ptr(in->xml_record),length,
		      this_node, in->layout->master_io_node);

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Done broadcasting user record info \"%s\"\n",
	   myname,this_node,QIO_string_ptr(in->xml_record));
  }

  /* Read the ILDG LFN if present and not already done */

  status = QIO_read_ILDG_LFN(in);
  if(status != QIO_SUCCESS)
    return status;

  /* Broadcast the ILDG LFN to all nodes */
  if(in->ildgLFN != NULL){
    /* First broadcast length */

    length = QIO_string_length(in->ildgLFN);
    DML_broadcast_bytes((char *)&length, sizeof(int), this_node, 
			in->layout->master_io_node);
    if(length > 0){
      /* Receiving nodes resize their strings */
      /* QIO_string_realloc is supposedly non-destructive. Can do it on
	 all nodes */
      QIO_string_realloc(in->ildgLFN,length);

      DML_broadcast_bytes(QIO_string_ptr(in->ildgLFN),length,
			  this_node, in->layout->master_io_node);
      
      if(QIO_verbosity() >= QIO_VERB_DEBUG){
	 printf("%s(%d): Done broadcasting ILDG LFN \"%s\"\n",
		myname,this_node,QIO_string_ptr(in->ildgLFN));
	 }
    }
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Finished\n",myname,this_node);
  }

  /* Copy user record info and ILDG LFN (for all nodes)*/
  if(xml_record != NULL)
    QIO_string_copy(xml_record,in->xml_record);

  return QIO_SUCCESS;
}

