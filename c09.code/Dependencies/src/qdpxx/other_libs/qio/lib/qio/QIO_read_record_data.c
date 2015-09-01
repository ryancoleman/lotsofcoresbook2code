/* QIO_read_record_data.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <qioxml.h>
#include <stdio.h>
#include <string.h>

/* Utility for printing the record info */
void QIO_print_record_info(QIO_RecordInfo *record_info){

  printf("Record header: datatype %s recordtype %d \n",
	 QIO_get_datatype(record_info),
	 QIO_get_recordtype(record_info));
  printf("precision %s colors %d spins %d count %d\n",
	 QIO_get_precision(record_info),
	 QIO_get_colors(record_info),
	 QIO_get_spins(record_info),
	 QIO_get_datacount(record_info));
  
  if(QIO_get_recordtype(record_info) == QIO_HYPER){
    int i;
    int n = QIO_get_hyper_spacetime(record_info);
    int *lower = QIO_get_hyperlower(record_info);
    int *upper = QIO_get_hyperupper(record_info);
    printf("Hypercube lower");
    for(i = 0; i < n; i++){
      printf(" %d",lower[i]);
    }
    printf("\n");
    printf("Hypercube upper");
    for(i = 0; i < n; i++){
      printf(" %d",upper[i]);
    }
    printf("\n");
  }
}


/* Read record data */
/* Can be called separately from QIO_read, but then QIO_read_record_info
   must be called first */

/* On failure calling program must signal abort to all nodes. */

int QIO_generic_read_record_data(QIO_Reader *in, 
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     size_t datum_size, int word_size, void *arg,
  	     DML_Checksum *checksum, uint64_t *nbytes)
{
  char myname[] = "QIO_generic_read_record_data";
  int count;
  LRL_RecordReader *lrl_record_in;
  size_t datum_size_info;
  DML_Layout *layout = in->layout;
  int this_node = layout->this_node; 
  int status;
  QIO_RecordInfo *record_info = &in->record_info;

  /* List of acceptable binary data LIME types */
  int ntypes = 2;
  /* Avoid a compiler bug */
  LIME_type lime_type0 = QIO_LIMETYPE_BINARY_DATA;
  LIME_type lime_type1 = QIO_LIMETYPE_ILDG_BINARY_DATA;
  LIME_type lime_type_list[2] = {lime_type0, lime_type1};
  /* LIME_type lime_type_list[2] = {
      QIO_LIMETYPE_BINARY_DATA,
      QIO_LIMETYPE_ILDG_BINARY_DATA
      }; */
  LIME_type lime_type;


  /* It is an error to call for the data before reading the info
     in a given record */
  if(in->read_state != QIO_RECORD_DATA_NEXT){
    printf("%s(%d): Bad read state %d\n",myname,this_node,in->read_state);
    return QIO_ERR_INFO_MISSED;
  }

  /* Require consistency between the byte count obtained from the
     private record metadata and the datum_size and count parameters
     (per site for field or hypercube data or total for global data) */
  count = QIO_get_datacount(record_info);
  if(datum_size !=  QIO_get_typesize(record_info) * count)
    {
      printf("%s(%d): requested byte count %lu disagrees with the record %d * %d\n",
	     myname,this_node,(unsigned long)datum_size,
	     QIO_get_typesize(record_info), count);
      
      if(this_node == layout->master_io_node)
	QIO_print_record_info(record_info);
      
      return QIO_ERR_BAD_READ_BYTES;
    }
  
  
  /* Verify byte count per site (for field or hypercube) or total (for
     global) */
  datum_size_info = QIO_get_typesize(record_info) * count;

  if(datum_size != datum_size_info){
    printf("%s(%d): byte count mismatch request %lu != actual %lu\n",
	   myname, this_node, (unsigned long)datum_size, 
	   (unsigned long)datum_size_info);

    return QIO_ERR_BAD_READ_BYTES;
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Calling QIO_open_read_field\n",myname,this_node);fflush(stdout);
  }

  /* Nodes read the field */

  /* Scan ahead and open the record with one of the listed LIME types */
  lrl_record_in = QIO_open_read_field(in, datum_size, 
         lime_type_list, ntypes, &lime_type, &status);
  if(status != QIO_SUCCESS)return status;

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Calling QIO_read_field_data\n",
	   myname,this_node);fflush(stdout);
  }

  /* Then read the data and close the record */
  status = QIO_read_field_data(in, lrl_record_in, 
			       put, count, datum_size, word_size, 
			       arg, checksum, nbytes);
  if(status != QIO_SUCCESS){
    printf("%s(%d): Error reading field data\n",myname,this_node);
    return status;
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): done reading field\n",myname,this_node);fflush(stdout);
  }

  /* Copy most recent node checksum into reader */
  memcpy(&(in->last_checksum), checksum, sizeof(DML_Checksum));

  if(this_node == layout->master_io_node){
    if(QIO_verbosity() >= QIO_VERB_REG){
      printf("%s(%d): Read field\n",myname,this_node);
      QIO_print_record_info(record_info);
    }
  }
      
  in->read_state = QIO_RECORD_CHECKSUM_NEXT;
  return QIO_SUCCESS;
}

QIO_ChecksumInfo *QIO_read_checksum(QIO_Reader *in)
{
  char myname[] = "QIO_read_checksum";
  QIO_String *xml_checksum;
  QIO_ChecksumInfo *checksum_info_expect = NULL;
  int this_node = in->layout->this_node;
  int status;
  LIME_type lime_type=NULL;
  
  if(in->read_state != QIO_RECORD_CHECKSUM_NEXT){
    printf("%s(%d): Bad read state %d\n",myname,this_node,in->read_state);
    return NULL;
  }
  
  /* Master node reads the checksum record */
  status = QIO_SUCCESS;  /* Changed if checksum does not match */
  if(this_node == in->layout->master_io_node){
    checksum_info_expect = QIO_create_checksum_info(0,0);
    /* No checksum record for non-native files */
    if(in->format == QIO_SCIDAC_NATIVE){
      xml_checksum = QIO_string_create();
      QIO_string_realloc(xml_checksum,QIO_STRINGALLOC);
      if((status=QIO_read_string(in, xml_checksum, &lime_type))
	 !=QIO_SUCCESS){
	printf("%s(%d): Error reading checksum\n",myname,this_node);
	return NULL;
      }
      if(QIO_verbosity() >= QIO_VERB_DEBUG){
	printf("%s(%d): checksum = %s\n",myname,this_node,
	       QIO_string_ptr(xml_checksum));
      }
      /* Extract checksum */
      if((status=QIO_decode_checksum_info(checksum_info_expect, xml_checksum))
	 !=0){
	printf("%s(%d): bad checksum record\n",myname,this_node);
	return NULL;
      }
      QIO_string_destroy(xml_checksum);
    }
  }
  
  in->read_state = QIO_RECORD_INFO_PRIVATE_NEXT;
  
  return checksum_info_expect;
}

int QIO_compare_checksum(int this_node, 
	 QIO_ChecksumInfo *checksum_info_expect, DML_Checksum *checksum)
{
  char myname[] = "QIO_compare_checksum";
  QIO_ChecksumInfo *checksum_info_found;
  int status;
  
  /* Compare checksums */
  checksum_info_found = QIO_create_checksum_info(checksum->suma, 
						 checksum->sumb);
  status = QIO_compare_checksum_info(checksum_info_found, 
				     checksum_info_expect,
				     myname,this_node);
  QIO_destroy_checksum_info(checksum_info_found);
  
  return status;
}

int QIO_read_record_data(QIO_Reader *in, 
	     void (*put)(char *buf, size_t index, int count, void *arg),
	     size_t datum_size, int word_size, void *arg){


  char myname[] = "QIO_read_record_data";
  DML_Checksum checksum;
  uint64_t nbytes, expect_bytes;
  QIO_ChecksumInfo *checksum_info_expect;
  size_t buf_size;
  size_t volume;
  int this_node = in->layout->this_node;
  int status;
  int recordtype = QIO_get_recordtype(&(in->record_info));

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): Calling QIO_generic_read_record_data\n",
	   myname,this_node);fflush(stdout);
  }

  status = QIO_generic_read_record_data(in, put, datum_size, word_size, arg,
					&checksum, &nbytes);
  if(status != QIO_SUCCESS)return status;

  /* Total number of sites in this record */
  volume = in->layout->subsetvolume;

  /* Compute the number of bytes read by all nodes */
  DML_sum_uint64_t(&nbytes);

  /* Check that the record size matches the expected size of the data */
  if(recordtype == QIO_GLOBAL){
    buf_size = datum_size; /* Global data */
    expect_bytes = buf_size;
  }
  else { /* Field data */
    expect_bytes = ((uint64_t)volume) * datum_size;
  }
  
  if(nbytes != expect_bytes){
    printf("%s(%d): bytes read %lu != expected rec_size %lu\n",
	   myname, in->layout->this_node,
	   (unsigned long)nbytes, (unsigned long)expect_bytes);
    return QIO_ERR_BAD_READ_BYTES;
  }
  
  /* Combine checksums over all nodes */
  DML_checksum_combine(&checksum);

  /* Copy most recent combined checksum into reader */
  memcpy(&(in->last_checksum), &checksum, sizeof(DML_Checksum));

  /* Everyone calls this. It changes the state machine */
  checksum_info_expect = QIO_read_checksum(in);
  if(this_node == in->layout->master_io_node)
    {
      if(in->format == QIO_SCIDAC_NATIVE){
	if(checksum_info_expect == NULL)return QIO_ERR_CHECKSUM_INFO;
	status = QIO_compare_checksum(this_node, 
				      checksum_info_expect, &checksum);
	if(status != QIO_SUCCESS)return status;
      }
    }

  /* Checksum info expect may be NULL on non master node.
     This is not an error. Freeing it is */
  if( checksum_info_expect != NULL ) { 
    QIO_destroy_checksum_info(checksum_info_expect);
  }

  if(QIO_verbosity() >= QIO_VERB_DEBUG){
    printf("%s(%d): done with QIO_read_record_data\n",
	   myname,in->layout->this_node);fflush(stdout);
  }
  
  return QIO_SUCCESS;
}

