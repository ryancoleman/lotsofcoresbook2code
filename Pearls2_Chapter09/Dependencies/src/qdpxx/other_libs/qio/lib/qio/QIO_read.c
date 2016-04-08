/* QIO_read.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <dml.h>
#include <qio_string.h>
#include <qioxml.h>
#include <stdio.h>
#include <string.h>

/* Reads a lattice field or global data from a record.  Includes XML
   and checksum */
/* Calls QIO_read_record_info and QIO_read_record_data */
/* Caller must allocate *record_info and *xml_record.
   Caller must signal abort to all nodes upon failure. */
int
QIO_read(QIO_Reader *in, QIO_RecordInfo *record_info,
	 QIO_String *xml_record, 
	 void (*put)(char *buf, size_t index, int count, void *arg),
	 size_t datum_size, int word_size, void *arg)
{
  int status;
  int this_node = in->layout->this_node;
  char myname[] = "QIO_read";

  /* Read info if not already done */
  status = QIO_read_record_info(in, record_info, xml_record);

  if(QIO_verbosity() >= QIO_VERB_DEBUG) {
    printf("%s(%d): QIO_read_record_info returned %d\n",
	   myname,this_node,status);
    fflush(stdout);
  }

  if(status!=QIO_SUCCESS) return status;

  /* Read data */
  status = QIO_read_record_data(in, put, datum_size, word_size, arg);

  if(QIO_verbosity() >= QIO_VERB_DEBUG) {
    printf("%s(%d): QIO_read_record_data returned %d\n",
	   myname,this_node,status);
  }

  return status;
}
