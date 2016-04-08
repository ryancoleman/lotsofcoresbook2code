/* QIO_next_record.c */

/* Skip to beginning of next field (LIME message) in file */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#include <stdio.h>

int QIO_next_record(QIO_Reader *in){
  if(LRL_next_message(in->lrl_file_in) != LRL_SUCCESS)
    return QIO_ERR_SKIP;
  in->read_state = QIO_RECORD_INFO_PRIVATE_NEXT;
  return QIO_SUCCESS;
}
