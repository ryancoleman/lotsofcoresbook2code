/* QIO_close_write.c */

#include <qio_config.h>
#include <qio.h>
#include <lrl.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

int QIO_close_write(QIO_Writer *out)
{
  /* Prevent premature file truncation in parallel writes */
  /* Note, the test will cause a hang if the oflag->serpar value is
     not the same for all nodes */
  if(out->serpar == QIO_PARALLEL) DML_sync();

  int status = LRL_close_write_file(out->lrl_file_out);
  if(out->layout) {
    free(out->layout->latsize);
    if(out->layout->hyperupper) free(out->layout->hyperupper);
    if(out->layout->hyperlower) free(out->layout->hyperlower);
    free(out->layout);
  }
  DML_free_sitelist(out->sites);
  if(out->ildgLFN != NULL) QIO_string_destroy(out->ildgLFN); 
  free(out);
  if(status != LRL_SUCCESS) return QIO_ERR_CLOSE;
  return QIO_SUCCESS;
}
