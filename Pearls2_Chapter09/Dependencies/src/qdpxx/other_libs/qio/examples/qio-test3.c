/* MULTIFILE test of QIO */

/* C. DeTar */
/* October 18, 2004 */

#include <qio.h>
#include "qio-test.h"

int main(int argc, char *argv[]){
  int output_volfmt = QIO_MULTIFILE;
  int output_serpar = QIO_SERIAL;
  int ildgstyle     = QIO_ILDGNO;
  int input_volfmt  = QIO_MULTIFILE;
  int input_serpar  = QIO_SERIAL;

  return qio_test(output_volfmt, output_serpar, ildgstyle,
		  input_volfmt, input_serpar, argc, argv);
}
