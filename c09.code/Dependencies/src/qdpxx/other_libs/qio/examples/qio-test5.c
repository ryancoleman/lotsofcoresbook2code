/* SINGLEFILE SERIAL out PARALLEL in test of QIO */

/* C. DeTar */
/* June 20, 2005 */

#include <qio.h>
#include "qio-test.h"

int main(int argc, char *argv[]){
  int output_volfmt = QIO_SINGLEFILE;
  int output_serpar = QIO_SERIAL;
  int ildgstyle     = QIO_ILDGNO;
  int input_volfmt  = QIO_UNKNOWN;
  int input_serpar  = QIO_PARALLEL;

  return qio_test(output_volfmt, output_serpar, ildgstyle,
		  input_volfmt, input_serpar, argc, argv);
}
