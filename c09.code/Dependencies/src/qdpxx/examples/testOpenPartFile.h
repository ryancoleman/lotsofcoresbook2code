#ifndef TEST_OPEN_PARTFILE_H
#define TEST_OPEN_PARTFILE_H

#include "qdp.h"
#include "unittest.h"

// Forward declarations
using namespace QDP;

class TestOpenPartFile : public TestCase {
private:
  const std::string file_name;
  const QDP_serialparallel_t serpar;
  const QDP_volfmt_t volfmt;
  XMLBufferWriter file_xml;

public:
  TestOpenPartFile(void) :
    file_name(std::string("t_pf_open.lime")),
    serpar(QDPIO_SERIAL),
    volfmt(QDPIO_PARTFILE)
  {
    push(file_xml, "openTest");
    write(file_xml, "dummy", "dummy string");
    pop(file_xml);
  }

  void run() {
    Assertions::assertNotEquals<QDPFileWriter*>(0x0, 
						new QDPFileWriter(file_xml,
								  file_name,
								  volfmt,
								  serpar,
								  QDPIO_OPEN) );
  }
};
#endif
