#include <xml_writer.h>
#include <xml_simpleschemawriter.h>

#include <xpath_reader.h>
#include <stdlib.h>
#include <iostream>
#include <string>
using namespace XMLWriterAPI;
using namespace XMLXPathReader;
using namespace std; 

int main(int argc, char *argv[]) 
{

  XMLSimpleFileWriter schema("testout.xsd");
  XMLSimpleFileWriter instance("testout.xml");

  XMLSimpleSchemaWriter writer(instance,schema);

  writer.openComplexElement("jim");
  writer.openComplexElement("foo");
  writer.openComplexElement("bar");
  writer.writeSimpleElement("LadyGodiva", (string)"Rides Naked through the night");
  writer.closeComplexElement();
  writer.closeComplexElement();
  writer.openComplexElement("fred");
  writer.openComplexElement("bar");
  writer.writeSimpleElement("Fester", (int)50);
  writer.closeComplexElement();
  writer.closeComplexElement();

  writer.closeComplexElement();



  return(EXIT_SUCCESS);
}

