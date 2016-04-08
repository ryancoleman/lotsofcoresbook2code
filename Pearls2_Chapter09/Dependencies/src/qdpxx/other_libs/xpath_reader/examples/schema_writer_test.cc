#include <xml_writer.h>
#include <xml_schemawriter.h>

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

  XMLSchemaWriter writer(instance,schema);

  writer.write("jimbob", (string)"foo");



  return(EXIT_SUCCESS);
}

