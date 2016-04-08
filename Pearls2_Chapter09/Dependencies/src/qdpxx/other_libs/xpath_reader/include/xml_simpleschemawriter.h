#ifndef XML_SIMPLESCHEMAWRITER_H
#define XML_SIMPLESCHEMAWRITER_H

#include <xml_attribute.h>
#include <ostream>
#include <string>
#include <stack>
#include <xml_writer.h>

namespace XMLWriterAPI {

  // Base class for the XMLSimpleWriter classes
  class XMLSimpleSchemaWriter {
  public:
    XMLSimpleSchemaWriter(XMLSimpleWriter& _instance, 
			  XMLSimpleWriter& _schema) : instance(_instance), 
      schema(_schema) { 


      // Construct attributes for schema element 
      AttributeList alist;
      alist.push_back(Attribute((std::string)"xmlns", (std::string)"xs", (std::string)"http://www.w3.org/2001/XMLSchema"));
      schema.openTag("xs", "schema", alist);
    }

    void writeSimpleElement(const std::string& name, const std::string& value) {
      instance.openTag(name);
      instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:string"));
      schema.emptyTag("xs", "element", schema_alist);
    }

    void writeSimpleElement(const std::string& name, const int& value) {
      instance.openTag(name);
      instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:int"));
      schema.emptyTag("xs", "element", schema_alist);
    }
     
   void writeSimpleElement(const std::string& name, const unsigned int& value) {
      instance.openTag(name);
      instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:unsignedInt"));
      schema.emptyTag("xs", "element", schema_alist);
   }

   void writeSimpleElement(const std::string& name, const short int& value) {
      instance.openTag(name);
      instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:short"));
      schema.emptyTag("xs", "element", schema_alist);
    }

   void writeSimpleElement(const std::string& name, const unsigned short int& value) {
      instance.openTag(name);
      instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:unsignedShort"));
      schema.emptyTag("xs", "element", schema_alist);
    }

   void writeSimpleElement(const std::string& name, const long int& value) {
      instance.openTag(name);
      instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:long"));
      schema.emptyTag("xs", "element", schema_alist);
   }

   void writeSimpleElement(const std::string& name, const unsigned long int& value) {
     instance.openTag(name);
     instance.write(value);
     instance.closeTag();
     AttributeList schema_alist;
     schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
     schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:unsignedLong"));
     schema.emptyTag("xs", "element", schema_alist);
   }

   void writeSimpleElement(const std::string& name, const float& value) {
     instance.openTag(name);
     instance.write(value);
     instance.closeTag();
     AttributeList schema_alist;
     schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
     schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:float"));
     schema.emptyTag("xs", "element", schema_alist);
   } 
   
   void writeSimpleElement(const std::string& name, const double& value) {
     instance.openTag(name);
     instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:double"));
      schema.emptyTag("xs", "element", schema_alist);
   } 

   void writeSimpleElement(const std::string& name, const bool& value) {
     instance.openTag(name);
     instance.write(value);
      instance.closeTag();
      AttributeList schema_alist;
      schema_alist.push_back(Attribute((std::string)"name", (std::string)name));
      schema_alist.push_back(Attribute((std::string)"type", (std::string)"xs:boolean"));
      schema.emptyTag("xs", "element", schema_alist);
   } 
   
   void openComplexElement(const std::string& name) {
     instance.openTag(name);
     
     AttributeList alist;
     alist.push_back(Attribute((std::string)"name", name));
     schema.openTag("xs", "element", alist);
     schema.openTag("xs", "complexType");
     schema.openTag("xs", "sequence");
   }
   
   void closeComplexElement(void) {
     instance.closeTag();
     schema.closeTag(); // sequence
     schema.closeTag(); // complexType
     schema.closeTag(); // element
   }

  protected:
    XMLSimpleWriter& instance;
    XMLSimpleWriter& schema;
  };
};


#endif



