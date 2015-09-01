
#include <xml_simplewriter.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <ios>

#define INDENT_SPACES ((unsigned int)2)

using namespace XMLWriterAPI;

void XMLSimpleWriter::openTag(const std::string& tagname)
{
  std::string empty_string;
  AttributeList empty_list;

  dumpTag(empty_string, tagname, empty_list, false);
}

void XMLSimpleWriter::openTag(const std::string& nsprefix, const std::string& tagname)
{
  AttributeList empty_list;

  dumpTag(nsprefix, tagname, empty_list, false);
}

void XMLSimpleWriter::openTag(const std::string& tagname, AttributeList& al)
{
  std::string empty_string;
  dumpTag(empty_string, tagname, al, false);
}

void XMLSimpleWriter::openTag(const std::string& nsprefix, 
			    const std::string& tagname, 
			    AttributeList& al)
{
  dumpTag(nsprefix, tagname, al, false);
}

void XMLSimpleWriter::emptyTag(const std::string& tagname)
{
  std::string empty_string;
  AttributeList empty_list;

  dumpTag(empty_string, tagname, empty_list, true);
}

void XMLSimpleWriter::emptyTag(const std::string& tagname,  AttributeList& al)
{
  std::string empty_string;
  dumpTag(empty_string, tagname, al, true);
}

void XMLSimpleWriter::emptyTag(const std::string& nsprefix, 
			     const std::string& tagname, 
			     AttributeList& al)
{
  dumpTag(nsprefix, tagname, al, true);
}


void XMLSimpleWriter::dumpTag(const std::string& nsprefix, 
			    const std::string& tagname, 
			    AttributeList& al,
			    bool is_empty)
{
  std::string qualified_tagname;
  std::ostream& os=getOstream();

  // Check whether we are trying to write a second 
  // root element. I thought this was allowed but apparently not
  if( doctag_written == true && namestack.empty() ) {
    std::ostringstream error_message;
    error_message << "Attempt to write second root tag -- this is not XML compliant: tagname = " 
		  << tagname << std::endl;
    throw error_message.str();
  }

  if( nsprefix.size() == 0 ) { 
    qualified_tagname = tagname;
  }
  else { 
    qualified_tagname = nsprefix + ":" + tagname;
  }
  
  std::string indent(indent_level*INDENT_SPACES, ' ');
  os << "\n" <<  indent << "<" << qualified_tagname;
  
  std::list<Attribute>::iterator the_iterator;
  for(the_iterator = al.begin(); the_iterator != al.end(); the_iterator++) {
    if( (*the_iterator).isEmpty() == false ) {
      os << "  " << the_iterator->getName() << "=\"" << the_iterator->getValue()
	 << "\"";
    }
  }
  
  if(is_empty == true) { 
    os << "/>";
  }
  else {
    os << ">" ;
    namestack.push(qualified_tagname);
    indent_level++;
  }
  
  if( doctag_written == false ) { 
    doctag_written = true; 
  }

  primitive_last = false;
}
  
void XMLSimpleWriter::closeTag(void)
{
  std::string qualified_tagname;
  std::ostream& os=getOstream();

  if( namestack.empty() == false ) { 
    qualified_tagname = namestack.top();
    namestack.pop();
    indent_level--;

    if(primitive_last == false) {
      std::string indent(indent_level*INDENT_SPACES, ' ');
      os << "\n" << indent;
    }
    os << "</" << qualified_tagname << ">" ;
  }
  else {
    std::ostringstream error_message;
    error_message << "Attempt to close non existent tag";
    throw error_message.str();
  }

  primitive_last = false;

}


void 
XMLSimpleWriter::write(const std::string& output)
{
  writePrimitive<std::string>(output);
}

void
XMLSimpleWriter::write(const int& output) 
{
  writePrimitive<int>(output);
}

void
XMLSimpleWriter::write(const unsigned int& output)
{
  writePrimitive<unsigned int>(output);
}

void
XMLSimpleWriter::write(const short int& output)
{
  writePrimitive<short int>(output);
}


void 
XMLSimpleWriter::write(const unsigned short int& output)
{
  writePrimitive<unsigned short int>(output);
}

void
XMLSimpleWriter::write(const long int& output)
{
  writePrimitive<long int>(output);
}

void 
XMLSimpleWriter::write(const unsigned long int& output)
{
  writePrimitive<unsigned long int>(output);
}

void 
XMLSimpleWriter::write(const float& output)
{
  writePrimitive(output);
}

void 
XMLSimpleWriter::write(const double& output)
{
  writePrimitive(output);
}

void
XMLSimpleWriter::write(const bool& output)
{
  writePrimitive<bool>(output);
}

using namespace std;
void
XMLSimpleWriter::writePrimitive(const float& output)
{
  std::ostream& os=getOstream();

  if( ! namestack.empty() ) { 
    streamsize initPrec = os.precision();
    os.precision(7);
    os << output;
    os.precision(initPrec);
  }
  else { 
    std::ostringstream error_string;
    error_string << "Attempt to write before opening root tag" << std::endl;
    throw error_string.str();
  }

  primitive_last = true;
}


void
XMLSimpleWriter::writePrimitive(const double& output)
{
  std::ostream& os=getOstream();

  if( ! namestack.empty() ) { 
    streamsize initPrec = os.precision();
    os.precision(15);
    os << output;
    os.precision(initPrec);
  }
  else { 
    std::ostringstream error_string;
    error_string << "Attempt to write before opening root tag" << std::endl;
    throw error_string.str();
  }

  primitive_last = true;
}

template < typename T > 
void
XMLSimpleWriter::writePrimitive(const T& output)
{
  std::ostream& os=getOstream();

  if( ! namestack.empty() ) { 

    os << std::boolalpha << output;
  }
  else { 
    std::ostringstream error_string;
    error_string << "Attempt to write before opening root tag" << std::endl;
    throw error_string.str();
  }

  primitive_last = true;
}


void 
XMLSimpleWriter::writeXML(const std::string& output)
{
  std::ostream& os=getOstream();

  // Should this have more checking, e.g. of primitive_last or doctag?
  os << output << "\n";
}

void 
XMLSimpleWriter::writePrologue(std::ostream& os) const
{
   os << "<?xml version=\"1.0\"?>" << "\n";
   os << "\n";
}







