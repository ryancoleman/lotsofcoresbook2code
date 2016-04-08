
INCFLAGS=-I@top_srcdir@/include

AM_CXXFLAGS = $(INCFLAGS) @LIBXML2_CXXFLAGS@

lib_LIBRARIES = libXPathReader.a libxmlWriter.a

libXPathReader_a_SOURCES = basic_xpath_reader.cc

libxmlWriter_a_SOURCES = xml_simplewriter.cc attribute.cc xml_struct_writer.cc
