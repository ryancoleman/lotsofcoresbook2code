/* ID: $Id: basic_xpath_reader_test3.cc,v 1.1 2004-04-27 11:22:59 bjoo Exp $
 *
 * file: basic_xpath_reader2.cc
 *
 * a second example of using BasicXPathReader -- now to bind 
 * primitive types 
 */

#include <basic_xpath_reader.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>


using namespace std;
using namespace XMLXPathReader;



int main(int argc, char *argv[])
{

  // Get the reader
  BasicXPathReader reader;
  // Dummy xml file
  std::string xml_input="<?xml version=\"1.0\"?><foo><bar>true</bar><bar2>true false true</bar2></foo>";
  // Turn into an instream
  std::istringstream input_is(xml_input);

  // Open the instream
  try { 
    reader.open(input_is);
  } 
  catch (string &error_mesg) {
    cerr << error_mesg << endl;
    cout << "FAILURE" << endl;
    throw;
  }

  bool fred;
  try {
    reader.get("/foo/bar", fred);
  }
  catch(string& error_mesg) { 
    cerr << error_mesg << endl;
    cout << "FAILURE" << endl;
    throw;
  }

  cout << "SUCCESS" << endl;

  return EXIT_SUCCESS;
}

