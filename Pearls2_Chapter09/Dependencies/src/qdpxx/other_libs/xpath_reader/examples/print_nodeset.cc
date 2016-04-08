#include <basic_xpath_reader.h>
#include <xpath_reader.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <xml_array.h>

using namespace std;
using namespace XMLXPathReader;

#define USE_BASIC

int main(int argc, char *argv[])
{
  // Instantiate a reader
#if defined(USE_BASIC)
  BasicXPathReader reader;
#else
  XPathReader reader;
#endif
 
  if (argc != 3)
  {
    cerr << "Usage: print_xpath  <xml file>  <xpath query>" << endl;
    exit(1);
  }

  // Try and open the reader
  try { 
    reader.open(argv[1]);
  } 
  catch (string &error_mesg) {
    cerr << error_mesg << endl;
    throw;
  }

  // Try and get a string 
  string sresult ="";  
  try {
#if defined(USE_BASIC)
    reader.evaluateXPath(string(argv[2]));
#else
    reader.getXPath<string>(argv[2], sresult);
    cout << sresult << endl;
#endif
  }
  catch (const string& e) {
    cout << e << endl;
    exit(1);
  }
    
#if defined(USE_BASIC)
  // If no error occurred, and also we didn't quit then print the query
  // results.
  reader.printQueryResult(cout);
#endif

  /* This bit kills the last refcount */
  reader.close();

  return EXIT_SUCCESS;
}

