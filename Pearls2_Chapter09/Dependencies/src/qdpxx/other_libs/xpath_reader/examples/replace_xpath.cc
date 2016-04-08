#include <basic_xpath_reader.h>
#include <xpath_reader.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <xml_array.h>

using namespace std;
using namespace XMLXPathReader;



int main(int argc, char *argv[])
{
  // Instantiate a reader
  BasicXPathReader reader;
 
  if (argc != 4)
  {
    cerr << "Usage: print_xpath  <xml file>  <xpath query> <replace_string>" << endl;
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


    reader.setPrimitiveString(string(argv[2]),string(argv[3]));

  }
  catch (const string& e) {
    cout << e << endl;
    exit(1);
  }
    
  reader.printRoot(cout);

  /* This bit kills the last refcount */
  reader.close();

  return EXIT_SUCCESS;
}

