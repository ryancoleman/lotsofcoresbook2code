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
  XPathReader reader;

 
  // Try and open the reader
  try { 
  reader.open("foo.xml");
  } 
  catch (string &error_mesg) {
    cerr << error_mesg << endl;
    throw;
  }
  cout << "Document Open Complete" << endl;


  // Register query namespace prefixes
  reader.registerNamespace("bj", "http://www.ph.ed.ac.uk/~bj/");
  reader.registerNamespace("fred", "http://www.ph.ed.ac.uk/~bj2");

  // Try and get a string 
  string sresult ="";  
  try {
    reader.getXPath<string>("/foo/funky", sresult);
    cout << "Query: /foo/funky returned string: " << sresult << endl;
  }
  catch (const string& e) {
    cout << e << endl;
  }
 
  // try and get an int
  int iresult;
  try {
    reader.getXPath<int>("//bj:armadillo", iresult);
    cout << "Query: //bj:armadillo returned int: " << iresult << endl;
  }
  catch (const string& e) {
    cout << e << endl;
  }

  // try and get an int that will fail
  // because armadillo2 is a string
  try {
    reader.getXPath<int>("//armadillo2", iresult);
    cout << "Query: //armadillo2 returned int: " << iresult << endl;
  }
  catch (const string& e) {
    cout << e << endl;
  }

  // try and get an int that will fail because the 
  // tag doesn't contain a single int
  try {
    reader.getXPath<int>("//funky", iresult);
    cout << "Query: //funky  returned int: " << iresult << endl;
  }
  catch (const string& e) {
    cout << e << endl;
  } 

  // try and get the third occurrance of <bar> within foo
  float fresult;
  try { 
    reader.getXPath<float>("/foo/bar[3]", fresult);
    cout << "Query: /foo/bar[3] returned float: " << fresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // try and get a float something unabmiguated by a namespace
  try { 
    reader.getXPath<float>("/foo/bar/fred:kappa", fresult);
    cout << "Query: /foo/bar/fred:kappa returned float: " << fresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // try and get a string value into a float 
  // armadillo2 contains a string -- this will fail.
  try { 
    reader.getXPath<float>("//armadillo2", fresult);
    cout << "Query: //armadillo2 returned float: " << fresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }
  
#if defined(XPATH_READER_TEST_DO_THIS_STUFF)
  // try and get a float out of a tag with more than one number
  // this will fail
  try { 
    reader.getXPath<float>("//funky", fresult);
    cout << "Query: //funky returned float: " << fresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // Try and get the third occurrance (document order) 
  // of <bar> within <foo> as a double
  double dresult;
  try { 
    reader.getXPath<double>("/foo/bar[3]", dresult);
    cout << "Query: /foo/bar[3] returned double: " << dresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // Try ang get a double where the tag is disambiguated by a 
  // namespace
  try { 
    reader.getXPath<double>("/foo/bar/fred:kappa", dresult);
    cout << "Query: /foo/bar/fred:kappa returned double: " << dresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // Try and read a string value into a double
  // this will fail.
  try { 
    reader.getXPath<double>("//armadillo2", dresult);
    cout << "Query: //armadillo2 returned double: " << dresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }
  
  // Try and get a double from a tag with two numbers 
  // this will fail
  try { 
    reader.getXPath<double>("//funky", dresult);
    cout << "Query: //funky returned double: " << dresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // try and get a boolean
  bool bresult;
  try { 
    reader.getXPath<bool>("//multi_quarkP", bresult);
    cout << "Query: //multi_quarkP returned boolean: " << boolalpha << bresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }
  
  // try ang get a boolean from a namespaced tag 
  try { 
    reader.getXPath<bool>("//bj:mybool", bresult);
    cout << "Query: //bj:mybool returned boolean: " << boolalpha << bresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // try and get a boolean from a string. This should work if the 
  // string in the tag is either "true" or "false"
  // which in this case it isn't so this will fail
  try { 
    reader.getXPath<bool>("//armadillo2", bresult);
    cout << "Query: //armadillo2 returned boolean: " << boolalpha << bresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // Try and get a bool from a tag that contains two bools instead of one
  // this will therefore fail.
  try { 
    reader.getXPath<bool>("//funkybool", bresult);
    cout << "Query: //funkybool returned boolean: " << boolalpha << bresult <<endl;
  }
  catch( const string& e) { 
    cout << e << endl;
  }

  // Try and get a complex with double data
  TComplex<double> rcomplex;
  
  try { 
    reader.getXPath("/foo/jim", rcomplex);
    cout << "Query: /foo/jim returned: Re: " << rcomplex.real() 
	 << " Im: " << rcomplex.imag()
	 << endl;
  }
  catch( const string& e) { 
    cerr << e << endl;
  }

  // Try and get a complex from /foo/bar.
  // /foo/bar is i) Not unique, ii) is not marked up as complex so this
  // should fail
  try {
       reader.getXPath("/foo/bar", rcomplex);
       cout << "Query: /foo/bar returned: Re: " << rcomplex.real()
            << " Im: " << rcomplex.imag()
            << endl;
  }
  catch( const string& e) {
         cerr << e << endl;
   }

  // Try and get a complex, made up of strings for the real and 
  // complex parts. Mathematically unsound but demonstrates
  // that we can read recursively...
  TComplex<string> cs;

  try {
    reader.getXPath("/foo/george", cs);
    cout << "Query: /foo/george returned: Re: " << cs.real()
	 << " Im: " << cs.imag()
	 << endl;
    cout.flush();
  }
  catch( const string& e) {
    cerr << e << endl;
  }

  // Try and get a complex made up of complexes
  // This demonstrates that our templating recursively works fine
  TComplex< TComplex<string> > ncs; 

  
  try {
    reader.getXPath("/foo/nested_george", ncs);
  }
  catch( const string& e) {
    cerr << e << endl;
  }
  
  cout << "Query: /foo/nested_george returned: Re: Re: " << ncs.real().real()
       << " Re: Im: " << ncs.real().imag() 
       << " Im: Re: " << ncs.imag().real() 
       << " Im: Im: " << ncs.imag().imag()
       << endl;
  cout.flush();
  
  // Try ang get some attributes
  string sattrib;
  int iattrib;
  
  // Get the index attribute from the first occurrance of /foo/bar
  try {
    reader.getXPathAttribute("/foo/bar[1]", "index", sattrib);
  }
  catch( const string& e) { 
    cout << e << endl;
  }
  cout << "String attribute index on node identified by" 
       << " /foo/bar[1]" << " has value " << sattrib << endl;

  // Get the intex attribute from the second occurrance of /foo/bar
  // this time it is an integer
  try {
    reader.getXPathAttribute("/foo/bar[2]", "index", iattrib);
  }
  catch( const string& e) { 
    cout << e << endl;
  }
  cout << "int attribute \"index\" on node identified by" 
       << " /foo/bar[2]" << " has value " << iattrib << endl;
   
  // Try ang get the attrite index from a query which 
  // finds multiple nodes... this should fail.
  try {
    reader.getXPathAttribute("/foo/bar", "index", iattrib);
  }
  catch( const string& e) { 
    cout << e << endl;
  } 
  
  // Count the number of tags /foo/bar returns 
  int occursFooBar;
  try { 
    occursFooBar = reader.countXPath("/foo/bar");
  }
  catch( const string& e) {
    cout << e << endl;
  }
  cout << "Query /foo/bar occurs " << occursFooBar << " times" << endl;

  // Read a string array
  Array<string> sa;
  try {
    reader.getXPath("/foo/arraytest_tag", sa);
  }
  catch(const string& e) {
    cout << e << endl;
  }

  int idx;
  for(idx = 0; idx < sa.size(); idx++) { 
    cout << "sa[" << idx << "] = " << sa[idx] << endl;
  }


  // Read an array of arrays whose elements are complexes templated on
  // strings -- this is the most difficult test. What is more the 
  // array elements are out of order
  Array< Array< TComplex< string > > > caa;
  
  try {
    reader.getXPath("//bj:rarraytest", caa);
  }
  catch(const string& e) { 
    cout << e << endl;
  }

  // Dump the array
  int alen1, alen2;
  int idx1, idx2;

  alen1 = caa.size();
  for(idx1 = 0; idx1 < alen1; idx1++) { 
    alen2 = caa[idx1].size();
    for(idx2 = 0; idx2 < alen2; idx2++) {
      cout << "caa[" << idx1<<"]["<<idx2<<"] = (" << caa[idx1][idx2].real()
	   <<" , " << caa[idx1][idx2].imag() << ")" << endl;
    }
  }

  XPathReader jimmi(reader, "//bj:pooh");
  try {
    //    jimmi.registerNamespace("bj", "http://www.ph.ed.ac.uk/~bj/");
    jimmi.getXPath("bj:bear", idx1);

    cout << "Reader jimmis query bj:bear returned " << idx1 << endl;
  }
  catch(const string& e) { 
    cout << e << endl;
  }
#endif

  // We're done.
  cout << "Closing reader." << endl;

  /* This bit kills the last refcount */
  reader.close();

  return EXIT_SUCCESS;
}

