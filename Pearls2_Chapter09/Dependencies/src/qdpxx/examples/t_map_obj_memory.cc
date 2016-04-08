#include "qdp.h"
#include <cstdlib>
#include <iostream>

#include "qdp_map_obj_memory.h"

namespace QDP
{
  //****************************************************************************
  //! Prop operator
  struct KeyPropColorVec_t
  {
    int        t_source;      /*!< Source time slice */
    int        colorvec_src;  /*!< Source colorvector index */
    int        spin_src;      /*!< Source spin index */
  };

  //----------------------------------------------------------------------------
  // Support for the keys of prop color vectors
  bool operator<(const KeyPropColorVec_t& a, const KeyPropColorVec_t& b)
  {
    multi1d<int> lgaa(3);
    lgaa[0] = a.t_source;
    lgaa[1] = a.colorvec_src;
    lgaa[2] = a.spin_src;

    multi1d<int> lgbb(3);
    lgbb[0] = b.t_source;
    lgbb[1] = b.colorvec_src;
    lgbb[2] = b.spin_src;

    return (lgaa < lgbb);
  }

  //----------------------------------------------------------------------------
  // KeyPropColorVec read
  void read(BinaryReader& bin, KeyPropColorVec_t& param)
  {
    read(bin, param.t_source);
    read(bin, param.colorvec_src);
    read(bin, param.spin_src);
  }

  // KeyPropColorVec write
  void write(BinaryWriter& bin, const KeyPropColorVec_t& param)
  {
    write(bin, param.t_source);
    write(bin, param.colorvec_src);
    write(bin, param.spin_src);
  }

}


using namespace QDP;

void fail(int line)
{
  QDPIO::cout << "FAIL: line= " << line << std::endl;
  QDP_finalize();
  exit(1);
}

void testMapObjInsertions(MapObjectMemory<char, float>& the_map)
{
 // Open the map for 'filling'
  QDPIO::cout << "Opening MapObjectMemory<char,float> for writing..."; 
  try { 
    QDPIO::cout << "OK" << std::endl;
  }
  catch(...) {
    fail(__LINE__);
  }
  
  QDPIO::cout << "Inserting (key,value): " << std::endl;
  // store a quadratic
  for(char i=0; i < 10; i++) { 
    float val = (float)i;
    val *= val;
    char key = 'a'+i; 

    try { 
      the_map.insert(key, val);
      QDPIO::cout << "  ( " << key <<", "<< val <<" )" << std::endl;
    }
    catch(...) { 
      fail(__LINE__);
    }
  }
  
  QDPIO::cout << "Closing MapObjectMemory<char,float> for writing..." ;
  try {
    QDPIO::cout << "... OK" << std::endl;
  }
  catch(...) {
    fail(__LINE__);
  }
}


void testMapObjLookups(MapObjectMemory<char, float>& the_map)
{
  /* Now reopen - random access */
  QDPIO::cout << "Opening MapObjectMemory<char,float> for reading..."; 
  try {  
    QDPIO::cout << "OK" << std::endl;
  }
  catch(...) { 
    fail(__LINE__);
  }

  QDPIO::cout << "Forward traversal test: " << std::endl;
  // Traverse in sequence
  
  QDPIO::cout << "Looking up key: " ;
  for(char i=0; i<10; i++) {
    char key = 'a'+i;
    float val;
    
    try{ 
      the_map.get(key, val);
      QDPIO::cout << " " << key;
    }
    catch(...) { 
      fail(__LINE__);
    }
    
    float refval = i*i;
    float diff = fabs(refval-val);
    if ( diff > 1.0e-8 ) {
      fail(__LINE__);
    }
    else { 
      QDPIO::cout << ".";
    }
  }

  QDPIO::cout << std::endl << "OK" << std::endl;


  QDPIO::cout << "Reverse traversal check" << std::endl;
  QDPIO::cout << "Looking up key: " ;
  // Traverse in reverse

  for(char i=9; i >= 0; i--) {
    char key='a'+i; 
    float val;
    try {
      the_map.get(key, val);
      QDPIO::cout << " " << key;
    }
    catch(...) { 
      fail(__LINE__);
    }

    float refval = i*i;
    float diff = fabs(refval-val);
    if ( diff > 1.0e-8 ) {
      fail(__LINE__);
    }
    else { 
      QDPIO::cout << "." << std::endl;
    }

  }
  QDPIO::cout << std::endl << "OK" << std::endl;


  // 20 'random' reads
  QDPIO::cout << "Random access... 20 reads" << std::endl;
  QDPIO::cout << "Looking up key: " ;
  for(int j=0; j < 20; j++){ 
    char key = 'a'+std::rand() % 10; // Pick randomly in 0..10
    float val;

    try { 
      the_map.get(key, val);
      QDPIO::cout << " " << key ;
    }
    catch(...) {
      fail(__LINE__);
    }

    float refval = (key-'a')*(key-'a');

    float diff = fabs(refval-val);
    if ( diff > 1.0e-8 ) {
      fail(__LINE__);
    }
    else { 
      QDPIO::cout <<".";
    }
  }
  QDPIO::cout << std::endl << "OK" << std::endl;
}



//**********************************************************************************************
void testMapKeyPropColorVecInsertions(MapObjectMemory<KeyPropColorVec_t, LatticeFermion>& pc_map, 
				      const multi1d<LatticeFermion>& lf_array)
{
  // Create the key-type
  KeyPropColorVec_t the_key = {0,0,0};

  // OpenMap for Writing
  QDPIO::cout << "Opening Map<KeyPropColorVec_t,LF> for writing..." << std::endl;
  try { 
    QDPIO::cout << "OK" << std::endl;
  }
  catch(...) {
    fail(__LINE__);
  }

  QDPIO::cout << "Inserting array element : ";
  for(int i=0; i < lf_array.size(); i++) { 
    the_key.colorvec_src = i;

    try { 
      pc_map.insert(the_key, lf_array[i]);
      QDPIO::cout << " "<< i << std::endl;
    }
    catch(...) {
      fail(__LINE__);
    }
  }

  QDPIO::cout << "Closing Map<KeyPropColorVec_t,LF> for writing..." << std::endl;
  try { 
    QDPIO::cout << "OK" << std::endl;
  }
  catch(...) { 
    fail(__LINE__);
  }
}


void testMapKeyPropColorVecLookups(MapObjectMemory<KeyPropColorVec_t, LatticeFermion>& pc_map, 
				   const multi1d<LatticeFermion>& lf_array)
{
  // Open map in read mode
  QDPIO::cout << "Opening Map<KeyPropColorVec_t,LF> for reading.." << std::endl;

  try { 
    QDPIO::cout << "OK" << std::endl;
  }
  catch(...) { 
    fail(__LINE__);
  }

  QDPIO::cout << "Increasing lookup test:" << std::endl;
  QDPIO::cout << "Looking up with colorvec_src = ";
  // Create the key-type
  KeyPropColorVec_t the_key = {0,0,0};

  for(int i=0; i < lf_array.size(); i++) {
    LatticeFermion lf_tmp;

    the_key.colorvec_src=i;
    try{
      pc_map.get(the_key, lf_tmp);
      QDPIO::cout << " " << i;
    }
    catch(...) { 
      fail(__LINE__);
    }

    // Compare with lf_array
    LatticeFermion diff;
    diff = lf_tmp - lf_array[i];
    Double diff_norm = sqrt(norm2(diff))/Double(Nc*Ns*Layout::vol());
    if(  toDouble(diff_norm) < 1.0e-6 )  { 
      QDPIO::cout << "." ;
    }
    else { 
      fail(__LINE__);
    }

  }
  QDPIO::cout << std::endl << "OK" << std::endl;

  QDPIO::cout << "Random access lookup test" << std::endl;
  QDPIO::cout << "Looking up with colorvec_src = " ;
  // Hey DJ! Spin that disk...
  for(int j=0; j < 100; j++) {
    int i = random() % lf_array.size();
    LatticeFermion lf_tmp;
    the_key.colorvec_src=i;
    try{
      pc_map.get(the_key, lf_tmp);
      QDPIO::cout << " " << i;
    }
    catch(...) {
      fail(__LINE__);
    }

    // Compare with lf_array
    LatticeFermion diff;
    diff = lf_tmp - lf_array[i];
    Double diff_norm = sqrt(norm2(diff))/Double(Nc*Ns*Layout::vol());
    if(  toDouble(diff_norm) < 1.0e-6 )  { 
      QDPIO::cout << ".";
    }
    else { 
      fail(__LINE__);
    }
  }
  QDPIO::cout << std::endl << "OK" << std::endl;
}



//**********************************************************************************************
int main(int argc, char *argv[])
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  // Setup the layout
  const int foo[] = {2,2,2,4};
  multi1d<int> nrow(Nd);
  nrow = foo;  // Use only Nd elements

  Layout::setLattSize(nrow);
  Layout::create();

  // Params to create a map object disk
  std::string map_obj_file("t_map_obj_disk.mod");

#if 1
  //
  // Test simple scalar
  //
  try {
    // Make a disk map object -- keys are ints, data floats
    MapObjectMemory<char,float> made_map;

    testMapObjInsertions(made_map);
    testMapObjLookups(made_map);
  }
  catch(const std::string& e) { 
    QDPIO::cout << "Caught: " << e << std::endl;
    fail(__LINE__);
  }
#endif


#if 1
  //
  // Test lattice objects
  //
  try {
    // Make an array of LF-s filled with noise
    multi1d<LatticeFermion> lf_array(10);
    for(int i=0; i < lf_array.size(); i++) { 
      gaussian(lf_array[i]);
    }
  
    MapObjectMemory<KeyPropColorVec_t, LatticeFermion> pc_map;
    
    testMapKeyPropColorVecInsertions(pc_map, lf_array);
    testMapKeyPropColorVecLookups(pc_map, lf_array);
    QDPIO::cout << std::endl << "OK" << std::endl;

    // Test an update 
    QDPIO::cout << "Doing update test ..." << std::endl;
    KeyPropColorVec_t the_key = {0,0,0};
    the_key.colorvec_src = 5;
    LatticeFermion f; gaussian(f);
    QDPIO::cout << "Updating..." ;
    pc_map.insert(the_key,f);
    QDPIO::cout << "OK" << std::endl;

    LatticeFermion f2;
    QDPIO::cout << "Re-Looking up...";
    pc_map.get(the_key,f2);
    QDPIO::cout << "OK" << std::endl;

    QDPIO::cout << "Comparing..." << std::endl;
    f2 -= f;
    if( toBool( sqrt(norm2(f2)) > toDouble(1.0e-6) ) ) {
      QDPIO::cout << "sqrt(norm2(f2))=" << sqrt(norm2(f2)) << std::endl;
      fail(__LINE__);
    }
    else { 
      QDPIO::cout << "OK" << std::endl ;
    }
    // Reinsert previous value
    pc_map.insert(the_key,lf_array[5]);
    testMapKeyPropColorVecLookups(pc_map, lf_array);
  }
  catch(const std::string& e) { 
    QDPIO::cout << "Caught: " << e << std::endl;
    fail(__LINE__);
  }
#endif


  QDP_finalize();
  return 0;
}

