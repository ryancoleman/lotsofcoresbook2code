// $Id: t_db.cc,v 1.10 2009-08-28 16:06:51 edwards Exp $
/*! \file
 *  \brief Test the database routines
 */

#include "qdp.h"
#include "qdp_db.h"

using namespace QDP;

namespace Chroma
{
  using namespace FILEDB;

  //---------------------------------------------------------------------
  //! Some struct to use
  struct TestDBKey_t
  {
    multi1d<int>  mom;
    int           gamma;
    std::string   mass_label;
  };

  //! Reader
  void read(BinaryReader& bin, TestDBKey_t& param)
  {
    read(bin, param.mom);
    read(bin, param.gamma);
    read(bin, param.mass_label, 100);
  }

  //! Writer
  void write(BinaryWriter& bin, const TestDBKey_t& param)
  {
    write(bin, param.mom);
    write(bin, param.gamma);
    write(bin, param.mass_label);
  }


  //---------------------------------------------------------------------
  //! Some struct to use
  struct TestDBData_t
  {
    multi1d<ComplexD>  op;
    int                type_of_data;
  };

  //! Reader
  void read(BinaryReader& bin, TestDBData_t& param)
  {
    read(bin, param.op);
    read(bin, param.type_of_data);
  }

  //! Writer
  void write(BinaryWriter& bin, const TestDBData_t& param)
  {
    write(bin, param.op);
    write(bin, param.type_of_data);
  }


  //---------------------------------------------------------------------
  // Simple concrete key class
  template<typename K>
  class TestDBKey : public DBKey
  {
  public:
    //! Default constructor
    TestDBKey() {} 

    //! Constructor from data
    TestDBKey(const K& k) : key_(k) {}

    //! Setter
    K& key() {return key_;}

    //! Getter
    const K& key() const {return key_;}

    // Part of Serializable
    const unsigned short serialID (void) const {return 456;}

    void writeObject (std::string& output) const throw (SerializeException) {
      BinaryBufferWriter bin;
      write(bin, key());
      output = bin.str();
    }

    void readObject (const std::string& input) throw (SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, key());
    }

    // Part of DBKey
    int hasHashFunc (void) const {return 0;}
    int hasCompareFunc (void) const {return 0;}

    /**
     * Empty hash and compare functions. We are using default functions.
     */
    static unsigned int hash (const void* bytes, unsigned int len) {return 0;}
    static int compare (const FFDB_DBT* k1, const FFDB_DBT* k2) {return 0;}
   
  private:
    K  key_;
  };


  //---------------------------------------------------------------------
  // Simple concrete data class
  template<typename D>
  class TestDBData : public DBData
  {
  public:
    //! Default constructor
    TestDBData() {} 

    //! Constructor from data
    TestDBData(const D& d) : data_(d) {}

    //! Setter
    D& data() {return data_;}

    //! Getter
    const D& data() const {return data_;}

    // Part of Serializable
    const unsigned short serialID (void) const {return 123;}

    void writeObject (std::string& output) const throw (SerializeException) {
      BinaryBufferWriter bin;
      write(bin, data());
      output = bin.str();
    }

    void readObject (const std::string& input) throw (SerializeException) {
      BinaryBufferReader bin(input);
      read(bin, data());
    }

  private:
    D  data_;
  };

} // namespace Chroma


using namespace Chroma;

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

  XMLFileWriter xml("t_db.xml");
  push(xml,"t_db");

//  push(xml,"lattice");
  write(xml,"Nd", Nd);
  write(xml,"Nc", Nc);
  write(xml,"Ns", Ns);
  write(xml,"nrow", nrow);
//  pop(xml);

  // Try out some simple DB stuff
  typedef BinaryStoreDB< TestDBKey<TestDBKey_t>, TestDBData<TestDBData_t> > DBType_t;

  try
  {
    TestDBKey<TestDBKey_t> testDBKey;
    {
      testDBKey.key().mom.resize(3);
      testDBKey.key().mom = 0;
      testDBKey.key().gamma = 15;
      testDBKey.key().mass_label = "-0.0840";
    }

    TestDBData<TestDBData_t> testDBData;
    {
      testDBData.data().op.resize(QDP::Layout::lattSize()[3]);
      testDBData.data().op = zero;
      testDBData.data().type_of_data = 2;
    }
    
    // Create it
    QDPIO::cout << "create" << std::endl;
    DBType_t db;

    // Create some meta data. Need to know length before the open call
    std::string meta_data("my name is fred");
    db.setMaxUserInfoLen(meta_data.size());

    // Open the db
    const std::string dbase("test.db");
    
    if (db.fileExists(dbase))
    {
      QDPIO::cout << "File " << dbase << " already exists - will truncate\n";
    }
    else
    {
      QDPIO::cout << "File " << dbase << " does not exist - will create\n";
    }

    db.open(dbase, O_RDWR | O_CREAT | O_TRUNC, 0664);

    // Insert meta data
    db.insertUserdata(meta_data);

    // Test it
    QDPIO::cout << "insert" << std::endl;
    db.insert(testDBKey, testDBData);

    // Flush it
    QDPIO::cout << "flush" << std::endl;
    db.flush();

    QDPIO::cout << "closing" << std::endl;
  }
  catch(std::exception &e) {
    QDPIO::cerr << "Std exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch (...) {
    QDPIO::cout << "Generic error in db routines" << std::endl;
    QDP_abort(1);
  }

  // See if we can lookup in the DB
  try
  {
    TestDBKey<TestDBKey_t> testDBKey;
    {
      testDBKey.key().mom.resize(3);
      testDBKey.key().mom = 0;
      testDBKey.key().gamma = 15;
      testDBKey.key().mass_label = "-0.0840";
    }

    TestDBData<TestDBData_t> testDBData;

    // Open it
    QDPIO::cout << "open" << std::endl;
    DBType_t db;
    db.open("test.db", O_RDONLY, 0644);

    // Test it
    QDPIO::cout << "get" << std::endl;
    if (db.get(testDBKey, testDBData) != 0)
    {
      QDPIO::cout << "Error: attempting a get, the key was not found" << std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "type_of_data= " << testDBData.data().type_of_data << std::endl;

    // Grab the user meta data
    std::string meta_data;
    db.getUserdata(meta_data);
    QDPIO::cout << "User meta_data = XX" << meta_data << "XX" << std::endl;

    QDPIO::cout << "closing" << std::endl;
  }
  catch(std::exception &e) {
    QDPIO::cerr << "Std exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch (...) {
    QDPIO::cout << "Generic error in db routines" << std::endl;
    QDP_abort(1);
  }


  // Find all the keys from the DB
  try
  {
    // Open it
    QDPIO::cout << "open" << std::endl;
    DBType_t db;
    db.open("test.db", O_RDONLY, 0644);

    // Test it
    QDPIO::cout << "find keys" << std::endl;
    std::vector< TestDBKey<TestDBKey_t> > keys;
    db.keys(keys);

    // Flush it
    QDPIO::cout << "printing" << std::endl;
    for(int i=0; i < keys.size(); ++i)
    {
      QDPIO::cout << "found a key: i= " << i 
		  << "  gamma= " << keys[i].key().gamma
		  << std::endl;
    }

    QDPIO::cout << "closing" << std::endl;
  }
  catch(std::exception &e) {
    QDPIO::cerr << "Std exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch (...) {
    QDPIO::cout << "Generic error in db routines" << std::endl;
    QDP_abort(1);
  }

  // Done
  pop(xml);

  // Time to bolt
  QDP_finalize();

  QDPIO::cout << "finished" << std::endl;
  exit(0);
}

