// -*- C++ -*-
/*! @file
 * @brief Stubs of wrappers over filedb
 */

#ifndef QDP_DB_STUB_H
#define QDP_DB_STUB_H

#include <vector>
#include <string>
#include <exception>
#include <fcntl.h>

//--------------------------------------------------------------------------------
namespace FILEDB
{
  //! Forward empty decl
  class FFDB_DBT {};



  //-------------------------------------------------------------------------
  /**
   * Dummy exception class: should emulate SerializeException 
   * from filedb/src/Serializable.h
   */
  class SerializeException : public std::exception {};


  //-------------------------------------------------------------------------
  /** Dummy serializable class
   * should emulated class from filedb/src/Serializable.h
   */
  class Serializable
  {
  public:
    /**
     * Destructor
     */
    virtual ~Serializable (void) {;}

    /**
     * Get the serial id of this class
     */
    virtual const unsigned short serialID (void) const = 0;

    /**
     * Return this object into a binary form
     */
    virtual void writeObject (std::string& output) const throw (SerializeException) = 0;


    /**
     * Convert input object retrieved from database or network into an object
     */
    virtual void readObject (const std::string& input) throw (SerializeException) = 0;


  protected:
    /**
     * Constructor
     */
    Serializable (void) {;}
  };


  //---------------------------------------------------------------------------
  /**  Dummy DBKey Base class: emulate class from filedb/src/DBKey.h
   */

  class DBKey : public Serializable
  {
  public:
    /**
     * Destructor
     */
    ~DBKey (void) {;}

    /**
     * Does this key provide its own hash function
     * If this class is going to provide the hash function, use the
     * above hash function definition to implement a static function
     * with name hash
     *
     * @return 1 this class provide hash function, 0 otherwise
     */
    virtual int hasHashFunc (void) const = 0;

    /**
     * Does this key provide its own btree key compare function
     * If this class is going to provide the compare function, use the
     * above compare function definition to implement a static function
     * with name compare
     *
     * @return 1 this class provide compare function, 0 otherwise
     */
    virtual int hasCompareFunc (void) const = 0;

  protected:
    /**
     * Constructor
     */
    DBKey (void) {;}
  };


  //--------------------------------------------------------------------------
  /**  Dummy DBKey Base class: emulate class from filedb/src/DBData.h 
   */
  class DBData : public Serializable
  {
  public:
    /**
     * Destructor
     */
    ~DBData (void) {;}

  protected:
    /**
     * Constructor
     */
    DBData (void) {;}

  };
} // end namespace FILEDB



//--------------------------------------------------------------------------------
namespace QDP
{
  //! Fake cache and page sizes
  const int db_cachesize = 0;
  const int db_pagesize = 0;

  //! Forward empty decl
//  class Db;
//  class Dbt;


  /*! @defgroup io IO
   *
   * FileDB support
   *
   * @{
   */

  //--------------------------------------------------------------------------------
  //!  DB Base class
  /*!
    This class is used for writing of user data (most usefully measurements)
    into a DB file with a key/value semantics. 
  */
  template<typename K, typename D>
  class BinaryStoreDB
  {
  public:
    /**
     * Empty constructor for a DB
     */
    BinaryStoreDB() {notImplemented();}

    /*!
      Destroy the object
    */
    virtual ~BinaryStoreDB() {notImplemented();}

    /**
     * How much data and keys should be kept in memory in bytes
o     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSize (const unsigned int size) {notImplemented();}
    
    /**
     * How much data and keys should be kept in memory in megabytes
     *
     * This should be called before the open is called
     * @param max_cache_size number of bytes of data and keys should be kept
     * in memory
     */
    virtual void setCacheSizeMB (const unsigned int size) {notImplemented();}

    /**
     * Page size used when a new data based is created
     * This only effects a new database
     *
     * @param pagesize the pagesize used in hash database. This value
     * should be power of 2. The minimum value is 512 and maximum value
     * is 262144
     *
     */
    virtual void setPageSize (const unsigned int size) {notImplemented();}


    /**
     * Set initial number of buckets
     *
     * This should be called before the open is called
     *
     * @param num the number of buckets should be used
     *
     */
    virtual void setNumberBuckets (const unsigned int num) {notImplemented();}


    /**
     * Set whether to move pages when close to save disk space
     *
     * This only effective on writable database
     *
     */
    virtual void enablePageMove (void) {notImplemented();}

    virtual void disablePageMove (void) {notImplemented();}


    /**
     * Set and get maximum user information length
     */
    virtual void setMaxUserInfoLen (unsigned int len) {notImplemented();}

    virtual unsigned int getMaxUserInfoLen (void) const {notImplemented();}

    /**
     * Set and get maximum number of configurations
     */
    virtual void setMaxNumberConfigs (unsigned int num) {notImplemented();}

    virtual unsigned int getMaxNumberConfigs (void) const {notImplemented();}

    /**
     * Check if a DB file exists before opening.
     */
    virtual bool fileExists (const std::string& file) const {notImplemented();}


    /**
     * Open
     * @param file filename holding all data and keys
     * @param open_flags: can be regular UNIX file open flags such as: O_RDONLY,
     * O_RDWR, O_TRUNC
     * @param mode regular unix file mode
     *
     * @return 0 on success, -1 on failure with proper errno set
     * 
     */
    virtual void open (const std::string& file, int open_flags, int mode) {notImplemented();}

    virtual void close (void) {notImplemented();}

    /**
     * Does this key exist in the store
     * @param key a key object
     * @return true if the answer is yes
     */
    bool exist(const K& key) {notImplemented();}

    /**
     * Insert a pair of data and key into the database
     * data is not ensemble, but a vector of complex.
     * @param key a key
     * @param data a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    void insert (const K& key, const D& data) {notImplemented();}

    /**
     * Insert a pair of data and key into the database in string format
     * @param key a key
     * @param data a user provided data
     *
     * @return 0 on successful write, -1 on failure with proper errno set
     */
    void insertBinary (const std::string& key, const std::string& data) {notImplemented();}

    /**
     * Get data for a given key
     * @param key user supplied key
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int get (K& key, D& data) {notImplemented();}

    /**
     * Get data for a given key in binary form
     *
     * @param key user supplied key in string format
     * @param data after the call data will be populated
     * @return 0 on success, otherwise the key not found
     */
    int getBinary (std::string& key, std::string& data) {notImplemented();}

    /**
     * Return all available keys to user
     * @param keys user suppled an empty vector which is populated
     * by keys after this call.
     */
    virtual void keys (std::vector<K>& keys_) {notImplemented();}

    /**
     * Return all pairs of keys and data
     * @param keys user supplied empty vector to hold all keys
     * @param data user supplied empty vector to hold data
     * @return keys and data in the vectors having the same size
     */
    virtual void keysAndData (std::vector<K>& keys_, std::vector<D>& values_) {notImplemented();}

    /**
     * Return all pairs of keys and data in binary string form
     * @param keys user supplied empty vector to hold all keys
     * @param data user supplied empty vector to hold data
     * @return keys and data in the vectors having the same size     
     */
    virtual void binaryKeysAndData (std::vector<std::string> keys_,
				    std::vector<std::string> values_)
      {notImplemented();}

    /**
     * Flush database in memory to disk
     */
    virtual void flush (void) {notImplemented();}

    /**
     * Name of database associated with this Data store
     *
     * @return database name
     */
    virtual std::string storageName (void) const {notImplemented();}

    
    /**
     * Insert user data into the  metadata database
     *
     * @param user_data user supplied data
     * @return returns 0 if success, else failure
     */
    virtual void insertUserdata (const std::string& user_data) {notImplemented();}
    
    /**
     * Get user user data from the metadata database
     *
     * @param user_data user supplied buffer to store user data
     * @return returns 0 if success. Otherwise failure.
     */
    virtual void getUserdata (std::string& user_data) {notImplemented();}

  private:
    void notImplemented() const
      {
	QDPIO::cerr << "BinaryStoreDB: not implemented - this is a stub version. You must --enable-filedb in qdp++" << std::endl;
	QDP_abort(1);
      }
  };
}  // namespace QDP

#endif
