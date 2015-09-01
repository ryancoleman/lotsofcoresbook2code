// -*- C++ -*-

/*! @file
 * @brief IO support
 */

#ifndef QDP_IO_H
#define QDP_IO_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "qdp_byteorder.h"
#include <vector>
#include <map>
#include <list>

#if 1
#include <complex>
#endif

namespace QDP
{
  /*! @defgroup io IO
   *
   * File input and output operations on QDP types
   *
   * @{
   */

  //--------------------------------------------------------------------------------
  //! Text input class
  /*!
    This class is used to read data from a text object. Input is done on the
    primary node and all nodes end up with the same data.

    The read methods are also wrapped by externally defined >> operators,
  */

  class TextReader
  {
  public:
    TextReader();

    /*!
      Closes the last file opened
    */
    virtual ~TextReader() {}

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in the previous IO operation
    */
    virtual bool fail();

    // Readers for builtin types
    virtual void read(std::string& result);
    virtual void read(char& result);
    virtual void read(int& result);
    virtual void read(unsigned int& result);
    virtual void read(short int& result);
    virtual void read(unsigned short int& result);
    virtual void read(long int& result);
    virtual void read(unsigned long int& result);
    virtual void read(long long int& result);
    virtual void read(float& result);
    virtual void read(double& result);
    virtual void read(bool& result);

  protected:
    //! The universal data-reader.
    /*!
      All the read functions call this.
      \param output The location to which the datum is read.
    */
    template< typename T>
    void
    readPrimitive(T& output);

    //! Get the internal input stream
    virtual std::istream& getIstream() = 0;
  };


  // Different bindings for same operators
  TextReader& operator>>(TextReader& txt, std::string& input);
  TextReader& operator>>(TextReader& txt, char& input);
  TextReader& operator>>(TextReader& txt, int& input);
  TextReader& operator>>(TextReader& txt, unsigned int& input);
  TextReader& operator>>(TextReader& txt, short int& input);
  TextReader& operator>>(TextReader& txt, unsigned short int& input);
  TextReader& operator>>(TextReader& txt, long int& input);
  TextReader& operator>>(TextReader& txt, unsigned long int& input);
  TextReader& operator>>(TextReader& txt, long long int& input);
  TextReader& operator>>(TextReader& txt, float& input);
  TextReader& operator>>(TextReader& txt, double& input);
  TextReader& operator>>(TextReader& txt, bool& input);


  //--------------------------------------------------------------------------------
  //! Text buffer input class
  /*!
    This class is used to read data from a text buffer. Input is done on the
    primary node and all nodes end up with the same data.

    The read methods are also wrapped by externally defined >> operators,
  */

  class TextBufferReader : public TextReader
  {
  public:
    TextBufferReader();

    //! Shutdown the buffer
    ~TextBufferReader(); 

    /*!
      Initialize a buffer for reading.
      \param p The string of input
    */
    explicit TextBufferReader(const std::string& p);

    //! Opens a buffer for reading.
    /*!
      \param p The string of input
    */
    void open(const std::string& p);

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Return entire buffer as a string
    std::string strPrimaryNode() const;
        
  protected:
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

  private:
    std::istringstream f;
  };


  //--------------------------------------------------------------------------------
  //! Text file input class
  /*!
    This class is used to read data from a text file. Input is done on the
    primary node and all nodes end up with the same data.

    The read methods are also wrapped by externally defined >> operators,
  */

  class TextFileReader : public TextReader
  {
  public:
    TextFileReader();

    /*!
      Closes the last file opened
    */
    ~TextFileReader(); 

    /*!
      Opens a file for reading.
      \param p The name of the file
    */
    explicit TextFileReader(const std::string& p);

    //! Opens a file for reading.
    /*!
      \param p The name of the file
    */
    void open(const std::string& p);

    //! Closes the last file opened
    void close();

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();

  protected:
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

  private:
    std::ifstream f;
  };


  //--------------------------------------------------------------------------------
  //! Text output base class
  /*!
    This class is used to write data to a text object.
    Output is done from the primary node only..
  
    The write methods are also wrapped by externally defined >> operators,
  */
  class TextWriter
  {
  public:
    TextWriter();

    /*!
      Closes the last file opened
    */
    virtual ~TextWriter() {}

    //! Flushes the object
    virtual void flush() = 0;

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in previous IO operation
    */
    virtual bool fail();
    
    // Overloaded Writer Functions
    virtual void write(const std::string& output);
    virtual void write(const char* output);
    virtual void write(const char& output);
    virtual void write(const int& output);
    virtual void write(const unsigned int& output);
    virtual void write(const short int& output);
    virtual void write(const unsigned short int& output);
    virtual void write(const long int& output);
    virtual void write(const unsigned long int& output);
    virtual void write(const long long int& output);
    virtual void write(const float& output);
    virtual void write(const double& output);
    virtual void write(const bool& output);

  protected:
    //! The universal data-writer.
    /*!
      All the write functions call this.
      \param output The location of the datum to be written.
    */
    template< typename T>
    void
    writePrimitive(const T& output);

    //! Get the internal output stream
    virtual std::ostream& getOstream() = 0;
  };


  // Different bindings for same operators
  TextWriter& operator<<(TextWriter& txt, const std::string& output);
  TextWriter& operator<<(TextWriter& txt, const char* output);
  TextWriter& operator<<(TextWriter& txt, char output);
  TextWriter& operator<<(TextWriter& txt, int output);
  TextWriter& operator<<(TextWriter& txt, unsigned int output);
  TextWriter& operator<<(TextWriter& txt, short int output);
  TextWriter& operator<<(TextWriter& txt, unsigned short int output);
  TextWriter& operator<<(TextWriter& txt, long int output);
  TextWriter& operator<<(TextWriter& txt, unsigned long int output);
  TextWriter& operator<<(TextWriter& txt, long long int output);
  TextWriter& operator<<(TextWriter& txt, float output);
  TextWriter& operator<<(TextWriter& txt, double output);
  TextWriter& operator<<(TextWriter& txt, bool output);


  //--------------------------------------------------------------------------------
  //! Text buffer output class
  /*!
    This class is used to write data to a text buffer.
    Output is done from the primary node only..
  
    The write methods are also wrapped by externally defined >> operators,
  */
  class TextBufferWriter : public TextWriter
  {
  public:
    TextBufferWriter();

    //! Shutdown the buffer
    ~TextBufferWriter();

    /*!
      Initialize a buffer
      \param p The string for initialization
    */
    explicit TextBufferWriter(const std::string& p);

    //! Construct from a string
    void open(const std::string& s);

    //! Flushes the buffer
    void flush() {}

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Return entire buffer as a string
    std::string strPrimaryNode() const;
        
  protected:
    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

  private:
    std::ostringstream f;
  };


  //--------------------------------------------------------------------------------
  //! Text output class
  /*!
    This class is used to write data to a text file.
    Output is done from the primary node only..
  
    The write methods are also wrapped by externally defined >> operators,
  */
  class TextFileWriter : public TextWriter
  {
  public:
    TextFileWriter();

    /*!
      Closes the last file opened
    */
    ~TextFileWriter();

    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    explicit TextFileWriter(const std::string& p);

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();

    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    void open(const std::string& p);

    //! Closes the last file opened
    void close();

    //! Flushes the buffer
    void flush();

  protected:
    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

  private:
    std::ofstream f;
  };


  //--------------------------------------------------------------------------------
  //!  Binary input base class
  /*!
    This class is used to read data from a binary object. The data in the file
    is assumed to be big-endian. If the host nachine is little-endian, the data
    is byte-swapped. All nodes end up with the same data
  
    The read methods are also wrapped by externally defined functions
    and >> operators,   
  */
  class BinaryReader
  {
  public:
    typedef std::istream::pos_type pos_type;  // position in buffer
    typedef std::istream::off_type off_type;  // offset in buffer

  public:
    /*!
      Shutdown the object
    */
    virtual ~BinaryReader() {}

    //! Read data on the primary node only
    /*!
      \param output The location to which data is read
      \param nbytes The size in bytes of each datum
      \param The number of data.
    */
    virtual void readArrayPrimaryNode(char* output, size_t nbytes, size_t nmemb);

    //! Read data on the primary node and broadcast to all nodes.
    /*!
      \param output The location to which data is read
      \param nbytes The size in bytes of each datum
      \param The number of data.
    */
    virtual void readArray(char* output, size_t nbytes, size_t nmemb);

    // Overloaded reader functions
    virtual void readDesc(std::string& result);

    //! Read some max number of characters - 1 upto and excluding a newline
    /*! This is the getline function for the underlying stream */
    virtual void read(std::string& result, size_t nbytes);

    virtual void read(char& result);
    virtual void read(int& result);
    virtual void read(unsigned int& result);
    virtual void read(short int& result);
    virtual void read(unsigned short int& result);
    virtual void read(long int& result);
    virtual void read(unsigned long int& result);
    virtual void read(long long int& result);
    virtual void read(float& result);
    virtual void read(double& result);
    virtual void read(bool& result);

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in the previous IO operation
    */
    virtual bool fail();

    //! Get the current checksum
    virtual QDPUtil::n_uint32_t getChecksum();
  
    //! Reset the current checksum
    virtual void resetChecksum();
  
    //! Current position
    virtual pos_type currentPosition();

    //! Set the current position
    /*! The checksum is reset */
    virtual void seek(pos_type off);

    //! Set the position relative from the start
    /*! The checksum is reset */
    virtual void seekBegin(off_type off);

    //! Set the position relative to the current position
    /*! The checksum is reset */
    virtual void seekRelative(off_type off);

    //! Set the position relative from the end
    /*! The checksum is reset */
    virtual void seekEnd(off_type off);

    //! Rewind object to the beginning
    /*! The checksum is reset */
    virtual void rewind();

  protected:
    //! The universal data-reader.
    /*!
      All the read functions call this.
      \param output The location to which the datum is read.
    */
    template< typename T>
    void
    readPrimitive(T& output);

    //! Get the current checksum to modify
    virtual QDPUtil::n_uint32_t& internalChecksum() = 0;
  
    //! Get the internal input stream
    virtual std::istream& getIstream() = 0;
  };


  // Telephone book of basic primitives
  void readDesc(BinaryReader& bin, std::string& input);
  void read(BinaryReader& bin, std::string& input, size_t maxBytes);
  void read(BinaryReader& bin, char& input);
  void read(BinaryReader& bin, int& input);
  void read(BinaryReader& bin, unsigned int& input);
  void read(BinaryReader& bin, short int& input);
  void read(BinaryReader& bin, unsigned short int& input);
  void read(BinaryReader& bin, long int& input);
  void read(BinaryReader& bin, unsigned long int& input);
  void read(BinaryReader& bin, long long int& input);
  void read(BinaryReader& bin, float& input);
  void read(BinaryReader& bin, double& input);
  void read(BinaryReader& bin, bool& input);

  // Different bindings for same operators
  BinaryReader& operator>>(BinaryReader& bin, char& input);
  BinaryReader& operator>>(BinaryReader& bin, int& input);
  BinaryReader& operator>>(BinaryReader& bin, unsigned int& input);
  BinaryReader& operator>>(BinaryReader& bin, short int& input);
  BinaryReader& operator>>(BinaryReader& bin, unsigned short int& input);
  BinaryReader& operator>>(BinaryReader& bin, long int& input);
  BinaryReader& operator>>(BinaryReader& bin, unsigned long int& input);
  BinaryReader& operator>>(BinaryReader& bin, long long int& input);
  BinaryReader& operator>>(BinaryReader& bin, float& input);
  BinaryReader& operator>>(BinaryReader& bin, double& input);
  BinaryReader& operator>>(BinaryReader& bin, bool& input);

#if 1
  //! Complex reader
  void read(BinaryReader& bin, std::complex<float>& param);
  void read(BinaryReader& bin, std::complex<double>& param);
#endif


  //! Read a binary multi1d object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The multi1d can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi1d<T>& d)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    d.resize(n);

    for(int i=0; i < d.size(); ++i)
      read(bin, d[i]);
  }


  //! Read a binary multi1d object
  /*!
    This assumes that the number of elements to be read is not written in
    the file, \e i.e. that the data was written with the corresponding write
    code. The number of elements must therefore be supplied by the caller
    \param bin The initialised binary reader
    \param d The data to be filled.
    \param num The number of elements.

    \pre The binary reader must have opened the file.
    \pre The multi1d must have space for at least \a num elements.  
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi1d<T>& d, int num)
  {
    for(int i=0; i < num; ++i)
      read(bin, d[i]);
  }

  //! Read a binary multi2d object
  /*!
    This assumes that the number of elements to be read is not written in
    the file, \e i.e. that the data was written with the corresponding write
    code. The number of elements must therefore be supplied by the caller
    \param bin The initialised binary reader
    \param d The data to be filled.
    \param num1 The first dimension of the array
    \param num2 The second dimension of the array..  

    \pre The binary reader must have opened the file.
    \pre The multi2d must have space for at least \a num elements.  
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi2d<T>& d, int num1, int num2)
  {
    for(int i=0; i < num2; ++i)
      for(int j=0; j < num1; ++j)
	read(bin, d[j][i]);

  }


  //! Read a binary multi2d element
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The multi2d can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi2d<T>& d)
  {
    int n1;
    int n2;
    read(bin, n1);    // the size is always written, even if 0
    read(bin, n2);    // the size is always written, even if 0
    d.resize(n1,n2);
  
    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
      {
	read(bin, d[j][i]);
      }
  }


  //! Read a binary multi3d element
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The multi2d can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi3d<T>& d)
  {
    int n1;
    int n2;
    int n3;
    read(bin, n3);    // the size is always written, even if 0
    read(bin, n2);    // the size is always written, even if 0
    read(bin, n1);    // the size is always written, even if 0

    // Destructively resize the array
    d.resize(n3,n2,n1);

    for(int i=0; i < d.size1(); ++i)
    {
      for(int j=0; j < d.size2(); ++j)
      {
	for(int k=0; k < d.size3(); ++k)
	{
	  read(bin, d[k][j][i]);
	}
      }
    }
  }


  //! Read a binary multi4d element
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The multi4d will be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi4d<T>& d)
  {
    int n1;
    int n2;
    int n3;
    int n4;
    read(bin, n4);    // the size is always written, even if 0
    read(bin, n3);    // the size is always written, even if 0
    read(bin, n2);    // the size is always written, even if 0
    read(bin, n1);    // the size is always written, even if 0

    // Destructively resize the array
    d.resize(n4,n3,n2,n1);

    for(int i=0; i < d.size1(); ++i)
    {
      for(int j=0; j < d.size2(); ++j)
      {
	for(int k=0; k < d.size3(); ++k)
	{
	  for(int l=0; l < d.size4(); ++l)
	  {
	    read(bin, d[l][k][j][i]);
	  }
	}
      }
    }
  }


  //! Read a binary multi5d element
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The multi5d will be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, multi5d<T>& d)
  {
    int n1;
    int n2;
    int n3;
    int n4;
    int n5;
    read(bin, n5);    // the size is always written, even if 0
    read(bin, n4);    // the size is always written, even if 0
    read(bin, n3);    // the size is always written, even if 0
    read(bin, n2);    // the size is always written, even if 0
    read(bin, n1);    // the size is always written, even if 0

    // Destructively resize the array
    d.resize(n5,n4,n3,n2,n1);

    for(int i=0; i < d.size1(); ++i)
    {
      for(int j=0; j < d.size2(); ++j)
      {
	for(int k=0; k < d.size3(); ++k)
	{
	  for(int l=0; l < d.size4(); ++l)
	  {
	    for(int m=0; m < d.size5(); ++m)
	    {
	      read(bin, d[m][l][k][j][i]);
	    }
	  }
	}
      }
    }
  }


  //! Read a binary multiNd element
  template<class T>
  inline
  void read(BinaryReader& bin, multiNd<T>& d)
  {
    multi1d<int> siz;
    read(bin, siz); // read the array of the sizes

    d.resize(siz);

    for(int i=0; i < d.numElem(); ++i)
      read(bin, d.getElem(i));
  }


  //! Read a binary std::vector object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The std::vector can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, std::vector<T>& d)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    d.resize(n);

    for(int i=0; i < d.size(); ++i)
      read(bin, d[i]);
  }


  //! Read a binary std::list object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The std::list can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, std::list<T>& d)
  {
    d.clear();

    int n;
    read(bin, n);    // the size is always written, even if 0

    for(int i=0; i < n; ++i)
    {
      T thingy;
      read(bin, thingy);

      d.push_back(thingy);
    }
  }


  //! Read a binary Array1dO object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The ADAT::Array1dO can be resized.
  */
  template<class T>
  inline
  void read(BinaryReader& bin, Array1dO<T>& d)
  {
    int n;
    read(bin, n);    // the size is always written, even if 0
    d.resize(n);

    for(int i=1; i <= d.size(); ++i)
      read(bin, d[i]);
  }


  //! Read a binary std::map object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The std::vector can be resized.
  */
  template<typename K, typename V>
  inline
  void read(BinaryReader& bin, std::map<K,V>& d)
  {
    d.clear();

    int n;
    read(bin, n);    // the size is always written, even if 0

    for(int i=0; i < n; ++i)
    {
      K key;
      read(bin, key);

      V val;
      read(bin, val);

      d.insert(std::make_pair(key,val));
    }
  }


  //! Read a binary std::pair object
  /*!
    This assumes that the number of elements to be read is also written in
    the file, \e i.e. that the data was written with the corresponding write
    code.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
    \post The std::vector can be resized.
  */
  template<typename T1, typename T2>
  inline
  void read(BinaryReader& bin, std::pair<T1,T2>& d)
  {
    T1 f;
    read(bin, f);

    T2 s;
    read(bin, s);

    d = std::make_pair(f,s);
  }


  //--------------------------------------------------------------------------------
  //!  Binary buffer input class
  /*!
    This class is used to read data from a binary buffer. The data in the file
    is assumed to be big-endian. If the host nachine is little-endian, the data
    is byte-swapped. All nodes end up with the same data
  
    The read methods are also wrapped by externally defined functions
    and >> operators,   
  */
  class BinaryBufferReader : public BinaryReader
  {
  public:
    BinaryBufferReader();

    //! Construct from a string
    explicit BinaryBufferReader(const std::string& s);

    //! Closes the buffer
    ~BinaryBufferReader();

    //! Construct from a string
    void open(const std::string& s);

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Return entire buffer as a string
    std::string strPrimaryNode() const;
        
    //! Clear the buffer
    void clear();

  protected:
    //! Get the current checksum to modify
    QDPUtil::n_uint32_t& internalChecksum() {return checksum;}
  
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

  private:
    //! Checksum
    QDPUtil::n_uint32_t checksum;
    std::istringstream f;
  };


  //--------------------------------------------------------------------------------
  //!  Binary file input class
  /*!
    This class is used to read data from a binary file. The data in the file
    is assumed to be big-endian. If the host machine is little-endian, the data
    is byte-swapped. All nodes end up with the same data
  
    The read methods are also wrapped by externally defined functions
    and >> operators,   
  */
  class BinaryFileReader : public BinaryReader
  {
  public:
    BinaryFileReader();

    /*!
      Closes the last file opened
    */
    ~BinaryFileReader();

    /*!
      Opens a file for reading.
      \param p The name of the file
    */
    explicit BinaryFileReader(const std::string& p);

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();

    //! Opens a file for reading.
    /*!
      \param p The name of the file
    */
    void open(const std::string& p);

    //! Closes the last file opened
    void close();

  protected:
    //! Get the current checksum to modify
    QDPUtil::n_uint32_t& internalChecksum() {return checksum;}
  
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

  private:
    //! Checksum
    QDPUtil::n_uint32_t checksum;
    std::ifstream f;
  };


  //--------------------------------------------------------------------------------
  //!  Binary writer base class
  /*!
    This class is used to write data to a binary object. The data in the object
    is big-endian. If the host machine is little-endian, the data
    is byte-swapped.   Output is done from the primary node only.

    Writer object need to be opened before any of the write methods are used.
    In the case of files, the files need to be opened first
  
    The write methods are also wrapped by externally defined functions
    and << operators,   
  */
  class BinaryWriter
  {
  public:
    typedef std::ostream::pos_type pos_type;  // position in buffer
    typedef std::ostream::off_type off_type;  // offset in buffer

  public:
    //! Default constructor
    BinaryWriter();

    /*!
      Closes the last file opened
    */
    virtual ~BinaryWriter() {}

    //! Flushes the buffer
    virtual void flush() = 0;

    //! Write data on the primary node only
    /*!
      \param output The location to which data is read
      \param nbytes The size in bytes of each datum
      \param The number of data.
    */
    virtual void writeArrayPrimaryNode(const char* output, size_t nbytes, size_t nmemb);

    //! Write data from the primary node.
    /*!
      \param output The data to write
      \param nbytes The size in bytes of each datum
      \param The number of data.
    */
    virtual void writeArray(const char* output, size_t nbytes, size_t nmemb);

    // Overloaded Writer Functions

    //! Writes a fixed length of characters like an array
    virtual void writeDesc(const std::string& output);

    /*!
      A newline is appended to the written string.
    */
    virtual void write(const std::string& output);
    /*!
      A newline is appended to the written string.
    */
    virtual void write(const char* output);
    virtual void write(const char& output);
    virtual void write(const int& output);
    virtual void write(const unsigned int& output);
    virtual void write(const short int& output);
    virtual void write(const unsigned short int& output);
    virtual void write(const long int& output);
    virtual void write(const unsigned long int& output);
    virtual void write(const long long int& output);
    virtual void write(const float& output);
    virtual void write(const double& output);
    virtual void write(const bool& output);

    //! Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in previous IO operation
    */
    virtual bool fail();

    //! Get the current checksum
    virtual QDPUtil::n_uint32_t getChecksum();
  
    //! Reset the current checksum
    virtual void resetChecksum();
  
    //! Current position
    virtual pos_type currentPosition();

    //! Set the current position
    /*! The checksum is reset */
    virtual void seek(pos_type off);

    //! Set the position relative from the start
    /*! The checksum is reset */
    virtual void seekBegin(off_type off);

    //! Set the position relative to the current position
    /*! The checksum is reset */
    virtual void seekRelative(off_type off);

    //! Set the position relative from the end
    /*! The checksum is reset */
    virtual void seekEnd(off_type off);

    //! Rewind object to the beginning
    /*! The checksum is reset */
    virtual void rewind();

  protected:

    //! The universal data-writer.
    /*!
      All the write functions call this.
      \param output The location of the datum to be written.
    */
    template< typename T>
    void
    writePrimitive(const T& output);

    //! Get the current checksum to modify
    virtual QDPUtil::n_uint32_t& internalChecksum() = 0;
  
    //! Get the internal output stream
    virtual std::ostream& getOstream() = 0;
  };


  // Telephone book of basic primitives
  void writeDesc(BinaryWriter& bin, const std::string& output);
  void write(BinaryWriter& bin, const std::string& output);
  void write(BinaryWriter& bin, const char* output);
  void write(BinaryWriter& bin, char output);
  void write(BinaryWriter& bin, int output);
  void write(BinaryWriter& bin, unsigned int output);
  void write(BinaryWriter& bin, short int output);
  void write(BinaryWriter& bin, unsigned short int output);
  void write(BinaryWriter& bin, long int output);
  void write(BinaryWriter& bin, unsigned long int output);
  void write(BinaryWriter& bin, long long int output);
  void write(BinaryWriter& bin, float output);
  void write(BinaryWriter& bin, double output);
  void write(BinaryWriter& bin, bool output);

  // Different bindings for same operators
  BinaryWriter& operator<<(BinaryWriter& bin, const std::string& output);
  BinaryWriter& operator<<(BinaryWriter& bin, const char* output);
  BinaryWriter& operator<<(BinaryWriter& bin, char output);
  BinaryWriter& operator<<(BinaryWriter& bin, int output);
  BinaryWriter& operator<<(BinaryWriter& bin, unsigned int output);
  BinaryWriter& operator<<(BinaryWriter& bin, short int output);
  BinaryWriter& operator<<(BinaryWriter& bin, unsigned short int output);
  BinaryWriter& operator<<(BinaryWriter& bin, long int output);
  BinaryWriter& operator<<(BinaryWriter& bin, unsigned long int output);
  BinaryWriter& operator<<(BinaryWriter& bin, long long int output);
  BinaryWriter& operator<<(BinaryWriter& bin, float output);
  BinaryWriter& operator<<(BinaryWriter& bin, double output);
  BinaryWriter& operator<<(BinaryWriter& bin, bool output);

#if 1
  //! Complex writer
  void write(BinaryWriter& bin, const std::complex<float>& param);
  void write(BinaryWriter& bin, const std::complex<double>& param);
#endif

  //! Write all of a binary multi1d object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const multi1d<T>& d)
  {
    write(bin, d.size());    // always write the size
    for(int i=0; i < d.size(); ++i)
      write(bin, d[i]);
  }

  //! Write some or all of a binary multi1d object
  /*!
    This does not write the number of elements to the file.
    \param bin The initialised binary writer
    \param d The data to be filled.
    \param num The number of elements to write.

    \pre The binary writer must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const multi1d<T>& d, int num)
  {
    for(int i=0; i < num; ++i)
      write(bin, d[i]);
  }


  //! Write a binary multi2d element
  template<class T>
  inline
  void write(BinaryWriter& bin, const multi2d<T>& d)
  {
    write(bin, d.size2());    // always write the size
    write(bin, d.size1());    // always write the size

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
      {
	write(bin, d[j][i]);
      }

  }



  //! Write a binary multi3d element
  template<class T>
  inline
  void write(BinaryWriter& bin, const multi3d<T>& d)
  {
    write(bin, d.size3());    // always write the size
    write(bin, d.size2());    // always write the size
    write(bin, d.size1());    // always write the size

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
	for(int k=0; k < d.size3(); ++k)
	  write(bin, d[k][j][i]);
   }


  //! Write a binary multi4d element
  template<class T>
  inline
  void write(BinaryWriter& bin, const multi4d<T>& d)
  {
    write(bin, d.size4());    // always write the size
    write(bin, d.size3());    // always write the size
    write(bin, d.size2());    // always write the size
    write(bin, d.size1());    // always write the size

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
	for(int k=0; k < d.size3(); ++k)
	  for(int l=0; l < d.size4(); ++l)
	    write(bin, d[l][k][j][i]);
   }


  //! Write a binary multi5d element
  template<class T>
  inline
  void write(BinaryWriter& bin, const multi5d<T>& d)
  {
    write(bin, d.size5());    // always write the size
    write(bin, d.size4());    // always write the size
    write(bin, d.size3());    // always write the size
    write(bin, d.size2());    // always write the size
    write(bin, d.size1());    // always write the size

    for(int i=0; i < d.size1(); ++i)
      for(int j=0; j < d.size2(); ++j)
	for(int k=0; k < d.size3(); ++k)
	  for(int l=0; l < d.size4(); ++l)
	    for(int m=0; m < d.size5(); ++m)
	      write(bin, d[m][l][k][j][i]);
   }


  //! Write a binary multiNd element
  template<class T>
  inline
  void write(BinaryWriter& bin, const multiNd<T>& d)
  {
    write(bin, d.size()); // write the array of the sizes

    for(int i=0; i < d.numElem(); ++i)
      write(bin, d.getElem(i));
  }


  //! Write all of a binary std::vector object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const std::vector<T>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size
    for(int i=0; i < d.size(); ++i)
      write(bin, d[i]);
  }


  //! Write all of a binary std::list object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const std::list<T>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size
    for(typename std::list<T>::const_iterator p = d.begin(); p != d.end(); ++p)
      write(bin, *p);
  }


  //! Write all of a binary Array object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<class T>
  inline
  void write(BinaryWriter& bin, const Array1dO<T>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size
    for(int i=1; i <= d.size(); ++i)
      write(bin, d[i]);
  }


  //! Write all of a binary std::map object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<typename K, typename V>
  inline
  void write(BinaryWriter& bin, const std::map<K,V>& d)
  {
    int n = d.size();
    write(bin, n);    // always write the size

    for(typename std::map<K,V>::const_iterator v = d.begin();
	v != d.end();
	++v)
    {
      write(bin, v->first);
      write(bin, v->second);
    }
  }


  //! Write all of a binary std::pair object
  /*!
    This also writes the number of elements to the file.
    \param bin The initialised binary reader
    \param d The data to be filled.

    \pre The binary reader must have opened the file.
  */
  template<typename T1, typename T2>
  inline
  void write(BinaryWriter& bin, const std::pair<T1,T2>& d)
  {
    write(bin, d.first);
    write(bin, d.second);
  }


  //--------------------------------------------------------------------------------
  //!  Binary buffer output class
  /*!
    This class is used to write data to a binary buffer. The data in the buffer
    is big-endian. If the host nachine is little-endian, the data
    is byte-swapped.   Output is done from the primary node only.

    Buffers need to be opened before any of the write methods are used  
  
    The write methods are also wrapped by externally defined functions
    and << operators,   
  */
  class BinaryBufferWriter : public BinaryWriter
  {
  public:
    BinaryBufferWriter();

    //! Construct from a string
    explicit BinaryBufferWriter(const std::string& s);

    //! Closes the buffer
    ~BinaryBufferWriter();

    //! Construct from a string
    void open(const std::string& s);

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Return entire buffer as a string
    std::string strPrimaryNode() const;
        
    //! Flushes the buffer
    void flush() {}

    //! Clear the buffer
    void clear();

  protected:
    //! Get the current checksum to modify
    QDPUtil::n_uint32_t& internalChecksum() {return checksum;}
  
    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

  private:
    //! Checksum
    QDPUtil::n_uint32_t checksum;
    std::ostringstream f;
  };

  //--------------------------------------------------------------------------------
  //!  Binary file output class
  /*!
    This class is used to write data to a binary file. The data in the file
    is big-endian. If the host nachine is little-endian, the data
    is byte-swapped.   Output is done from the primary node only.

    Files need to be opened before any of the write methods are used  
  
    The write methods are also wrapped by externally defined functions
    and << operators,   
  */
  class BinaryFileWriter : public BinaryWriter
  {
  public:
    explicit BinaryFileWriter();

    /*!
      Closes the last file opened
    */
    ~BinaryFileWriter();

    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    explicit BinaryFileWriter(const std::string& p);

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();
    
    /*!
      Opens a file for writing.
      \param p The name of the file
    */
    void open(const std::string& p);

    //! Closes the last file opened   
    void close();

    //! Flushes the buffer
    void flush();

  protected:
    //! Get the current checksum to modify
    QDPUtil::n_uint32_t& internalChecksum() {return checksum;}
  
    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

  private:
    //! Checksum
    QDPUtil::n_uint32_t checksum;
    std::ofstream f;
  };


  //--------------------------------------------------------------------------------
  //!  Binary input/output base class
  /*!
    This class is used to read/write data from a binary object. The data in the file
    is assumed to be big-endian. If the host nachine is little-endian, the data
    is byte-swapped. 
  
    The read methods are also wrapped by externally defined functions
    and >> operators,   
  */
  class BinaryReaderWriter : public BinaryReader, public BinaryWriter
  {
  public:
    typedef std::iostream::pos_type pos_type;  // position in buffer
    typedef std::iostream::off_type off_type;  // offset in buffer

  public:
    /*!
      Shutdown the object
    */
    virtual ~BinaryReaderWriter() {}

    //!Checks status of the previous IO operation.
    /*!
      \return true if some failure occurred in the previous IO operation
    */
    virtual bool fail();

    //! Get the current checksum
    virtual QDPUtil::n_uint32_t getChecksum();
  
    //! Reset the current checksum
    virtual void resetChecksum();
  
    //! Current position
    virtual pos_type currentPosition();

    //! Set the current position
    /*! The checksum is reset */
    virtual void seek(pos_type off);

    //! Set the position relative from the start
    /*! The checksum is reset */
    virtual void seekBegin(off_type off);

    //! Set the position relative to the current position
    /*! The checksum is reset */
    virtual void seekRelative(off_type off);

    //! Set the position relative from the end
    /*! The checksum is reset */
    virtual void seekEnd(off_type off);

    //! Rewind object to the beginning
    /*! The checksum is reset */
    virtual void rewind();


  protected:
    //! Get the current checksum to modify
    virtual QDPUtil::n_uint32_t& internalChecksum() = 0;
  
    //! Get the internal input stream
    virtual std::istream& getIstream() = 0;
  
    //! Get the internal output stream
    virtual std::ostream& getOstream() = 0;
  
    //! Get the internal input/output stream
    virtual std::iostream& getIOstream() = 0;
  };


  //--------------------------------------------------------------------------------
  //!  Binary buffer input/ouput class
  /*!
    This class is used to read data from a binary buffer.
  
    The read methods are also wrapped by externally defined functions
    and >> operators,   
  */
  class BinaryBufferReaderWriter : public BinaryReaderWriter
  {
  public:
    BinaryBufferReaderWriter();

    //! Construct from a string
    explicit BinaryBufferReaderWriter(const std::string& s);

    //! Closes the buffer
    ~BinaryBufferReaderWriter();

    //! Construct from a string
    void open(const std::string& s);

    //! Return entire buffer as a string
    std::string str() const;
        
    //! Return entire buffer as a string
    std::string strPrimaryNode() const;
        
    //! Flushes the buffer
    void flush() {}

    //! Clear the buffer
    void clear();

  protected:
    //! Get the current checksum to modify
    QDPUtil::n_uint32_t& internalChecksum() {return checksum;}
  
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

    //! Get the internal inpu/output stream
    std::iostream& getIOstream() {return f;}

  private:
    //! Checksum
    QDPUtil::n_uint32_t checksum;
    std::stringstream f;
  };


  //--------------------------------------------------------------------------------
  //!  Binary file input/output class
  /*!
    This class is used to read data from a binary file. The data in the file
    is assumed to be big-endian. If the host machine is little-endian, the data
    is byte-swapped. All nodes end up with the same data
  
    The read methods are also wrapped by externally defined functions
    and >> operators,   
  */
  class BinaryFileReaderWriter : public BinaryReaderWriter
  {
  public:
    BinaryFileReaderWriter();

    /*!
      Closes the last file opened
    */
    ~BinaryFileReaderWriter();

    /*!
      Opens a file for reading.
      \param p The name of the file
    */
    explicit BinaryFileReaderWriter(const std::string& p, std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out);

    //! Queries whether the file is open
    /*!
      \return true if the file is open; false otherwise.
    */
    bool is_open();

    //! Opens a file for reading.
    /*!
      \param p The name of the file
    */
    void open(const std::string& p, std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out);

    //! Closes the last file opened
    void close();

    //! Flushes the buffer
    void flush();

  protected:
    //! Get the current checksum to modify
    QDPUtil::n_uint32_t& internalChecksum() {return checksum;}
  
    //! Get the internal input stream
    std::istream& getIstream() {return f;}

    //! Get the internal output stream
    std::ostream& getOstream() {return f;}

    //! Get the internal inpu/output stream
    std::iostream& getIOstream() {return f;}

  private:
    //! Checksum
    QDPUtil::n_uint32_t checksum;
    std::fstream f;
  };


  /*! @} */   // end of group io
} // namespace QDP

#endif
