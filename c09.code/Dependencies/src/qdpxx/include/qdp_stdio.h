// -*- C++ -*-

/*! @file
 * @brief Parallel version of stdio
 *
 * These are the QDP parallel versions of the STL cout,cerr,cin.
 * The user can control which nodes are involved in output/input.
 */

#ifndef QDP_STDIO_H
#define QDP_STDIO_H

#include <istream>
#include <ostream>
//#include <streambuf>
#include <string>

namespace QDP {


/*! @defgroup stdio STD IO
 *
 * Standard-IO-like input and output operations on QDP types
 *
 * @{
 */

//--------------------------------------------------------------------------------
//! StandardInputStream class
/*! Parallel version of standard input
 */
class StandardInputStream
{
public:
  //! Constructor
  StandardInputStream();

  //! destructor
  ~StandardInputStream();

  //! Constructor from a stream
  void init(std::istream *is);

//  //! Redirect input stream
//  void rdbuf(std::streambuf* b);

  //! Return true if some failure occurred in previous IO operation
  bool fail();

  // Overloaded input functions
  StandardInputStream& operator>>(std::string& input);
  StandardInputStream& operator>>(char& input);
  StandardInputStream& operator>>(int& input);
  StandardInputStream& operator>>(unsigned int& input);
  StandardInputStream& operator>>(short int& input);
  StandardInputStream& operator>>(unsigned short int& input);
  StandardInputStream& operator>>(long int& input);
  StandardInputStream& operator>>(unsigned long int& input);
  StandardInputStream& operator>>(long long int& input);
  StandardInputStream& operator>>(float& input);
  StandardInputStream& operator>>(double& input);
  StandardInputStream& operator>>(long double& input);
  StandardInputStream& operator>>(bool& input);


private:
  //! Hide copy constructor and =
  StandardInputStream(const StandardInputStream&) {}
  void operator=(const StandardInputStream&) {}

protected:
  // The universal data-reader. All the read functions call this
  template<typename T>
  StandardInputStream& readPrimitive(T& output);

  // Get the internal istream
  std::istream& getIstream() {return *is;}

  //! Is the stream open?
  bool is_open() {return open;}

private:
  std::istream* is;
  bool open;
};


//! Read an array
template<class T>
inline
StandardInputStream& operator>>(StandardInputStream& s, multi1d<T>& d)
{
  for(int i=0; i < d.size(); ++i)
    s >> d[i];

  return s;
}



//--------------------------------------------------------------------------------
//! StandardOutputStream class
/*! Parallel version of standard output
 */
class StandardOutputStream
{
public:
  //! Constructor
  StandardOutputStream();

  ~StandardOutputStream();

  //! Constructor from a stream
  void init(std::ostream *os);

//  //! Redirect output stream
//  void rdbuf(std::streambuf* b);

  //! Flush the buffer
  void flush();

  //! Return true if some failure occurred in previous IO operation
  bool fail();

  //! Hook for manipulators
  StandardOutputStream& operator<<(std::ostream& (*op)(std::ostream&));

  // Overloaded output functions
  StandardOutputStream& operator<<(const std::string& output);
  StandardOutputStream& operator<<(const char* output);
  StandardOutputStream& operator<<(char output);
  StandardOutputStream& operator<<(int output);
  StandardOutputStream& operator<<(unsigned int output);
  StandardOutputStream& operator<<(short int output);
  StandardOutputStream& operator<<(unsigned short int output);
  StandardOutputStream& operator<<(long int output);
  StandardOutputStream& operator<<(unsigned long int output);
  StandardOutputStream& operator<<(long long int output);
  StandardOutputStream& operator<<(float output);
  StandardOutputStream& operator<<(double output);
  StandardOutputStream& operator<<(long double output);
  StandardOutputStream& operator<<(bool output);


private:
  //! Hide copy constructor and =
  StandardOutputStream(const StandardOutputStream&) {}
  void operator=(const StandardOutputStream&) {}

protected:
  // The universal data-writer. All the write functions call this
  template<typename T>
  StandardOutputStream& writePrimitive(T output);

  // Get the internal ostream
  std::ostream& getOstream() {return *os;}

  //! Is the stream open?
  bool is_open() {return open;}

private:
  std::ostream* os;
  bool open;
};


// NOTE: A write of an array is *NOT* provided since these are usual special purpose
//template<class T>
//inline
//StandardOutputStream& operator<<(StandardOutputStream& s, const multi1d<T>& d);



// Make global (parallel) versions of stdin, stdout, stderr
// Put this in a different namespace to avoid collisions
namespace QDPIO
{
  extern StandardInputStream  cin;
  extern StandardOutputStream cout;
  extern StandardOutputStream cerr;
}


/*! @} */   // end of group stdio

} // namespace QDP

#endif
