// -*- C++ -*-
/*! \file
 *  \brief Useful for QDP to read/write slices of a lattice object
 */


#ifndef __qdp_disk_map_slice_h__
#define __qdp_disk_map_slice_h__

#include "qdp.h"
#include "qdp_util.h"


namespace QDP
{
  //! Time-slice wrapper
  /*! 
   * A reference object holding a lattice object. 
   * Used to read/write 3d objects from a 4d chroma. 
   */
  template <typename T>
  class TimeSliceIO
  {
    T& data;
    int current_time;
    int start_lexico;
    int stop_lexico;

  public:

    TimeSliceIO(T& lattice, int time_slice);
    ~TimeSliceIO() {}

    void setTimeSlice(int time_slice);
    void getTimeSlice() const {return current_time;}
    T&   getObject() const {return data;}

    void binaryRead(BinaryReader& bin);
    void binaryWrite(BinaryWriter& bin) const;

  };


  template <typename T>
  TimeSliceIO<T>::TimeSliceIO(T& lattice, int time_slice) : data(lattice)
  {
    setTimeSlice(time_slice);
  }

  template <typename T>
  void TimeSliceIO<T>::setTimeSlice(int time_slice)
  {
    current_time = time_slice;
    int tDir = Nd-1;

    if ((current_time < 0) || (current_time >= Layout::lattSize()[tDir]))
    {
      throw(std::string("Invalid time in Lattice TimeSliceIO"));
    }

    multi1d<int> coord(Nd);
    coord = 0;
    coord[tDir] = current_time;
    start_lexico = QDP::local_site(coord,Layout::lattSize());
    if (current_time == Layout::lattSize()[tDir]-1) 
      stop_lexico = Layout::vol();
    else
    {
      coord[tDir]++; 
      stop_lexico = QDP::local_site(coord,Layout::lattSize());
    }
  }

  template <typename T>
  void TimeSliceIO<T>::binaryRead(BinaryReader& bin)
  {
    LatticeTimeSliceIO::readSlice(bin,data,start_lexico,stop_lexico);
  }

  template <typename T>
  void TimeSliceIO<T>::binaryWrite(BinaryWriter& bin) const
  {
    LatticeTimeSliceIO::writeSlice(bin,data,start_lexico,stop_lexico);
  }


  // overloaded read/write functions (used by MapObjDisk)
  template <typename T>
  void read(BinaryReader& bin, TimeSliceIO<T>& val)
  {
    val.binaryRead(bin);
  }

  template <typename T>
  void write(BinaryWriter& bin, const TimeSliceIO<T>& val)
  {
    val.binaryWrite(bin);
  }

} 
#endif
