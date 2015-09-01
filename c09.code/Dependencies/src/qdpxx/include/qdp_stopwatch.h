// -*- C++ -*-
/*! @file
 * @brief Timer support
 *
 * A stopwatch like timer.
 */

#ifndef QDP_STOPWATCH_H
#define QDP_STOPWATCH_H

#include<sys/time.h>

namespace QDP {

/*! @defgroup timer Timer
 *
 * @ingroup qdp
 *
 * @{
 */
class StopWatch 
{
public:
  //! Constructor
  StopWatch();

  //! Destructor
  ~StopWatch();

  //! Reset the timer
  void reset();

  //! Start the timer
  void start();

  //! Stop the timer
  void stop();

  //! Get time in microseconds
  double getTimeInMicroseconds();

  //! Get time in seconds
  double getTimeInSeconds();

private:
  long sec;
  long usec;
  bool startedP;
  bool stoppedP;

  struct timeval t_start;
  struct timeval t_end;
  
  void calcDuration(long&,long&);
};

} // namespace QDP  

#endif
