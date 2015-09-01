/*! @file
 * @brief Timer support
 *
 * A stopwatch like timer.
 */


#include "qdp.h"
#include<sys/time.h>

namespace QDP {

  StopWatch::StopWatch() 
  {
    stoppedP=false;
    startedP=false;
    sec=0;
    usec=0;
  }

  StopWatch::~StopWatch() {}
	
  void StopWatch::calcDuration(long& secs, long& usecs){
    secs=0;
    usecs=0;
    if( startedP && stoppedP ) 
    { 
      if( t_end.tv_sec < t_start.tv_sec ) 
      { 
	QDPIO::cerr << __func__ << ": critical timer rollover" << std::endl;
	usecs = 0;
      }
      else 
      { 
	secs = t_end.tv_sec - t_start.tv_sec;

	if( t_end.tv_usec < t_start.tv_usec ) 
	{
	  secs -= 1;
	  usecs = 1000000;
	}
	usecs += t_end.tv_usec - t_start.tv_usec;
      }
    }
    else 
    {
      QDPIO::cerr << __func__ << ": either stopwatch not started, or not stopped. I will return 0 instead!" << std::endl;
      secs=0;
      usecs=0;
    }
  }

  void StopWatch::reset() 
  {
    startedP = false;
    stoppedP = false;
    sec=0;
    usec=0;
  }

  void StopWatch::start() 
  {
    int ret_val;
    ret_val = gettimeofday(&t_start, NULL);
    if( ret_val != 0 ) 
    {
      QDPIO::cerr << __func__ << ": gettimeofday failed in StopWatch::start()" << std::endl;
      QDP_abort(1);
    }
    startedP = true;
    stoppedP = false;
  }

  void StopWatch::stop() 
  {
    if( !startedP ) 
    { 
      QDPIO::cerr << __func__ << ": attempting to stop a non running stopwatch in StopWatch::stop()" << std::endl;
    }
    else{
      int ret_val;
      ret_val = gettimeofday(&t_end, NULL);
      if( ret_val != 0 ) 
      {
	QDPIO::cerr << __func__ << ": gettimeofday failed in StopWatch::end()" << std::endl;
	QDP_abort(1);
      }
      stoppedP = true;
		
      long usecs, secs;
      calcDuration(secs,usecs);
			
      usec+=usecs;
      sec+=secs;
		
      startedP=false;
      stoppedP=false;
    }
  }

  double StopWatch::getTimeInMicroseconds() 
  {
    return static_cast<double>(usec)+static_cast<double>(sec*1.e6);
  }
    
  double StopWatch::getTimeInSeconds()  
  {
    return static_cast<double>(sec)+static_cast<double>(usec)/1.e6;
  }


} // namespace QDP;
