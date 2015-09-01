// -*- C++ -*-

/*! @file
 * @brief Profiling info
 *
 * Diagnostics to print profiling info.
 */

#include "qdp.h"
#include <time.h>

#if defined(QDP_USE_PROFILING)   
#include <stack>
#include <queue>
#endif


namespace QDP {

static int prof_level = 0;
static int prog_prof_level = 0;
#ifdef QDP_USE_PROFILING
static bool prof_init = false;
#endif
int getProfileLevel() {return prof_level;}
int getProgramProfileLevel() {return prog_prof_level;}

QDPTime_t
getClockTime()
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return((unsigned long)tp.tv_sec * 1000000 + (unsigned long)tp.tv_usec);
  //struct tms buf;
  //times(&buf);
  //return buf.tms_utime + buf.tms_stime;
  //  return (double)clock()/CLOCKS_PER_SEC;
}


//--------------------------------------------------------------------------------------
// Selectively turn on profiling
//--------------------------------------------------------------------------------------

#if ! defined(QDP_USE_PROFILING)   

//--------------------------------------------------------------------------------------
// No profiling
//--------------------------------------------------------------------------------------

void initProfile(const std::string& file, const std::string& caller, int line) {}
void closeProfile() {}
void printProfile() {}
int setProfileLevel(int n) {return prof_level;}
int setProgramProfileLevel(int n) {return prog_prof_level;}

void pushProfileInfo(int level, const std::string& file, const std::string& caller, int line) {}
void popProfileInfo() {}

#else   // Profiling enabled

//--------------------------------------------------------------------------------------
// Profiling enabled
//--------------------------------------------------------------------------------------

// A stack to hold profile info
std::stack<QDPProfileInfo_t> infostack;

// A queue for the profile data
std::queue<QDPProfileHead_t> profqueue;



void
QDPProfile_t::init()
{
  time = 0;
  expr = "";
  count = 0;
  next = 0;
}

void 
QDPProfile_t::print()
{
  if ((getProfileLevel() & 2) > 0)
  {
    QDPIO::cout << expr
		<< "\t[" << time << "]" << std::endl;
  }
}

void 
initProfile(const std::string& file, const std::string& caller, int line)
{
  if (prof_init)
    return;

  pushProfileInfo(0, file, caller, line);
  prof_init = true;
}

void 
closeProfile()
{
  if (! prof_init)
    return;

  popProfileInfo();
}

void 
printProfile()
{
  if (profqueue.size() == 0)
    return;

  while(! profqueue.empty())
  {
    QDPProfileHead_t& head = profqueue.front();

    if (head.start != 0)
    {
      QDPProfile_t *qp = head.start;

      QDPIO::cout << std::endl
		  << "func = " << head.info.caller 
		  << "   line = " << head.info.line 
		  << "   file = " << head.info.file
		  << std::endl;

      while(qp)
      {
	if (qp->count > 0)
	{
//	  QDPIO::cout << "   " << qp->count << "  [" << qp->time << "]  " << qp->expr << std::endl;

	  // Gag, I still have not really grokked the fricken SL iostream field
	  // width manipulators. Drop down to low tech just to get out
	  // the bleaping thing. Note: the real problem is the QDPIO::cout class
	  // is not really a stream instantiation, hence the needed member functions
	  // are not there. Sigh.
	  char lin[80];  // more than adequate
	  sprintf(lin, "  %7d   [%8d]  ", qp->count, qp->time);
	  QDPIO::cout << lin << qp->expr << std::endl;
	}

	qp = qp->next;
      }
    }

    profqueue.pop();
  }
}

int
setProfileLevel(int n)
{
  int old = prof_level;
  prof_level = n;
  return old;
}

int
setProgramProfileLevel(int n)
{
  int old = prog_prof_level;
  prog_prof_level = n;
  return old;
}

void pushProfileInfo(int level, const std::string& file, const std::string& caller, int line)
{
  QDPProfileInfo_t info(level, file, caller, line);
  infostack.push(info);

  setProfileLevel(level);

  QDPProfileHead_t  head(infostack.top());
  profqueue.push(head);

  if ((level & 2) > 0)
  {
    QDPIO::cout << std::endl
		<< "func = " << head.info.caller 
		<< "   line = " << head.info.line 
		<< "   file = " << head.info.file
		<< std::endl;
  }
}

void popProfileInfo()
{
  if (infostack.empty())
  {
    QDPIO::cerr << "popProfileInfo: invalid pop" << std::endl;
    QDP_abort(1);
  }

  QDPProfileInfo_t& info = infostack.top();
  setProfileLevel(info.level);

  infostack.pop();
}

void
registerProfile(QDPProfile_t* qp)
{
  if (profqueue.empty())
  {
    QDPIO::cerr << "registerProfile: profile queue empty" << std::endl;
    QDP_abort(1);
  }

  QDPProfileHead_t& head = profqueue.back();
  if (head.start == 0)
  {
    head.start = head.end = qp;
  }
  else
  {
    head.end->next = qp;
    head.end = qp;
  }
}

#endif  // ! defined(QDP_USE_PROFILING)

} // namespace QDP;
