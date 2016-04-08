/***********************************************************************
 *   DL_MESO       Version 2.6                                         *
 *   Authors   :   R. S. Qin, M. A. Seaton                             *
 *   Copyright :   STFC Daresbury Laboratory                           *
 *             :   07/08/2014                                          *
 ***********************************************************************
 *   This extract used with permission of the copyright owners         *
 ***********************************************************************
 *   The full DL_MESO package is supplied to individuals under an      *
 *   academic licence, which is free of cost to academic scientists    *
 *   pursuing scientific research of a non-commercial nature           *
 *   To register please visit www.ccp5.ac.uk/DL_MESO                   *
 *   Commercial organisations interested in acquiring the package      *
 *   should approach Dr. M. A. Seaton at Daresbury Laboratory in the   *
 *   first instance. Daresbury Laboratory is the sole centre for       * 
 *   distribution of the package                                       *
 ***********************************************************************/


// Functions in this package are general purpose and are not directly 
// related to Mesoscale Methods. These may be replaced with any suitable 
// functions in C++ standard libraries.

template <class T>
T fCppAbs(T a){return (a<0)?-a:a;}


template <class T>
T fReciprocal(T a){return (a!=0)?1/a:0;}


template <class T>
T fEvapLimit(T a){return (a<lbevaplim)?0:a;}


template <class T> void fSwapPair ( T& a, T& b)
{
  T c(a); a=b; b=c;
}


template <typename T> T fStringToNumber ( const string &text )
  {
     istringstream ss(text);
     T result;
     return ss >> result ? result : 0;
  }


template <typename T> T fCppMax ( T& a, T& b)
{
    return (a>b?a:b);
}


template <typename T> T fCppMin ( T& a, T& b)
{
    return (a<b?a:b);
}










int fCppMod(int a, int b)
{

  // ensure integer number a is in a range between 0 and b-1

  if( a < 0 )
    return (a+b);

  else if(a >= b)
    return (a-b);

  else
    return a;
}




int fPrintLine()
{
  cout << string(76, '-') << endl;
  return 0;
}


int fPrintDoubleLine()
{
  cout << string(76, '=') << endl;
  return 0;
}

double fRandom()
{

  //  create a random number between -1 and 1

  int seed = 50200 + lbdm.rank;
  static bool notFirstTime1;
  double randomNumber;
  static int theNumber;
  if(!notFirstTime1) theNumber=seed;   
  theNumber = 2045*theNumber + 1;
  theNumber = theNumber - (theNumber/1048576)*1048576;
  randomNumber=2.0*(theNumber+1)/1048577.0-1.0;
  notFirstTime1=1;
  return randomNumber;  
}


int fBigEndian()
{

  // indicate endianness of machine: 1 = big endian, 0 = little endian

  short int n = 1;
  char *ep = (char *)&n;

  return (*ep == 0);
}


#ifdef __linux__
double fCheckTimeSerial()
{

  // check the elapsed time (in seconds)

  static int timeini;
  static double starttime;
  double timeelapsed;
  struct timeval Tp;

  gettimeofday(&Tp, NULL);

  if(timeini == 0)
    {  
      timeini = 1;
      starttime = Tp.tv_sec + 0.000001 * (double) (Tp.tv_usec);
      timeelapsed = 0.0;
    }
  else
    timeelapsed = Tp.tv_sec + 0.000001 * (double) (Tp.tv_usec) - starttime;
  
  return timeelapsed;
}
#else


double fCheckTimeSerial()
{
  // check the elapsed time (in seconds)
  static int timeini;
  static double starttime;
  double timeelapsed;
  LARGE_INTEGER ticks;
  LARGE_INTEGER frequency;
  QueryPerformanceCounter(&ticks);
  QueryPerformanceFrequency(&frequency);
  
  if(timeini == 0)
    {  
      timeini = 1;
      starttime = (double)(ticks.QuadPart/(double)frequency.QuadPart);
      timeelapsed = 0.0;
    }
  else
    timeelapsed = (double)(ticks.QuadPart/(double)frequency.QuadPart) - starttime;
  
  return timeelapsed;
}
#endif

string fReadString(string line, int i)
{
    // return i-th word of input string (line)

    string buf;
    vector<string> tokens;
    stringstream ss(line);

    while (ss >> buf)
        tokens.push_back(buf);

    if(i>=tokens.size()) {
      buf = string();
    }
    else {
      buf = tokens.at(i);
    }

    return(buf);
}


