
/*! @file
 * @brief Parscalarvec init routines
 * 
 * Initialization routines for parscalarvec implementation
 */

#include <stdlib.h>
#include <unistd.h>

#include "qdp.h"
#include "qmp.h"

namespace QDP {

//! Private flag for status
static bool isInit = false;

//! Turn on the machine
void QDP_initialize(int *argc, char ***argv)
{
  if (isInit)
  {
    QDPIO::cerr << "QDP already inited" << std::endl;
    QDP_abort(1);
  }

  //
  // Process command line
  //

  // Look for help
  bool help_flag = false;
  for (int i=0; i<*argc; i++) 
  {
    if (strcmp((*argv)[i], "-h")==0)
      help_flag = true;
  }

  bool setGeomP = false;
  multi1d<int> logical_geom(Nd);   // apriori logical geometry of the machine
  logical_geom = 0;

  int rtiP = 0;
  int QMP_verboseP = 0;
  const int maxlen = 256;
  char rtinode[maxlen];
  strncpy(rtinode, "your_local_food_store", maxlen);

  // Usage
  if (Layout::primaryNode()) 
    if (help_flag) 
    {
      fprintf(stderr,"Usage:    %s options\n",(*argv)[0]);
      fprintf(stderr,"options:\n");
      fprintf(stderr,"    -h        help\n");
      fprintf(stderr,"    -V        %%d [%d] verbose mode for QMP\n", 
	      QMP_verboseP);
#if defined(QDP_USE_PROFILING)   
      fprintf(stderr,"    -p        %%d [%d] profile level\n", 
	      getProfileLevel());
#endif

      // logical geometry info
      fprintf(stderr,"    -geom     %%d");
      for(int i=1; i < Nd; i++) 
	fprintf(stderr," %%d");

      fprintf(stderr," [-1");
      for(int i=1; i < Nd; i++) 
	fprintf(stderr,",-1");
      fprintf(stderr,"] logical machine geometry\n");

#ifdef USE_REMOTE_QIO
      fprintf(stderr,"    -cd       %%s [.] set working dir for QIO interface\n");
      fprintf(stderr,"    -rti      %%d [%d] use run-time interface\n", 
	      rtiP);
      fprintf(stderr,"    -rtinode  %%s [%s] run-time interface fileserver node\n", 
	      rtinode);
#endif
      
      QDP_abort(1);
    }

  for (int i=1; i<*argc; i++) 
  {
    if (strcmp((*argv)[i], "-V")==0) 
    {
      QMP_verboseP = 1;
    }
#if defined(QDP_USE_PROFILING)   
    else if (strcmp((*argv)[i], "-p")==0) 
    {
      int lev;
      sscanf((*argv)[++i], "%d", &lev);
      setProgramProfileLevel(lev);
    }
#endif
    else if (strcmp((*argv)[i], "-geom")==0) 
    {
      setGeomP = true;
      for(int j=0; j < Nd; j++) 
      {
	int uu;
	sscanf((*argv)[++i], "%d", &uu);
	logical_geom[j] = uu;
      }
    }
#ifdef USE_REMOTE_QIO
    else if (strcmp((*argv)[i], "-cd")==0) 
    {
      /* push the dir into the environment vars so qio.c can pick it up */
      setenv("QHOSTDIR", (*argv)[++i], 0);
    }
    else if (strcmp((*argv)[i], "-rti")==0) 
    {
      sscanf((*argv)[++i], "%d", &rtiP);
    }
    else if (strcmp((*argv)[i], "-rtinode")==0) 
    {
      int n = strlen((*argv)[++i]);
      if (n >= maxlen)
      {
	QDPIO::cerr << "rtinode name too long" << std::endl;
	QDP_abort(1);
      }
      sscanf((*argv)[i], "%s", rtinode);
    }
#endif
#if 0
    else 
    {
      QDPIO::cerr << __func__ << ": Unknown argument = " << (*argv)[i] << std::endl;
      QDP_abort(1);
    }
#endif
 
    if (i >= *argc) 
    {
      QDPIO::cerr << __func__ << ": missing argument at the end" << std::endl;
      QDP_abort(1);
    }
  }


  QMP_verbose (QMP_verboseP);

#if QDP_DEBUG >= 1
  // Print command line args
  for (int i=0; i<*argc; i++) 
    QDP_info("QDP_init: arg[%d] = XX%sXX",i,(*argv)[i]);
#endif

#if QDP_DEBUG >= 1
  QDP_info("Now initialize QMP");
#endif

  if (QMP_is_initialized() == QMP_FALSE)
  {
    QMP_thread_level_t prv;
    if (QMP_init_msg_passing(argc, argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS)
    {
      QDPIO::cerr << "QMP_init_msg_passing failed" << std::endl;
      QDP_abort(1);
    }
  }

#if QDP_DEBUG >= 1
  QDP_info("QMP inititalized");
#endif

  if (setGeomP)
    if (QMP_declare_logical_topology(logical_geom.slice(), Nd) != QMP_SUCCESS)
    {
      QDPIO::cerr << "QMP_declare_logical_topology failed" << std::endl;
      QDP_abort(1);
    }

#if QDP_DEBUG >= 1
  QDP_info("Some layout init");
#endif

  Layout::init();   // setup extremely basic functionality in Layout

  isInit = true;

#if QDP_DEBUG >= 1
  QDP_info("Init qio");
#endif

  // initialize the global streams
  QDPIO::cin.init(&std::cin);
  QDPIO::cout.init(&std::cout);
  QDPIO::cerr.init(&std::cerr);

  initProfile(__FILE__, __func__, __LINE__);

  QDPIO::cout << "Initialize done" << std::endl;
}

//! Is the machine initialized?
bool QDP_isInitialized() {return isInit;}

//! Turn off the machine
void QDP_finalize()
{
  if ( ! QDP_isInitialized() )
  {
    QDPIO::cerr << "QDP is not inited" << std::endl;
    QDP_abort(1);
  }

  printProfile();

#if defined(QDP_USE_HDF5)
  H5close();
#endif

  QMP_finalize_msg_passing();

  isInit = false;
}

//! Panic button
void QDP_abort(int status)
{
  QMP_abort(status); 
}

//! Resumes QDP communications
void QDP_resume() {}

//! Suspends QDP communications
void QDP_suspend() {}


} // namespace QDP;
