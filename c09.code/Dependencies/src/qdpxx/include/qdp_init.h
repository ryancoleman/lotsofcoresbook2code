// -*- C++ -*-

/*! \file
 * \brief Routines for top level QDP management
 *
 * Fundamental routines for turning on/off and inserting/extracting
 * variables.
 */

// Info/error routines
namespace QDP {

//! Turn on the machine
void QDP_initialize (int *argc, char ***argv);

//! Is the machine initialized?
bool QDP_isInitialized ();

//! Turn off the machine
void QDP_finalize ();

//! Panic button
void QDP_abort (int status);

//! Simple information display routine
int  QDP_info (const char* format, ...);

//! Simple error display routine
int  QDP_error (const char* format, ...);

//! Simple error display and abort routine
void QDP_error_exit (const char *format, ...);

//! Resumes QDP communications
void QDP_resume();

//! Suspends QDP communications
void QDP_suspend();

}

