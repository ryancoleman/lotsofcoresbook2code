/*----------------------------------------------------------------------------
 * Copyright (c) 2007      Jefferson Science Associates, LLC               
 *                         Under U.S. DOE Contract No. DE-AC05-06OR23177  
 *                                                                      
 *                         Thomas Jefferson National Accelerator Facility
 *
 *                         Jefferson Lab 
 *                         Scientific Computing Group,
 *                         12000 Jefferson Ave.,      
 *                         Newport News, VA 23606 
 *----------------------------------------------------------------------------
 *  
 * Description:
 *     Template Functions for Database Storage and Retrieve
 *
 *     Functions are dealing with classes derived from DBKey and DBData
 * 
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: DBCursor.cpp,v $
 *   Revision 1.2  2009-03-04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *   Revision 1.1.1.1  2008/05/29 19:31:39  chen
 *   Initial CVS import of FFDB
 *
 *
 */
#include "DBCursor.h"

using namespace std;

namespace FILEDB
{
  /***********************************************************************
   * Implementation of wrapper class for Database cursor with reference  *
   ***********************************************************************/
  DBCursorRep::DBCursorRep (void)
    :count_ (1), cursor_ (0)
  {
    // empty
  }

  DBCursorRep::DBCursorRep (ffdb_cursor_t* cursor)
    :count_(1), cursor_ (cursor)
  {
    // empty
  }


  DBCursorRep::DBCursorRep (const DBCursorRep& c)
    :count_ (1)
  {
    if (cursor_)
      cursor_->close (cursor_);
    cursor_ = c.cursor_;
  }

  DBCursorRep&
  DBCursorRep::operator = (const DBCursorRep& c)
  {
    if (this != &c) {
      if (cursor_ != c.cursor_) {
	if (cursor_)
	  cursor_->close(cursor_);
      }
      cursor_ = c.cursor_;
      count_ = 1;
    }
    return *this;
  }

  DBCursorRep::~DBCursorRep(void)
  {
    if (cursor_) {
      cursor_->close (cursor_);
    }
  }


  /*******************************************************************
   *   Implmentation of DBCursor Class on top of DBCursorRep         *
   *******************************************************************/
  DBCursor::DBCursor (void)
  {
    rep_ = new DBCursorRep ();
  }
  
  DBCursor::DBCursor (ffdb_cursor_t* cursor)
  {
    rep_ = new DBCursorRep (cursor);
  }

  DBCursor::DBCursor (const DBCursor& c)
  {
    rep_ = c.rep_;
    rep_->count_++;
  }

  DBCursor&
  DBCursor::operator = (const DBCursor& c)
  {
    if (this != &c) {
      c.rep_->count_++;
      if (--rep_->count_ <= 0) delete rep_;
      rep_ = c.rep_;
    }
    return *this;
  }

  DBCursor::~DBCursor (void)
  {
    if (--rep_->count_ <= 0)
      delete rep_;
  }
}
