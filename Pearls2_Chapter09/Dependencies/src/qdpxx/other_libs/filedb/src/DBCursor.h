// -*- C++ -*-
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
 *     
 *     Database Cursor Class
 *     
 *
 * Author:
 *     Jie Chen
 *     Scientific Computing Group
 *     Jefferson Lab
 *      
 * Revision History:
 *   $Log: DBCursor.h,v $
 *   Revision 1.4  2009-03-05 00:40:05  edwards
 *   Changed include path of filehash files back to no relative path.
 *
 *   Revision 1.3  2009/03/04 19:13:05  edwards
 *   Changed some include guards and paths to filehash to be relative.
 *
 *   Revision 1.2  2009/03/04 15:55:25  chen
 *   Change Namespace from FFDB to FILEDB
 *
 *   Revision 1.1  2009/02/20 20:44:48  chen
 *   initial import
 *
 *   Revision 1.4  2008/08/12 17:37:00  chen
 *   Remove unnecessary code
 *
 *   Revision 1.3  2008/08/02 20:27:51  edwards
 *   Changed get_next_key to new semantics.
 *
 *   Revision 1.2  2008/06/16 04:19:21  edwards
 *   Added emacs controls to throw into C++ mode.
 *
 *   Revision 1.1.1.1  2008/05/29 19:31:39  chen
 *   Initial CVS import of FFDB
 *
 *
 */
#ifndef _FILEDB_DBCURSOR_H
#define _FILEDB_DBCURSOR_H

#include "DBKey.h"
#include "DBData.h"

#include "ffdb_db.h"

namespace FILEDB 
{
  /**
   * A database cursor class with reference counting
   */
  class DBCursor;

  class DBCursorRep
  {
    DBCursorRep (void);
    DBCursorRep (ffdb_cursor_t* cursor);
    DBCursorRep (const DBCursorRep& c);
    DBCursorRep& operator = (const DBCursorRep& c);
    ~DBCursorRep (void);

    int  count_;
    ffdb_cursor_t* cursor_;
    friend class DBCursor;
  };

  /**
   * Database cursor class used by applications
   * This class should be used as dynamic variable
   */
  class DBCursor
  {
  public:
    /**
     * Default constructor
     * The internal cursor is empty one
     */
    DBCursor (void);
    
    /**
     * Constructor
     * @param cursor real Database cursor
     */
    DBCursor (ffdb_cursor_t* cursor);

    /**
     * Copy constructor
     */
    DBCursor (const DBCursor& c);
    
    /**
     * Assignment operator
     */
    DBCursor& operator = (const DBCursor& c);

    /**
     * Destructor
     */
    ~DBCursor (void);

    /**
     * Return internal database cursor pointer
     * This pointer could be null.
     */
    ffdb_cursor_t* cursor (void) const {return rep_->cursor_;}
    
  private:
    DBCursorRep* rep_;
  };


  /**
   * Template class of DB iterator
   * Template class T is key for the database. T should be derived
   * from DBKey Class. T should have =  and == operators defined
   *
   * Template class D is data for the database. D should be derived
   * from DBData Class
   *
   * The parameter T is what returned by iterator
   *
   */
  template <typename T, typename D>
  class DBKeyIterator
  {
  private:
    // internal cursor
    DBCursor cursor_;

    // current key value
    T keyvalue_;

    // end flag
    int endflag_;

    void
    get_next_key (ffdb_cursor_t* cursor)
    {
      T    mykey;
      FFDB_DBT  key, data;
      int  ret;
      
      // get everything from meta data
      key.data = data.data = 0;
      key.size = data.size = 0;
      ret = cursor->get (cursor, &key, &data, FFDB_NEXT);

      if (ret == 0) {
	// get network byte order and convert to host order
	try {
	  // convert into key object
	  std::string keyObj;
	  keyObj.assign((char*)key.data, key.size);
	  mykey.readObject (keyObj);
	}
	catch (SerializeException& e) {
	  std::cerr << e.what () << std::endl;
	  cursor->close (cursor);
	  cursor = 0;
	  ::exit (143);
	}
	free (key.data); free (data.data);
	// depends on operator = defined in T
	keyvalue_ = mykey;
      }
      else if (ret == FFDB_NOT_FOUND) {
	// this is beyond the last one
	endflag_ = 1;
      }
      else {
	std::cerr << "DB Cursor has fatal error" << std::endl;
	cursor->close (cursor);
	cursor = 0;
	::exit (1);
      }
    }

  public:
    typedef DBKeyIterator _self;
    
    /**
     * Constructor
     */
    DBKeyIterator (void)
      :cursor_(), keyvalue_(), endflag_ (0) { }

    /**
     * Explicit constructor
     */
    explicit DBKeyIterator (const DBCursor& cursor,
			    const T& arg, int endflag)
      :cursor_ (cursor), keyvalue_(arg), endflag_ (endflag) { }

    /**
     * Copy constructor
     */
    DBKeyIterator (const DBKeyIterator& _x)
      :cursor_ (_x.cursor_), keyvalue_ (_x.keyvalue_), endflag_ (_x.endflag_){}

    /**
     * Destructor
     */
    ~DBKeyIterator (void) { }

    /**
     * Get key value pointed by this iterator
     */
    T operator * (void) const {return keyvalue_;}

    /**
     * To Get a pointer of the key value of this iterator
     */
    const T* operator -> (void) const {return &keyvalue_;}

    /**
     * Prefix incremental operator
     */
    _self& operator ++ (void) 
    {
      ffdb_cursor_t* c= cursor_.cursor(); 
      if (c) {
	get_next_key (c);
      }
      return *this;
    }
    
    /**
     * Postfix incremental operator
     */
    void operator ++ (int) 
    {
      ffdb_cursor_t* c= cursor_.cursor(); 
      if (c) {
	get_next_key (c);
      }
    }

    /**
     * Comparison operators
     */
    int operator == (const _self& __x)
    {
      ffdb_cursor_t* c1 = cursor_.cursor();
      ffdb_cursor_t* c2 = __x.cursor_.cursor();
      if ((c1 && c2 ) && (endflag_ == __x.endflag_)) {
	// depends on T == operator
	if (keyvalue_ == __x.keyvalue_)
	  return 1;
      }
      else if (endflag_ == 1 && __x.endflag_ == 1)
	return 1;
    
      return 0;
    }

    int operator != (const _self& __x)
    {
      return (!(operator == (__x)));
    }


    /**
     * Assignment operator
     */
    _self& operator = (const _self& __x)
    {
      if (this != &__x) {
	cursor_ = __x.cursor_;
	// depends T = operator
	keyvalue_ = __x.keyvalue_;
	endflag_ = __x.endflag_;
      }
      return *this;
    }

    /**
     * Initialize a database key iterator
     */
    static DBKeyIterator<T, D> init (FFDB_DB* dbh)
    {
      T mykey;
      D mydata;
      FFDB_DBT  key, data;
      ffdb_cursor_t  *crp;
      int  ret;
      
      try {
	dbh->cursor (dbh, &crp, FFDB_KEY_CURSOR);
	
	key.data = data.data = 0;
	key.size = data.size = 0;
	ret = crp->get (crp, &key, &data, FFDB_NEXT);

	if (ret == 0) {
	  // convert into key object
	  std::string keyObj;
	  keyObj.assign((char*)key.data, key.size);
	  mykey.readObject (keyObj);

	  free (key.data); free(data.data);

	  return DBKeyIterator<T, D>(DBCursor(crp), mykey, 0);
	}
	else if (ret != FFDB_NOT_FOUND) {
	  std::cerr <<"Fatal: DB cursor  error " << std::endl;
	  ::exit (143);
	}
	else
	  return DBKeyIterator<T, D>(DBCursor(crp), mykey, 1);
      }
      catch (SerializeException& e) {
	std::cerr << "Fatal: DB cursor error: " << e.what () << std::endl;
      }
      return DBKeyIterator<T, D>(DBCursor(), mykey, 1);
    }
  };
}
#endif
