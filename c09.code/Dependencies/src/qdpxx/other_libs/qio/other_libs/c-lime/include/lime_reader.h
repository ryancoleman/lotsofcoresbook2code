#ifndef LIME_READER_H
#define LIME_READER_H

#include <stdio.h>
#include <lime_defs.h>
#include <lime_header.h>
#include <sys/types.h>


/** A structure to hold the state of the reader */
typedef struct { 
  int first_read;        /**< Have we read the first record header yet? */
  int is_last;           /**< Is  the last record */
  int header_nextP;      /**< Are we supposed to be reading a header? */
  FILE *fp;              /**< The input file stream */
  LimeRecordHeader *curr_header; /**< The current record header */

  n_uint64_t bytes_left;      /**< Data bytes still unread in the current record */
  n_uint64_t bytes_total;  /**< Total data bytes in the current record */
  n_uint64_t rec_ptr;         /**< Next byte to be read relative to 
			      the start of the record payload.
			      ranges 0 to bytes_total - 1 */
  n_uint64_t rec_start;       /**< File pointer at start of record payload */
  size_t bytes_pad;      /**< Padding bytes at end of current record */
} LimeReader;

/** \brief Create a LIME reader 
 *  
 *  \param fp is the file pointer to the open file to read from
 *  \returns a LimerReader*, on error NULL is returned
 */
LimeReader* limeCreateReader(FILE *fp);

/** \brief Set file pointer to a new position (assumed to start a LIME record)
 *  
 *  \param r is the LIME reader
 *  \param offset is the new file pointer.
 *  \returns a LimerReader*, on error NULL is returned
 */
int limeSetReaderPointer(LimeReader *r, off_t offset);

/** \brief Get the file pointer of the next LIME record
 *  
 *  \param  r is the LIME reader
 */
off_t limeGetReaderPointer(LimeReader *r);

/** \brief Destroy a LIME reader 
 *
 *  \param r is a pointer to the LimeReader to be destroyed 
 */
void limeDestroyReader(LimeReader *r);

/** \brief Skip to Next Record
 *
 * If there is a next record, we will skip to it.
 * By this we mean:
 *   i) If there is data left in the current record, we will read it
 *      and throw it away. (Unrecoverably)
 *
 *  ii) We read the header of the next record (if there is one)
 *      and check it for draft internet standard compliance.
 *
 *  iii) We position the stream pointer at the start of the data
 *       portion of the next record (if there is one)
 *
 * \params r is a pointer to a LimeReader 
 * \return a status code to indicate success or failure.
 * 
 */
int limeReaderNextRecord(LimeReader *r);

/** \brief Accessor for MB Flag in current header
 *  \params r points to a LimeReader
 *  \returns -1 if r is null, otherwise the MB_flag
 */
int limeReaderMBFlag(LimeReader *r);

/** \brief Accessor for ME Flag in current header
 *  \params r points to a LimeReader
 *  \returns -1 if r is null, otherwise the MB_flag
 */
int limeReaderMEFlag(LimeReader *r);

/** \brief Accessor for type string in current header
 *  \params r points to a LimeReader
 *  \returns NULL if r is null, otherwise the MB_flag
 */
char *limeReaderType(LimeReader *r);

/** \brief Accessor for number of total bytes in current record
 *  \params r points to a LimeReader
 *  \returns 0 if r is null, otherwise the byte count
 */
n_uint64_t limeReaderBytes(LimeReader *r);

/** \brief Accessor for number of pad bytes in current record
 *  \params r points to a LimeReader
 *  \returns 0 if r is null, otherwise the byte count
 */
size_t limeReaderPadBytes(LimeReader *r);

/** \brief Read data from a record 
 *  \params dest is the destination buffer for the bytes to read.
 *  \params nbytes is a pointer to the integer containing the no of 
 *          bytes to read.
 *  \params r is the LimeReader to read from 
 *
 *  \returns a status code, and sets nbytes to the actual number of 
 *           bytes written.
 */
int limeReaderReadData(void *dest, n_uint64_t *nbytes, LimeReader *r);

/** \brief Advance to end of current record
 *  \params r points to a LimeReader
 */

int limeReaderCloseRecord(LimeReader *r);

/** \brief Seek bytes within current open record 
 *  \params r is a pointer to a LimeReader
 *  \params offset counts bytes from the current position ("WHENCE")
 *
 *  \returns a status code
 */

int limeReaderSeek(LimeReader *r, off_t offset, int whence);

/** \brief Sets the LIME reader state to a prescribed state
 *         and positions the file accordingly.
 *  \params rdest our private LIME reader
 *  \params rsrc LIME reader whose state we want to duplicate.
 *  \returns a status code
 */
int limeReaderSetState(LimeReader *rdest, LimeReader *rsrc );

int limeEOM(LimeReader *r);


#endif
