#ifndef LIME_WRITER_H
#define LIME_WRITER_H

#include <lime_defs.h>
#include <lime_header.h>
#include <sys/types.h>
#include <stdio.h>

/** \brief This describes the state of LIME */
typedef struct { 
  int first_record;        /**< Are we to write the first record */
  int last_written;        /**< has the last record been written */
  FILE* fp;                /**< What file are we writing to */
  int header_nextP;        /**< Are we to write a header next or data */
  n_uint64_t bytes_total;       /**< Total no of bytes in this record */
  n_uint64_t bytes_left;        /**< The number of bytes left to write in
                                 the current record */
  n_uint64_t rec_ptr;           /**< Next byte to be written relative to 
			       the start of the record payload.
			       ranges 0 to bytes_total - 1 */
  n_uint64_t rec_start;         /**< File pointer at start of record payload */
  size_t bytes_pad;        /**< The number of bytes to pad the record */
  int isLastP;             /**< Is this the last record in the message? */
} LimeWriter;

/** \brief Initialise the lime state. Currently with an open file
    -- this may change
    \param fp the FILE to which we want to write 
    \returns a LimeWriter* which you need for further
    transactions
*/
LimeWriter* limeCreateWriter(FILE *fp);

/** \brief Close a lime generator */
int limeDestroyWriter(LimeWriter *s);

/** Write a lime record header
 *
 *\param d is the handle of the LIME stream
 *\param props contains data for the record header 
 *             -- note we fill in some fields depending on 
 *                circumstances
 *\param lastP  is a boolean indicating whether this is to be the 
 *              last record
 */
int limeWriteRecordHeader( LimeRecordHeader *props, LimeWriter* d);

/** Write a lime record
 *
 *\param source is the buffer or writing
 *\param nbytes number of bytes to write. Is modified with number 
 *              of bytes actually written
 *\param d  is the LimeWriter/generator to use
 */
int limeWriteRecordData( void *source, n_uint64_t *nbytes,  LimeWriter* d);

/** \brief Advance to end of current record
 *  \params w points to a LimeWriter
 */

int limeWriterCloseRecord(LimeWriter *w);

/** \brief Seek bytes within current open record 
 *  \params r is a pointer to a LimeWriter
 *  \params offset counts bytes from the current position ("WHENCE")
 *    or relative to the start of the binary data.
 *  \params whence = SEEK_CUR or SEEK_SET
 *  \returns a status code
 */
int limeWriterSeek(LimeWriter *r, off_t offset, int whence);

/** \brief Sets the LIME writer state to a prescribed state
 *         and positions the file accordingly.
 *  \params wdest our private LIME writer
 *  \params wsrc LIME writer whose state we want to duplicate.
 *  \returns a status code
 */
int limeWriterSetState(LimeWriter *wdest, LimeWriter *wsrc);

#endif
