#ifndef LIME_HEADER_H
#define LIME_HEADER_H

#include <lime_defs.h>
#include <sys/types.h>
#include <lime_fixed_types.h>

/** \brief -- The header structure for a lime record */
typedef struct { 
  unsigned int lime_version;   /**< LIME Version */
  int MB_flag;                 /**< Message begin Flag */
  int ME_flag;                 /**< Message end Flag */
  char *type;                  /**< The LimeType string */
  n_uint64_t data_length;        /**< The length of the data to follow */
} LimeRecordHeader;


/** \brief Creates a partially filled out header. 
    The header does not contain the MB flag, ME flag as this
    is filled out when it comes time to write the record */
LimeRecordHeader *limeCreateHeader(int MB_flag,
				   int ME_flag,
				   char *type, 
				   n_uint64_t reclen);

void limeDestroyHeader(LimeRecordHeader *h);

				     
#endif
