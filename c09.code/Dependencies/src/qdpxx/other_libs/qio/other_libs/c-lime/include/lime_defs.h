#ifndef LIME_DEFS_H
#define LIME_DEFS_H

/**< The magic number */
#define LIME_MAGIC_NO     1164413355

/**< The lime version we are implementing (v.1) */
#define LIME_VERSION              ((unsigned int)0x01)

/**< The maximum internal buffer size I can allocate (this is 2Gb). This
     value appears to be unused in the code so I'll leave it at this enormous
     amount (BJ: 28/02/05) */
#define LIME_MAX_BUFSIZE          ((unsigned int)2147483647)

#define LIME_SUCCESS              ((int) 0) /**< Success status code */
#define LIME_ERR_LAST_NOT_WRITTEN ((int)-1) 
#define LIME_ERR_PARAM            ((int)-2) /**< Something with the inputs is 
					       wrong */
#define LIME_ERR_HEADER_NEXT      ((int)-3) 
#define LIME_LAST_REC_WRITTEN     ((int)-4)
#define LIME_ERR_WRITE            ((int)-5) /**< Some sort of write error occurred */
#define LIME_EOR                  ((int)-6) /**< End of Record (EOR) */
#define LIME_EOF                  ((int)-7) /**< End of File (EOF) */
#define LIME_ERR_READ             ((int)-8) /**< A read error has occurred */
#define LIME_ERR_SEEK             ((int)-9) /**< Error seeking */
#define LIME_ERR_MBME             ((int)-10) /**< MB/ME flags incorrect */
#define LIME_ERR_CLOSE            ((int)-11) /**< Error closing file */

#endif
