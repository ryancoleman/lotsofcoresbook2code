/*
 * MyFile.h
 *
 *  Created on: Jan 25, 2012
 *      Author: yongchao
 */

#ifndef MYFILE_H_
#define MYFILE_H_
#include "Macros.h"

#ifdef COMPRESSED_INPUT
typedef gzFile MyFilePt;

/*open the file*/
#define myfopen(fileName, mode)	gzopen(fileName, mode)

/*open the stdin*/
#define myopenstdin(mode)	gzopen("stdin", mode)

/*close the file*/
#define myfclose(file) gzclose(file)

/*check end of file*/
#define myfeof(file) gzeof(file)

/*read data from file*/
#define myfread(buffer, size, nmemb, file) gzread(file, buffer, (size) * (nmemb))

#else	/*COMPRESSED_INPUT*/
typedef FILE* MyFilePt;

/*open the file*/
#define myfopen(fileName, mode) fopen(fileName, mode)

/*open the stdin*/
#define myopenstdin(mode) 		stdin

/*close the file*/
#define myfclose(file) fclose(file)

/*check the end-of-file*/
#define myfeof(file) feof(file)

/*read data from file*/
#define myfread(buffer, size, nmemb, file) fread(buffer, size, nmemb, file)

#endif	/*COMPRESSED_INPUT*/

#endif /* MYFILE_H_ */
