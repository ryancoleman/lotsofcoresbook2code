#ifndef QIOXML_PRIVATE_H
#define QIOXML_PRIVATE_H

#define QIO_FILEFORMATVERSION "1.1"
#define QIO_RECORDFORMATVERSION "1.1"
#define QIO_CHECKSUMFORMATVERSION "1.0"

#ifdef __cplusplus
extern "C"
{
#endif

/*******************************************************************/
/* Top level wrapper for private record XML

   tag           member           description          
   ------------------------------------------------------------
   scidacRecord  recordinfo_tags  string of private record tags (see below)

*/

typedef struct {
  QIO_TagCharValue     recordinfo_tags;
} QIO_RecordInfoWrapper;

#define QIO_RECORD_INFO_WRAPPER {\
  {"scidacRecord", "", "" , 0}       \
}


/*******************************************************************/
/* Contents of private record XML

   tag        member       description                  e.g. gauge config
   ------------------------------------------------------------
   version    version      record format version number    1.0
   date       date         creation date in UT     Wed Oct 22 14:58:08 UTC 2003
   recordtype recordtype   0 if field 1 if global 2 if hypercube subset 0
   spacetime  spacetime    lattice dimensions              4
   hyperlower hyperlower   lower bounds of hypercube       0 0 0 5
   hyperupper hyperupper   upper bounds of hypercube       32 32 32 5
   datatype   datatype     QLA type                        ColorMatrix 
   precision  precision	   I, F, D, or S (random no state) F
   colors     colors       number of colors	           3
   spins      spins        number of spins	 	   --
   typesize   typesize     byte length of datum	           72           
   datacount  datacount    number of data per site	   4

*/

typedef struct {
  QIO_TagCharValue     version    ;
  QIO_TagCharValue     date       ;
  QIO_TagIntValue      recordtype ;
  QIO_TagIntValue      spacetime  ;
  QIO_TagIntArrayValue hyperlower ;     
  QIO_TagIntArrayValue hyperupper ;     
  QIO_TagCharValue     datatype   ;
  QIO_TagCharValue     precision  ;
  QIO_TagIntValue      colors     ;
  QIO_TagIntValue      spins      ;
  QIO_TagIntValue      typesize   ;
  QIO_TagIntValue      datacount  ;
} QIO_RecordInfo;

#define QIO_RECORD_INFO_TEMPLATE {  \
  {"version",   "", "", 0},         \
  {"date",      "", "", 0},         \
  {"recordtype","", 0 , 0},         \
  {"spacetime", "", 0 , 0},         \
  {"hyperlower","", {0} , 1, 0},    \
  {"hyperupper","", {0} , 1, 0},    \
  {"datatype",  "", "", 0},         \
  {"precision", "", "", 0},         \
  {"colors",    "", 0 , 0},         \
  {"spins",     "", 0 , 0},         \
  {"typesize",  "", 0 , 0},         \
  {"datacount", "", 0 , 0}          \
}

/* Obsolete version 1.0 format */
typedef struct {
  QIO_TagCharValue version    ;
  QIO_TagCharValue date       ;
  QIO_TagIntValue  globaldata ;
  QIO_TagCharValue datatype   ;
  QIO_TagCharValue precision  ;
  QIO_TagIntValue  colors     ;
  QIO_TagIntValue  spins      ;
  QIO_TagIntValue  typesize   ;
  QIO_TagIntValue  datacount  ;
} QIO_RecordInfo_v1p0;

#define QIO_RECORD_INFO_TEMPLATE_v1p0 { \
  {"version",   "", "", 0},         \
  {"date",      "", "", 0},         \
  {"globaldata","", 0 , 0},         \
  {"datatype",  "", "", 0},         \
  {"precision", "", "", 0},         \
  {"colors",    "", 0 , 0},         \
  {"spins",     "", 0 , 0},         \
  {"typesize",  "", 0 , 0},         \
  {"datacount", "", 0 , 0}          \
}

/*******************************************************************/
/* Top level wrapper for private file XML

   tag           member           description          
   ------------------------------------------------------------
   scidacFile  fileinfo_tags  string of private file tags (see below)

*/

typedef struct {
  QIO_TagCharValue     fileinfo_tags;
} QIO_FileInfoWrapper;

#define QIO_FILE_INFO_WRAPPER {\
  {"scidacFile", "", "" , 0}       \
}


/*******************************************************************/
/* Contents of private file XML

   tag        member       description                  e.g. gauge config
   ------------------------------------------------------------
   version    version      file format version number      1.0
   spacetime  spacetime    dimensions of space plus time   4
   dims       dims         lattice dimensions              20 20 20 64
   volfmt     volfmt       QIO_SINGLEFILE, QIO_PARTFILE, 
                           QIO_MULTIFILE                   0
*/

typedef struct {
  QIO_TagCharValue      version;
  QIO_TagIntValue       spacetime;
  QIO_TagIntArrayValue  dims;
  QIO_TagIntValue       volfmt;
} QIO_FileInfo;

#define QIO_FILE_INFO_TEMPLATE  {\
  {"version",  "",  "", 0},       \
  {"spacetime","",  0,  0},       \
  {"dims",     "",  {0}, 1, 0},   \
  {"volfmt",   "", 0 , 0}       \
}

/* Obsolete version 1.0 format */

typedef struct {
  QIO_TagCharValue      version;
  QIO_TagIntValue       spacetime;
  QIO_TagIntArrayValue  dims;
  QIO_TagIntValue       multifile;
} QIO_FileInfo_v1p0;


#define QIO_FILE_INFO_TEMPLATE_v1p0  {\
  {"version",   "", "", 0},       \
  {"spacetime", "", 0,  0},       \
  {"dims",      "", {0}, 1, 0},   \
  {"multifile", "",  0 , 0}       \
}


/*******************************************************************/
/* Top level wrapper for private checksum XML

   tag           member           description          
   ------------------------------------------------------------
   scidacChecksum  checksuminfo_tags  string of checksum tags (see below)

*/

typedef struct {
  QIO_TagCharValue     checksuminfo_tags;
} QIO_ChecksumInfoWrapper;

#define QIO_CHECKSUM_INFO_WRAPPER {\
  {"scidacChecksum", "", "" , 0}       \
}


/*******************************************************************/
/* Contents of record checksum XML

   tag        member       description
   ------------------------------------------------------------
   version    version      checksum version number      1.0
   suma       suma         
   sumb       sumb         
*/

typedef struct {
  QIO_TagCharValue      version;
  QIO_TagHex32Value    suma;
  QIO_TagHex32Value    sumb;
} QIO_ChecksumInfo;

#define QIO_CHECKSUM_INFO_TEMPLATE  {\
  {"version",   "", "", 0},       \
  {"suma",      "", 0,  0},       \
  {"sumb",      "", 0,  0}        \
}

/*********************************************************************/
int QIO_decode_file_info(QIO_FileInfo *file_info, 
			 QIO_String *file_string);
void QIO_encode_file_info(QIO_String *file_string, 
			   QIO_FileInfo *file_info);
int QIO_decode_record_info(QIO_RecordInfo *record_info, 
			    QIO_String *record_string);
void QIO_encode_record_info(QIO_String *record_string, 
			      QIO_RecordInfo *record_info);
int QIO_decode_checksum_info(QIO_ChecksumInfo *checksum, 
			     QIO_String *file_string);
void QIO_encode_checksum_info(QIO_String *file_string, 
			      QIO_ChecksumInfo *checksum);


int QIO_insert_file_tag_string(QIO_FileInfoWrapper *wrapper, 
			       char *fileinfo_tags);
int QIO_insert_spacetime_dims(QIO_FileInfo *file_info, 
			      int spacetime, int *dims);
int QIO_insert_volfmt(QIO_FileInfo *file_info, int volfmt);

int QIO_insert_record_tag_string(QIO_RecordInfoWrapper *wrapper, 
				 char *recordinfo_tags);
int QIO_insert_record_date(QIO_RecordInfo *record_info, char* date);
int QIO_insert_recordtype(QIO_RecordInfo *record_info, int recordtype);
int QIO_insert_hypercube_bounds(QIO_RecordInfo *record_info, 
				int *lower, int *upper, int n);
int QIO_insert_datatype(QIO_RecordInfo *record_info, char* datatype);
int QIO_insert_precision(QIO_RecordInfo *record_info, char* precision);
int QIO_insert_colors(QIO_RecordInfo *record_info, int colors);
int QIO_insert_spins(QIO_RecordInfo *record_info, int spins);
int QIO_insert_typesize(QIO_RecordInfo *record_info, int typesize);
int QIO_insert_datacount(QIO_RecordInfo *record_info, int datacount);

int QIO_insert_checksum_tag_string(QIO_ChecksumInfoWrapper *wrapper, 
				   char *checksuminfo_tags);
int QIO_insert_suma_sumb(QIO_ChecksumInfo *checksum_info, 
			 uint32_t suma, uint32_t sumb);

char *QIO_get_file_info_tag_string(QIO_FileInfoWrapper *wrapper);
char *QIO_get_file_version(QIO_FileInfo *file_info);
int QIO_get_spacetime(QIO_FileInfo *file_info);
int *QIO_get_dims(QIO_FileInfo *file_info);
int QIO_get_volfmt(QIO_FileInfo *file_info);
int QIO_get_multifile(QIO_FileInfo_v1p0 *file_info);
int QIO_defined_spacetime(QIO_FileInfo *file_info);
int QIO_defined_dims(QIO_FileInfo *file_info);
int QIO_defined_volfmt(QIO_FileInfo *file_info);

char *QIO_get_record_info_tag_string(QIO_RecordInfoWrapper *wrapper);
char *QIO_get_record_info_version(QIO_RecordInfo *record_info);
char *QIO_get_record_date(QIO_RecordInfo *record_info);
int QIO_get_globaldata(QIO_RecordInfo_v1p0 *record_info);
int QIO_get_recordtype(QIO_RecordInfo *record_info);
int QIO_get_hyper_spacetime(QIO_RecordInfo *record_info);
int *QIO_get_hyperlower(QIO_RecordInfo *record_info);
int *QIO_get_hyperupper(QIO_RecordInfo *record_info);
char *QIO_get_datatype(QIO_RecordInfo *record_info);
char *QIO_get_precision(QIO_RecordInfo *record_info);
int QIO_get_colors(QIO_RecordInfo *record_info);
int QIO_get_spins(QIO_RecordInfo *record_info);
int QIO_get_typesize(QIO_RecordInfo *record_info);
int QIO_get_datacount(QIO_RecordInfo *record_info);

void QIO_set_recordtype(QIO_RecordInfo *record_info, int recordtype);
void QIO_set_datatype(QIO_RecordInfo *record_info, char *datatype);
void QIO_set_precision(QIO_RecordInfo *record_info, char *precision);
void QIO_set_record_date(QIO_RecordInfo *record_info, char *date);
void QIO_set_colors(QIO_RecordInfo *record_info, int colors);
void QIO_set_spins(QIO_RecordInfo *record_info, int spins);
void QIO_set_typesize(QIO_RecordInfo *record_info, int typesize);
void QIO_set_datacount(QIO_RecordInfo *record_info, int datacount);

char *QIO_get_checksum_info_tag_string(QIO_ChecksumInfoWrapper *wrapper);
uint32_t QIO_get_suma(QIO_ChecksumInfo *checksum_info);
uint32_t QIO_get_sumb(QIO_ChecksumInfo *checksum_info);
int QIO_defined_suma(QIO_ChecksumInfo *checksum_info);
int QIO_defined_sumb(QIO_ChecksumInfo *checksum_info);

int QIO_defined_record_date(QIO_RecordInfo *record_info);
int QIO_defined_recordtype(QIO_RecordInfo *record_info);
int QIO_defined_datatype(QIO_RecordInfo *record_info);
int QIO_defined_precision(QIO_RecordInfo *record_info);
int QIO_defined_colors(QIO_RecordInfo *record_info);
int QIO_defined_spins(QIO_RecordInfo *record_info);
int QIO_defined_typesize(QIO_RecordInfo *record_info);
int QIO_defined_datacount(QIO_RecordInfo *record_info);

QIO_FileInfo *QIO_create_file_info(int spacetime, int *dims, int volfmt);
void QIO_destroy_file_info(QIO_FileInfo *file_info);
int QIO_compare_file_info(QIO_FileInfo *found, QIO_FileInfo *expect,
			  char *myname, int this_node);
QIO_RecordInfo *QIO_create_record_info(int recordtype, int lower[],
				       int upper[], int n,
				       char *datatype, char *precision, 
				       int colors, int spins, int typesize, 
				       int datacount);
void QIO_destroy_record_info(QIO_RecordInfo *record_info);
int QIO_compare_record_info(QIO_RecordInfo *r, QIO_RecordInfo *s);

QIO_ChecksumInfo *QIO_create_checksum_info(uint32_t suma, uint32_t sumb);
void QIO_destroy_checksum_info(QIO_ChecksumInfo *checksum_info);
int QIO_compare_checksum_info(QIO_ChecksumInfo *found, 
			      QIO_ChecksumInfo *expect, 
			      char *myname, int this_node);

#ifdef __cplusplus
}
#endif

#endif /* QIOXML_PRIVATE_H */
