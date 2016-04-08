#ifndef QIOXML_USQCD_KSPROP_H
#define QIOXML_USQCD_KSPROP_H

#ifdef __cplusplus
extern "C"
{
#endif

/*********************************************************************/
/* Top level wrapper for USQCD KS propagator file XML 

   tag           member           description          
   ------------------------------------------------------------
 
*/

typedef struct {
  QIO_TagCharValue usqcdkspropfileinfo_tags;
} QIO_USQCDKSPropFileInfoWrapper;

#define QIO_USQCD_KSPROPFILE_INFO_WRAPPER {\
  {"usqcdKSPropFile", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD kspropfile XML

   tag           member           description          
   ------------------------------------------------------------
   version       version          file XML version
   type          type             file format type string
   info          info             collaboration discretion
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagCharValue type;
  QIO_TagCharValue info;
} QIO_USQCDKSPropFileInfo;


#define QIO_USQCDKSPROPFILEFORMATVERSION "1.0"

#define QIO_USQCD_KSPROPFILE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"type"   , "", "", 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/
/* Top level wrapper for USQCD KS propagator source record XML 

   tag           member           description          
   ------------------------------------------------------------
   
*/

typedef struct {
  QIO_TagCharValue usqcdkspropsourceinfo_tags;
} QIO_USQCDKSPropSourceInfoWrapper;

#define QIO_USQCD_KSPROPSOURCE_INFO_WRAPPER {\
  {"usqcdSourceInfo", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD kspropsource record XML

   tag           member           description          
   ------------------------------------------------------------
   version       version        kspropsource record version number   1.0
   color         color          
   info          info           XML string     collaboration option
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagIntValue color;
  QIO_TagCharValue info;
} QIO_USQCDKSPropSourceInfo;


#define QIO_USQCDKSPROPSOURCEFORMATVERSION "1.0"

#define QIO_USQCD_KSPROPSOURCE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"color" ,  "", 0, 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/
/* Top level wrapper for USQCD ksprop record XML 

   tag           member           description          
   ------------------------------------------------------------
   
*/

typedef struct {
  QIO_TagCharValue usqcdksproprecordinfo_tags;
} QIO_USQCDKSPropRecordInfoWrapper;

#define QIO_USQCD_KSPROPRECORD_INFO_WRAPPER {\
  {"usqcdKSPropInfo", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD ksprop record XML

   tag           member           description          
   ------------------------------------------------------------
   version       version        ksprop record version number   1.0
   color         color          
   info          info           XML string     collaboration option
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagIntValue color;
  QIO_TagCharValue info;
} QIO_USQCDKSPropRecordInfo;


#define QIO_USQCDKSPROPRECORDFORMATVERSION "1.0"

#define QIO_USQCD_KSPROPRECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"color" ,  "", 0, 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/

void QIO_encode_usqcd_kspropfile_info(QIO_String *file_string, 
				     QIO_USQCDKSPropFileInfo *file_info);
int QIO_decode_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info,
				    QIO_String *file_string);

void QIO_encode_usqcd_kspropsource_info(QIO_String *record_string, 
				     QIO_USQCDKSPropSourceInfo *record_info);
int QIO_decode_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *record_info,
				    QIO_String *record_string);

void QIO_encode_usqcd_ksproprecord_info(QIO_String *record_string, 
				     QIO_USQCDKSPropRecordInfo *record_info);
int QIO_decode_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info,
				    QIO_String *record_string);

/* KS Propagator file types */

#define QIO_USQCDKSPROPFILETYPESTRING_C1V3  \
 "USQCD_ColorVector_ScalarSource_ThreeSink"
#define QIO_USQCDKSPROPFILETYPESTRING_VV_PAIRS \
 "USQCD_ColorVector_Source_Sink_Pairs"
#define QIO_USQCDKSPROPFILETYPESTRING_CV_PAIRS \
 "USQCD_ColorVector_ScalarSource_Sink_Pairs"

#define QIO_USQCDKSPROPFILETYPE_C1V3       0
#define QIO_USQCDKSPROPFILETYPE_VV_PAIRS   1
#define QIO_USQCDKSPROPFILETYPE_CV_PAIRS   2

int QIO_insert_usqcdkspropfile_version(QIO_USQCDKSPropFileInfo *file_info, char *version);
int QIO_insert_usqcdkspropfile_type( QIO_USQCDKSPropFileInfo *file_info, int type);
int QIO_insert_usqcdkspropfile_info( QIO_USQCDKSPropFileInfo *file_info, char *info);
int QIO_insert_usqcdkspropfile_tag_string(QIO_USQCDKSPropFileInfoWrapper *wrapper,
					char *fileinfo_tags);

int QIO_insert_usqcdkspropsource_version(QIO_USQCDKSPropSourceInfo *source_info, char *version);
int QIO_insert_usqcdkspropsource_info( QIO_USQCDKSPropSourceInfo *source_info, char *info);
int QIO_insert_usqcdkspropsource_tag_string(QIO_USQCDKSPropSourceInfoWrapper *wrapper,
					char *sourceinfo_tags);
int QIO_insert_usqcd_kspropsource_color( QIO_USQCDKSPropSourceInfo *source_info, int color);

int QIO_insert_usqcd_ksproprecord_version(QIO_USQCDKSPropRecordInfo *record_info, char *version);
int QIO_insert_usqcd_ksproprecord_info( QIO_USQCDKSPropRecordInfo *record_info, char *info);
int QIO_insert_usqcd_ksproprecord_tag_string(QIO_USQCDKSPropRecordInfoWrapper *wrapper,
					char *recordinfo_tags);
int QIO_insert_usqcd_ksproprecord_color( QIO_USQCDKSPropRecordInfo *record_info, int color);

char *QIO_get_usqcd_kspropfile_info_tag_string(QIO_USQCDKSPropFileInfoWrapper *wrapper);
int QIO_get_usqcd_kspropfile_type(QIO_USQCDKSPropFileInfo *file_info);
char *QIO_get_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info);
int QIO_defined_usqcd_kspropfile_type(QIO_USQCDKSPropFileInfo *file_info);
int QIO_defined_usqcd_kspropfileinfo(QIO_USQCDKSPropFileInfo *file_info);

char *QIO_get_usqcd_kspropsource_info_tag_string(QIO_USQCDKSPropSourceInfoWrapper *wrapper);
char *QIO_get_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *source_info);
int QIO_get_usqcd_kspropsource_color(QIO_USQCDKSPropSourceInfo *source_info);
int QIO_defined_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *source_info);
int QIO_defined_usqcd_kspropsource_color(QIO_USQCDKSPropSourceInfo *source_info);

char *QIO_get_usqcd_ksproprecord_info_tag_string(QIO_USQCDKSPropRecordInfoWrapper *wrapper);
char *QIO_get_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info);
int QIO_get_usqcd_ksproprecord_color(QIO_USQCDKSPropRecordInfo *record_info);
int QIO_defined_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info);
int QIO_defined_usqcd_ksproprecord_color(QIO_USQCDKSPropRecordInfo *record_info);

QIO_USQCDKSPropFileInfo *QIO_create_usqcd_kspropfile_info(int type, char *info);
void QIO_destroy_usqcd_kspropfile_info(QIO_USQCDKSPropFileInfo *file_info);

QIO_USQCDKSPropSourceInfo *QIO_create_usqcd_kspropsource_info(char *info);
QIO_USQCDKSPropSourceInfo *QIO_create_usqcd_kspropsource_c_info(
	     int color, char *info);
void QIO_destroy_usqcd_kspropsource_info(QIO_USQCDKSPropSourceInfo *record_info);

QIO_USQCDKSPropRecordInfo *QIO_create_usqcd_ksproprecord_info(char *info);
QIO_USQCDKSPropRecordInfo *QIO_create_usqcd_ksproprecord_c_info(
	     int color, char *info);
void QIO_destroy_usqcd_ksproprecord_info(QIO_USQCDKSPropRecordInfo *record_info);


#ifdef __cplusplus
}
#endif

#endif /* QIOXML_USQCD_KSPROP_H */
