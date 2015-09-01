#ifndef QIOXML_USQCD_PROP_H
#define QIOXML_USQCD_PROP_H

#ifdef __cplusplus
extern "C"
{
#endif

/*********************************************************************/
/* Top level wrapper for USQCD propagator file XML 

   tag           member           description          
   ------------------------------------------------------------
 
*/

typedef struct {
  QIO_TagCharValue usqcdpropfileinfo_tags;
} QIO_USQCDPropFileInfoWrapper;

#define QIO_USQCD_PROPFILE_INFO_WRAPPER {\
  {"usqcdPropFile", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD propfile XML

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
} QIO_USQCDPropFileInfo;


#define QIO_USQCDPROPFILEFORMATVERSION "1.0"

#define QIO_USQCD_PROPFILE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"type"   , "", "", 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/
/* Top level wrapper for USQCD propagator source record XML 

   tag           member           description          
   ------------------------------------------------------------
   
*/

typedef struct {
  QIO_TagCharValue usqcdpropsourceinfo_tags;
} QIO_USQCDPropSourceInfoWrapper;

#define QIO_USQCD_PROPSOURCE_INFO_WRAPPER {\
  {"usqcdSourceInfo", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD propsource record XML

   tag           member           description          
   ------------------------------------------------------------
   version       version        propsource record version number   1.0
   spin          spin           
   color         color          
   info          info           XML string     collaboration option
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagIntValue spin;
  QIO_TagIntValue color;
  QIO_TagCharValue info;
} QIO_USQCDPropSourceInfo;


#define QIO_USQCDPROPSOURCEFORMATVERSION "1.0"

#define QIO_USQCD_PROPSOURCE_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"spin"   , "", 0, 0}, \
   {"color" ,  "", 0, 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/
/* Top level wrapper for USQCD prop record XML 

   tag           member           description          
   ------------------------------------------------------------
   
*/

typedef struct {
  QIO_TagCharValue usqcdproprecordinfo_tags;
} QIO_USQCDPropRecordInfoWrapper;

#define QIO_USQCD_PROPRECORD_INFO_WRAPPER {\
  {"usqcdPropInfo", "", "" , 0}       \
}

  /* Backward compatibility feature */
#define QIO_USQCD_PROPRECORD_INFO_WRAPPER_LEGACY {\
  {"usqcdInfo", "", "" , 0}       \
}

/*******************************************************************/
/* Contents of USQCD prop record XML

   tag           member           description          
   ------------------------------------------------------------
   version       version        prop record version number   1.0
   spin          spin           
   color         color          
   info          info           XML string     collaboration option
*/


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagIntValue spin;
  QIO_TagIntValue color;
  QIO_TagCharValue info;
} QIO_USQCDPropRecordInfo;


#define QIO_USQCDPROPRECORDFORMATVERSION "1.0"

#define QIO_USQCD_PROPRECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"spin"   , "", 0, 0}, \
   {"color" ,  "", 0, 0}, \
   {"info"   , "", "", 0}  \
}

/*********************************************************************/

void QIO_encode_usqcd_propfile_info(QIO_String *file_string, 
				     QIO_USQCDPropFileInfo *file_info);
int QIO_decode_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info,
				    QIO_String *file_string);

void QIO_encode_usqcd_propsource_info(QIO_String *record_string, 
				     QIO_USQCDPropSourceInfo *record_info);
int QIO_decode_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info,
				    QIO_String *record_string);

void QIO_encode_usqcd_proprecord_info(QIO_String *record_string, 
				     QIO_USQCDPropRecordInfo *record_info);
int QIO_decode_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info,
				    QIO_String *record_string);

/* Propagator file types */

#define QIO_USQCDPROPFILETYPESTRING_C1D12  \
 "USQCD_DiracFermion_ScalarSource_TwelveSink"
#define QIO_USQCDPROPFILETYPESTRING_DD_PAIRS \
 "USQCD_DiracFermion_Source_Sink_Pairs"
#define QIO_USQCDPROPFILETYPESTRING_CD_PAIRS \
 "USQCD_DiracFermion_ScalarSource_Sink_Pairs"
#define QIO_USQCDPROPFILETYPESTRING_LHPC \
 "LHPC_DiracPropagator"

#define QIO_USQCDPROPFILETYPE_C1D12      0
#define QIO_USQCDPROPFILETYPE_DD_PAIRS   1
#define QIO_USQCDPROPFILETYPE_CD_PAIRS   2
#define QIO_USQCDPROPFILETYPE_LHPC       3


int QIO_insert_usqcd_propfile_version(QIO_USQCDPropFileInfo *file_info, char *version);
int QIO_insert_usqcd_propfile_type( QIO_USQCDPropFileInfo *file_info, int type);
int QIO_insert_usqcd_propfile_info( QIO_USQCDPropFileInfo *file_info, char *info);
int QIO_insert_usqcd_propfile_tag_string(QIO_USQCDPropFileInfoWrapper *wrapper,
					char *fileinfo_tags);

int QIO_insert_usqcd_propsource_version(QIO_USQCDPropSourceInfo *record_info, char *version);
int QIO_insert_usqcd_propsource_info( QIO_USQCDPropSourceInfo *record_info, char *info);
int QIO_insert_usqcd_propsource_tag_string(QIO_USQCDPropSourceInfoWrapper *wrapper,
					   char *recordinfo_tags);
int QIO_insert_usqcd_propsource_spin( QIO_USQCDPropSourceInfo *record_info, int spin);
int QIO_insert_usqcd_propsource_color( QIO_USQCDPropSourceInfo *record_info, int color);


int QIO_insert_usqcd_proprecord_version(QIO_USQCDPropRecordInfo *record_info, char *version);
int QIO_insert_usqcd_proprecord_info( QIO_USQCDPropRecordInfo *record_info, char *info);
int QIO_insert_usqcd_proprecord_tag_string(QIO_USQCDPropRecordInfoWrapper *wrapper,
					char *recordinfo_tags);
/* added EES (these two functions were present in QIO_info_usqcd_prop.c but not in this header file) */
int QIO_insert_usqcd_proprecord_spin( QIO_USQCDPropRecordInfo *record_info, int spin);
int QIO_insert_usqcd_proprecord_color( QIO_USQCDPropRecordInfo *record_info, int color);

char *QIO_get_usqcd_propfile_info_tag_string(QIO_USQCDPropFileInfoWrapper *wrapper);
int QIO_get_usqcd_propfile_type(QIO_USQCDPropFileInfo *file_info);
char *QIO_get_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info);
int QIO_defined_usqcd_propfile_type(QIO_USQCDPropFileInfo *file_info);
int QIO_defined_usqcd_propfileinfo(QIO_USQCDPropFileInfo *file_info);

char *QIO_get_usqcd_propsource_info_tag_string(QIO_USQCDPropSourceInfoWrapper *wrapper);
char *QIO_get_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info);
int QIO_get_usqcd_propsource_spin(QIO_USQCDPropSourceInfo *record_info);
int QIO_get_usqcd_propsource_color(QIO_USQCDPropSourceInfo *record_info);
int QIO_defined_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info);
int QIO_defined_usqcd_propsource_spin(QIO_USQCDPropSourceInfo *record_info);
int QIO_defined_usqcd_propsource_color(QIO_USQCDPropSourceInfo *record_info);

char *QIO_get_usqcd_proprecord_info_tag_string(QIO_USQCDPropRecordInfoWrapper *wrapper);
char *QIO_get_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info);
int QIO_get_usqcd_proprecord_spin(QIO_USQCDPropRecordInfo *record_info);
int QIO_get_usqcd_proprecord_color(QIO_USQCDPropRecordInfo *record_info);
int QIO_defined_usqcd_proprecordinfo(QIO_USQCDPropRecordInfo *record_info);
/* added EES (two new functions in QIO_info_usqcd_prop.c)*/
int QIO_defined_usqcd_proprecord_spin(QIO_USQCDPropRecordInfo *record_info);
int QIO_defined_usqcd_proprecord_color(QIO_USQCDPropRecordInfo *record_info);

QIO_USQCDPropFileInfo *QIO_create_usqcd_propfile_info(int type, char *info);
void QIO_destroy_usqcd_propfile_info(QIO_USQCDPropFileInfo *file_info);

QIO_USQCDPropSourceInfo *QIO_create_usqcd_propsource_info(char *info);
QIO_USQCDPropSourceInfo *QIO_create_usqcd_propsource_sc_info(int spin, 
	     int color, char *info);
void QIO_destroy_usqcd_propsource_info(QIO_USQCDPropSourceInfo *record_info);

QIO_USQCDPropRecordInfo *QIO_create_usqcd_proprecord_info(char *info);
QIO_USQCDPropRecordInfo *QIO_create_usqcd_proprecord_sc_info(int spin, 
	     int color, char *info);
void QIO_destroy_usqcd_proprecord_info(QIO_USQCDPropRecordInfo *record_info);


#ifdef __cplusplus
}
#endif

#endif /* QIOXML_USQCD_PROP_H */
