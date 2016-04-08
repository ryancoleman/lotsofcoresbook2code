#ifndef QIOXML_ILDG_RECORD_H
#define QIOXML_ILDG_RECORD_H

#define QIO_ILDGFORMATVERSION "1.0"

/*******************************************************************/
/* Top level wrapper for ILDG Lattice XML

   tag           member           description          
   ------------------------------------------------------------
   ildgFormat  ildgformat_tags  string of ILDG format tags (see below)

*/

typedef struct {
  QIO_TagCharValue     ildgformatinfo_tags;
} QIO_ILDGFormatInfoWrapper;

#define QIO_ILDGFORMATSCHEMA "xmlns=\"http://www.lqcd.org/ildg\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.lqcd.org/ildg/filefmt.xsd\""

#define QIO_ILDG_FORMAT_INFO_WRAPPER {\
  {"ildgFormat", QIO_ILDGFORMATSCHEMA, "" , 0}       \
}


/*******************************************************************/
/* Contents of ILDG format XML

   tag        member       description                  e.g. gauge config
   ------------------------------------------------------------
   version    version      ILDG format version number      1.0
   precision  precision	   32 or 64
   field      field        "su3_gauge"
   lx         lx           x dimension
   ly         ly           y dimension
   lz         lz           z dimension
   lt         lt           t dimension
*/

typedef struct {
  QIO_TagCharValue version    ;
  QIO_TagCharValue field      ;
  QIO_TagIntValue  precision  ;
  QIO_TagIntValue  lx         ;
  QIO_TagIntValue  ly         ;
  QIO_TagIntValue  lz         ;
  QIO_TagIntValue  lt         ;
} QIO_ILDGFormatInfo;

#define QIO_ILDG_FORMAT_INFO_TEMPLATE { \
  {"version",   "", "", 0},         \
  {"field",     "", "", 0},         \
  {"precision", "", 0 , 0},         \
  {"lx",        "", 0 , 0},         \
  {"ly",        "", 0 , 0},         \
  {"lz",        "", 0 , 0},         \
  {"lt",        "", 0 , 0}          \
}

int QIO_decode_ILDG_format_info(QIO_ILDGFormatInfo *ildg_info, 
				QIO_String *ildg_string);
void QIO_encode_ILDG_format_info(QIO_String *ildg_string, 
				 QIO_ILDGFormatInfo *ildg_info);
int QIO_insert_ildgformat_tag_string(QIO_ILDGFormatInfoWrapper *wrapper, 
				     char *ildgformatinfo_tags);
int QIO_insert_ildgformat_version(QIO_ILDGFormatInfo *ildg_info, 
				  char *version);
int QIO_insert_ildgformat_field(QIO_ILDGFormatInfo *ildg_info, 
				char *field_string);
int QIO_insert_ildgformat_precision(QIO_ILDGFormatInfo *ildg_info, 
				    int precision);
int QIO_insert_ildgformat_lx(QIO_ILDGFormatInfo *ildg_info, int lx);
int QIO_insert_ildgformat_ly(QIO_ILDGFormatInfo *ildg_info, int ly);
int QIO_insert_ildgformat_lz(QIO_ILDGFormatInfo *ildg_info, int lz);
int QIO_insert_ildgformat_lt(QIO_ILDGFormatInfo *ildg_info, int lt);

char *QIO_get_ildgformat_info_tag_string(QIO_ILDGFormatInfoWrapper *wrapper);
char *QIO_get_ildgformat_field(QIO_ILDGFormatInfo *ildg_info);
int QIO_get_ildgformat_precision(QIO_ILDGFormatInfo *ildg_info);
int QIO_get_ildgformat_lx(QIO_ILDGFormatInfo *ildg_info);
int QIO_get_ildgformat_ly(QIO_ILDGFormatInfo *ildg_info);
int QIO_get_ildgformat_lz(QIO_ILDGFormatInfo *ildg_info);
int QIO_get_ildgformat_lt(QIO_ILDGFormatInfo *ildg_info);

QIO_ILDGFormatInfo *QIO_create_ildg_format_info(int precision, int *dims);
void QIO_destroy_ildg_format_info(QIO_ILDGFormatInfo *ildg_info);

#endif /* QIOXML_ILDG_RECORD_H */
