#ifndef QIOXML_INTERNAL_H
#define QIOXML_INTERNAL_H

#define QIO_MAXTAG 64
#define QIO_MAXATTR 512
#define QIO_MAXVALUESTRING 32768
#define QIO_MAXINTARRAY 8
#define QIO_STRINGALLOC 1024

#ifdef __cplusplus
extern "C"
{
#endif

/*******************************************************************/
/* Datatype structures */

typedef struct {
  char tag[QIO_MAXTAG];
  char attr[QIO_MAXATTR];
  char value[QIO_MAXVALUESTRING];
  short occur;
} QIO_TagCharValue;

typedef struct {
  char tag[QIO_MAXTAG];
  char attr[QIO_MAXATTR];
  uint32_t value;     /* Must be 32-bit type */
  short occur;
} QIO_TagHex32Value;

typedef struct {
  char tag[QIO_MAXTAG];
  char attr[QIO_MAXATTR];
  int  value;
  short occur;
} QIO_TagIntValue;

typedef struct {
  char tag[QIO_MAXTAG];
  char attr[QIO_MAXATTR];
  int  value[QIO_MAXINTARRAY];
  int  n;
  short occur;
} QIO_TagIntArrayValue;


/*********************************************************************/
/* Internal utilities */

char *QIO_next_token(char *parse_pt, char *tokens);
char *QIO_next_nontoken(char *parse_pt, char *tokens);
char *QIO_strncat(char *s1, char *s2, int *n);
char *QIO_next_tag(char *parse_pt, char *tag, char **left_angle);
char *QIO_get_tag_value(char *parse_pt, char *tag, char *value_string);
void QIO_decode_as_string(char *tag, char *value_string, 
			  QIO_TagCharValue *tag_value);
void QIO_decode_as_int(char *tag, char *value_string, 
		       QIO_TagIntValue *tag_value);
void QIO_decode_as_intlist(char *tag, char *value_string, 
			   QIO_TagIntArrayValue *tag_value);
void QIO_decode_as_hex32(char *tag, char *value_string, 
			 QIO_TagHex32Value *tag_value);
char *QIO_write_tag(char *buf, char *tag, char *attr, int *remainder);
char *QIO_write_endtag(char *buf,char *tag,int *remainder);
char *QIO_encode_as_string(char *buf, QIO_TagCharValue *tag_value,
			   int *remainder);
char *QIO_encode_as_int(char *buf, QIO_TagIntValue *tag_value,
			int *remainder);
char *QIO_encode_as_hex32(char *buf, QIO_TagHex32Value *tag_value,
			  int *remainder);
char *QIO_encode_as_intlist(char *buf, 
			      QIO_TagIntArrayValue *tag_value, 
			    int n, int *remainder);
int QIO_check_string_occur(QIO_TagCharValue *tag_value);
int QIO_check_int_occur(QIO_TagIntValue *tag_value);
int QIO_check_intarray_occur(QIO_TagIntArrayValue *tag_value);
int QIO_check_hex32_occur(QIO_TagHex32Value *tag_value);

#ifdef __cplusplus
}
#endif

#endif /* QIOXML_INTERNAL_H */
