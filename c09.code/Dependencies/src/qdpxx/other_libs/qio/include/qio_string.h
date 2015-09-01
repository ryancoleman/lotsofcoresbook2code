/* String support for a safe string type */

#ifndef QIO_STRING_H
#define QIO_STRING_H

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct {
  char *string;
  size_t length;
} QIO_String;

QIO_String *QIO_string_create(void);
void QIO_string_destroy(QIO_String *qs);
void QIO_string_set(QIO_String *qs, const char *const string);
size_t QIO_string_length(const QIO_String *const qs);
char * QIO_string_ptr(QIO_String *qs);
void QIO_string_copy(QIO_String *dest, QIO_String *src);
void QIO_string_realloc(QIO_String *qs, int length);
void QIO_string_append(QIO_String *qs, const char *const string);

#ifdef __cplusplus
}
#endif

#endif /* QIO_STRING_H */
