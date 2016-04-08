AM_CFLAGS = -D_FILE_OFFSET_BITS=64 @FILEHASH_DEBUG@ -I@top_srcdir@/filehash/include -I@top_builddir@/filehash/include

lib_LIBRARIES = libfilehash.a

libfilehash_a_SOURCES = \
	ffdb_cq.h \
	ffdb_db.h \
	ffdb_db.c \
	ffdb_hash.h \
	ffdb_hash.c \
	ffdb_hash_func.h \
	ffdb_hash_func.c \
	ffdb_page.h \
	ffdb_page.c \
	ffdb_pagepool.h \
	ffdb_pagepool.c

include_HEADERS = \
	ffdb_db.h

check_PROGRAMS = tcreate twrite tdualwrite treplace tread tread-th tkeyseq tcheck createStrings replaceStrings

LDADD = libfilehash.a -lpthread
clean-local:
	-rm -rf *~ core*
