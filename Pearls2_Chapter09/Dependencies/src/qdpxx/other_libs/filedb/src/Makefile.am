AM_CXXFLAGS = @FILEDB_BYTE_ORDER_DEF@ @FILEDB_TEMP_COMP_TEST@ -I@top_srcdir@/filehash

lib_LIBRARIES = libfiledb.a

libfiledb_a_SOURCES = \
	FileDB.cpp \
	FileDB.h \
	Serializable.h \
	DBKey.h \
	DBData.h \
	DBCursor.h \
	DBCursor.cpp \
	DBFunc.cpp \
	DBFunc.h \
	DBString.cpp \
	DBString.h \
	ConfDataStoreDB.cpp \
	ConfDataStoreDB.h \
	ConfDataStoreMultipleDB.h \
	AllConfStoreDB.h \
	AllConfStoreDB.cpp \
	AllConfStoreMultipleDB.h \
	ConfDataDBMerger.h \
	ConfDataDBMerger.cpp \
	ConfigInfo.h \
	ConfigInfo.cpp \
	HashFunc.h \
	HashFunc.cpp \
	VFloatData.h \
	VFloatData.cpp

include_HEADERS = \
	filedb_config.h \
	FileDB.h \
	Serializable.h \
	DBKey.h \
	DBData.h \
	DBCursor.h \
	DBFunc.h \
	DBString.h \
	ConfDataStoreDB.h \
	ConfDataStoreMultipleDB.h \
	AllConfStoreDB.h \
	AllConfStoreMultipleDB.h \
	ConfDataDBMerger.h \
	ConfigInfo.h \
	HashFunc.h

check_PROGRAMS = ConfCreateTest SConfCreateTest SConfReadTest ConfReadTest AllConfCreateTest AllConfReadTest AllConfUpdateTest DBMergeTest

ConfCreateTest_SOURCES=ConfCreateTest.cpp
SConfCreateTest_SOURCES=SConfCreateTest.cpp
SConfReadTest_SOURCES=SConfReadTest.cpp
ConfReadTest_SOURCES=ConfReadTest.cpp
AllConfCreateTest_SOURCES=AllConfCreateTest.cpp
AllConfReadTest_SOURCES=AllConfReadTest.cpp
AllConfUpdateTest_SOURCES=AllConfUpdateTest.cpp
DBMergeTest_SOURCES=DBMergeTest.cpp

LDADD = libfiledb.a ../filehash/libfilehash.a -lpthread

clean-local:
	-rm -rf *~ core*

nodist_include_HEADERS = filedb_config_internal.h
