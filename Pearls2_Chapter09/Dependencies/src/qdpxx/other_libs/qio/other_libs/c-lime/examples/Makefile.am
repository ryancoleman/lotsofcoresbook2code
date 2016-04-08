
AM_CFLAGS  = -I@top_srcdir@/include -I@top_builddir@/include
AM_LDFLAGS = -L@top_builddir@/lib 
LDADD = @LIBS@ -llime

#
# The programs to build
#
bin_PROGRAMS = lime_contents lime_extract_record lime_extract_type lime_pack lime_unpack
check_PROGRAMS = lime_writer_test1

#
# The program and its dependencies
#
lime_contents_SOURCES = lime_contents.c
lime_writer_test1_SOURCES = lime_writer_test1.c
lime_extract_record_SOURCES = lime_extract_record.c
lime_extract_type_SOURCES = lime_extract_type.c
lime_pack_SOURCES = lime_pack.c
lime_unpack_SOURCES = lime_unpack.c


.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
