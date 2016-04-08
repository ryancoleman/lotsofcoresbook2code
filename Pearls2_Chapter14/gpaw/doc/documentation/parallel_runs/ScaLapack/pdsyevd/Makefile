include Makefile.fjorm

COMPILE = $(CC) $(CFLAGS) $(INCLUDEMPI) -c

%.o: %.c
	$(COMPILE) $*.c -o $@

SOURCES = main.c
OBJECTS = $(SOURCES:%.c=%.o)

PDSYEV = test.exe

all: $(PDSYEV)
pdsyev: $(PDSYEV)

LINK = $(F77)
LIBS = $(LIBSCALAPACK) $(LIBBLACS) $(LIBBLAS) $(LIBMPI) $(LDFLAGS)

test.exe: main.o
	$(LINK) main.o $(LIBS) -o $@

clean:
	rm -f $(PDSYEV) main.o
