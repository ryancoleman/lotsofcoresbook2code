# NIM Makefile

include macros.make

SRCDIR      = $(CURDIR)

all: info build

info:
	$(info Inside top level Makefile)
	$(info ARCH   = $(ARCH))
	$(info OPENMP = $(OPENMP))
	$(info CHUNK  = $(CHUNK))
	$(info NZ     = $(NZ))
	$(info TRANSLATEI = $(TRANSLATEI))
	$(info FPMP   = $(FPMP))
	$(info SRCDIR = $(SRCDIR))

build:
	$(MAKE) $(GMAKEMINUSJ) SRCDIR=$(SRCDIR) -C kernel

clean:
	$(MAKE) -C kernel clean
	$(RM) bin/* lib/*
