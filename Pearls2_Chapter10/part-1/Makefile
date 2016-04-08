
ifndef offload
CFLAGS += -no-offload
endif

ifndef version
version = 1
version_message = "version not specified, assuming version 1"
endif

# different version of the code

ifeq ($(version), 1)
# Options
PROGRAMS = modal_v1
# dirs
SRCDIR := Code-V1
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif

ifeq ($(version), 2)
# Options
PROGRAMS = modal_v2
# dirs
SRCDIR := Code-V2
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif

ifeq ($(version), 3)
# Options
PROGRAMS = modal_v3
# dirs
SRCDIR := Code-V3
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif
ifeq ($(version), 4)
# Options
PROGRAMS = modal_v4
# dirs
SRCDIR := Code-V4
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif

ifeq ($(version), 5)
# Options
PROGRAMS = modal_v5
# dirs
SRCDIR := Code-V5
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif

ifeq ($(version), 6)
# Options
PROGRAMS = modal_v6
# dirs
SRCDIR := Code-V6
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif
ifeq ($(version), 7)
# Options
PROGRAMS = modal_v7
# dirs
SRCDIR := Code-V7
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif
ifeq ($(version), 8)
# Options
PROGRAMS = modal_v8
# dirs
SRCDIR := Code-V8
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif
ifeq ($(version), 9)
# Options
PROGRAMS = modal_v9
# dirs
SRCDIR := Code-V9
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif
ifeq ($(version), 10)
# Options
PROGRAMS = modal_v10
# dirs
SRCDIR := Code-V10
OBJDIR := Objects
BINDIR := .
# files
CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o
endif
#PROGRAMS = rmpi

rmpi = main

# rules

all : $(PROGRAMS)

$(PROGRAMS): $(CFILESR:%=$(OBJDIR)/%.o) $(rmpi:%=$(OBJDIR)/%.o)
	$(CC) $(CFLAGS) -I$(SRCDIR) $^ -o $(BINDIR)/$@ $(LIBS) $(LIBSMPI)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I$(SRCDIR) -c $< -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/global.h

.PHONY: clean
.PHONY: veryclean

clean:
	rm -f $(OBJDIR)/*.o

veryclean:
	rm -f $(OBJDIR)/*.o modal_v1 modal_v2 modal_v3 modal_v4 modal_v5 modal_v6 modal_v7 modal_v8 modal_v9 modal_v10 
