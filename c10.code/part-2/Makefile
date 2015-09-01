
# Options
ifndef nesting
nesting = manual
nesting_message = "nesting type not defined, assuming manual"
endif
ifeq ($(nesting), manual)
CFLAGS += -D MANUAL_NESTED
PROGRAMS = modal_man
endif
ifeq ($(nesting), omp)
CFLAGS += -D OMP_NESTED
PROGRAMS = modal_omp
endif
ifeq ($(nesting), teams)
CFLAGS += -D OMP_TEAMS
PROGRAMS = modal_teams
endif

# dirs

SRCDIR := Code
OBJDIR := Objects
BINDIR := .

# files

CFILESR =  local_barrier calculate_gamma gamma_mpi calculate_geometric calculate_int orthogonalisation_late basis_functions parameter arrays mtime io fixed_data
OBJFILESR = $(CFILESR:%=%.o) $(OBJDIR)/affinity.o

#PROGRAMS = rmpi

rmpi = main

# rules

all : $(PROGRAMS)

$(PROGRAMS): $(CFILESR:%=$(OBJDIR)/%.o) $(rmpi:%=$(OBJDIR)/%.o) $(OBJDIR)/affinity.o
	@echo ${nesting_message}
	$(CC) $(CFLAGS) -I$(SRCDIR) $^ -o $(BINDIR)/$@ $(LIBS) $(LIBSMPI)

$(OBJDIR)/affinity.o: $(SRCDIR)/affinity.cpp
	mpiicpc $(CFLAGS) -I$(SRCDIR) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I$(SRCDIR) -c $< -o $@ 

$(OBJDIR)/%.o: $(SRCDIR)/global.h

clean:
	rm -f $(OBJDIR)/*.o 
veryclean:
	rm -f $(OBJDIR)/*.o modal_man modal_nested modal_teams
