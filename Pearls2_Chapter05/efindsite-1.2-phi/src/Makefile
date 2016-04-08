#===============================================================================
#         ___________.__            .____________.__  __          
#     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
#   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
#   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
#    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
#        \/    \/            \/      \/       \/               \/ 
#
#                                                  
#   eFindSite - ligand-binding site prediction from meta-threading
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#   Report bugs to michal@brylinski.org wfeinstein@lsu.edu
#
#   Copyright 2013 Michal Brylinski, Wei Feinstein
#
#   This file is part of eFindSite.
#
#   eFindSite is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eFindSite is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


EXE = efindsite_omp efindsite_mic efindsite_scr efindsite_map efindsite_phi

OBJ_OMP = cluster.o cmps.o coords.o distance.o data.o efindsite_omp.o frtmalign_omp.o list.o nwalign.o pocket.o refine.o runsvm.o tanimoto.o target.o template.o walltime.o

OBJ_MIC = cluster.o cmps.o coords.o distance.o data.o efindsite_mic.o frtmalign_mic.o list.o nwalign.o pocket.o refine.o runsvm.o tanimoto.o target.o template.o walltime.o

OBJ_SCR = efindsite_scr.o tanimoto.o walltime.o runsvm.o

SH = sh

CXX = icpc

FC = ifort

CPPFLAGS = -O2 -openmp -Wall -Wno-write-strings -I../lib -I. -I/usr/include

FFLAGS = -O2 -openmp

MICFLAGS = -offload-attribute-target=mic

OMPLDFLAGS = -lgzstream -lsvm -lz -lgfortran -openmp -lm -ldl -L../lib -L. -L/usr/lib

MICLDFLAGS = -lgzstream -lsvm -lz -lgfortran -openmp -lm -ldl -L../lib -L. -L/usr/lib -offload-option,mic,ld,"-z defs"

SCRLDFLAGS = -lgzstream -lsvm -lz -openmp -lm -L../lib -L. -L/usr/lib

default: $(EXE)

all: $(EXE)

efindsite_omp: $(OBJ_OMP)
	$(CXX) -o $@ $(OBJ_OMP) $(OMPLDFLAGS)
	@mv efindsite_omp ../bin/

efindsite_mic: $(OBJ_MIC)
	$(CXX) -o $@ $(OBJ_MIC) $(MICLDFLAGS)
	@mv efindsite_mic ../bin/

efindsite_scr: $(OBJ_SCR)
	$(CXX) -o $@ $(OBJ_SCR) $(SCRLDFLAGS)
	@mv efindsite_scr ../bin/

efindsite_map:
	$(SH) efindsite_map.shar
	@chmod +x efindsite_map
	@mv efindsite_map ../bin/

efindsite_phi:
	$(SH) efindsite_phi.shar
	@chmod +x efindsite_phi
	@mv efindsite_phi ../bin/

cluster.o: cluster.C
	$(CXX) $(CPPFLAGS) -c -o cluster.o cluster.C

cmps.o: cmps.C
	$(CXX) $(CPPFLAGS) -c -o cmps.o cmps.C

coords.o: coords.C
	$(CXX) $(CPPFLAGS) -c -o coords.o coords.C

distance.o: distance.C
	$(CXX) $(CPPFLAGS) -c -o distance.o distance.C

data.o: data.C
	$(CXX) $(CPPFLAGS) -c -o data.o data.C

efindsite_scr.o: efindsite_scr.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_scr.o efindsite_scr.C

list.o: list.C
	$(CXX) $(CPPFLAGS) -c -o list.o list.C

nwalign.o: nwalign.f
	$(FC) $(FFLAGS) -c -o nwalign.o nwalign.f

pocket.o: pocket.C
	$(CXX) $(CPPFLAGS) -c -o pocket.o pocket.C

refine.o: refine.C
	$(CXX) $(CPPFLAGS) -c -o refine.o refine.C

runsvm.o: runsvm.C
	$(CXX) $(CPPFLAGS) -c -o runsvm.o runsvm.C

target.o: target.C
	$(CXX) $(CPPFLAGS) -c -o target.o target.C

tanimoto.o: tanimoto.C
	$(CXX) $(CPPFLAGS) -c -o tanimoto.o tanimoto.C

template.o: template.C
	$(CXX) $(CPPFLAGS) -c -o template.o template.C

walltime.o: walltime.C
	$(CXX) $(CPPFLAGS) -c -o walltime.o walltime.C

efindsite_omp.o: efindsite_omp.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_omp.o efindsite_omp.C

efindsite_mic.o: efindsite_mic.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_mic.o efindsite_mic.C

frtmalign_omp.o: frtmalign_omp.f
	$(FC) $(FFLAGS) -c -o frtmalign_omp.o frtmalign_omp.f

frtmalign_mic.o: frtmalign_mic.f
	$(FC) $(FFLAGS) ${MICFLAGS} -c -o frtmalign_mic.o frtmalign_mic.f

clean:
	@(rm -f ${EXE} *MIC.o ../bin/efindsite_omp ../bin/efindsite_mic ../bin/efindsite_scr ../bin/efindsite_map ../bin/efindsite_phi ${OBJ_OMP} ${OBJ_MIC} ${OBJ_SCR})
