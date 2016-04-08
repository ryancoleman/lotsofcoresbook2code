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
#   Report bugs to michal@brylinski.org
#
#   Copyright 2013 Michal Brylinski
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


EXE = efindsite efindsite_scr efindsite_map

OBJ_F = cluster.o cmps.o coords.o distance.o data.o efindsite.o frtmalign.o list.o nwalign.o pocket.o refine.o runsvm.o tanimoto.o target.o template.o walltime.o

OBJ_S = efindsite_scr.o tanimoto.o walltime.o runsvm.o

SH = sh

CXX = g++

FC = gfortran

CPPFLAGS = -O2 -Wall -Wno-write-strings -I../lib -I.

FFLAGS = -O2

LDFLAGS = -lgzstream -lsvm -lz -lgfortran -lm -ldl -L../lib -L.

default: $(EXE)

all: $(EXE)

efindsite: $(OBJ_F)
	$(CXX) -o $@ $(OBJ_F) $(LDFLAGS)
	@mv efindsite ../bin/

efindsite_scr: $(OBJ_S)
	$(CXX) -o $@ $(OBJ_S) $(LDFLAGS)
	@mv efindsite_scr ../bin/

efindsite_map:
	$(SH) efindsite_map.shar
	@chmod +x efindsite_map
	@mv efindsite_map ../bin/

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

efindsite.o: efindsite.C
	$(CXX) $(CPPFLAGS) -c -o efindsite.o efindsite.C

efindsite_scr.o: efindsite_scr.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_scr.o efindsite_scr.C

frtmalign.o: frtmalign.f
	$(FC) $(FFLAGS) -c -o frtmalign.o frtmalign.f

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

clean:
	@(rm -f ${EXE} ../bin/efindsite ../bin/efindsite_scr ../bin/efindsite_map ${OBJ_F} ${OBJ_S})
