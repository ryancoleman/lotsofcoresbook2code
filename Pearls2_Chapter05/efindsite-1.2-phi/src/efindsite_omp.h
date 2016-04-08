/*
===============================================================================
         ___________.__            .____________.__  __          
     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
        \/    \/            \/      \/       \/               \/ 

                                                  
   eFindSite - ligand-binding site prediction from meta-threading

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org wfeinstein@lsu.edu

   Copyright 2013 Michal Brylinski, Wei Feinstein

   This file is part of eFindSite.

   eFindSite is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eFindSite is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#ifndef __EFINDSITE_OMP_H_
#define __EFINDSITE_OMP_H_

#include<cstring>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<list>
#include<map>
#include<algorithm>
#include<omp.h>

#ifdef __linux
#include<dlfcn.h>
#endif

#include "size.h"
#include "walltime.h"
#include "runsvm.h"
#include "coords.h"
#include "data.h"
#include "target.h"
#include "cmps.h"
#include "list.h"
#include "template.h"
#include "nwalign.h"
#include "distance.h"
#include "cluster.h"
#include "refine.h"
#include "pocket.h"
#include "tanimoto.h"
#include "frtmalign_omp.h"

#endif
