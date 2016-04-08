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


#ifndef __EFINDSITE_SCR_H_
#define __EFINDSITE_SCR_H_

#include<fstream>
#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cstdlib>
#include<vector>
#include<list>
#include<map>
#include<set>
#include<bitset>
#include<ctime>
#include<cmath>
#include<algorithm>
#include<gzstream.h>

#include "size.h"
#include "tanimoto.h"
#include "walltime.h"
#include "runsvm.h"

using namespace std;

struct cmp1
{
 double         weight;
 bitset<MAXSMI> fingerprint;
};

struct cmp2
{
 double         weight;
 bitset<MAXMAC> fingerprint;
};

struct cmp3
{
 double mw;
 double logp;
 double psa;
 
 int hbd;
 int hba;
 
 bitset<MAXSMI> fpt1;
 bitset<MAXMAC> fpt2;
 
 double sco_tst;
 double sco_tsa;
 double sco_tsc;
 
 double sco_tmt;
 double sco_tma;
 double sco_tmc;
 
 double sco_mw;
 double sco_logp;
 double sco_psa;
 
 double sco_hbd;
 double sco_hba;
 
 double ran_tst;
 double ran_tsa;
 double ran_tsc;
 
 double ran_tmt;
 double ran_tma;
 double ran_tmc;
 
 double dat_sum;
 double dat_max;
 double dat_min;
 
 double sco_svm;
};

#endif
