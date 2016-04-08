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


#ifndef __RUNSVM_H_
#define __RUNSVM_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<sstream>
#include<cstdio>
#include<cctype>
#include<string>

#include "size.h"
#include "svm.h"

using namespace std;

class ModelSVM {
        
  private:
    
    bool _comp_loaded; /* aa composition flag */
    bool _prof_loaded; /* aa profiles flag */
    bool _bres_loaded; /* binding residues flag */
    bool _alig_loaded; /* aux ligand flag */
    bool _ran7_loaded; /* pocket ranking flag */
    bool _ran8_loaded; /* pocket ranking flag */
    bool _scr1_loaded; /* screening conf 1% */
    bool _scr2_loaded; /* screening conf 10% */
    bool _scor_loaded; /* svm scoring for screen */
    
    int _efs_attr;
    
    struct svm_model * _efs_model_comp;
    struct svm_node * _efs_node_comp;
    double _efs_scale_comp[MAXSV1][2];
    
    struct svm_model * _efs_model_prof;
    struct svm_node * _efs_node_prof;
    double _efs_scale_prof[MAXSV2][2];
    
    struct svm_model * _efs_model_bres;
    struct svm_node * _efs_node_bres;
    double _efs_scale_bres[MAXSV3][2];
    
    struct svm_model * _efs_model_alig;
    struct svm_node * _efs_node_alig;
    double _efs_scale_alig[MAXSV4][2];
    
    struct svm_model * _efs_model_ran7;
    struct svm_node * _efs_node_ran7;
    double _efs_scale_ran7[MAXSV5][2];
    
    struct svm_model * _efs_model_ran8;
    struct svm_node * _efs_node_ran8;
    double _efs_scale_ran8[MAXSV6][2];
    
    struct svm_model * _efs_model_scr1;
    struct svm_node * _efs_node_scr1;
    double _efs_scale_scr1[MAXSV7][2];
    
    struct svm_model * _efs_model_scr2;
    struct svm_node * _efs_node_scr2;
    double _efs_scale_scr2[MAXSV7][2];
    
    struct svm_model * _efs_model_scor;
    struct svm_node * _efs_node_scor;
    double _efs_scale_scor[MAXSV8][2];
    
  public:
    
    ModelSVM( bool, bool, bool, bool, bool, bool, bool, bool, bool );
    
    ModelSVM( void );
    
    ~ModelSVM();
    
    void loadModel( int, std::string );
    
    void loadScale( int, std::string );
    
    double SVMpredict( int, double [] );
};

#endif
