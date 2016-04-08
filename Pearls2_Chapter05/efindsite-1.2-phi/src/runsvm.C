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


#include "runsvm.h"

using namespace std;

ModelSVM::ModelSVM( bool aa, bool ab, bool ac, bool ad, bool ae, bool af, bool ag, bool ah, bool ai )
{
 _comp_loaded = aa;
 _prof_loaded = ab;
 _bres_loaded = ac;
 _alig_loaded = ad;
 _ran7_loaded = ae;
 _ran8_loaded = af;
 _scr1_loaded = ag;
 _scr2_loaded = ah;
 _scor_loaded = ai;
 
 _efs_attr = 64;
 
 _efs_node_comp = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_prof = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_bres = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_alig = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_ran7 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_ran8 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_scr1 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_scr2 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_scor = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
}

ModelSVM::ModelSVM( void )
{
 _comp_loaded = false;
 _prof_loaded = false;
 _bres_loaded = false;
 _alig_loaded = false;
 _ran7_loaded = false;
 _ran8_loaded = false;
 _scr1_loaded = false;
 _scr2_loaded = false;
 _scor_loaded = false;
 
 _efs_attr = 64;
 
 _efs_node_comp = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_prof = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_bres = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_alig = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_ran7 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_ran8 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_scr1 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_scr2 = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
 _efs_node_scor = ( struct svm_node * ) malloc( _efs_attr * sizeof( struct svm_node ) );
}

ModelSVM::~ModelSVM() {}


// ==================================================================================   loadModel

void ModelSVM::loadModel( int m1_num, std::string m1_name )
{
 switch ( m1_num )
 {
  case 1: _efs_model_comp = svm_load_model( m1_name.c_str() ); break;
  case 2: _efs_model_prof = svm_load_model( m1_name.c_str() ); break;
  case 3: _efs_model_bres = svm_load_model( m1_name.c_str() ); break;
  case 4: _efs_model_alig = svm_load_model( m1_name.c_str() ); break;
  case 5: _efs_model_ran7 = svm_load_model( m1_name.c_str() ); break;
  case 6: _efs_model_ran8 = svm_load_model( m1_name.c_str() ); break;
  case 7: _efs_model_scr1 = svm_load_model( m1_name.c_str() ); break;
  case 8: _efs_model_scr2 = svm_load_model( m1_name.c_str() ); break;
  case 9: _efs_model_scor = svm_load_model( m1_name.c_str() ); break;
 }
}


// ==================================================================================   loadScale

void ModelSVM::loadScale( int s1_num, std::string s1_name )
{
 string line1, line2;
 
 ifstream s1_file( s1_name.c_str() );
 
 int i1 = 0;
 
 while (getline(s1_file,line1))
  if ( i1++ > 1 )
  {
   istringstream iss(line1, istringstream::in);
   
   int i2 = 0;
   int i3 = 0;
   double i4 = 0.0;
   double i5 = 0.0;
   
   while( iss >> line2 )     
   {
         if ( i2 == 0 ) { i3 = atoi(line2.c_str()) - 1; }
    else if ( i2 == 1 ) { i4 = atof(line2.c_str()); }
    else if ( i2 == 2 ) { i5 = atof(line2.c_str()); }
    
    i2++;
   }
   
   switch ( s1_num )
   {
    case 1: _efs_scale_comp[i3][0] = i4; _efs_scale_comp[i3][1] = i5; break;
    case 2: _efs_scale_prof[i3][0] = i4; _efs_scale_prof[i3][1] = i5; break;
    case 3: _efs_scale_bres[i3][0] = i4; _efs_scale_bres[i3][1] = i5; break;
    case 4: _efs_scale_alig[i3][0] = i4; _efs_scale_alig[i3][1] = i5; break;
    case 5: _efs_scale_ran7[i3][0] = i4; _efs_scale_ran7[i3][1] = i5; break;
    case 6: _efs_scale_ran8[i3][0] = i4; _efs_scale_ran8[i3][1] = i5; break;
    case 7: _efs_scale_scr1[i3][0] = i4; _efs_scale_scr1[i3][1] = i5; break;
    case 8: _efs_scale_scr2[i3][0] = i4; _efs_scale_scr2[i3][1] = i5; break;
    case 9: _efs_scale_scor[i3][0] = i4; _efs_scale_scor[i3][1] = i5; break;
   }
  }
 
 s1_file.close();
}


// ==================================================================================   SVMpredict

double ModelSVM::SVMpredict( int m1_num, double m1_fet [] )
{
 int m2; m2 = 0;
 
 switch ( m1_num )
 {
  case 1: m2 = MAXSV1; break;
  case 2: m2 = MAXSV2; break;
  case 3: m2 = MAXSV3; break;
  case 4: m2 = MAXSV4; break;
  case 5: m2 = MAXSV5; break;
  case 6: m2 = MAXSV6; break;
  case 7: m2 = MAXSV7; break;
  case 8: m2 = MAXSV7; break;
  case 9: m2 = MAXSV8; break;
 }
 
 for ( int m3 = 0; m3 < m2; m3++ )
 {
  double lb; lb = 0.0;
  double ub; ub = 0.0;
  
  switch ( m1_num )
  {
   case 1: lb = _efs_scale_comp[m3][0]; ub = _efs_scale_comp[m3][1]; break;
   case 2: lb = _efs_scale_prof[m3][0]; ub = _efs_scale_prof[m3][1]; break;
   case 3: lb = _efs_scale_bres[m3][0]; ub = _efs_scale_bres[m3][1]; break;
   case 4: lb = _efs_scale_alig[m3][0]; ub = _efs_scale_alig[m3][1]; break;
   case 5: lb = _efs_scale_ran7[m3][0]; ub = _efs_scale_ran7[m3][1]; break;
   case 6: lb = _efs_scale_ran8[m3][0]; ub = _efs_scale_ran8[m3][1]; break;
   case 7: lb = _efs_scale_scr1[m3][0]; ub = _efs_scale_scr1[m3][1]; break;
   case 8: lb = _efs_scale_scr2[m3][0]; ub = _efs_scale_scr2[m3][1]; break;
   case 9: lb = _efs_scale_scor[m3][0]; ub = _efs_scale_scor[m3][1]; break;
  }
  
  m1_fet[m3] = ( ( m1_fet[m3] - lb ) / ( ub - lb ) ) * 2.0 - 1.0;
  
  if ( m1_fet[m3] < -1.0 )
   m1_fet[m3] = -1.0;
  if ( m1_fet[m3] > 1.0 )
   m1_fet[m3] = 1.0;
 }
 
 int nr_class = 0;
 
 double *prob_estimates = NULL;
 
 int *labels = ( int * ) malloc( nr_class * sizeof( int ) );
 
 switch ( m1_num )
 {
  case 1: nr_class = svm_get_nr_class(_efs_model_comp); svm_get_labels(_efs_model_comp,labels); break;
  case 2: nr_class = svm_get_nr_class(_efs_model_prof); svm_get_labels(_efs_model_prof,labels); break;
  case 3: nr_class = svm_get_nr_class(_efs_model_bres); svm_get_labels(_efs_model_bres,labels); break;
  case 4: nr_class = svm_get_nr_class(_efs_model_alig); svm_get_labels(_efs_model_alig,labels); break;
  case 5: nr_class = svm_get_nr_class(_efs_model_ran7); svm_get_labels(_efs_model_ran7,labels); break;
  case 6: nr_class = svm_get_nr_class(_efs_model_ran8); svm_get_labels(_efs_model_ran8,labels); break;
  case 7: nr_class = svm_get_nr_class(_efs_model_scr1); svm_get_labels(_efs_model_scr1,labels); break;
  case 8: nr_class = svm_get_nr_class(_efs_model_scr2); svm_get_labels(_efs_model_scr2,labels); break;
  case 9: nr_class = svm_get_nr_class(_efs_model_scor); svm_get_labels(_efs_model_scor,labels); break;
 }
 
 prob_estimates = ( double * ) malloc( nr_class * sizeof( double ) );
 
 int r1 = 0;
 
 if ( labels[1] )
  r1 = 1;
 
 free(labels);
 
 for ( int i = 0; i < m2; i++ )
 {
  switch ( m1_num )
  {
   case 1: _efs_node_comp[i].index = i + 1; _efs_node_comp[i].value = m1_fet[i]; break;
   case 2: _efs_node_prof[i].index = i + 1; _efs_node_prof[i].value = m1_fet[i]; break;
   case 3: _efs_node_bres[i].index = i + 1; _efs_node_bres[i].value = m1_fet[i]; break;
   case 4: _efs_node_alig[i].index = i + 1; _efs_node_alig[i].value = m1_fet[i]; break;
   case 5: _efs_node_ran7[i].index = i + 1; _efs_node_ran7[i].value = m1_fet[i]; break;
   case 6: _efs_node_ran8[i].index = i + 1; _efs_node_ran8[i].value = m1_fet[i]; break;
   case 7: _efs_node_scr1[i].index = i + 1; _efs_node_scr1[i].value = m1_fet[i]; break;
   case 8: _efs_node_scr2[i].index = i + 1; _efs_node_scr2[i].value = m1_fet[i]; break;
   case 9: _efs_node_scor[i].index = i + 1; _efs_node_scor[i].value = m1_fet[i]; break;
  }
 }
 
 double p1 = 0.0;
 double p2 = 0.0;
 
 switch ( m1_num )
 {
  case 1: _efs_node_comp[m2].index = -1; p1 = svm_predict_probability(_efs_model_comp,_efs_node_comp,prob_estimates); p2 = prob_estimates[r1]; break;
  case 2: _efs_node_prof[m2].index = -1; p1 = svm_predict_probability(_efs_model_prof,_efs_node_prof,prob_estimates); p2 = prob_estimates[r1]; break;
  case 3: _efs_node_bres[m2].index = -1; p1 = svm_predict_probability(_efs_model_bres,_efs_node_bres,prob_estimates); p2 = prob_estimates[r1]; break;
  case 4: _efs_node_alig[m2].index = -1; p1 = svm_predict_probability(_efs_model_alig,_efs_node_alig,prob_estimates); p2 = prob_estimates[r1]; break;
  case 5: _efs_node_ran7[m2].index = -1; p1 = svm_predict_probability(_efs_model_ran7,_efs_node_ran7,prob_estimates); p2 = prob_estimates[r1]; break;
  case 6: _efs_node_ran8[m2].index = -1; p1 = svm_predict_probability(_efs_model_ran8,_efs_node_ran8,prob_estimates); p2 = prob_estimates[r1]; break;
  case 7: _efs_node_scr1[m2].index = -1; p1 = svm_predict_probability(_efs_model_scr1,_efs_node_scr1,prob_estimates); p2 = prob_estimates[r1]; break;
  case 8: _efs_node_scr2[m2].index = -1; p1 = svm_predict_probability(_efs_model_scr2,_efs_node_scr2,prob_estimates); p2 = prob_estimates[r1]; break;
  case 9: _efs_node_scor[m2].index = -1; p1 = svm_predict_probability(_efs_model_scor,_efs_node_scor,prob_estimates); p2 = prob_estimates[r1]; break;
 }
 
 free(prob_estimates);
 
 if ( p2 < 0.0 )
  p2 = 0.0;
 if ( p2 > 1.0 )
  p2 = 1.0;
 
 return p2;
}
